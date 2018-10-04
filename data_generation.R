#' ARB Data Generation Functions
#' 
#' Functions should be named with genXXXX, where XXXX is the version number.
#' 
#' Document the parameters and changes each makes to the structure, and the 
#' model numbers the generators work with.
#' 
#' Generic functions should not be versioned.


load_mira_phyloseq <- function(css_normalize=FALSE, rarefy_depth=1000) {
  suppressMessages({
    mira.all <- load_mira_data(
      seqtab_fp = here("../shared/data/seqtab.rds"),
      taxa_fp = here("../shared/data/taxa.rds"),
      meds_fp = here("../shared/data/MIRA_Medications_Table.csv"),
      specimen_fp = here("../shared/data/MIRA_Specimen_Table.csv")
    )
  })
  
  seqs <- mira.all$seqs
  meds <- mira.all$meds
  
  mira <- mira.all$ps
  # Remove non-MIRA subjects from dataset (incl. blanks and controls)
  .nonmira <- is.na(sample_data(mira)$subject_id)
  print(sprintf("Removing %d non-MIRA samples...", sum(.nonmira)))
  mira <- prune_samples(!.nonmira, mira)
  # Remove culture samples (they break the subject/type/date unique constraint)
  .culture <- grepl("Culture", sample_data(mira)$specimen_type)
  print(sprintf("Removing %d culture samples...", sum(.culture)))
  mira <- prune_samples(!.culture, mira)
  # Remove empty samples
  .empty <- sample_sums(mira) == 0
  print(sprintf("Removing %d empty samples...", sum(.empty)))
  mira <- prune_samples(!.empty, mira)
  
  # Identify "duplicated" specimens (same subject, specimen type, and study day)
  sample_data(mira) <- sample_data(mira) %>% 
    group_by(subject_id, specimen_type, study_day) %>%
    # "specimen_id3" has the same value for duplicated specimens so phyloseq can 
    # use it as a grouping level
    mutate(specimen_id3 = case_when(
      n() > 1 ~ paste0(first(as.character(specimen_id2)), "_D"),
      TRUE ~ as.character(specimen_id2)
    )) %>%
    ungroup() %>% as.data.frame() %>%
    set_rownames(.$specimen_id2)
  
  # Sum abundances and merge sample table for duplicates
  mira <- phyloseq::merge_samples(mira, "specimen_id3", fun=sum)
  # Re-add the relevant metadata since merge_samples mangles factors and dates
  sample_data(mira) <- sample_data(mira) %>% 
    mutate(specimen_id2 = rownames(.)) %>%
    select(specimen_id2, specimen_id3) %>%
    left_join(mira.all$samples) %>%
    ungroup() %>% as.data.frame() %>%
    set_rownames(.$specimen_id2)
  
  # Restrict to only samples for which we have abx data
  .abx_specimens <- as.character(inner_join(sample_data(mira), meds)$specimen_id2)
  mira.abx <- prune_samples(
    sample_data(mira)$specimen_id2 %in% .abx_specimens, mira)
  
  # Normalize if requested
  if (css_normalize) {
    gt1_feature_samples <- rowSums(otu_table(mira) > 0) > 1
    message(sprintf("removing %d samples for having < 2 features", length(gt1_feature_samples)-sum(gt1_feature_samples)))
    mira.abx.gt1 <- prune_samples(names(which(gt1_feature_samples)), mira.abx)
    normalized <- t(metagenomeSeq::cumNormMat(t(mira.abx.gt1@otu_table@.Data)))
    otu_table(mira.abx.gt1) <- otu_table(normalized, taxa_are_rows = F)
    mira.abx <- mira.abx.gt1
  }
  
  # Rarefy if requested
  if (rarefy_depth > 0) {
    mira.abx <- rarefy_even_depth(mira.abx, sample.size = rarefy_depth, rngseed = 1988)
  }
  
  return(list(mira.abx=mira.abx, meds=meds))
}

prepare_data <- function(mira.abx, meds, max_lag) {
  
  # Convert to melted form
  agg <- phyloseq_to_agglomerated(mira.abx, "specimen_id2", "otu_id", "read_count") %>%
    group_by(specimen_id2) %>%
    mutate(total_reads = sum(read_count)) %>%
    ungroup()
  
  # Clean up subject/timepoint data
  subjects <- sample_data(mira.abx) %>%
    group_by(subject_id) %>% 
    mutate(exit_date = max(collection_date)) %>%
    distinct(subject_id, enroll_date, exit_date) %>%
    right_join(meds) %>%
    group_by(subject_id) %>%
    mutate(study_day = as.integer(collection_date - enroll_date)) %>%
    mutate(exit_day = study_day[collection_date == exit_date]) %>%
    # Limit to only a week before enrollment date and nothing after
    filter(study_day > -7, collection_date <= exit_date)
  
  # Manually split MIRA_024 into two sub-subjects
  subjects <- subjects %>%
    mutate(subject_id2 = case_when(
      subject_id != "MIRA_024" ~ subject_id,
      study_day <= 33 ~ "MIRA_024a",
      study_day >= 73 ~ "MIRA_024b",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(subject_id2)) %>%
    mutate(study_day = case_when(
      subject_id2 == "MIRA_024b" ~ study_day - 73,
      TRUE ~ as.double(study_day)
    )) %>%
    mutate(exit_day = ifelse(subject_id2 == "MIRA_024b", exit_day-73, exit_day)) %>%
    ungroup() %>%
    mutate(subject_id=subject_id2) %>%
    select(-subject_id2)
  
  # Clean up antibiotics data
  subjects <- subjects %>%
    ungroup() %>%
    separate_rows(abx_b) %>%
    # Handle coformulations
    mutate(abx_b2 = case_when(
      abx_b == "clavulanate" ~ "amoxicillin.clavulanate",
      abx_b == "sulbactam" ~ "ampicillin.sulbactam",
      abx_b == "tazobactam" ~ "piperacillin.tazobactam",
      abx_b == "trimethoprim" ~ "trimethoprim.sulfamethoxazole",
      TRUE ~ abx_b)) %>%
    group_by(subject_id, collection_date) %>%
    # Remove remaining part of coformulation if coformulation is present
    filter(
      (abx_b2 != "amoxicillin" & ("amoxicillin.clavulanate" %in% abx_b2)) | 
        !("amoxicillin.clavulanate" %in% abx_b2) |
        is.na(abx_b2)
    ) %>%
    filter(
      (abx_b2 != "ampicillin" & ("ampicillin.sulbactam" %in% abx_b2)) |
        !("ampicillin.sulbactam" %in% abx_b2) |
        is.na(abx_b2)
    ) %>%
    filter(
      (abx_b2 != "sulfamethoxazole" & ("trimethoprim.sulfamethoxazole" %in% abx_b2)) |
        !("trimethoprim.sulfamethoxazole" %in% abx_b2) |
        is.na(abx_b2)
    ) %>%
    ungroup() %>%
    filter((abx_b2 != "piperacillin") | is.na(abx_b2)) %>%
    mutate(abx_b = abx_b2) %>%
    select(-c(abx_b2, enroll_date, exit_date, exit_day, collection_date, abx_b_number)) %>%
    group_by(subject_id, study_day) %>%
    nest(.key="abx")
  
  # Add lagged values of antibiotics up to defined `max_lag` days in the past
  subjects <- subjects %>%
    select(subject_id, study_day, abx) %>%  # Redundant for now but the following is column-order sensitive
    arrange(subject_id, study_day) %>%
    group_by(subject_id) %>%
    nest() %>%
    # Here we're mapping a tibble created by adply that contains a list-column for each day to lag
    # back to the 'lagged' column. This operates subject-by-subject (see group_by above)
    mutate(lagged = map(data, function(df) {
      as_tibble(plyr::alply(1:max_lag, 1, function(n) lag(df$abx,  default = list(tibble(abx_b=NA)), n=n)))
    })) %>% unnest() %>%
    # Columns besides subject_id and study_day reflect the number of days lagged
    set_colnames(c("subject_id", "study_day", "0", as.character(1:max_lag))) %>%
    gather(key="day_lag", value="abx_lagged", -c(subject_id, study_day)) %>%
    mutate(day_lag = as.numeric(day_lag)) %>% 
    unnest()

  subjects <- subjects %>%
    group_by(subject_id, study_day, abx_b) %>%
    filter(!is.na(abx_b)) %>%
    # 'gap' is a T/F column where T if more than a 1-day gap since previous
    mutate(gap = c(FALSE, diff(day_lag) > 1)) %>% 
    # We remove any rows at or above the gap since we're only interested in days adjacent to this day
    filter(day_lag < min(day_lag[gap])) %>%
    # 'continuous' is the number of non-gap days. We zero it out if it does not include today (day_lag==0)
    summarize(continuous = (sum(!gap)) * 0 %in% day_lag) %>%
    # any values of 'continuous' that are zero are irrelevant and are removed
    filter(continuous > 0) %>%
    right_join(subjects) %>%
    group_by(subject_id, study_day) %>%
    nest(.key="abx_lagged")
  
  left_join(subjects, agg) 
}

#' Frameshifts antibiotics by multiple days.
gen3.01 <- function(agg, .specimen_type, min_times_abx=2, min_median_abx=2, max_lag=3, ...) {

  d <- select(
    agg, otu_id, specimen_id2, specimen_type, subject_id, study_day, read_count, total_reads, abx_lagged) %>%
    filter(specimen_type == .specimen_type)
  d2 <- d %>% 
    group_by(specimen_id2) %>%
    filter(read_count > 0) %>%
    summarize(D1 = exp(-sum((read_count/total_reads) * log(read_count/total_reads))))
  d3 <- d %>%
    select(-c(otu_id, read_count, total_reads)) %>% 
    distinct(specimen_id2, .keep_all=T) %>%
    left_join(d2)
  
  # Add lag columns for previous day's proportions and antibiotics
  d4 <- d3 %>%
    arrange(subject_id, study_day) %>%
    group_by(subject_id) %>%
    mutate(D1_lag = lag(D1)) %>%
    mutate(D1_delta = D1/D1_lag) %>%
    filter(!is.na(D1_lag)) %>% 
    unnest() %>%
    filter(day_lag <= max_lag) %>%
    rename(specimens=specimen_id2, subjects=subject_id, abx_lag=abx_b) %>%
    mutate(abx_yn = ifelse(is.na(abx_lag), 0, 1)) %>% 
    mutate(abx_lag = ifelse(is.na(abx_lag), "none", abx_lag)) 
  
  abx_to_keep <- (d4 %>% group_by(abx_lag, subjects) %>%
                    filter(abx_lag != "none", day_lag==0) %>%
                    summarize(n=n()) %>% 
                    mutate(n_subjects=n()) %>%
                    filter(n_subjects >= min_times_abx, median(n) >= min_median_abx) %>%
                    distinct(abx_lag))$abx_lag
  
  subjects_to_keep <- (filter(d4, abx_lag %in% abx_to_keep) %>% distinct(subjects))$subjects 
  
  stopifnot(length(subjects_to_keep) > 0)
  # browser()
  d5 <- d4 %>% ungroup() %>%
    filter(subjects %in% subjects_to_keep)
  d_abx <- d5 %>%
    reshape2::dcast(specimens + subjects + D1_delta ~ abx_lag, value.var="abx_yn")
  abx <- reshape2::acast(d5, specimens ~ abx_lag ~ day_lag, value.var="abx_yn", fill=0)
  abx <- abx[match(d_abx$specimens, rownames(abx)), , ]
  stopifnot(all(rownames(abx) == d_abx$specimens))
  result <- tidybayes::compose_data(d_abx[, c(1:3)])
  # Remove abxs seen less than the minimum, and the NA column
  abx <- abx[, abx_to_keep, , drop=F]
  result$abx_b <- as.numeric(as.factor(colnames(abx)))
  result$n_abx <- ncol(abx)
  result$abx <- abx
  result$d5 <- d5
  result$d_abx <- d_abx[, c(1:3)]
  result$n_lag <- dim(abx)[3]  # include day 0
  result
}

#' Frameshifts antibiotics by multiple days.
#' Returns D1 and D1_lag in output list.
gen3.02 <- function(agg, .specimen_type, min_times_abx=2, min_median_abx=2, max_lag=3, ...) {
  
  d <- select(
    agg, otu_id, specimen_id2, specimen_type, subject_id, study_day, read_count, total_reads, abx_lagged) %>%
    filter(specimen_type == .specimen_type)
  d2 <- d %>% 
    group_by(specimen_id2) %>%
    filter(read_count > 0) %>%
    summarize(D1 = exp(-sum((read_count/total_reads) * log(read_count/total_reads))))
  d3 <- d %>%
    select(-c(otu_id, read_count, total_reads)) %>% 
    distinct(specimen_id2, .keep_all=T) %>%
    left_join(d2)
  
  # Add lag columns for previous day's proportions and antibiotics
  d4 <- d3 %>%
    arrange(subject_id, study_day) %>%
    group_by(subject_id) %>%
    mutate(D1_lag = lag(D1)) %>%
    mutate(D1_delta = D1/D1_lag) %>%
    filter(!is.na(D1_lag)) %>% 
    unnest() %>%
    filter(day_lag <= max_lag) %>%
    rename(specimens=specimen_id2, subjects=subject_id, abx_lag=abx_b) %>%
    mutate(abx_yn = ifelse(is.na(abx_lag), 0, 1)) %>% 
    mutate(abx_lag = ifelse(is.na(abx_lag), "none", abx_lag)) 
  
  abx_to_keep <- (d4 %>% group_by(abx_lag, subjects) %>%
                    filter(abx_lag != "none", day_lag==0) %>%
                    summarize(n=n()) %>% 
                    mutate(n_subjects=n()) %>%
                    filter(n_subjects >= min_times_abx, median(n) >= min_median_abx) %>%
                    distinct(abx_lag))$abx_lag
  
  subjects_to_keep <- (filter(d4, abx_lag %in% abx_to_keep) %>% distinct(subjects))$subjects 
  
  stopifnot(length(subjects_to_keep) > 0)
  # browser()
  d5 <- d4 %>% ungroup() %>%
    filter(subjects %in% subjects_to_keep)
  d_abx <- d5 %>%
    reshape2::dcast(specimens + subjects + D1 + D1_lag ~ abx_lag, value.var="abx_yn")
  abx <- reshape2::acast(d5, specimens ~ abx_lag ~ day_lag, value.var="abx_yn", fill=0)
  abx <- abx[match(d_abx$specimens, rownames(abx)), , ]
  stopifnot(all(rownames(abx) == d_abx$specimens))
  result <- tidybayes::compose_data(d_abx[, c(1:4)])
  # Remove abxs seen less than the minimum, and the NA column
  abx <- abx[, abx_to_keep, , drop=F]
  result$abx_b <- as.numeric(as.factor(colnames(abx)))
  result$n_abx <- ncol(abx)
  result$abx <- abx
  result$d5 <- d5
  result$d_abx <- d_abx[, c(1:4)]
  result$n_lag <- dim(abx)[3]  # include day 0
  result
}

#' Frameshifts antibiotics by multiple days.
#' Returns D1 and D1_lag in output list.
#' Returns days of continuous administration rather than abx_lag
gen3.03 <- function(agg, .specimen_type, min_times_abx=2, min_median_abx=2, max_continuous=5, ...) {
  d <- select(
    agg, otu_id, specimen_id2, specimen_type, subject_id, study_day, read_count, total_reads, abx_lagged) %>%
    filter(specimen_type == .specimen_type)
  d2 <- d %>% 
    group_by(specimen_id2) %>%
    filter(read_count > 0) %>%
    summarize(D1 = exp(-sum((read_count/total_reads) * log(read_count/total_reads))))
  d3 <- d %>%
    select(-c(otu_id, read_count, total_reads)) %>% 
    distinct(specimen_id2, .keep_all=T) %>%
    left_join(d2)
  
  # Add lag columns for previous day's proportions and antibiotics
  d4 <- d3 %>%
    arrange(subject_id, study_day) %>%
    group_by(subject_id) %>%
    mutate(D1_lag = lag(D1)) %>%
    mutate(D1_delta = D1/D1_lag) %>%
    filter(!is.na(D1_lag)) %>% 
    unnest() %>%
    filter(day_lag == 0) %>%
    rename(specimens=specimen_id2, subjects=subject_id, abx_lag=abx_b) %>%
    mutate(abx_yn = ifelse(is.na(abx_lag), 0, 1)) %>% 
    mutate(abx_lag = ifelse(is.na(abx_lag), "none", abx_lag)) %>%
    # Truncate continuous administration data if longer than max_continuous
    mutate(continuous = ifelse(continuous > max_continuous, max_continuous, continuous))
  
  abx_to_keep <- (d4 %>% group_by(abx_lag, subjects) %>%
                    filter(abx_lag != "none") %>%
                    summarize(n=n()) %>% 
                    mutate(n_subjects=n()) %>%
                    filter(n_subjects >= min_times_abx, median(n) >= min_median_abx) %>%
                    distinct(abx_lag))$abx_lag
  
  subjects_to_keep <- (filter(d4, abx_lag %in% abx_to_keep) %>% distinct(subjects))$subjects 
  
  stopifnot(length(subjects_to_keep) > 0)
  # browser()
  d5 <- d4 %>% ungroup() %>%
    filter(subjects %in% subjects_to_keep)
  d_abx <- d5 %>%
    reshape2::dcast(specimens + subjects + D1 + D1_lag ~ abx_lag, value.var="abx_yn")
  abx <- reshape2::acast(d5, specimens ~ abx_lag ~ continuous, value.var="abx_yn", fill=0)
  abx <- abx[match(d_abx$specimens, rownames(abx)), , ]
  stopifnot(all(rownames(abx) == d_abx$specimens))
  result <- tidybayes::compose_data(d_abx[, c(1:4)])
  # Remove abxs seen less than the minimum, and the NA column
  abx <- abx[, abx_to_keep, , drop=F]
  result$abx_b <- as.numeric(as.factor(colnames(abx)))
  result$n_abx <- ncol(abx)
  result$abx <- abx
  result$d5 <- d5
  result$d_abx <- d_abx[, c(1:4)]
  result$n_continuous <- dim(abx)[3]  # include day 0
  result
}

#' Frameshifts antibiotics by multiple days.
#' Returns D1 and D1_lag in output list.
#' Returns days of continuous administration rather than abx_lag
gen3.03 <- function(agg, .specimen_type, min_times_abx=2, min_median_abx=2, max_continuous=5, ...) {
  d <- select(
    agg, otu_id, specimen_id2, specimen_type, subject_id, study_day, read_count, total_reads, abx_lagged) %>%
    filter(specimen_type == .specimen_type)
  d2 <- d %>% 
    group_by(specimen_id2) %>%
    filter(read_count > 0) %>%
    summarize(D1 = exp(-sum((read_count/total_reads) * log(read_count/total_reads))))
  d3 <- d %>%
    select(-c(otu_id, read_count, total_reads)) %>% 
    distinct(specimen_id2, .keep_all=T) %>%
    left_join(d2)
  
  # Add lag columns for previous day's proportions and antibiotics
  d4 <- d3 %>%
    arrange(subject_id, study_day) %>%
    group_by(subject_id) %>%
    mutate(D1_lag = lag(D1)) %>%
    mutate(D1_delta = D1/D1_lag) %>%
    filter(!is.na(D1_lag)) %>% 
    unnest() %>%
    filter(day_lag == 0) %>%
    rename(specimens=specimen_id2, subjects=subject_id, abx_lag=abx_b) %>%
    mutate(abx_yn = ifelse(is.na(abx_lag), 0, 1)) %>% 
    mutate(abx_lag = ifelse(is.na(abx_lag), "none", abx_lag)) %>%
    # Truncate continuous administration data if longer than max_continuous
    mutate(continuous = ifelse(continuous > max_continuous, max_continuous, continuous))
  
  abx_to_keep <- (d4 %>% group_by(abx_lag, subjects) %>%
                    filter(abx_lag != "none") %>%
                    summarize(n=n()) %>% 
                    mutate(n_subjects=n()) %>%
                    filter(n_subjects >= min_times_abx, median(n) >= min_median_abx) %>%
                    distinct(abx_lag))$abx_lag
  
  subjects_to_keep <- (filter(d4, abx_lag %in% abx_to_keep) %>% distinct(subjects))$subjects 
  
  stopifnot(length(subjects_to_keep) > 0)
  # browser()
  d5 <- d4 %>% ungroup() %>%
    filter(subjects %in% subjects_to_keep)
  d_abx <- d5 %>%
    reshape2::dcast(specimens + subjects + D1 + D1_lag ~ abx_lag, value.var="abx_yn")
  abx <- reshape2::acast(d5, specimens ~ abx_lag ~ continuous, value.var="abx_yn", fill=0)
  abx <- abx[match(d_abx$specimens, rownames(abx)), , ]
  stopifnot(all(rownames(abx) == d_abx$specimens))
  result <- tidybayes::compose_data(d_abx[, c(1:4)])
  # Remove abxs seen less than the minimum, and the NA column
  abx <- abx[, abx_to_keep, , drop=F]
  result$abx_b <- as.numeric(as.factor(colnames(abx)))
  result$n_abx <- ncol(abx)
  result$abx <- abx
  result$d5 <- d5
  result$d_abx <- d_abx[, c(1:4)]
  result$n_continuous <- dim(abx)[3]  # include day 0
  result
}