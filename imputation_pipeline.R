library(mice)
library(dplyr)
library(tidyr)

set.seed(42)

input_dir  <- "C:/Users/Martin/Desktop/Github/MICPR/CDC_data"
output_dir <- "C:/Users/Martin/Desktop/Github/MICPR/Result"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
setwd(output_dir)

years <- 2014:2019


# read the yearly files for one cause and stack them
read_yearly <- function(prefix) {
  out <- list()
  for (yr in years) {
    f <- paste0(prefix, yr, ".csv")
    cat("Loading:", f, "\n")
    df <- read.csv(file.path(input_dir, f), stringsAsFactors = FALSE, check.names = FALSE)
    df[] <- lapply(df, function(x) as.numeric(as.character(x)))  # "Suppressed" etc. -> NA
    df <- df[!is.na(df$`County Code`), ]
    df$Year <- yr
    out[[length(out) + 1]] <- df
  }
  bind_rows(out) %>% arrange(`County Code`, Year)
}


# stage 1: use the neighbouring years to fill suppressed cells
fill_from_neighbours <- function(d) {
  for (i in which(is.na(d$Deaths))) {
    cc <- d$`County Code`[i]
    yr <- d$Year[i]
    before <- d$Deaths[d$`County Code` == cc & d$Year == yr - 1][1]
    after  <- d$Deaths[d$`County Code` == cc & d$Year == yr + 1][1]
    
    if (!is.na(before) && before >= 10 && !is.na(after) && after >= 10) {
      d$Deaths[i] <- 9
    } else if ((!is.na(before) && before >= 10 && (is.na(after) || after == 0)) ||
               (!is.na(after) && after >= 10 && (is.na(before) || before == 0))) {
      known <- ifelse(!is.na(before) && before >= 10, before, after)
      d$Deaths[i] <- round(known / 2)
    }
  }
  d
}


# stage 2: impute the 6 year file with mice (pmm on log scale), rescale to 1-9
impute_6y <- function(file_name) {
  d6 <- read.csv(file.path(input_dir, file_name), stringsAsFactors = FALSE, check.names = FALSE)
  d6$Deaths        <- as.numeric(as.character(d6$Deaths))
  d6$Population    <- as.numeric(as.character(d6$Population))
  d6$`County Code` <- as.numeric(as.character(d6$`County Code`))
  d6 <- d6[!is.na(d6$`County Code`), ]
  
  dat  <- d6 %>% select(Deaths, Population, `County Code`)
  miss <- which(is.na(dat$Deaths))
  dat$Deaths_log <- log(dat$Deaths + 1)
  
  meth <- make.method(dat)
  meth["Deaths_log"] <- "pmm"
  pred <- make.predictorMatrix(dat)
  pred["Deaths_log", "Population"]  <- 1
  pred["Deaths_log", "County Code"] <- 0
  pred["Deaths_log", "Deaths"]      <- 0
  
  fit <- mice(dat, m = 5, maxit = 5, method = meth, predictorMatrix = pred, printFlag = FALSE)
  
  # back to counts and squeeze into 1-9
  for (j in 1:fit$m) {
    raw <- fit$imp$Deaths_log[, j]
    cnt <- exp(raw) - 1
    cnt[cnt < 0] <- 0
    ok <- cnt[is.finite(cnt)]
    
    if (length(ok) == 0) {
      sc <- rep(1, length(raw))
    } else {
      lo <- min(ok); hi <- max(ok)
      sc <- if (hi == lo) rep(1, length(raw)) else 1 + (cnt - lo) * (9 - 1) / (hi - lo)
    }
    fit$imp$Deaths_log[, j] <- round(pmax(1, pmin(9, sc)))
  }
  
  comp <- complete(fit, 1)
  d6$Deaths[miss] <- comp$Deaths_log[miss]
  d6
}


# one cause of death from start to finish
process_cause <- function(prefix, file_6y) {
  long <- read_yearly(prefix) %>% fill_from_neighbours()
  
  cat("Processing 6Y aggregate:", file_6y, "\n")
  d6 <- impute_6y(file_6y)
  long <- long %>%
    left_join(d6 %>% select(`County Code`, Deaths_6yr = Deaths), by = "County Code")
  
  # stage 3: whatever is still missing gets half of the 6 year value
  still_na <- which(is.na(long$Deaths))
  if (length(still_na) > 0) {
    half <- long$Deaths_6yr[still_na] / 2
    half[half < 0] <- 0
    long$Deaths[still_na] <- round(pmax(1, half))
  }
  long
}


cat("\n=== Mental & Behavioural Disorders ===\n")
mental <- process_cause("Mentalandbehaviouraldisorders3Y", "Mentalandbehaviouraldisorders6Y.csv")

cat("\n=== Suicides ===\n")
suicide <- process_cause("Suicide3Y", "Suicides6Y.csv")


# merge both outcomes and go wide
master_long <- full_join(
  mental  %>% select(`County Code`, Year, Deaths_Mental = Deaths, Population),
  suicide %>% select(`County Code`, Year, Deaths_Suicide = Deaths),
  by = c("County Code", "Year")
)

master_wide <- master_long %>%
  pivot_wider(
    id_cols     = `County Code`,
    names_from  = Year,
    values_from = c(Deaths_Mental, Deaths_Suicide, Population),
    names_sep   = "_"
  )

output_file <- file.path(output_dir, "Master_Imputed_Outcomes_Wide.csv")
write.csv(master_wide, file = output_file, row.names = FALSE)

cat("\nSaved:", output_file, "\n")
cat("Dimensions:", dim(master_wide), "\n")
