library(readxl)

dat <- read_excel(
  "water_quality_compared.xlsx",
  sheet = "Table4_Comparison"
)

parameters <- c("EC", "pH", "DO", "Salinity",
                "Ammonia", "TDS", "Temperature")

# Predefined acceptable absolute differences
margin <- c(
  EC = 100,
  pH = 0.20,
  DO = 0.50,
  Salinity = 2.00,
  Ammonia = 1.00,
  TDS = 50.00,
  Temperature = 1.00
)

analyze_parameter <- function(parameter) {

  system <- as.numeric(dat[[paste0(parameter, " System")]])
  reference <- as.numeric(dat[[paste0(parameter, " Reference")]])

  valid <- complete.cases(system, reference)
  system <- system[valid]
  reference <- reference[valid]

  difference <- system - reference
  n <- length(difference)
  mean_difference <- mean(difference)
  sd_difference <- sd(difference)
  se <- sd_difference / sqrt(n)
  delta <- margin[parameter]

  # Paired t-test
  if (sd_difference == 0) {
    paired_p <- ifelse(mean_difference == 0, 1, 0)
    ci95 <- c(mean_difference, mean_difference)
  } else {
    test <- t.test(system, reference, paired = TRUE)
    paired_p <- test$p.value
    ci95 <- test$conf.int
  }

  # TOST equivalence test
  if (sd_difference == 0) {
    tost_p1 <- ifelse(mean_difference > -delta, 0, 1)
    tost_p2 <- ifelse(mean_difference <  delta, 0, 1)
    ci90 <- c(mean_difference, mean_difference)
  } else {
    tost_p1 <- 1 - pt((mean_difference + delta) / se, df = n - 1)
    tost_p2 <- pt((mean_difference - delta) / se, df = n - 1)
    ci90 <- mean_difference +
      c(-1, 1) * qt(0.95, df = n - 1) * se
  }

  tost_p <- max(tost_p1, tost_p2)

  # Bland–Altman limits of agreement
  ba_limits <- mean_difference + c(-1, 1) * 1.96 * sd_difference

  data.frame(
    Parameter = parameter,
    N = n,
    Mean_difference = mean_difference,
    Paired_p = paired_p,
    CI95_low = ci95[1],
    CI95_high = ci95[2],
    TOST_p = tost_p,
    TOST_CI90_low = ci90[1],
    TOST_CI90_high = ci90[2],
    BA_low = ba_limits[1],
    BA_high = ba_limits[2],
    No_significant_difference = paired_p >= 0.05,
    TOST_equivalent = tost_p < 0.05 &
      ci90[1] > -delta & ci90[2] < delta,
    BA_within_margin =
      abs(ba_limits[1]) <= delta &
      abs(ba_limits[2]) <= delta
  )
}

results <- do.call(
  rbind,
  lapply(parameters, analyze_parameter)
)

results$Paired_p_Holm <- p.adjust(results$Paired_p, method = "holm")
results$TOST_p_Holm <- p.adjust(results$TOST_p, method = "holm")

print(results)
write.csv(results, "water_quality_test_results.csv", row.names = FALSE)