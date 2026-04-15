#R code for hypothesis test-based sample size calculation.

# Illustrative exact one-sample binomial design table
# H0: p <= p0 vs H1: p > p1
# One-sided alpha = 0.05, target power = 0.80, p1 = 0.30

alpha_target <- 0.05
power_target <- 0.80
p1 <- 0.30
p0_list <- c(0.05, 0.10, 0.15, 0.20)

# Exact type I error for rejection rule X >= r under null p0
alpha_exact <- function(n, p0, r) {
  pbinom(r - 1, size = n, prob = p0, lower.tail = FALSE)
}

# Exact power for rejection rule X >= r under alternative p1
power_exact <- function(n, p1, r) {
  pbinom(r - 1, size = n, prob = p1, lower.tail = FALSE)
}

# For a given n, find the rejection rule r that satisfies alpha <= alpha_target
# and gives the largest power among such rules
find_best_r <- function(n, p0, p1, alpha_target = 0.05) {
  candidates <- lapply(0:n, function(r) {
    a <- alpha_exact(n, p0, r)
    pw <- power_exact(n, p1, r)
    data.frame(r = r, alpha = a, power = pw)
  })
  candidates <- do.call(rbind, candidates)
  ok <- subset(candidates, alpha <= alpha_target)
  if (nrow(ok) == 0) return(NULL)
  ok[which.max(ok$power), ]
}

# Find the minimum n that satisfies both alpha and power constraints
find_min_n <- function(p0, p1, alpha_target = 0.05, power_target = 0.80,
                       n_max = 500) {
  for (n in 1:n_max) {
    best <- find_best_r(n, p0, p1, alpha_target)
    if (!is.null(best) && best$power >= power_target) {
      return(data.frame(
        p0 = p0,
        p1 = p1,
        rejection_rule = paste0("Reject if \u2265 ", best$r, " pCRs"),
        min_n = n,
        type_I_error = best$alpha,
        power = best$power,
        n26_sufficient = ifelse(26 >= n, "Yes", "No")
      ))
    }
  }
  return(NULL)
}

# Build the table
result_list <- lapply(p0_list, function(p0) {
  find_min_n(p0 = p0, p1 = p1,
             alpha_target = alpha_target,
             power_target = power_target)
})

result <- do.call(rbind, result_list)

# Nicely formatted version
result$type_I_error <- round(result$type_I_error, 3)
result$power <- round(result$power, 3)

print(result, row.names = FALSE)
