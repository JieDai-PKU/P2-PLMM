#R code for error bound based sample size calculation.

# Target: smallest n such that the 95% exact binomial CI width
# is <= 0.40 when the anticipated response rate is 30%

target_p <- 0.30
target_width <- 0.40
conf_level <- 0.95

# Exact Clopper-Pearson CI width
cp_width <- function(n, p = target_p, conf.level = conf_level) {
  x <- round(n * p) # anticipated number of responses
  ci <- binom.test(x, n, conf.level = conf.level)$conf.int
  width <- ci[2] - ci[1]
  data.frame(
    n = n,
    x = x,
    p_hat = x / n,
    lower = ci[1],
    upper = ci[2],
    width = width
  )
}

# Search over candidate sample sizes
results <- do.call(rbind, lapply(15:40, cp_width))

# Show all candidates
print(results)

# Smallest n meeting the width criterion
subset(results, width <= target_width)[1, ]
