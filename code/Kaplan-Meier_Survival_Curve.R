library(tidyverse)
library(openxlsx)
library(factoextra)
library(cluster)
library(survival)
library(survminer)

Data <- openxlsx::read.xlsx("input.xlsx", sheet = 1, rowNames = TRUE)

#HR = high / low
Data$Group <- factor(Data$Group, levels = c("low", "high"))

#RFS
fitd <- survdiff(Surv(RFS_Month, RFS_Event) ~ Group, data = Data, na.action = na.exclude)

# Log-rank test
p.val <- 1 - pchisq(fitd$chisq, df = length(fitd$n) - 1)

# calculate Hazard Ratio (HR = high / low)
HR_val <- (fitd$obs[2] / fitd$exp[2]) / (fitd$obs[1] / fitd$exp[1])

# calculate 95% CI
se_logHR <- sqrt(1 / fitd$exp[2] + 1 / fitd$exp[1])
logHR <- log(HR_val)
up95  <- exp(logHR + qnorm(0.975) * se_logHR)
low95 <- exp(logHR - qnorm(0.975) * se_logHR)

HR_text <- paste("Hazard Ratio = ", round(HR_val, 2), sep = "")
CI_text <- paste("95% CI: ", round(low95, 2), " - ", round(up95, 2), sep = "")

# Kaplan-Meier Survival Curve
fit <- survfit(Surv(RFS_Month, RFS_Event) ~ Group, data = Data, conf.type = "log-log")

palette_colors <- c("#104E8B", "#8B0000") 
legend_labels <- c("Low", "High")

p <- ggsurvplot(
  fit,
  pval = paste(
    ifelse(p.val < 0.001, "Log-rank p < 0.001", paste("Log-rank p =", round(p.val, 3))),
    HR_text,
    CI_text,
    sep = "\n"
  ),
  linetype = "solid",
  palette = palette_colors,
  censor.size = 4.5,
  censor.shape = 124,
  xlab = "Time (months)",
  xlim = c(0, 50),
  ylab = "Recurrence-Free Survival",
  break.y.by = 0.1,
  break.x.by = 12,
  axes.offset = FALSE,
  surv.median.line = "hv",
  legend.title = "",
  legend = c(0.8, 0.95),
  legend.labs = legend_labels,
  risk.table.title = "No. at risk",
  font.legend = 8,
  conf.int = FALSE,
  risk.table = TRUE,
  ggtheme = theme_classic(),
  tables.theme = theme_cleantable()
)