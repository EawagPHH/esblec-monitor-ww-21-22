
library(ggplot2)
library(MASS)

# Your data
data <- c(100,
          100,
          100,
          100,
          2.66,
          7.57995758,
          2,
          23,
          2.142857143,
          50.74076338,
          0.49999125,
          26.11111111,
          10.41666667,
          8.518518519,
          1.6,
          14.99999991,
          0.183606557,
          46.36363636,
          0.116666667,
          0.25,
          0.928571429,
          0.176470588,
          0.54999695,
          1.000001572,
          0.272727273,
          4.174666667,
          0.17500065,
          57.14285595,
          41.84640523,
          54.71698113,
          0.023275864,
          0.18571428,
          2.314286,
          0.935897436,
          100,
          4.166666667,
          48.51485149,
          74.15730337,
          0.593023256,
          100,
          1.368932039,
          0.130769231,
          2.469879518,
          100,
          0.513513514,
          8.307692308,
          0.444833103,
          57.46336996,
          2.3125,
          0.001894737,
          7.631578947,
          0.06085977,
          0.033316049,
          3.557692308,
          0.004901961,
          0.435897436,
          0.083108091,
          5.097270229,
          0.120772947,
          0.032837472,
          1.757316641,
          10,
          0.081434263,
          0.003781484,
          0.083249093,
          0.016056984,
          0.869189907
  
)
# Create a dataframe from the data
df <- data.frame(value = data)

# Create empirical CDF data
df$ecdf <- ecdf(df$value)(df$value)

# Fit lognormal distribution to data
fit <- fitdistr(df$value, "lognormal")

# Generate fitted CDF values based on the lognormal parameters
df$lognorm_cdf <- plnorm(df$value, meanlog = fit$estimate[1], sdlog = fit$estimate[2])

# Calculate empirical arithmetic mean and median
mean_val <- mean(df$value)
median_val <- median(df$value)

ggplot(df, aes(x=value)) + 
  geom_line(aes(y=ecdf, color="Empirical CDF"), size=1) +
  geom_line(aes(y=lognorm_cdf, color="Lognormal CDF"), linetype="solid", size=0.5) +
  scale_x_log10(labels = scales::comma) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  ) +
  labs(
    x = expression(paste("Percentage of ESBL-", italic("E. coli"), " out of total ", italic("E. coli"), " in the gut of children in Bangladesh")),
    y = "CDF",
    title = paste("Lognormal Parameters: Meanlog =", round(fit$estimate[1], 3), 
                  "Sdlog =", round(fit$estimate[2], 3))) +
  geom_textvline(label = "Median = 2%", vjust = 1.3, hjust=0.01,
                 aes(xintercept = median_val), linetype = "dashed", color="black", alpha =1)+
  geom_textvline(label = "Mean = 19.3%", vjust = 1.3, hjust=0.01,
                 aes(xintercept = mean_val), linetype = "dotted", color="black", alpha =1)+
  scale_color_manual("Legend", 
                     values = c("Empirical CDF" = "blue", "Lognormal CDF" = "red"))
