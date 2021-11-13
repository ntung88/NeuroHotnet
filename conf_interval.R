confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}
neurconf = confidence_interval(finalneur,0.95)
higconf = confidence_interval(finalhig,0.95)
naiconf = confidence_interval(finalnai,0.95)

data <- data.frame(
  name=c('Neuro-Hotnet','SiGGM with Diffusion','SC Naive'),
  value=c(mean(finalneur),mean(finalhig),mean(finalnai)),
  mins=c(neurconf[[1]],higconf[[1]],naiconf[[1]]),
  maxs=c(neurconf[[2]],higconf[[2]],naiconf[[2]])
)

# Most basic error bar
pdf('recrates.pdf')
ggplot(data) +
  geom_bar( aes(x=name, y=value), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=name, ymin=mins, ymax=maxs), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  ylim(0,1) +
  labs(title=expression(paste("G(", gamma, ') Recovery Rates')),x='',y='Average Recovery Rate') +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# barplot(c(mean(finalneur),mean(finalhig),mean(finalnai)),
#         names.arg=c('Neuro-Hotnet','SiGGM','SC Naive'),
#         ylim=0:1,
#         ylab='Average Recovery Rate',
#         main='Recovery Rates For All Methods')