library(tidyverse)
library(data.table)
library(forecast)
library(fpp2)
library(zoo)
library(tseries)
library(cowplot)
library(kableExtra)
library(xtable)

# Data ----
riders = read_csv('data/portland-oregon-average-monthly.csv')
colnames(riders) = c('month', 'bus_ridership_avg')
riders = riders %>% 
  slice(1:114) %>% 
  mutate(month = as.Date(as.yearmon(month)),
         bus_ridership_avg = as.numeric(bus_ridership_avg))

riders_ts = ts(data=riders[, 'bus_ridership_avg'],
               start=c(1960, 1),
               end=c(1969, 6),
               frequency = 12)
x = window(riders_ts, start = c(1960, 1), end = c(1967, 12))
n=length(x)
xx = window(riders_ts, start = c(1968, 1), end = c(1969, 6))

# Time series plot ----
x %>%
  autoplot() +
  labs(x='Year', y='Observed Value') +
  theme_bw()
ggsave("images/serie.pdf", width = 158, height = 93, units = "mm")

# STL decomposition
x %>%
  mstl() %>%
  autoplot() +
  labs(x='Year') +
  theme_bw()
ggsave("images/decomposition.pdf", width = 158, height = 93, units = "mm")

# ARIMA ----
# Diffs
diffs <- c()
if (x %>% ndiffs() == 1) {
  w <- diff(x)
  diffs <- c(diffs, 1)
  if (w %>% nsdiffs() == 1) {
    w <- diff(w, lag=12)
    diffs <- c(diffs, 12)
  }
}

cbind(Original=x, Differentiated=w) %>%
  autoplot(facets=T) +
  labs(x='Year', y='Observed Value') +
  theme_bw()
ggsave("images/original-diff.pdf", width = 158, height = 93, units = "mm")

w %>% adf.test()

# Manual selection
w %>%
  ggAcf(lag.max=12*5) +
  labs(title='') +
  theme_bw()
ggsave("images/serie-acf.pdf", width = 158, height = 93, units = "mm")

w %>%
  ggPacf(lag.max=12*5) +
  labs(title='') +
  theme_bw()
ggsave("images/serie-pacf.pdf", width = 158, height = 93, units = "mm")

p = 0
q = 0
P = 0:2
Q = 0:2

# Without Box-Cox
best_aicc <- Inf
for (Pi in P) {
  for (Qi in Q) {
    mod <- x %>% Arima(order=c(p, 1, q), seasonal=c(Pi, 1, Qi),
                       include.mean=FALSE, lambda=NULL)
    if (mod$aicc < best_aicc) {
      best_mod <- mod
      best_aicc <- mod$aicc
    }
  }
} 

mod1 <- best_mod
mod1 %>% summary()
# With Box-Cox
best_aicc <- Inf
for (Pi in P) {
  for (Qi in Q) {
    mod <- x %>% Arima(order=c(p, 1, q), seasonal=c(Pi, 1, Qi),
                       include.mean=FALSE, lambda='auto')
    if (mod$aicc < best_aicc) {
      best_mod <- mod
      best_aicc <- mod$aicc
    }
  }
}

mod2 <- best_mod
mod2 %>% summary()

# Residual analysis without Box-Cox
res1 <- mod1 %>% residuals()
est1 <- res1 %>%
  autoplot() +
  labs(x='Year', y='Residuals') +
  theme_bw()
qq1 <- res1 %>%
  data.frame() %>%
  ggplot(aes(sample=res1)) +
  stat_qq() +
  stat_qq_line() +
  labs(x='Theoretical Quantiles', y='Sample Quantiles') +
  theme_bw()
acf1 <- res1 %>%
  ggAcf(lag.max=12*5) +
  labs(title='') +
  theme_bw()
pacf1 <- res1 %>%
  ggPacf(lag.max=12*5) +
  labs(title='') +
  theme_bw()
plot_grid(est1, qq1, acf1, pacf1, nrow=2)
ggsave("images/arima-residuals.pdf", width = 158, height = 93, units = "mm")
res1 %>% adf.test() 
res1 %>% Box.test(lag=20, type='Ljung-Box', fitdf=2)
res1 %>% shapiro.test()
mod1 %>% checkresiduals()

# Residual analysis with Box-Cox
res2 <- mod2 %>% residuals()
est2 <- res2 %>%
  autoplot() +
  labs(x='Year', y='Residuals') +
  theme_bw()
qq2 <- res2 %>%
  data.frame() %>%
  ggplot(aes(sample=res2)) +
  stat_qq() +
  stat_qq_line() +
  labs(x='Theoretical Quantiles', y='Sample Quantiles') +
  theme_bw()
acf2 <- res2 %>%
  ggAcf(lag.max=12*5) +
  labs(title='') +
  theme_bw()
pacf2 <- res2 %>%
  ggPacf(lag.max=12*5) +
  labs(title='') +
  theme_bw()
plot_grid(est2, qq2, acf2, pacf2, nrow=2)
ggsave("images/arima-box-cox-residuals.pdf", width = 158, height = 93, units = "mm")
mod2
res2 %>% adf.test()
res2 %>% Box.test(lag=20, type='Ljung-Box', fitdf=2)
res2 %>% shapiro.test()
mod2 %>% checkresiduals()

# ETS---
# Auto selection without Box-Cox
mod3 <- x %>% ets()
mod3 %>% summary()
# Auto selection with Box-Cox
mod4 <- x %>% ets(lambda='auto')
mod4 %>% summary()
# Residual analysis without Box-Cox
res3 <- mod3 %>% residuals()
est3 <- res3 %>%
  autoplot() +
  labs(x='Year', y='Residuals') +
  theme_bw()
qq3 <- res3 %>%
  data.frame() %>%
  ggplot(aes(sample=res3)) +
  stat_qq() +
  stat_qq_line() +
  labs(x='Theoretical Quantiles', y='Sample Quantiles') +
  theme_bw()
acf3 <- res3 %>%
  ggAcf(lag.max=12*5) +
  labs(title='') +
  theme_bw()
pacf3 <- res3 %>%
  ggPacf(lag.max=12*5) +
  labs(title='') +
  theme_bw()
plot_grid(est3, qq3, acf3, pacf3, nrow=2)
ggsave("images/ets-residuals.pdf", width = 158, height = 93, units = "mm")
res3 %>% adf.test()
res3 %>% Box.test(lag=20, type='Ljung-Box', fitdf=17)
res3 %>% shapiro.test()
mod3 %>% checkresiduals()
# Residual analysis with Box-Cox
res4 <- mod4 %>% residuals()
est4 <- res4 %>%
  autoplot() +
  labs(x='Year', y='Residuals') +
  theme_bw()
qq4 <- res4 %>%
  data.frame() %>%
  ggplot(aes(sample=res4)) +
  stat_qq() +
  stat_qq_line() +
  labs(x='Theoretical Quantiles', y='Sample Quantiles') +
  theme_bw()
acf4 <- res4 %>%
  ggAcf(lag.max=12*5) +
  labs(title='') +
  theme_bw()
pacf4 <- res4 %>%
  ggPacf(lag.max=12*5) +
  labs(title='') +
  theme_bw()
plot_grid(est4, qq4, acf4, pacf4, nrow=2)
ggsave("images/ets-box-cox-residuals.pdf", width = 158, height = 93, units = "mm")
res4 %>% adf.test()
res4 %>% Box.test(lag=20, type='Ljung-Box', fitdf=17)
res4 %>% shapiro.test()
mod4 %>% checkresiduals()

# Sliding window validation
f_arima1 <- function(y, h){
  fit <- Arima(y, order=c(0, 1, 0), seasonal=c(1, 1, 1), include.mean=F, lambda=NULL)
  forecast(fit, h)
}
f_arima2 <- function(y, h){
  fit <- Arima(y, order=c(0, 1, 0), seasonal=c(1, 1, 1), include.mean=F, lambda='auto')
  forecast(fit, h)
}
f_ets1 <- function(y, h){
  fit <- ets(y, model='MAM', damped = TRUE)
  forecast(fit, h)
}
f_ets2 <- function(y, h){
  fit <- ets(y, model='AAA', damped = TRUE)
  forecast(fit, h)
}
h=5
CV_arima1 <- x %>% tsCV(forecastfunction=f_arima1, h=h, initial=n-14)
CV_arima2 <- x %>% tsCV(forecastfunction=f_arima2, h=h, initial=n-14)
CV_ets1 <- x %>% tsCV(forecastfunction=f_ets1, h=h, initial=n-14)
CV_ets2 <- x %>% tsCV(forecastfunction=f_ets2, h=h, initial=n-14)
MAE_arima1 <- CV_arima1 %>% abs() %>% colMeans(na.rm=T)
MAE_arima2 <- CV_arima2 %>% abs() %>% colMeans(na.rm=T)
MAE_ets1 <- CV_ets1 %>% abs() %>% colMeans(na.rm=T)
MAE_ets2 <- CV_ets2 %>% abs() %>% colMeans(na.rm=T)
tab <- cbind(MAE_arima1, MAE_arima2, MAE_ets1, MAE_ets2)
tab %>%
  kable(
    col.names=c('ARIMA', 'ARIMA + Box-Cox', 'ETS', 'ETS + Box-Cox'),
    caption='MAE por horizonte de predição.',
    digits=0,
    format.args=list(decimal.mark=',', scientific=F),
    align='c'
  ) %>%
  kable_styling(
    position='center',
    bootstrap_options=c('striped', 'hover', 'condensed', 'responsive')
  )

tab_plot <- tab %>%
  as.data.frame() %>%
  mutate(Horizon=1:h) %>%
  gather(key='Model', value='MAE', -Horizon)
tab_plot %>%
  ggplot(aes(x=Horizon, y=MAE)) +
  geom_line(aes(color=Model)) + 
  scale_color_manual(
    values=c('black', 'red', '#0000AA', 'darkgreen'),
    breaks=c('MAE_arima1', 'MAE_arima2', 'MAE_ets1', 'MAE_ets2'),
    labels=c('ARIMA', 'ARIMA + Box-Cox', 'ETS', 'ETS + Box-Cox')
  ) +
  theme_bw()

ggsave("images/horizonte-predicao.pdf", width = 158, height = 93, units = "mm")

# Forecast ----
# tables
h=18
preds1 <- forecast(mod1, h=h, level=95)
preds2 <- forecast(mod2, h=h, level=95)
preds3 <- forecast(mod3, h=h, level=95)
preds4 <- forecast(mod4, h=h, level=95)
pontual <- t(cbind(xx, preds1$mean, preds2$mean, preds3$mean, preds4$mean))
colnames(pontual) <- 1:h
row.names(pontual) <- c('Observado', 'ARIMA', 'ARIMA + Box-Cox', 'ETS', 'ETS + Box-Cox')
pontual %>%
  kable(
    caption='Previsões pontuais por horizonte de predição.',
    digits=0,
    format.args=list(decimal.mark=',', scientific=F),
    align='c'
  ) %>%
  kable_styling(
    position='center',
    bootstrap_options=c('striped', 'hover', 'condensed', 'responsive')
  ) 

t(pontual) %>% xtable(digits=0)

intervalares <- t(cbind(xx, preds1$lower, preds1$upper, preds2$lower, preds2$upper,
                        preds3$lower, preds3$upper, preds4$lower, preds4$upper))
colnames(intervalares) <- 1:h
row.names(intervalares) <- c('Observado', 'ARIMA Inf', 'ARIMA Sup', 'ARIMA + Box-Cox Inf',
                             'ARIMA + Box-Cox Sup', 'ETS Inf', 'ETS Sup', 'ETS + Box-Cox Inf',
                             'ETS + Box-Cox Sup')
intervalares %>%
  kable(
    caption='Previsões intervalares de 95% de confiança por horizonte de predição.',
    digits=0,
    format.args=list(decimal.mark=',', scientific=F),
    align='c'
  ) %>%
  kable_styling(
    position='center',
    bootstrap_options=c('striped', 'hover', 'condensed', 'responsive')
  )

t(intervalares) %>% xtable(digits = 0)
# plots
plot_preds <- function(mod, nome='') {
  vec <- c(nome, 'Observed')
  cores <- c('#0000AA', 'red')
  names(cores) <- vec
  preds <- forecast(mod, h=h, level=95)
  plot_obj <- x %>%
    autoplot() + xlab('Year') + ylab('Observed Value') + theme_bw() +
    autolayer(preds, series=nome) +
    autolayer(xx, series='Observed') +
    scale_colour_manual(
      values=cores,
      breaks=vec,
      name='')
  return(plot_obj)
}
plot_preds(mod1, 'ARIMA')
ggsave("images/arima-forecast.pdf", width = 158, height = 93, units = "mm")
plot_preds(mod2, 'ARIMA + Box-Cox')
ggsave("images/arima-box-cox-forecast.pdf", width = 158, height = 93, units = "mm")
plot_preds(mod3, 'ETS')
ggsave("images/ets-forecast.pdf", width = 158, height = 93, units = "mm")
plot_preds(mod4, 'ETS + Box-Cox')
ggsave("images/ets-box-cox-forecast.pdf", width = 158, height = 93, units = "mm")

# Benchmark comparison
preds <- list(
  'ARIMA' = forecast(mod1, h=h),
  'ARIMA + Box-Cox' = forecast(mod2, h=h),
  'ETS' = forecast(mod3, h=h),
  'ETS + Box-Cox' = forecast(mod4, h=h),
  'auto.arima' = forecast(auto.arima(x), h=h),
  'SES' = ses(x, h=h),
  'Holt' = holt(x, h=h),
  'sltf' = stlf(x, h=h),
  'BATS' = forecast(bats(x), h=h),
  'TBATS' = forecast(tbats(x), h=h)
)
mae <- unlist(lapply(preds, function(m) return(mean(abs(xx - m$mean)))))
final <- data.frame(MAE=mae)
final %>%
  kable(
    caption='MAE nos dados de teste.',
    digits=0,
    format.args=list(decimal.mark=',', scientific=F),
    align='c'
  ) %>%
  kable_styling(
    position='center',
    bootstrap_options=c('striped', 'hover', 'condensed', 'responsive')
  )

final %>% 
  arrange(MAE) %>% 
  xtable()