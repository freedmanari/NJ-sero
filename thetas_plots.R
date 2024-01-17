### to be run after thetas.R
require(ggeasy)
require(bayestestR)
require(tidyverse)
require(colorspace)
require(rstan)
require(cowplot)
require(funtimes)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# make nice-looking scientific notation for plots
fancy_scientific <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("0e\\+00","0",l)
  l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("e", "%*%10^", l)
  l <- gsub("\\+","",l)
  parse(text=l)
}


plot_theta <- function(theta_best, theta_lower, theta_upper, zero_line=T, ylab="log odds ratio", sci_labels=T, time_range=c(1,W)) {
  data <- data.frame(date=get_week_start_from_test_week(time_range[1]:time_range[2]),
                     mean=theta_best,
                     CI_lower=theta_lower,
                     CI_upper=theta_upper)
  if (any(is.na(data))) {
    data <- data[(tail(which(apply(data,1,function(r) any(is.na(r)))),1)+1):W,]
  }
  plot <-
    data %>%
    pivot_longer(cols=-date, names_to="quantile", values_to="theta") %>% 
    mutate(type=ifelse(quantile=="mean","mean","CI")) %>% 
    ggplot() +
    geom_line(aes(date,theta,group=quantile,linetype=type,color=type)) +
    scale_linetype_manual(values=c("mean"=1,"CI"=2)) +
    scale_color_manual(values=c("mean"="black","CI"="steelblue2")) +
    scale_x_date(date_labels="%b '%y",date_breaks="1 month",expand=c(0,0),limits=get_week_start_from_test_week(time_range)) +
    theme_bw() +
    theme(panel.grid=element_blank(),
          panel.border=element_blank(),
          axis.text.x = element_text(angle = 90),
          legend.position = "none",
          axis.line = element_line(colour = "black", linewidth=.3))
    if (zero_line) {
      plot +
        geom_hline(yintercept=0, col="red", linewidth=.3) +
        scale_y_continuous(name=ylab)
    } else if (sci_labels) {
      plot +
        scale_y_continuous(name=ylab, labels=fancy_scientific, expand=expansion(c(0,.03))) +
        coord_cartesian(ylim = c(0,NA))
    } else {
      plot +
        scale_y_continuous(name=ylab, expand=expansion(c(0,.03))) +
        coord_cartesian(ylim = c(0,NA))
    }
}


plot_theta(I_best, I_lower, I_upper, zero_line=F, ylab="weekly incidence")



plot_theta(theta_I_P_best[1:(W-1)], theta_I_P_lower[1:(W-1)], theta_I_P_upper[1:(W-1)],
           time_range=c(1,W-1))

plot_theta(theta_I_S_best[1:(W-1)], theta_I_S_lower[1:(W-1)], theta_I_S_upper[1:(W-1)],
           time_range=c(1,W-1))

plot_theta(theta_Pp_P_best[1:(W-1)], theta_Pp_P_lower[1:(W-1)], theta_Pp_P_upper[1:(W-1)],
           time_range=c(1,W-1))

plot_theta(theta_Pp_S_best[1:(W-1)], theta_Pp_S_lower[1:(W-1)], theta_Pp_S_upper[1:(W-1)],
           time_range=c(1,W-1))

# only if theta_Pp_Sp was defined,
# either calculated with access to sero_tests.csv and all_sero, or loaded from thetas_simulations
# plot_theta(theta_Pp_Sp_best[6:(W-1)], theta_Pp_Sp_lower[6:(W-1)], theta_Pp_Sp_upper[6:(W-1)], time_range=c(6,W-1))




p1 <- plot_theta(I_best[1:(W-1)], I_lower[1:(W-1)], I_upper[1:(W-1)],
                 zero_line=F, ylab="weekly incidence", time_range=c(1,W-1))
p2 <- plot_theta(theta_I_P_best[1:(W-1)], theta_I_P_lower[1:(W-1)], theta_I_P_upper[1:(W-1)],
                 time_range=c(1,W-1))
plot_grid(p1, p2, ncol=1, align="v")


p2 <- plot_theta(theta_Pp_S_best[1:(W-1)], theta_Pp_S_lower[1:(W-1)], theta_Pp_S_upper[1:(W-1)],
                 time_range=c(1,W-1))
plot_grid(p1, p2, ncol=1, align="v")



## plots for correlations between I and thetas
plot_correlation <- function(I, theta, time_range=c(1,W)) {
  fit <- lm(theta[time_range[1]:time_range[2]] ~ log(I[time_range[1]:time_range[2]]))
  data.frame(I=I[time_range[1]:time_range[2]], theta=theta[time_range[1]:time_range[2]],
             week=get_week_start_from_test_week(time_range[1]:time_range[2])) %>%
    ggplot() +
    geom_point(aes(log(I), theta, color=week)) +
    geom_abline(intercept=fit$coefficients[1],
                slope=fit$coefficients[2],
                lty="dashed") +
    annotate("text", label=paste("R^2==",round(summary(fit)$r.squared,2)), parse=TRUE,
             x = min(log(I[time_range[1]:time_range[2]])) + .9 * (max(log(I[time_range[1]:time_range[2]]))-min(log(I[time_range[1]:time_range[2]]))),
             y = min(theta[time_range[1]:time_range[2]]) + .95 * (max(theta[time_range[1]:time_range[2]])-min(theta[time_range[1]:time_range[2]]))) +
    theme_classic() +
    xlab("log weekly incidence") +
    ylab("log odds ratio") +
    scale_color_date(name="week start",date_labels="%b '%y",date_breaks="4 month",
                     limits=as.Date(c("2020-05-01","2021-9-21")))
}

plot_correlation(I_best, theta_I_P_best, time_range=c(6,(W-1)))
plot_correlation(I_best, theta_I_S_best, time_range=c(6,(W-1)))
plot_correlation(I_best, theta_Pp_P_best, time_range=c(6,(W-1)))
plot_correlation(I_best, theta_Pp_S_best, time_range=c(6,(W-1)))

ccf_boot(I_best[6:(W-1)],theta_I_P_best[6:(W-1)],plot="Spearman",B=10000,lag.max = 11)
ccf_boot(I_best[6:(W-1)],theta_I_S_best[6:(W-1)],plot="Spearman",B=10000,lag.max = 11)
ccf_boot(I_best[6:(W-1)],theta_Pp_P_best[6:(W-1)],plot="Spearman",B=10000,lag.max = 11)
ccf_boot(I_best[6:(W-1)],theta_Pp_S_best[6:(W-1)],plot="Spearman",B=10000,lag.max = 11)



## plot of fatality rate used in estimating incidence
plot_theta(IFR_mean*CFR/mean(CFR), IFR_lower*CFR/mean(CFR), IFR_upper*CFR/mean(CFR),
           ylab="fatality rate used for incidence estimation", zero_line=F, sci_labels=F, time_range=c(1,W_deaths))


## plot of underreporting multiplier
plot_theta(I_best[4:(W-1)]/Pp[4:(W-1)], I_lower[4:(W-1)]/Pp[4:(W-1)], I_upper[4:(W-1)]/Pp[4:(W-1)],
           ylab="underreporting multiplier", zero_line=F, sci_labels=F, time_range=c(4,W-1))

