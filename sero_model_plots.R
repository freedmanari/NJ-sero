### to be run after sero_model.R
require(tidyverse)
require(rstan)
require(lubridate)
require(cowplot)
require(ggtext)
require(stringr)
require(scales)
require(bayestestR)
require(colorspace)
require(viridis)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# for saving the results of the sero model simulations

# save(sero_model_summs_normal_r_SV,
#      summ_prev_PCR_normal_r_SV,
#      sero_model_normal_r_SV,
#      sero_model_prev_PCR_normal_r_SV,
#      model_sero_normal_r_SV,
#      sero_model_summs_low_r_SV,
#      summ_prev_PCR_low_r_SV,
#      sero_model_low_r_SV,
#      sero_model_prev_PCR_low_r_SV,
#      model_sero_low_r_SV,
#      sero_model_summs_high_r_SV,
#      summ_prev_PCR_high_r_SV,
#      sero_model_high_r_SV,
#      sero_model_prev_PCR_high_r_SV,
#      model_sero_high_r_SV,
#      file="data/sero_model_simulations.RData")

# for loading the results of the sero model simulations

# load("data/sero_model_simulations.RData")







## choose which set of sero model outputs to use based on r_SV value

# for main text

r_SV <- 1
sero_model_summs <- sero_model_summs_normal_r_SV
summ_prev_PCR <- summ_prev_PCR_normal_r_SV
sero_model <- sero_model_normal_r_SV
sero_model_prev_PCR <- sero_model_prev_PCR_normal_r_SV
model_sero <- model_sero_normal_r_SV
model_sero_prev_PCR <- model_sero %>% filter(prev_PCR)
model_sero_no_prev_PCR <- model_sero %>% filter(!prev_PCR)


# options for supplementary sensitivity analysis

# r_SV <- .5
# sero_model_summs <- sero_model_summs_low_r_SV
# summ_prev_PCR <- summ_prev_PCR_low_r_SV
# sero_model <- sero_model_low_r_SV
# sero_model_prev_PCR <- sero_model_prev_PCR_low_r_SV
# model_sero <- model_sero_low_r_SV
# model_sero_prev_PCR <- model_sero %>% filter(prev_PCR)
# model_sero_no_prev_PCR <- model_sero %>% filter(!prev_PCR)

# r_SV <- 2
# sero_model_summs <- sero_model_summs_high_r_SV
# summ_prev_PCR <- summ_prev_PCR_high_r_SV
# sero_model <- sero_model_high_r_SV
# sero_model_prev_PCR <- sero_model_prev_PCR_high_r_SV
# model_sero <- model_sero_high_r_SV
# model_sero_prev_PCR <- model_sero %>% filter(prev_PCR)
# model_sero_no_prev_PCR <- model_sero %>% filter(!prev_PCR)



# make nice-looking scientific notation for plots
fancy_scientific <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("0e\\+00","0",l)
  l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("e", "%*%10^", l)
  l <- gsub("\\+","",l)
  parse(text=l)
}



## titers over time out of just positive serology tests

model_sero %>%
  filter(result=="P") %>% 
  group_by(test_week) %>%
  summarise(week_start=min(test_date),
            y_data=mean(log(numeric))) %>%
  ungroup() %>%
  filter(test_week>min(test_week)) %>% 
  left_join(data.frame(test_week=1:W, prop_vaccinated=cumsum(V)/N)) %>% 
  ggplot() +
  geom_point(aes(week_start, y_data)) +
  geom_vline(xintercept=get_week_start_from_test_week(w_vac), color="red", lty="dashed") +
  scale_x_date(date_labels="%b '%y",date_breaks="1 month") +
  scale_y_continuous(breaks=seq(3.75,5,.25),
                     sec.axis=sec_axis(~ .-4,
                                       name = 'proportion vaccinated in NJ',
                                       breaks=seq(0,1,.25))) +
  expand_limits(y=5) +
  geom_line(aes(week_start, prop_vaccinated+4), col="blue") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.line = element_line(colour = "black", linewidth=.3),
        axis.line.y.right = element_line(color = "blue"), 
        axis.ticks.y.right = element_line(color = "blue"),
        axis.text.y.right = element_text(color="blue"),
        axis.title.y.right = element_text(color="blue"),) +
  ylab("mean log titer value") +
  xlab("test date")


# broken up by past PCR positive or not

model_sero %>% 
  filter(test_date>="2020-05-01",test_date<"2021-09-01",result=="P") %>% 
  group_by(year=year(test_date)+1/12*month(test_date), prev_PCR=factor(ifelse(prev_PCR,"yes","no"),levels=c("yes","no"))) %>%
  summarise(y=mean(y_data),
            lower=t.test(y_data)[[4]][1],
            upper=t.test(y_data)[[4]][2]) %>%
  ungroup() %>% 
  mutate(date=as.Date(sapply(year, function(d) paste(floor(d-1/12),"-",round((d-floor(d-1/12))*12), "-1",sep="")))) %>% 
  ggplot() +
  geom_point(aes(date,y,group=prev_PCR,color=prev_PCR),position=position_dodge(10)) +
  geom_errorbar(aes(date,ymin=lower,ymax=upper,color=prev_PCR),width=8,position=position_dodge(10)) +
  ylab("mean log titer value") +
  xlab("test date") +
  scale_color_manual(name="past PCR positive", values=c("#00BFC4", "#F8766D")) +
  scale_x_date(date_labels="%b '%y",date_breaks="1 month") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.line = element_line(colour = "black", linewidth=.3),
        legend.position = c(.8,.18),
        legend.title=element_text(size=9)) +
  geom_vline(xintercept=get_week_start_from_test_week(w_vac), color="red", lty="dashed")







## plots for titers by delay

y_pred_total <- lapply(0:W, function(delay_int)
  c(sero_model_prev_PCR@sim$samples[[1]][[paste("y_pred_age1[",delay_int+1,"]",sep="")]][1001:2000],
    sero_model_prev_PCR@sim$samples[[2]][[paste("y_pred_age1[",delay_int+1,"]",sep="")]][1001:2000],
    sero_model_prev_PCR@sim$samples[[1]][[paste("y_pred_age2[",delay_int+1,"]",sep="")]][1001:2000],
    sero_model_prev_PCR@sim$samples[[2]][[paste("y_pred_age2[",delay_int+1,"]",sep="")]][1001:2000],
    sero_model_prev_PCR@sim$samples[[1]][[paste("y_pred_age3[",delay_int+1,"]",sep="")]][1001:2000],
    sero_model_prev_PCR@sim$samples[[2]][[paste("y_pred_age3[",delay_int+1,"]",sep="")]][1001:2000],
    sero_model_prev_PCR@sim$samples[[1]][[paste("y_pred_age4[",delay_int+1,"]",sep="")]][1001:2000],
    sero_model_prev_PCR@sim$samples[[2]][[paste("y_pred_age4[",delay_int+1,"]",sep="")]][1001:2000],
    sero_model_prev_PCR@sim$samples[[1]][[paste("y_pred_age5[",delay_int+1,"]",sep="")]][1001:2000],
    sero_model_prev_PCR@sim$samples[[2]][[paste("y_pred_age5[",delay_int+1,"]",sep="")]][1001:2000]))

y_pred_age_data <- 
  data.frame(delay_int=0:W,
             mean_total=y_pred_total %>% lapply(mean) %>% unlist,
             low_total=y_pred_total %>% lapply(function(l) hdi(l,ci=.95)$CI_low) %>% unlist,
             high_total=y_pred_total %>% lapply(function(l) hdi(l,ci=.95)$CI_high) %>% unlist,
             mean_age1=summary(sero_model_prev_PCR)$summary[grep("y_pred_age1",names(summ_prev_PCR)), "50%"],
             low_age1=summary(sero_model_prev_PCR)$summary[grep("y_pred_age1",names(summ_prev_PCR)), "2.5%"],
             high_age1=summary(sero_model_prev_PCR)$summary[grep("y_pred_age1",names(summ_prev_PCR)), "97.5%"],
             mean_age2=summary(sero_model_prev_PCR)$summary[grep("y_pred_age2",names(summ_prev_PCR)), "50%"],
             low_age2=summary(sero_model_prev_PCR)$summary[grep("y_pred_age2",names(summ_prev_PCR)), "2.5%"],
             high_age2=summary(sero_model_prev_PCR)$summary[grep("y_pred_age2",names(summ_prev_PCR)), "97.5%"],
             mean_age3=summary(sero_model_prev_PCR)$summary[grep("y_pred_age3",names(summ_prev_PCR)), "50%"],
             low_age3=summary(sero_model_prev_PCR)$summary[grep("y_pred_age3",names(summ_prev_PCR)), "2.5%"],
             high_age3=summary(sero_model_prev_PCR)$summary[grep("y_pred_age3",names(summ_prev_PCR)), "97.5%"],
             mean_age4=summary(sero_model_prev_PCR)$summary[grep("y_pred_age4",names(summ_prev_PCR)), "50%"],
             low_age4=summary(sero_model_prev_PCR)$summary[grep("y_pred_age4",names(summ_prev_PCR)), "2.5%"],
             high_age4=summary(sero_model_prev_PCR)$summary[grep("y_pred_age4",names(summ_prev_PCR)), "97.5%"],
             mean_age5=summary(sero_model_prev_PCR)$summary[grep("y_pred_age5",names(summ_prev_PCR)), "50%"],
             low_age5=summary(sero_model_prev_PCR)$summary[grep("y_pred_age5",names(summ_prev_PCR)), "2.5%"],
             high_age5=summary(sero_model_prev_PCR)$summary[grep("y_pred_age5",names(summ_prev_PCR)), "97.5%"]) %>%
  filter(delay_int <= max(model_sero_prev_PCR$delay_int)) %>% 
  pivot_longer(cols=-delay_int,
               names_to="age_group",
               names_prefix="^[^_]*_",
               values_to="y_pred") %>%
  cbind(type=c("mean","low","high")) %>%
  pivot_wider(id_cols=c(delay_int, age_group),
              names_from=type,
              values_from=y_pred)


# titers by delay with confidence bars encompassing all age classes

model_sero_prev_PCR %>%
  group_by(delay_int) %>%
  summarise(y=mean(log(numeric))) %>%
  ggplot() +
  geom_ribbon(data=filter(y_pred_age_data,age_group=="total"), aes(delay_int,ymin=low,ymax=high), fill="grey", alpha=.3) +
  geom_line(data=filter(y_pred_age_data,age_group=="total"), aes(delay_int,y=mean), col="red") +
  geom_point(aes(delay_int,y)) +
  geom_hline(aes(yintercept=filter(y_pred_age_data,age_group=="total",delay_int==0)$mean),linewidth=.3) +
  xlab("delay from PCR positive to serology test (weeks)") +
  ylab("mean log titer value") +
  scale_x_continuous(expand=expansion(c(0,0)), breaks=c(0,20,160/7,40,60), labels=c(0,20,"x",40,60)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(color="black",linewidth=.3))



# titers by delay with age class specific confidence intervals

age_colors <- c("orange","magenta","green","red","blue")

model_sero_prev_PCR %>%
  group_by(delay_int) %>%
  summarise(y=mean(log(numeric))) %>%
  ggplot() +
  geom_ribbon(data=filter(y_pred_age_data, age_group!="total"),
              aes(x=delay_int,ymin=low,ymax=high,group=age_group,fill=age_group,color=age_group),
              alpha=.1) +
  xlab("delay from PCR positive to serology test (weeks)") +
  ylab("mean log titer value") +
  scale_x_continuous(expand=expansion(c(0,0))) +
  scale_y_continuous(expand=expansion(c(0,.05))) +
  scale_fill_manual(values=age_colors,
                    name="age group",
                    labels=c("0-15","16-29","30-49","50-64","65+")) +
  scale_color_manual(values=adjust_transparency(age_colors,alpha=.65),
                     name="age group",
                     labels=c("0-15","16-29","30-49","50-64","65+")) +
  geom_point(aes(delay_int,y)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(linewidth=.3),
        plot.margin = margin(t = 5, r = 5, b = 10, l = 5, unit = "pt"))





## plot the interval weights used for simulating delays for serology tests with no past PCR positive


for (i in 1:length(interval_mins)) {
  print(
    interval_probs %>%
    filter(interval==i) %>%
    ggplot() +
    geom_point(aes(delay,p)) +
    geom_function(fun = ~dggamma(.x, interval_weight_pars[i,1], interval_weight_pars[i,2], interval_weight_pars[i,3]) / interval_weight_pars[i,4],
                  col="red",
                  xlim=c(.0001,max(interval_probs$delay)),
                  n=50000) +
    scale_x_continuous(name="delay from PCR positive to serology test (weeks)",
                       expand=expansion(c(.01,.05))) +
    scale_y_continuous(name="interval weight",
                       expand=expansion(c(.01,.05)),
                       limits=c(0,1.1*max(interval_probs$p[interval_probs$interval==i]))) +
    ggtitle(paste("weights for log-titers in interval [",i,",",i+1,")",sep="")) +
    theme_bw() +
    theme(panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid=element_blank(),
          axis.line=element_line(linewidth=.3),
          plot.title = element_text(hjust = 0.5)))
}







## sero model fit

summ <- sero_model_summs[[sero_model_reps]]

model_sero %>%
  cbind(y_pred_stan=summ[grep("y_pred[1-9]", names(summ))]) %>% 
  filter(test_date>="2020-05-01" & test_date<"2021-09-01") %>% 
  pivot_longer(cols=c(y_data, y_pred_stan), names_to="y_type", values_to="y") %>%
  group_by(year=year(test_date)+1/12*month(test_date), prev_PCR=factor(ifelse(prev_PCR,"yes","no"),levels=c("yes","no")),y_type) %>%
  summarise(lwr=t.test(y)[[4]][1],
            upr=t.test(y)[[4]][2],
            y=mean(y)) %>%
  ungroup() %>% 
  mutate(date=as.Date(sapply(year, function(d) paste(floor(d-1/12),"-",round((d-floor(d-1/12))*12), "-1",sep=""))) + ifelse(prev_PCR=="yes",-2.5,2.5)) %>% 
  ggplot() +
  geom_line(aes(date,y,group=interaction(prev_PCR,y_type),color=prev_PCR,linetype=factor(y_type))) +
  geom_errorbar(aes(date,ymin=lwr,ymax=upr,group=interaction(prev_PCR,y_type),color=prev_PCR,alpha=ifelse(y_type=="y_data",1,0))) +
  ylab("mean log titer value") +
  xlab("test date") +
  ggtitle(bquote("serology model fit for" ~ italic(r)[italic(SV)] == .(r_SV))) +
  scale_color_manual(name="past PCR\npositive", values=c("#00BFC4", "#F8766D")) +
  scale_linetype_discrete(name="titers from...", labels=c("data","predictions")) +
  scale_x_date(date_labels="%b '%y",date_breaks="1 month") +
  scale_y_continuous(breaks=seq(2.5,5,.5)) +
  scale_alpha(guide = 'none', range=c(0,1)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.line = element_line(colour = "black", linewidth=.3),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        plot.title=element_text(hjust = 0.5)) +
  geom_vline(xintercept=get_week_start_from_test_week(w_vac), color="red", lty="dashed")


## sero model fits for just tests with past PCR positive, across all simulations

stan_preds_prev_PCR <- data.frame(r=numeric(), date=Date(), y=numeric())
for (r in 1:sero_model_reps) {
  summ <- sero_model_summs[[r]]
  temp <-
    model_sero %>%
    cbind(y_pred_stan=summ[grep("y_pred[1-9]", names(summ))]) %>% 
    filter(test_date>="2020-05-01" & test_date<"2021-09-01", prev_PCR) %>%
    group_by(year=year(test_date)+1/12*month(test_date)) %>%
    summarise(y=mean(y_pred_stan)) %>%
    ungroup() %>% 
    mutate(date=as.Date(sapply(year, function(d) paste(floor(d-1/12),"-",round((d-floor(d-1/12))*12), "-1",sep=""))))
  
  stan_preds_prev_PCR <- rbind(stan_preds_prev_PCR, data.frame(r=r, date=temp$date, y=temp$y))
}

y_data_prev_PCR <-
  model_sero %>%
  filter(test_date>="2020-05-01" & test_date<"2021-09-01", prev_PCR) %>% 
  group_by(year=year(test_date)+1/12*month(test_date)) %>%
  summarise(lwr=t.test(y_data)[[4]][1],
            upr=t.test(y_data)[[4]][2],
            y=mean(y_data)) %>%
  mutate(date=as.Date(sapply(year, function(d) paste(floor(d-1/12),"-",round((d-floor(d-1/12))*12), "-1",sep=""))))


stan_preds_no_prev_PCR <- data.frame(r=numeric(), date=Date(), y=numeric())
for (r in 1:sero_model_reps) {
  summ <- sero_model_summs[[r]]
  temp <-
    model_sero %>%
    cbind(y_pred_stan=summ[grep("y_pred[1-9]", names(summ))]) %>% 
    filter(test_date>="2020-05-01" & test_date<"2021-09-01", !prev_PCR) %>%
    group_by(year=year(test_date)+1/12*month(test_date)) %>%
    summarise(y=mean(y_pred_stan)) %>%
    ungroup() %>% 
    mutate(date=as.Date(sapply(year, function(d) paste(floor(d-1/12),"-",round((d-floor(d-1/12))*12), "-1",sep=""))))
  
  stan_preds_no_prev_PCR <- rbind(stan_preds_no_prev_PCR, data.frame(r=r, date=temp$date, y=temp$y))
}

y_data_no_prev_PCR <-
  model_sero %>%
  filter(test_date>="2020-05-01" & test_date<"2021-09-01", !prev_PCR) %>% 
  group_by(year=year(test_date)+1/12*month(test_date)) %>%
  summarise(lwr=t.test(y_data)[[4]][1],
            upr=t.test(y_data)[[4]][2],
            y=mean(y_data)) %>%
  mutate(date=as.Date(sapply(year, function(d) paste(floor(d-1/12),"-",round((d-floor(d-1/12))*12), "-1",sep=""))))


stan_preds_prev_PCR %>%
  ggplot() +
  geom_errorbar(data=y_data_prev_PCR, aes(x=date,ymin=lwr, ymax=upr), color=hue_pal()(2)[2], width=8) +
  geom_point(data=y_data_prev_PCR, aes(x=date,y=y), color=hue_pal()(2)[2]) +
  geom_vline(xintercept=get_week_start_from_test_week(w_vac), color="red", lty="dashed") +
  geom_line(aes(date,y,color=factor(r),group=r),alpha=.3) +
  ylab("mean log titer value") +
  xlab("test date") +
  scale_x_date(date_labels="%b '%y",date_breaks="1 month") +
  scale_y_continuous(breaks=seq(2.5,5,.5),
                     limits=c(min(c(y_data_prev_PCR$lwr,y_data_no_prev_PCR$lwr)),max(c(y_data_prev_PCR$upr,y_data_no_prev_PCR$upr)))) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.line = element_line(colour = "black", linewidth=.3),
        legend.position="none",
        plot.title=element_text(hjust = 0.5))
  

## sero model fits for just tests without past PCR positive, across all simulations

stan_preds_no_prev_PCR %>%
  ggplot() +
  geom_point(data=y_data_no_prev_PCR, aes(x=date,y=y), color=hue_pal()(2)[1]) +
  geom_errorbar(data=y_data_no_prev_PCR, aes(x=date,ymin=lwr, ymax=upr), color=hue_pal()(2)[1], width=8) +
  geom_vline(xintercept=get_week_start_from_test_week(w_vac), color="red", lty="dashed") +
  geom_line(aes(date,y,color=factor(r),group=r),alpha=.3) +
  ylab("mean log titer value") +
  xlab("test date") +
  scale_x_date(date_labels="%b '%y",date_breaks="1 month") +
  scale_y_continuous(breaks=seq(2.5,5,.5),
                     limits=c(min(c(y_data_prev_PCR$lwr,y_data_no_prev_PCR$lwr)),max(c(y_data_prev_PCR$upr,y_data_no_prev_PCR$upr)))) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.line = element_line(colour = "black", linewidth=.3),
        legend.position="none",
        plot.title=element_text(hjust = 0.5))




## parameter fit plots

# r_VP fits
r_VPs <- data.frame(r=numeric(), x=numeric(), y=numeric())
for (r in 1:sero_model_reps) {
  summ <- sero_model_summs[[r]]
  r_VPs <- rbind(r_VPs, data.frame(r=r, x=c(w_vac,summ["w_switch"],W), y=summ[c("r_VP_init","r_VP_init","r_VP_end")]))
}
r_VPs <- r_VPs %>% mutate(date=get_week_start_from_test_week(x))

ggplot(r_VPs) +
  geom_line(aes(date,y,color=factor(r),group=r),alpha=.4) +
  ylab(expression(paste(italic(r[VP]),"(",italic(w),")"))) +
  scale_x_date(name="test date",date_labels="%b '%y",date_breaks="1 month") +
  scale_y_continuous(name=bquote(italic(r)[italic(VP)](italic(w))),
                     limits=c(.07,.73),
                     expand=expansion(c(0,0)),
                     breaks=seq(.1,.7,.1)) +
  ggtitle(bquote(italic(r)[italic(VP)](italic(w)) ~ "fits for" ~ italic(r)[italic(SV)] == .(r_SV))) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(linewidth=.3),
        legend.position="none",
        axis.text.x=element_text(angle = 90),
        plot.margin=margin(5,10,5,5),
        plot.title=element_text(hjust=.5))



# other important parameters, for supplement

p1 <- data.frame(rep=factor(1:sero_model_reps),
                 I_total=sapply(I_indices, function(i) sum(unlist(lapply(Is, function(w) w[i]))))) %>%
  pivot_longer(cols=-rep, names_to="parameter", values_to="value") %>% 
  ggplot() +
  geom_point(aes(x=parameter,y=value,col=rep),position=position_jitter(w=0.1, h=0)) +
  scale_x_discrete(labels="<i>I</i><sub>total</sub>") +
  scale_y_continuous(limits=c(2.5e6-(5e6-2.5e6)*.05,5e6+(5e6-2.5e6)*.05),
                     expand=expansion(c(0,0)),
                     labels=fancy_scientific) +
  xlab("") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(linewidth=.3),
        legend.position="none",
        axis.text.x=element_markdown())

p2 <- sero_model_summs %>%
  lapply(function(s) s[c("sigma","psi","r_SI","r_VP_init","r_VP_end")]) %>%
  do.call(what=rbind) %>%
  as.data.frame() %>% 
  cbind(rep=factor(1:sero_model_reps)) %>% 
  pivot_longer(cols=-rep, names_to="parameter", values_to="value") %>%
  mutate(parameter=factor(parameter, levels=c("sigma","psi","r_SI","r_VP_init","r_VP_end"))) %>% 
  ggplot() +
  geom_point(aes(x=parameter,y=value,col=rep),position=position_jitter(w=0.1, h=0)) +
  scale_x_discrete(labels=c("<i>\u03C3</i>",
                            "<i>\u03C8</i>",
                            "<i>r</i><sub><i>SI</i></sub>",
                            "<i>r</i><sub><i>VP</i>,init</sub>",
                            "<i>r</i><sub><i>VP</i>,end</sub>")) +
  scale_y_continuous(limits=c(0,2.6),breaks=c(seq(0,2.5,.5))) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(linewidth=.3),
        legend.position="none",
        axis.title.y=element_blank(),
        axis.text.x=element_markdown())

p3 <- data.frame(rep=factor(1:sero_model_reps),
                 w_switch=get_week_start_from_test_week(unlist(lapply(sero_model_summs, function(s) s["w_switch"])))) %>% 
  pivot_longer(cols=-rep, names_to="parameter", values_to="value") %>% 
  ggplot() +
  geom_point(aes(x=parameter,y=value,col=rep),position=position_jitter(w=0.1, h=0)) +
  scale_x_discrete(labels="<i>w</i><sub>switch</sub>") +
  scale_y_date(limits=c(as.Date("2021-06-12")-as.numeric(as.Date("2021-07-07")-as.Date("2021-06-12"))*.05,
                        as.Date("2021-07-07")+as.numeric(as.Date("2021-07-07")-as.Date("2021-06-12"))*.05),
               expand=expansion(c(0,0)),
               breaks=seq(as.Date("2021-06-12"),as.Date("2021-07-07"),5),
               labels=function(d) str_remove_all(format(d,"%D"),"0(?=[1-9])"),
               position="right") +
  xlab("") +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(linewidth=.3),
        legend.position="none",
        axis.title.y=element_blank(),
        axis.text.x=element_markdown())

pg <- plot_grid(p1,p2,p3,nrow=1,rel_widths=c(.5,1,.4))

title <-
  ggdraw() +
  draw_label(bquote("parameter fits for" ~ italic(r)[italic(SV)] == .(r_SV)),
             fontface="bold")

# uncomment the next two commented-out lines in order to save pdf with the special characters coming out right
# quartz(type = "pdf", file = "figures/supplement/parameter_fits.pdf", width=par()$fin[1], height=par()$fin[2])
plot_grid(title,pg,ncol=1,rel_heights=c(.1,1))
# dev.off()






## predicated vaccination probability by past PCR positive over time

summ <- sero_model_summs[[sero_model_reps]]
model_sero %>% 
  cbind(vac_pred=summ[grep("^vac_pred", names(summ))]) %>% 
  filter(test_week>=w_vac-2) %>% 
  group_by(test_week, prev_PCR=factor(ifelse(prev_PCR,"yes","no"),levels=c("yes","no"))) %>%
  summarise(y=mean(vac_pred)) %>%
  ungroup() %>% 
  mutate(date=get_week_start_from_test_week(test_week)) %>% 
  ggplot() +
  geom_line(data=data.frame(date=get_week_start_from_test_week(w_vac:W),
                            prop_vac=cumsum(V[w_vac:W])/N),
            aes(x=date,y=prop_vac), lty="dashed") +
  geom_line(aes(date,y,group=prev_PCR,color=prev_PCR)) +
  ggtitle(bquote("predicted vaccination probabilities for" ~ italic(r)[italic(SV)] == .(r_SV))) +
  scale_color_manual(name="past PCR positive", values=c("#00BFC4", "#F8766D")) +
  scale_x_date(name="date",date_labels="%b '%y",date_breaks="1 month",expand=expansion(c(0,0))) +
  scale_y_continuous(name="probability of vaccination",limits=c(0,1),expand=expansion(c(0,0))) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(linewidth=.3),
        axis.text.x = element_text(angle = 90),
        plot.title=element_text(hjust=.5))

# same thing but with heatmap of r_VP below
p1 <-
  model_sero %>% 
  cbind(vac_pred=summ[grep("^vac_pred", names(summ))]) %>% 
  filter(test_week>=w_vac) %>% 
  group_by(test_week, prev_PCR=factor(ifelse(prev_PCR,"yes","no"),levels=c("yes","no"))) %>%
  summarise(y=mean(vac_pred)) %>%
  ungroup() %>% 
  mutate(date=get_week_start_from_test_week(test_week)) %>% 
  ggplot() +
  geom_line(data=data.frame(date=get_week_start_from_test_week(w_vac:W),
                            prop_vac=cumsum(V[w_vac:W])/N),
            aes(x=date,y=prop_vac), lty="dashed") +
  geom_line(aes(date,y,group=prev_PCR,color=prev_PCR)) +
  scale_color_manual(name="past PCR\npositive", values=c("#00BFC4", "#F8766D")) +
  scale_x_date(date_labels="%b '%y",date_breaks="1 month",expand=expansion(c(0,0))) +
  scale_y_continuous(name="probability of vaccination",limits=c(0,1),expand=expansion(c(0,0))) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(linewidth=.3),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title=element_text(hjust=.5))
r_VP_init <- summ["r_VP_init"]
r_VP_end <- summ["r_VP_end"]
w_switch <- summ["w_switch"]
p2 <-
  data.frame(date = get_week_start_from_test_week(seq(w_vac,W,.01)),
             r_VP = r_VP_init + pmax(0, (seq(w_vac,W,.01)-w_switch) / (W-w_switch) * (r_VP_end-r_VP_init))) %>%
  ggplot() +
  geom_tile(aes(x=date,y=1,fill=r_VP,color=r_VP)) +
  ylab(expression(paste(italic(r[VP]),"(",italic(w),")"))) +
  scale_x_date(date_labels="%b '%y",date_breaks="1 month",expand=expansion(c(0,0))) +
  scale_y_discrete(expand=c(0,0)) +
  scale_fill_viridis(option = "plasma",
                     name = expression(paste(italic(r[VP]),"(",italic(w),")")),
                     breaks = c(.2,.4,.6),
                     guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_color_viridis(option = "plasma",
                     name = expression(paste(italic(r[VP]),"(",italic(w),")")),
                     breaks = c(.2,.4,.6),
                     guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line=element_line(linewidth=.3),
        axis.text.x = element_text(angle = 90),
        axis.title.y = element_text(angle = 0),
        plot.margin = unit(c(0,5.5,-2.3,5.5), "pt"),
        legend.key.height = unit(11, "pt"),
        legend.key.width = unit(17, "pt"))

plot_grid(p1, p2, ncol=1, align="v", rel_heights=c(2.6,1))





## proportions vaccinated by vaccination age groups, used in serology model

prop_vaccinated_by_age %>%
  ggplot() +
  geom_line(aes(get_week_start_from_test_week(test_week),
                prop_vaccinated,group=age_group_min,col=factor(age_group_min)),
            linewidth=.8) +
  scale_x_date(name="date",date_labels="%b '%y",date_breaks="1 month",expand=expansion(c(0,0))) +
  scale_y_continuous(name="proportion vaccinated",limits=c(0,1),expand=expansion(c(.003,0))) +
  scale_color_discrete(name="age group",
                       labels=c("0-4","5-11","12-15","16-17","18-29","30-49","50-64","65-79","80+")) +
  geom_hline(yintercept=prop_vaccinated_cap, lty="dashed",alpha=.5) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(linewidth=.3),
        axis.text.x = element_text(angle = 90))



## proportions of people who have ever been infected by infection age groups, used in serology model

# using the mean estimate for I (I_best) and
# assuming infections are at all times split equally between different age groups
# (see next two plots for justification of this)
prop_infected_by_age(I_best) %>%
  ggplot() +
  geom_line(aes(get_week_start_from_test_week(test_week),
                prop_infected,group=age_group_min,col=factor(age_group_min)),
            linewidth=.8) +
  scale_x_date(name="date",date_labels="%b '%y",date_breaks="1 month",expand=expansion(c(0,0))) +
  scale_y_continuous(name="estimated cumulative proportion infected in NJ",limits=c(0,1),expand=expansion(c(.003,0))) +
  scale_color_discrete(name="age group",
                       labels=c("0-17","18-49","50-64","65+")) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(linewidth=.3),
        axis.text.x = element_text(angle = 90))



# measured NJ seroprevalence by age
NJ_seroprevalence_by_age_raw %>%
  ggplot() +
  geom_line(aes(date, prevalence, col=factor(age_group_min)),
            linewidth=.8) +
  scale_x_date(name="date",date_labels="%b '%y",date_breaks="1 month",expand=expansion(c(0,0))) +
  scale_y_continuous(name="measured seroprevalence in NJ",limits=c(0,1),expand=expansion(c(0,0))) +
  scale_color_discrete(name="age group",
                       labels=c("0-17","18-49","50-64","65+")) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(linewidth=.3),
        axis.text.x = element_text(angle = 90))

# NJ seroprevalence distributions between different age groups over time
# Because it these are relatively constant, we average these proportions for estimating
# cumulative infections per age group
NJ_seroprevalence_by_age_raw %>%
  left_join(group_census_data(prop_infected_age_group_mins)) %>%
  mutate(num_infected=prevalence*pop) %>%
  group_by(date) %>%
  reframe(age_group_min=age_group_min,
          num_infected=num_infected,
          pop=pop,
          frac_of_infecteds=num_infected/sum(num_infected)) %>%
  ggplot() +
  geom_bar(aes(date, frac_of_infecteds, fill=factor(age_group_min)),
           col="black", linewidth=.1, stat="identity", width=16) +
  scale_x_date(name="date",date_labels="%b '%y",date_breaks="1 month",expand=expansion(c(.02,.02))) +
  scale_y_continuous(name="measured seroprevalence age distributions in NJ",limits=c(0,1),expand=expansion(c(0,.001))) +
  scale_fill_discrete(name="age group",
                      labels=c("0-17","18-49","50-64","65+")) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(linewidth=.3),
        axis.text.x = element_text(angle = 90))




## proportions of people in model_sero who have had past PCR positive, by whether they are younger or older than 40

model_sero %>%
  group_by(age_group=ifelse(age>=40, "older than 40", "younger than 40"),
           month=floor_date(min_date %m+% months(test_month-1), "month")) %>%
  summarise(prop_prev_PCR = sum(prev_PCR) / n()) %>%
  ggplot() +
  geom_line(aes(month, prop_prev_PCR, col=age_group),linewidth=.7) +
  scale_color_manual(name="age group", values=c("limegreen","brown2")) +
  theme_classic() +
  scale_x_date(date_labels="%b '%y",date_breaks="1 month") +
  scale_y_continuous(expand=c(0,0), limits=c(0,.3001)) +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("test month") +
  ylab("proportion of serology tests\nwith past PCR positive")





