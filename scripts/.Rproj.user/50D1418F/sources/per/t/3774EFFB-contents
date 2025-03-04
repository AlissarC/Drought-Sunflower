#########################################
# packages necessary
#########################################

library(lme4)
library(car)
library(emmeans)
library(ggeffects)
library(RColorBrewer)
library(multcompView)
library(nlme)
library(marginaleffects)
library(rstantools)
library(multcomp)
library(treemapify)
library(relaimpo)
library(r2glmm)
library(patchwork)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(MuMIn)
library(boot)
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(ggfortify)
library(lavaan)
library(patchwork)
library(dplyr)
library(piecewiseSEM)

############## load data #################################################################### 
stat_sunflower <- read.csv("../output/sunflower.csv")
names(stat_sunflower)

hist(log(stat_sunflower$vcmax_tleaf))# normal distribution
hist(log(stat_sunflower$jmax_tleaf))# quasi normal distribution
hist(log(stat_sunflower$LMA)) 
hist(log(stat_sunflower$Nmass))
hist(stat_sunflower$Nmass)
hist(log(stat_sunflower$Narea))
hist(log(stat_sunflower$nrubisco))
hist(log(stat_sunflower$nbioe))
hist(stat_sunflower$chl_gm2)
hist(stat_sunflower$chi)
hist(stat_sunflower$WUE)


########################## Analysis for chi #####################################################
names(stat_sunflower)
##################### select predrought as baseline values #################### 
baseline_chi_values <- stat_sunflower[stat_sunflower$time == "t1", c("n_trt", "ID_trt", "plant", "chi")]
colnames(baseline_chi_values)[4] <- "baseline_chi"
stat_sunflower <- merge(stat_sunflower, baseline_chi_values, 
            by = c("n_trt", "ID_trt", "plant"), suffixes = c("", "_baseline"))
names(stat_sunflower)
stat_sunflower$centered_chi <- ((stat_sunflower$chi - stat_sunflower$baseline_chi)/stat_sunflower$baseline_chi)*100
hist(stat_sunflower$centered_chi)

####################### Lmer for Chi #####################################################################
stat_sunflower_t234 = subset(stat_sunflower, time!="t1")
hist(stat_sunflower_t234$centered_chi) # normal distribution
ggplot(data=stat_sunflower_t234, aes(x=Cumul_N, y=centered_chi, color= drought, shape = n_trt)) + 
  geom_point()+
  theme_minimal()

chi_lmer <- lmer(centered_chi~ Cumul_N*Sdt_DD  + 
                   (1|plant),
                 data = stat_sunflower_t234)

plot(chi_lmer)
outlierTest(chi_lmer) ## no outliers
shapiro.test(residuals(chi_lmer)) #
Anova(chi_lmer) ##  DD, N effects
AIC(chi_lmer) 
plot(resid(chi_lmer) ~ fitted(chi_lmer)) ## good
r.squaredGLMM(chi_lmer) 
summary(chi_lmer) 
qqnorm(residuals(chi_lmer))
qqline(residuals(chi_lmer))
densityPlot(residuals(chi_lmer))
residuals <- resid(chi_lmer)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
Sdt_DD_range = seq(min(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
             max(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
             0.001)
slope_350_chi <- test(emtrends(chi_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_chi <- test(emtrends(chi_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_chi <- test(emtrends(chi_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_chi <- test(emtrends(chi_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_chi <- summary(emmeans(chi_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_chi <- summary(emmeans(chi_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_chi <- summary(emmeans(chi_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_chi <- summary(emmeans(chi_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_chi_350 <- intercept_350_chi + slope_350_chi * Sdt_DD_range
reg_chi_560 <- intercept_560_chi + slope_560_chi * Sdt_DD_range
reg_chi_1050 <- intercept_1050_chi + slope_1050_chi * Sdt_DD_range
reg_chi_1680 <- intercept_1680_chi + slope_1680_chi * Sdt_DD_range

regline_chi_350 <- data.frame(Sdt_DD = Sdt_DD_range, centered_chi = reg_chi_350, Cumul_N = 350)
regline_chi_560 <- data.frame(Sdt_DD = Sdt_DD_range, centered_chi = reg_chi_560, Cumul_N = 560)
regline_chi_1050 <- data.frame(Sdt_DD = Sdt_DD_range, centered_chi = reg_chi_1050, Cumul_N = 1050)
regline_chi_1680<- data.frame(Sdt_DD = Sdt_DD_range, centered_chi = reg_chi_1680, Cumul_N = 1680)

slope_results_chi <- emtrends(chi_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_chi <- contrast(slope_results_chi, method = "pairwise")
slope_summary_chi <- summary(slope_results_chi, infer = c(TRUE, TRUE))
slope_summary_chi

slopes_chi_test <- cld(emtrends(chi_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

######################## Figure ####################################################################################
Centered_Chi_plot <- ggplot(data = stat_sunflower_t234, aes(x = Sdt_DD, y = centered_chi, color=Cumul_N)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_chi_350, aes(x = Sdt_DD  , y = centered_chi  ), col = 'cyan2', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_chi_560, aes(x = Sdt_DD  , y = centered_chi  ), col = 'cyan3', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_chi_1050, aes(x = Sdt_DD  , y = centered_chi  ), col = 'darkcyan', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_chi_1680, aes(x = Sdt_DD  , y = centered_chi  ), col = 'lightcyan4', lwd = 2, alpha = 0.8) +
  
  theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +

  ylab(expression(Delta * " " * chi * " (" * "%" * ")"))+
  xlab(expression (italic('Drought Severity Index')))
Centered_Chi_plot

names(stat_sunflower)
##################### select predrought as baseline values for gs #################### 
baseline_gs_values <- stat_sunflower[stat_sunflower$time == "t1", c("n_trt", "ID_trt", "plant", "gs_420")]
colnames(baseline_gs_values)[4] <- "baseline_gs"
stat_sunflower <- merge(stat_sunflower, baseline_gs_values, 
                        by = c("n_trt", "ID_trt", "plant"), suffixes = c("", "_baseline"))
names(stat_sunflower)
stat_sunflower$centered_gs <- ((stat_sunflower$gs_420 - stat_sunflower$baseline_gs)/stat_sunflower$baseline_gs)*100
hist(stat_sunflower$centered_gs)

####################### Lmer for Chi #####################################################################
stat_sunflower_t234 = subset(stat_sunflower, time!="t1")
hist(stat_sunflower_t234$centered_gs) # normal distribution
stat_sunflower_t234_gs= subset(stat_sunflower_t234,centered_gs<500)
nrow(stat_sunflower_t234_gs)

ggplot(data=stat_sunflower_t234, aes(x=Cumul_N, y=centered_gs, color= drought, shape = n_trt)) + 
  geom_point()+
  theme_minimal()

gs_lmer <- lmer(centered_gs~ Cumul_N*Sdt_DD  + 
                   (1|plant),
                 data = stat_sunflower_t234_gs)

plot(gs_lmer)
outlierTest(gs_lmer) ## no outliers
shapiro.test(residuals(gs_lmer)) 
Anova(gs_lmer) 
AIC(gs_lmer) 
plot(resid(gs_lmer) ~ fitted(gs_lmer)) ## good
r.squaredGLMM(gs_lmer) 
summary(gs_lmer) 
qqnorm(residuals(gs_lmer))
qqline(residuals(gs_lmer))
densityPlot(residuals(gs_lmer))
residuals <- resid(gs_lmer)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
Sdt_DD_range = seq(min(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   max(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_gs <- test(emtrends(gs_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_gs <- test(emtrends(gs_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_gs <- test(emtrends(gs_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_gs <- test(emtrends(gs_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_gs <- summary(emmeans(gs_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_gs <- summary(emmeans(gs_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_gs <- summary(emmeans(gs_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_gs <- summary(emmeans(gs_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_gs_350 <- intercept_350_gs + slope_350_gs * Sdt_DD_range
reg_gs_560 <- intercept_560_gs + slope_560_gs * Sdt_DD_range
reg_gs_1050 <- intercept_1050_gs + slope_1050_gs * Sdt_DD_range
reg_gs_1680 <- intercept_1680_gs + slope_1680_gs * Sdt_DD_range

regline_gs_350 <- data.frame(Sdt_DD = Sdt_DD_range, centered_gs = reg_gs_350, Cumul_N = 350)
regline_gs_560 <- data.frame(Sdt_DD = Sdt_DD_range, centered_gs = reg_gs_560, Cumul_N = 560)
regline_gs_1050 <- data.frame(Sdt_DD = Sdt_DD_range, centered_gs = reg_gs_1050, Cumul_N = 1050)
regline_gs_1680<- data.frame(Sdt_DD = Sdt_DD_range, centered_gs = reg_gs_1680, Cumul_N = 1680)

slope_results_gs <- emtrends(gs_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_gs <- contrast(slope_results_gs, method = "pairwise")
slope_summary_gs <- summary(slope_results_gs, infer = c(TRUE, TRUE))
slope_summary_gs

slopes_gs_test <- cld(emtrends(gs_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

######################## Figure ####################################################################################
Centered_gs_plot <- ggplot(data = stat_sunflower_t234_gs, aes(x = Sdt_DD, y = centered_gs, color=Cumul_N)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_gs_350, aes(x = Sdt_DD  , y = centered_gs), col = 'cyan2', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_gs_560, aes(x = Sdt_DD  , y = centered_gs), col = 'cyan3', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_gs_1050, aes(x = Sdt_DD  , y = centered_gs), col = 'darkcyan', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_gs_1680, aes(x = Sdt_DD  , y = centered_gs), col = 'lightcyan4', lwd = 2, alpha = 0.8,linetype ="dashed") +
  
  theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  
  ylab(expression(Delta * " " * gs * " (" * "%" * ")"))+
  xlab(expression (italic('Drought Severity Index')))
Centered_gs_plot

########################## Analysis for leaf beta leaf #####################################################
names(stat_sunflower)

baseline_betaleaf_values <- stat_sunflower[stat_sunflower$time == "t1", c("n_trt", "ID_trt", "plant", "betaleaf")]
colnames(baseline_betaleaf_values)[4] <- "baseline_betaleaf"
stat_sunflower <- merge(stat_sunflower, baseline_betaleaf_values, 
                        by = c("n_trt", "ID_trt", "plant"), suffixes = c("", "_baseline"))
names(stat_sunflower)
stat_sunflower$centered_betaleaf <- ((stat_sunflower$betaleaf - stat_sunflower$baseline_betaleaf)/stat_sunflower$baseline_betaleaf)*100
hist(stat_sunflower$centered_betaleaf)

####################### Lmer for beta leaf #################
stat_sunflower_t234 = subset(stat_sunflower, time!="t1")
hist(log(stat_sunflower_t234$centered_betaleaf))## Normal dist
nrow(stat_sunflower_t234)

stat_sunflower_beta = subset(stat_sunflower_t234, centered_betaleaf<200)
nrow(stat_sunflower_beta)

Betaleaf_lmer <- lmer(centered_betaleaf~ Cumul_N*Sdt_DD +
                        (1|plant),
                      data = stat_sunflower_beta)
outlierTest(Betaleaf_lmer)  ##########
plot(Betaleaf_lmer)
Anova(Betaleaf_lmer) 
AIC(Betaleaf_lmer) 
vif(Betaleaf_lmer) 
plot(resid(Betaleaf_lmer) ~ fitted(Betaleaf_lmer)) ## good
r.squaredGLMM(Betaleaf_lmer) 

summary(Betaleaf_lmer) 
qqnorm(residuals(Betaleaf_lmer))
qqline(residuals(Betaleaf_lmer))
densityPlot(residuals(Betaleaf_lmer))
shapiro.test(residuals(Betaleaf_lmer))

cumul_N_levels_Betaleaf<- c(350, 560, 1050, 1680)

Sdt_DD_range = seq(min(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   max(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   0.001)

slope_350_Betaleaf <- test(emtrends(Betaleaf_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_Betaleaf <- test(emtrends(Betaleaf_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_Betaleaf <- test(emtrends(Betaleaf_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_Betaleaf <- test(emtrends(Betaleaf_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_Betaleaf <- summary(emmeans(Betaleaf_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_Betaleaf <- summary(emmeans(Betaleaf_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_Betaleaf <- summary(emmeans(Betaleaf_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_Betaleaf <- summary(emmeans(Betaleaf_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_Betaleaf_350 <- intercept_350_Betaleaf + slope_350_Betaleaf * Sdt_DD_range
reg_Betaleaf_560 <- intercept_560_Betaleaf + slope_560_Betaleaf* Sdt_DD_range
reg_Betaleaf_1050 <- intercept_1050_Betaleaf + slope_1050_Betaleaf * Sdt_DD_range
reg_Betaleaf_1680 <- intercept_1680_Betaleaf + slope_1680_Betaleaf * Sdt_DD_range

regline_Betaleaf_350 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Betaleaf= reg_Betaleaf_350, Cumul_N = 350)
regline_Betaleaf_560 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Betaleaf = reg_Betaleaf_560, Cumul_N = 525)
regline_Betaleaf_1050 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Betaleaf = reg_Betaleaf_1050, Cumul_N = 1050)
regline_Betaleaf_1680 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Betaleaf = reg_Betaleaf_1680, Cumul_N = 1680)

slope_results_Betaleaf <- emtrends(Betaleaf_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_Betaleaf <- contrast(slope_results_Betaleaf, method = "pairwise")
slope_summary_Betaleaf <- summary(slope_results_Betaleaf, infer = c(TRUE, TRUE))
slope_summary_Betaleaf

cld(emtrends(Betaleaf_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))
################  plot #####################################################################
Centered_Betaleaf_plot <- ggplot(data = stat_sunflower_beta, aes(x = Sdt_DD, y = centered_betaleaf)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_Betaleaf_350, aes(x = Sdt_DD  , y = centered_Betaleaf), col = 'cyan2', lwd = 2, alpha = 0.8, , linetype="dashed") +
  geom_line(data = regline_Betaleaf_560, aes(x = Sdt_DD  , y = centered_Betaleaf), col = 'cyan3', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Betaleaf_1050, aes(x = Sdt_DD  , y = centered_Betaleaf), col = 'darkcyan', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Betaleaf_1680, aes(x = Sdt_DD  , y = centered_Betaleaf), col = 'lightcyan4', lwd = 2, alpha = 0.8) +
  
  
  theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  
  ylab(expression('∆ '* italic('β') ['leaf'] * ' (%)')) + 
  xlab(expression (italic('Drought Severity Index')))
Centered_Betaleaf_plot

########################## Analysis for Nmass #####################################################
##################### select predrought as baseline values #################### 
baseline_Nmass_values <- stat_sunflower[stat_sunflower$time == "t1", c("n_trt", "ID_trt", "plant", "Nmass")]
colnames(baseline_Nmass_values)[4] <- "baseline_Nmass"
stat_sunflower <- merge(stat_sunflower, baseline_Nmass_values, 
                        by = c("n_trt", "ID_trt", "plant"), suffixes = c("", "_baseline"))
names(stat_sunflower)
stat_sunflower$centered_Nmass <- ((stat_sunflower$Nmass - stat_sunflower$baseline_Nmass)/stat_sunflower$baseline_Nmass)*100
hist(stat_sunflower$centered_Nmass)

####################### Lmer for  centerd Nmass
stat_sunflower_t234 = subset(stat_sunflower, time!="t1")
hist(stat_sunflower_t234$centered_Nmass) # Normal distribution
nrow(stat_sunflower_t234)

Nmass_lmer <- lmer(centered_Nmass~ Cumul_N*Sdt_DD + 
                     (1|plant),
                 data = stat_sunflower_t234)
plot(Nmass_lmer)
outlierTest(Nmass_lmer)  
Anova(Nmass_lmer) 
shapiro.test(residuals(Nmass_lmer)) 
AIC(Nmass_lmer) 
vif(Nmass_lmer) 
plot(resid(Nmass_lmer) ~ fitted(Nmass_lmer)) 
r.squaredGLMM(Nmass_lmer) 

summary(Nmass_lmer) 
qqnorm(residuals(Nmass_lmer))
qqline(residuals(Nmass_lmer))
densityPlot(residuals(Nmass_lmer))

cumul_N_levels <- c(350, 560, 1050, 1560)
Sdt_DD_range = seq(min(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   max(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_Nmass <- test(emtrends(Nmass_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_Nmass  <- test(emtrends(Nmass_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_Nmass  <- test(emtrends(Nmass_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_Nmass  <- test(emtrends(Nmass_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_Nmass  <- summary(emmeans(Nmass_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_Nmass  <- summary(emmeans(Nmass_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_Nmass  <- summary(emmeans(Nmass_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_Nmass  <- summary(emmeans(Nmass_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_Nmass_350 <- intercept_350_Nmass  + slope_350_Nmass* Sdt_DD_range
reg_Nmass_560  <- intercept_560_Nmass  + slope_560_Nmass* Sdt_DD_range
reg_Nmass_1050 <- intercept_1050_Nmass  + slope_1050_Nmass*Sdt_DD_range
reg_Nmass_1680  <- intercept_1680_Nmass  + slope_1680_Nmass*Sdt_DD_range

regline_Nmass_350 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Nmass = reg_Nmass_350, Cumul_N = 350)
regline_Nmass_560 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Nmass = reg_Nmass_560, Cumul_N = 560)
regline_Nmass_1050 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Nmass = reg_Nmass_1050, Cumul_N = 1050)
regline_Nmass_1680 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Nmass = reg_Nmass_1680, Cumul_N = 1680)

slope_results_Nmass  <- emtrends(Nmass_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_Nmass  <- contrast(slope_results_Nmass , method = "pairwise")
slope_summary_Nmass  <- summary(slope_results_Nmass , infer = c(TRUE, TRUE))
slope_summary_Nmass 

slopes_Nmass_test <- cld(emtrends(Nmass_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))
################ Nmass plot #####################################################################
Centred_Nmass_plot <- ggplot(data = stat_sunflower_t234, aes(x = Sdt_DD, y = centered_Nmass)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  
  geom_line(data = regline_Nmass_350, aes(x = Sdt_DD  , y = centered_Nmass  ), col = 'cyan2', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Nmass_560, aes(x = Sdt_DD  , y = centered_Nmass  ), col = 'cyan3', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Nmass_1050, aes(x = Sdt_DD  , y = centered_Nmass  ), col = 'darkcyan', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Nmass_1680, aes(x = Sdt_DD  , y = centered_Nmass  ), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype ="dashed") +
  
  theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  
  ylab(expression('∆ '* italic('N')['mass'] * ' (%)')) +
  xlab(expression (italic('Drought Severity Index')))
Centred_Nmass_plot

######################################################################################################################
########################## Analysis for Narea ########################################################################
baseline_Narea_values <- stat_sunflower[stat_sunflower$time == "t1", c("n_trt", "ID_trt", "plant", "Narea")]
colnames(baseline_Narea_values)[4] <- "baseline_Narea"
stat_sunflower <- merge(stat_sunflower, baseline_Narea_values, 
                        by = c("n_trt", "ID_trt", "plant"), suffixes = c("", "_baseline"))
names(stat_sunflower)
stat_sunflower$centered_Narea <- ((stat_sunflower$Narea - stat_sunflower$baseline_Narea)/stat_sunflower$baseline_Narea)*100
hist(stat_sunflower$centered_Narea)
####################### Lmer for  centerd Narea #################
stat_sunflower_t234 = subset(stat_sunflower, time!="t1")
hist(stat_sunflower_t234$centered_Narea)
Narea_lmer <- lmer(centered_Narea~ Cumul_N*Sdt_DD + 
                     (1|plant),
                   data = stat_sunflower_t234)
outlierTest(Narea_lmer) # One outlier to remove 
plot(Narea_lmer)
shapiro.test(residuals(Narea_lmer))
residuals_data <- data.frame(
  row_index = 1:nrow(stat_sunflower_t234),
  residuals = rstudent(Narea_lmer)
)
hist(residuals_data$residuals)
threshold <- 4
outliers <- abs(residuals_data$residuals) >threshold
# Filter the dataset to remove outliers
stat_sunflower_t234_clean <- stat_sunflower_t234[!outliers, ]
nrow(stat_sunflower_t234_clean)

hist(stat_sunflower_t234_clean$centered_Narea) # Normal dist
Narea_lmer <- lmer(centered_Narea~ Cumul_N*Sdt_DD + 
                     (1|plant),
                   data = stat_sunflower_t234_clean)
outlierTest(Narea_lmer) # No outlier
plot(Narea_lmer)
shapiro.test(residuals(Narea_lmer))
Anova(Narea_lmer)
AIC(Narea_lmer) 
vif(Narea_lmer) 
plot(resid(Narea_lmer) ~ fitted(Narea_lmer)) ## good
r.squaredGLMM(Narea_lmer) 

summary(Narea_lmer) 
qqnorm(residuals(Narea_lmer))
qqline(residuals(Narea_lmer))
densityPlot(residuals(Narea_lmer))
shapiro.test(residuals(Narea_lmer)) 

cumul_N_levels <- c(350, 560, 1050, 1680)

Sdt_DD_range = seq(min(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   max(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_Narea  <- test(emtrends(Narea_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_Narea  <- test(emtrends(Narea_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_Narea  <- test(emtrends(Narea_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_Narea  <- test(emtrends(Narea_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_Narea  <- summary(emmeans(Narea_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_Narea  <- summary(emmeans(Narea_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_Narea  <- summary(emmeans(Narea_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_Narea  <- summary(emmeans(Narea_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_Narea_350 <- intercept_350_Narea  + slope_350_Narea  * Sdt_DD_range
reg_Narea_560 <- intercept_560_Narea  + slope_560_Narea  * Sdt_DD_range
reg_Narea_1050 <- intercept_1050_Narea  + slope_1050_Narea  * Sdt_DD_range
reg_Narea_1680 <- intercept_1680_Narea  + slope_1680_Narea  * Sdt_DD_range

regline_Narea_350 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Narea = reg_Narea_350, Cumul_N = 350)
regline_Narea_560 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Narea = reg_Narea_560, Cumul_N = 525)
regline_Narea_1050 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Narea = reg_Narea_1050, Cumul_N = 1050)
regline_Narea_1680 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Narea = reg_Narea_1680, Cumul_N = 1575)

slope_results_Narea  <- emtrends(Narea_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_Narea  <- contrast(slope_results_Narea , method = "pairwise")
slope_summary_Narea  <- summary(slope_results_Narea , infer = c(TRUE, TRUE))
slope_summary_Narea 

slopes_Narea_test <- cld(emtrends(Narea_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))
################ Narea plot #####################################################################
Centered_Narea_plot <- ggplot(data = stat_sunflower_t234_clean, aes(x = Sdt_DD, y = centered_Narea)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_Narea_350, aes(x = Sdt_DD  , y = centered_Narea), col = 'cyan2', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_Narea_560, aes(x = Sdt_DD  , y = centered_Narea), col = 'cyan3', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_Narea_1050, aes(x = Sdt_DD  , y = centered_Narea), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_Narea_1680, aes(x = Sdt_DD  , y = centered_Narea), col = 'lightcyan4', lwd = 2, alpha = 0.8) +
  
  theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  
  ylab(expression('∆ '* italic('N')['area'] * ' (%)')) +
  xlab(expression (italic('Drought Severity Index')))
Centered_Narea_plot

########################## Analysis for LMA #####################################################
baseline_LMA_values <- stat_sunflower[stat_sunflower$time == "t1", c("n_trt", "ID_trt", "plant", "LMA")]
colnames(baseline_LMA_values)[4] <- "baseline_LMA"
stat_sunflower <- merge(stat_sunflower, baseline_LMA_values, 
                        by = c("n_trt", "ID_trt", "plant"), suffixes = c("", "_baseline"))
names(stat_sunflower)
stat_sunflower$centered_LMA <- ((stat_sunflower$LMA - stat_sunflower$baseline_LMA)/stat_sunflower$baseline_LMA)*100
hist(stat_sunflower$centered_LMA)
####################### Lmer for  centerd Narea #################
stat_sunflower_t234 = subset(stat_sunflower, time!="t1")
hist(stat_sunflower_t234$centered_LMA) # Normal distribution

LMA_lmer <- lmer(centered_LMA~ Cumul_N*Sdt_DD +
                     (1|plant),
                   data = stat_sunflower_t234)
plot(LMA_lmer) 
outlierTest(LMA_lmer) 
shapiro.test(residuals(LMA_lmer)) 
Anova(LMA_lmer) #
AIC(LMA_lmer) 
vif(LMA_lmer) 
plot(resid(LMA_lmer) ~ fitted(LMA_lmer)) 
r.squaredGLMM(LMA_lmer) 
summary(LMA_lmer) 
qqnorm(residuals(LMA_lmer))
qqline(residuals(LMA_lmer))
densityPlot(residuals(LMA_lmer))

cumul_N_levels <- c(350, 560, 1050, 1680)

Sdt_DD_range = seq(min(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   max(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_LMA <- test(emtrends(LMA_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_LMA <- test(emtrends(LMA_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_LMA <- test(emtrends(LMA_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_LMA <- test(emtrends(LMA_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_LMA <- summary(emmeans(LMA_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_LMA <- summary(emmeans(LMA_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_LMA <- summary(emmeans(LMA_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_LMA <- summary(emmeans(LMA_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_LMA_350 <- intercept_350_LMA + slope_350_LMA * Sdt_DD_range
reg_LMA_560 <- intercept_560_LMA + slope_560_LMA * Sdt_DD_range
reg_LMA_1050 <- intercept_1050_LMA + slope_1050_LMA * Sdt_DD_range
reg_LMA_1680 <- intercept_1680_LMA + slope_1680_LMA * Sdt_DD_range

regline_LMA_350 <- data.frame(Sdt_DD = Sdt_DD_range, centered_LMA = reg_LMA_350, Cumul_N = 350)
regline_LMA_560 <- data.frame(Sdt_DD = Sdt_DD_range, centered_LMA = reg_LMA_560, Cumul_N = 560)
regline_LMA_1050 <- data.frame(Sdt_DD = Sdt_DD_range, centered_LMA = reg_LMA_1050, Cumul_N = 1050)
regline_LMA_1680 <- data.frame(Sdt_DD = Sdt_DD_range, centered_LMA = reg_LMA_1680, Cumul_N = 1680)

slope_results_LMA <- emtrends(LMA_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_LMA <- contrast(slope_results_LMA, method = "pairwise")
slope_summary_LMA <- summary(slope_results_LMA, infer = c(TRUE, TRUE))
slope_summary_LMA

slopes_LMA_test <- cld(emtrends(LMA_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))
################ LMA plot #####################################################################
Centered_LMA_plot <- ggplot(data = stat_sunflower_t234, aes(x = Sdt_DD, y = centered_LMA)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_LMA_350, aes(x = Sdt_DD  , y = centered_LMA), col = 'cyan2', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_LMA_560, aes(x = Sdt_DD  , y = centered_LMA), col = 'cyan3', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_LMA_1050, aes(x = Sdt_DD  , y = centered_LMA), col = 'darkcyan', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_LMA_1680, aes(x = Sdt_DD  , y = centered_LMA), col = 'lightcyan4', lwd = 2, alpha = 0.8) +

  theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  
  ylab(expression('∆ '* italic('LMA')[''] * ' (%)')) +
  xlab(expression (italic('Drought Severity Index')))
Centered_LMA_plot

###################################################################################################################
########################## Analysis for Vcmax #####################################################################
baseline_Vcmax_values <- stat_sunflower[stat_sunflower$time == "t1", c("n_trt", "ID_trt", "plant", "vcmax_tleaf")]
colnames(baseline_Vcmax_values)[4] <- "baseline_Vcmax"
stat_sunflower <- merge(stat_sunflower, baseline_Vcmax_values, 
                        by = c("n_trt", "ID_trt", "plant"), suffixes = c("", "_baseline"))
names(stat_sunflower)
stat_sunflower$centered_Vcmax <- ((stat_sunflower$vcmax_tleaf - stat_sunflower$baseline_Vcmax)/stat_sunflower$baseline_Vcmax)*100
####################### Lmer for  centerd Vcmax #################
stat_sunflower_t234 = subset(stat_sunflower, time!="t1")
nrow(stat_sunflower_t234)
hist(stat_sunflower_t234$centered_Vcmax) ## some extreme values 

Vcmax_lmer <- lmer(centered_Vcmax~ Cumul_N*Sdt_DD + 
                     (1|plant),
                 data = stat_sunflower_t234)
outlierTest(Vcmax_lmer)  ## two outliers to remove
plot(Vcmax_lmer)
shapiro.test(residuals(Vcmax_lmer))
residuals_data <- data.frame(
  row_index = 1:nrow(stat_sunflower_t234),
  residuals = rstudent(Vcmax_lmer)
)
hist(residuals_data$residuals)
threshold <- 3
outliers <- abs(residuals_data$residuals) >threshold
# Filter the dataset to remove outliers
stat_sunflower_vcmax_clean <- stat_sunflower_t234[!outliers, ]
nrow(stat_sunflower_vcmax_clean) # 142 , three removed
hist(stat_sunflower_vcmax_clean$centered_Vcmax) # better dist
######## refit
Vcmax_lmer <- lmer(centered_Vcmax~ Cumul_N*Sdt_DD + 
                     (1|plant),
                   data = stat_sunflower_vcmax_clean)
outlierTest(Vcmax_lmer) ## P = 0.082, perfect
plot(Vcmax_lmer)
shapiro.test(residuals(Vcmax_lmer))
Anova(Vcmax_lmer) 
AIC(Vcmax_lmer) 
vif(Vcmax_lmer) 
plot(resid(Vcmax_lmer) ~ fitted(Vcmax_lmer)) 
r.squaredGLMM(Vcmax_lmer) 

summary(Vcmax_lmer) 
qqnorm(residuals(Vcmax_lmer))
qqline(residuals(Vcmax_lmer))
densityPlot(residuals(Vcmax_lmer))

cumul_N_levels <- c(350, 560, 1050, 1680)
Sdt_DD_range = seq(min(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   max(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_Vcmax <- test(emtrends(Vcmax_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_Vcmax <- test(emtrends(Vcmax_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_Vcmax <- test(emtrends(Vcmax_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_Vcmax <- test(emtrends(Vcmax_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_Vcmax <- summary(emmeans(Vcmax_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_Vcmax <- summary(emmeans(Vcmax_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_Vcmax <- summary(emmeans(Vcmax_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_Vcmax <- summary(emmeans(Vcmax_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_Vcmax_350 <- intercept_350_Vcmax + slope_350_Vcmax * Sdt_DD_range
reg_Vcmax_560 <- intercept_560_Vcmax + slope_560_Vcmax * Sdt_DD_range
reg_Vcmax_1050 <- intercept_1050_Vcmax + slope_1050_Vcmax * Sdt_DD_range
reg_Vcmax_1680 <- intercept_1680_Vcmax + slope_1680_Vcmax * Sdt_DD_range

regline_Vcmax_350 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Vcmax = reg_Vcmax_350, Cumul_N = 350)
regline_Vcmax_560 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Vcmax = reg_Vcmax_560, Cumul_N = 560)
regline_Vcmax_1050 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Vcmax = reg_Vcmax_1050, Cumul_N = 1050)
regline_Vcmax_1680 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Vcmax = reg_Vcmax_1680, Cumul_N = 1680)

slope_results_Vcmax <- emtrends(Vcmax_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_Vcmax<- contrast(slope_results_Vcmax, method = "pairwise")
slope_summary_Vcmax <- summary(slope_results_Vcmax, infer = c(TRUE, TRUE))
slope_summary_Vcmax

slopes_Vcmax_test <- cld(emtrends(Vcmax_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))
################ Vcmax plot #####################################################################
Centred_Vcmax_plot <- ggplot(data = stat_sunflower_vcmax_clean, aes(x = Sdt_DD, y = centered_Vcmax)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_Vcmax_350, aes(x = Sdt_DD  , y = centered_Vcmax), col = 'cyan2', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Vcmax_560, aes(x = Sdt_DD  , y = centered_Vcmax), col = 'cyan3', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Vcmax_1050, aes(x = Sdt_DD  , y = centered_Vcmax), col = 'darkcyan', lwd = 2, alpha = 0.8,linetype ="dashed" ) +
  geom_line(data = regline_Vcmax_1680, aes(x = Sdt_DD  , y = centered_Vcmax), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype ="dashed") +
 
   theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  
  ylab(expression('∆ '* italic('V')['cmax'] * ' (%)')) +
  xlab(expression (italic('Drought Severity Index')))
Centred_Vcmax_plot
#################################################################################################################
################################################################################################################

########################## Analysis for Jmax ########################################################################
baseline_Jmax_values <- stat_sunflower[stat_sunflower$time == "t1", c("n_trt", "ID_trt", "plant", "jmax_tleaf")]
colnames(baseline_Jmax_values)[4] <- "baseline_Jmax"
stat_sunflower <- merge(stat_sunflower, baseline_Jmax_values, 
                        by = c("n_trt", "ID_trt", "plant"), suffixes = c("", "_baseline"))
names(stat_sunflower)
stat_sunflower$centered_Jmax <- ((stat_sunflower$jmax_tleaf - stat_sunflower$baseline_Jmax)/stat_sunflower$baseline_Jmax)*100
hist(stat_sunflower$centered_Jmax)
####################### Lmer for  centerd Jmax #################
stat_sunflower_t234 = subset(stat_sunflower, time!="t1")
hist(stat_sunflower_t234$centered_Jmax) ## some extreme values 

Jmax_lmer <- lmer(centered_Jmax~ Cumul_N*Sdt_DD +
                     (1|plant),
                   data = stat_sunflower_t234)
outlierTest(Jmax_lmer)  ### two to remove

residuals_data <- data.frame(
  row_index = 1:nrow(stat_sunflower_t234),
  residuals = rstudent(Jmax_lmer)
)
hist(residuals_data$residuals)
threshold <- 2
outliers <- abs(residuals_data$residuals) > threshold
# Filter the dataset to remove outliers
stat_sunflower_Jmax_clean <- stat_sunflower_t234[!outliers, ]
nrow(stat_sunflower_Jmax_clean) # 139 , three removed
hist(stat_sunflower_Jmax_clean$centered_Jmax) # better dist

############# Refit ###########################
Jmax_lmer <- lmer(centered_Jmax~ Cumul_N*Sdt_DD +
                    (1|plant),
                  data = stat_sunflower_Jmax_clean)
outlierTest(Jmax_lmer) 
shapiro.test(residuals(Jmax_lmer))
plot(Jmax_lmer)
Anova(Jmax_lmer) 
AIC(Jmax_lmer) 
vif(Jmax_lmer) # 
plot(resid(Jmax_lmer) ~ fitted(Jmax_lmer)) 
r.squaredGLMM(Jmax_lmer)
summary(Jmax_lmer) 
qqnorm(residuals(Jmax_lmer))
qqline(residuals(Jmax_lmer))
densityPlot(residuals(Jmax_lmer))

residuals <- resid(Jmax_lmer)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350, 560, 1050, 1680)

Sdt_DD_range = seq(min(stat_sunflower_Jmax_clean$Sdt_DD, na.rm = TRUE), 
                   max(stat_sunflower_Jmax_clean$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_jmax <- test(emtrends(Jmax_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_jmax <- test(emtrends(Jmax_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_jmax <- test(emtrends(Jmax_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_jmax <- test(emtrends(Jmax_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_jmax <- summary(emmeans(Jmax_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_jmax <- summary(emmeans(Jmax_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_jmax <- summary(emmeans(Jmax_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_jmax <- summary(emmeans(Jmax_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_Jmax_350 <- intercept_350_jmax + slope_350_jmax * Sdt_DD_range
reg_Jmax_560 <- intercept_560_jmax + slope_560_jmax * Sdt_DD_range
reg_Jmax_1050 <- intercept_1050_jmax + slope_1050_jmax * Sdt_DD_range
reg_Jmax_1680 <- intercept_1680_jmax + slope_1680_jmax * Sdt_DD_range

regline_Jmax_350 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Jmax = reg_Jmax_350, Cumul_N = 350)
regline_Jmax_560 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Jmax = reg_Jmax_560, Cumul_N = 560)
regline_Jmax_1050 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Jmax = reg_Jmax_1050, Cumul_N = 1050)
regline_Jmax_1680 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Jmax = reg_Jmax_1680, Cumul_N = 1680)

slope_results_jmax <- emtrends(Jmax_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_jmax <- contrast(slope_results_jmax, method = "pairwise")
slope_summary_jmax <- summary(slope_results_jmax, infer = c(TRUE, TRUE))
slope_summary_jmax

slopes_jmax_test <- cld(emtrends(Jmax_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))
################ Jmax plot #####################################################################
Centered_Jmax_plot <- ggplot(data = stat_sunflower_Jmax_clean, aes(x = Sdt_DD, y = centered_Jmax)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_Jmax_350, aes(x = Sdt_DD  , y = centered_Jmax  ), col = 'cyan2', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Jmax_560, aes(x = Sdt_DD  , y = centered_Jmax  ), col = 'cyan3', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Jmax_1050, aes(x = Sdt_DD  , y = centered_Jmax  ), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_Jmax_1680, aes(x = Sdt_DD  , y = centered_Jmax  ), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype ="dashed") +
  
  theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  
  ylab(expression('∆'* italic('J')['max'] * ' (%)')) +
  xlab(expression (italic('Drought Severity Index')))

Centered_Jmax_plot
##################################################################################################################
##################################################################################################################
##################### Anet ############################################################################################## 
baseline_Anet_values <- stat_sunflower[stat_sunflower$time == "t1", c("n_trt", "ID_trt", "plant", "anet_420")]
colnames(baseline_Anet_values)[4] <- "baseline_Anet"
stat_sunflower <- merge(stat_sunflower, baseline_Anet_values, 
                        by = c("n_trt", "ID_trt", "plant"), suffixes = c("", "_baseline"))
names(stat_sunflower)
stat_sunflower$centered_Anet <- ((stat_sunflower$anet_420 - stat_sunflower$baseline_Anet)/stat_sunflower$baseline_Anet)*100
hist(stat_sunflower$centered_Anet)
####################### Lmer for  centerd Anet #################
stat_sunflower_t234 = subset(stat_sunflower,time != "t1")
hist(stat_sunflower_t234$centered_Anet)

Anet_lmer <- lmer(centered_Anet~ Cumul_N*Sdt_DD + 
                    (1|plant),
                  data = stat_sunflower_t234)
outlierTest(Anet_lmer)  ## One outlier to remove
plot(Anet_lmer)

residuals_data <- data.frame(
  row_index = 1:nrow(stat_sunflower_t234),
  residuals = rstudent(Anet_lmer)
)
hist(residuals_data$residuals)
threshold <- 3
outliers <- abs(residuals_data$residuals) > threshold
# Filter the dataset to remove outliers
stat_sunflower_anet_clean <- stat_sunflower_t234[!outliers, ]
nrow(stat_sunflower_anet_clean) # 141 , three removed
hist(stat_sunflower_anet_clean$centered_Anet) # better dist

############### Refit ##########################################
Anet_lmer <- lmer(centered_Anet~ Cumul_N*Sdt_DD + 
                    (1|plant),
                  data = stat_sunflower_anet_clean)
outlierTest(Anet_lmer)  ## p = 0.07
plot(Anet_lmer) 
Anova(Anet_lmer) 
AIC(Anet_lmer) 
vif(Anet_lmer) 
plot(resid(Anet_lmer) ~ fitted(Anet_lmer)) ## good
r.squaredGLMM(Anet_lmer) 

summary(Anet_lmer) 
qqnorm(residuals(Anet_lmer))
qqline(residuals(Anet_lmer))
densityPlot(residuals(Anet_lmer))
shapiro.test(residuals(Anet_lmer))

residuals <- resid(Anet_lmer)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels_anet <- c(350, 560, 1050, 1680)

Sdt_DD_range = seq(min(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   max(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   0.001)

slope_350_anet <- test(emtrends(Anet_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_anet <- test(emtrends(Anet_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_anet <- test(emtrends(Anet_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_anet <- test(emtrends(Anet_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_anet <- summary(emmeans(Anet_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_anet <- summary(emmeans(Anet_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_anet <- summary(emmeans(Anet_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_anet <- summary(emmeans(Anet_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_Anet_350 <- intercept_350_anet + slope_350_anet * Sdt_DD_range
reg_Anet_560 <- intercept_560_anet + slope_560_anet* Sdt_DD_range
reg_Anet_1050 <- intercept_1050_anet + slope_1050_anet * Sdt_DD_range
reg_Anet_1680 <- intercept_1680_anet + slope_1680_anet * Sdt_DD_range

regline_Anet_350 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Anet= reg_Anet_350, Cumul_N = 350)
regline_Anet_560 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Anet = reg_Anet_560, Cumul_N = 560)
regline_Anet_1050 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Anet = reg_Anet_1050, Cumul_N = 1050)
regline_Anet_1680 <- data.frame(Sdt_DD = Sdt_DD_range, centered_Anet = reg_Anet_1680, Cumul_N = 1680)

slope_results_anet <- emtrends(Anet_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_anet <- contrast(slope_results_anet, method = "pairwise")
slope_summary_anet <- summary(slope_results_anet, infer = c(TRUE, TRUE))
slope_summary_anet

slopes_Anet_test <- cld(emtrends(Anet_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))
################ Anet plot #####################################################################
Centered_Anet_plot <- ggplot(data = stat_sunflower_anet_clean, aes(x = Sdt_DD, y = centered_Anet)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_Anet_350, aes(x = Sdt_DD  , y = centered_Anet), col = 'cyan2', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Anet_560, aes(x = Sdt_DD  , y = centered_Anet), col = 'cyan3', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Anet_1050, aes(x = Sdt_DD  , y = centered_Anet), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_Anet_1680, aes(x = Sdt_DD  , y = centered_Anet), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype="dashed") +
  
  theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  
  ylab(expression('∆ '* italic('A')['net'] * ' (%)')) +
  xlab(expression (italic('Drought Severity Severity')))
Centered_Anet_plot


##############################################################################################################
########################## Analysis for nrubisco #####################################################
baseline_nrubisco_values <- stat_sunflower[stat_sunflower$time == "t1", c("n_trt", "ID_trt", "plant", "nrubisco")]
colnames(baseline_nrubisco_values)[4] <- "baseline_nrubisco"
stat_sunflower <- merge(stat_sunflower, baseline_nrubisco_values, 
                        by = c("n_trt", "ID_trt", "plant"), suffixes = c("", "_baseline"))
names(stat_sunflower)
stat_sunflower$centered_nubisco <- ((stat_sunflower$nrubisco - stat_sunflower$baseline_nrubisco)/stat_sunflower$baseline_nrubisco)*100
hist(stat_sunflower$centered_nubisco)

####################### Lmer for  centerd Nrubisco #############################
stat_sunflower_t234 = subset(stat_sunflower, time!="t1")
hist(stat_sunflower_t234$centered_nubisco) # Normal dist

nrubisco_lmer <- lmer(centered_nubisco~ Cumul_N*Sdt_DD +
                    (1|plant),
                  data = stat_sunflower_t234)
outlierTest(nrubisco_lmer)  ### No outliers 
plot(nrubisco_lmer)
Anova(nrubisco_lmer) ## Only N effect
AIC(nrubisco_lmer) # 1490.795
vif(nrubisco_lmer) # no colinerarity
plot(resid(nrubisco_lmer) ~ fitted(nrubisco_lmer)) ## good
r.squaredGLMM(nrubisco_lmer) ## marginal : 0.11, conditional: 0.13

summary(nrubisco_lmer) # minimal variability in the intercepts and slopes across plants: variability in chi is likely being captured by the fixed effects 
qqnorm(residuals(nrubisco_lmer))
qqline(residuals(nrubisco_lmer))
densityPlot(residuals(nrubisco_lmer))
shapiro.test(residuals(nrubisco_lmer))

cumul_N_levels <- c(350, 560, 1050, 1680)
trends_results_nrubisco <- emtrends(nrubisco_lmer, ~ Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
print(trends_results_nrubisco)
trends_results_df_nrubisco <- as.data.frame(trends_results_nrubisco)

Sdt_DD_range = seq(min(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   max(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_nrubisco <- test(emtrends(nrubisco_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_nrubisco <- test(emtrends(nrubisco_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_nrubisco <- test(emtrends(nrubisco_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_nrubisco <- test(emtrends(nrubisco_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_nrubisco <- summary(emmeans(nrubisco_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_nrubisco <- summary(emmeans(nrubisco_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_nrubisco <- summary(emmeans(nrubisco_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_nrubisco <- summary(emmeans(nrubisco_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_nrubisco_350 <- intercept_350_nrubisco + slope_350_nrubisco * Sdt_DD_range
reg_nrubisco_560 <- intercept_560_nrubisco + slope_560_nrubisco * Sdt_DD_range
reg_nrubisco_1050 <- intercept_1050_nrubisco + slope_1050_nrubisco * Sdt_DD_range
reg_nrubisco_1680 <- intercept_1680_nrubisco + slope_1680_nrubisco * Sdt_DD_range

regline_nrubisco_350 <- data.frame(Sdt_DD = Sdt_DD_range, centered_nrubisco = reg_nrubisco_350, Cumul_N = 350)
regline_nrubisco_560 <- data.frame(Sdt_DD = Sdt_DD_range, centered_nrubisco = reg_nrubisco_560, Cumul_N = 560)
regline_nrubisco_1050 <- data.frame(Sdt_DD = Sdt_DD_range, centered_nrubisco = reg_nrubisco_1050, Cumul_N = 1050)
regline_nrubisco_1680 <- data.frame(Sdt_DD = Sdt_DD_range, centered_nrubisco = reg_nrubisco_1680, Cumul_N = 1680)

slope_results_nrubisco <- emtrends(nrubisco_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_nrubisco <- contrast(slope_results_nrubisco, method = "pairwise")
slope_summary_nrubisco <- summary(slope_results_nrubisco, infer = c(TRUE, TRUE))
slope_summary_nrubisco
cld(emtrends(nrubisco_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))
################ Nrubisco  plot #####################################################################
Centered_Nrubisco_plot <- ggplot(data = stat_sunflower_t234, aes(x = Sdt_DD, y = centered_nubisco)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +

  geom_line(data = regline_nrubisco_350, aes(x = Sdt_DD  , y = centered_nrubisco  ), col = 'cyan2', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_nrubisco_560, aes(x = Sdt_DD  , y = centered_nrubisco  ), col = 'cyan3', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_nrubisco_1050, aes(x = Sdt_DD  , y = centered_nrubisco  ), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_nrubisco_1680, aes(x = Sdt_DD  , y = centered_nrubisco  ), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype ="dashed") +
  
  theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  
  ylab(expression('∆'* italic('ρ')['rubisco'] * ' (%)')) +
  xlab(expression (italic('Drought Severity Index')))
Centered_Nrubisco_plot
###################################################################################################################
##################################################################################################################

########################## Analysis for nbioe #####################################################
baseline_nbioe_values <- stat_sunflower[stat_sunflower$time == "t1", c("n_trt", "ID_trt", "plant", "nbioe")]
colnames(baseline_nbioe_values)[4] <- "baseline_nbioe"
stat_sunflower <- merge(stat_sunflower, baseline_nbioe_values, 
                        by = c("n_trt", "ID_trt", "plant"), suffixes = c("", "_baseline"))
names(stat_sunflower)
stat_sunflower$centered_nbioe <- ((stat_sunflower$nbioe - stat_sunflower$baseline_nbioe)/stat_sunflower$baseline_nbioe)*100
hist(stat_sunflower$centered_nbioe)
####################### Lmer for  centerd Nrubisco #############################
stat_sunflower_t234 = subset(stat_sunflower, time!="t1")
hist(stat_sunflower_t234$centered_nbioe)

nbioe_lmer <- lmer(centered_nbioe~ Cumul_N*Sdt_DD +
                        (1|plant),
                      data = stat_sunflower_t234)
outlierTest(nbioe_lmer)  
plot(nbioe_lmer) ## two outliers to remove 

residuals_data <- data.frame(
  row_index = 1:nrow(stat_sunflower_t234),
  residuals = rstudent(nbioe_lmer)
)
hist(residuals_data$residuals)
threshold <- 2
outliers <- abs(residuals_data$residuals) > threshold
# Filter the dataset to remove outliers
stat_sunflower_nbioe_clean <- stat_sunflower_t234[!outliers, ]
nrow(stat_sunflower_nbioe_clean) # 139 , three removed
hist(stat_sunflower_nbioe_clean$centered_nbioe) # better dist

############### Refit ##########################################

nbioe_lmer <- lmer(centered_nbioe~ Cumul_N*Sdt_DD +
                     (1|plant),
                   data = stat_sunflower_nbioe_clean)
outlierTest(nbioe_lmer)  # P=0.169 perfect!
plot(nbioe_lmer) 
Anova(nbioe_lmer) ## Only N effect
AIC(nbioe_lmer) # 1401.708
vif(nbioe_lmer) # no colinerarity
plot(resid(nbioe_lmer) ~ fitted(nbioe_lmer)) ## good
r.squaredGLMM(nbioe_lmer) ## marginal : 0.04, conditional: 0.08

summary(nbioe_lmer) 
qqline(residuals(nbioe_lmer))
densityPlot(residuals(nbioe_lmer))
shapiro.test(residuals(nbioe_lmer))
residuals <- resid(nbioe_lmer)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350, 560, 1050, 1680)
Sdt_DD_range = seq(min(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   max(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_nbioe  <- test(emtrends(nbioe_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_nbioe  <- test(emtrends(nbioe_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_nbioe  <- test(emtrends(nbioe_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_nbioe  <- test(emtrends(nbioe_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_nbioe  <- summary(emmeans(nbioe_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_nbioe  <- summary(emmeans(nbioe_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_nbioe  <- summary(emmeans(nbioe_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_nbioe  <- summary(emmeans(nbioe_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_nbioe_350 <- intercept_350_nbioe  + slope_350_nbioe  * Sdt_DD_range
reg_nbioe_560 <- intercept_560_nbioe  + slope_560_nbioe  * Sdt_DD_range
reg_nbioe_1050 <- intercept_1050_nbioe  + slope_1050_nbioe  * Sdt_DD_range
reg_nbioe_1680 <- intercept_1680_nbioe  + slope_1680_nbioe  * Sdt_DD_range

regline_nbioe_350 <- data.frame(Sdt_DD = Sdt_DD_range, centered_nbioe = reg_nbioe_350, Cumul_N = 350)
regline_nbioe_560 <- data.frame(Sdt_DD = Sdt_DD_range, centered_nbioe = reg_nbioe_560, Cumul_N = 560)
regline_nbioe_1050 <- data.frame(Sdt_DD = Sdt_DD_range, centered_nbioe = reg_nbioe_1050, Cumul_N = 1050)
regline_nbioe_1680 <- data.frame(Sdt_DD = Sdt_DD_range, centered_nbioe = reg_nbioe_1680, Cumul_N = 1680)

slope_results_nbioe  <- emtrends(nbioe_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_nbioe  <- contrast(slope_results_nbioe , method = "pairwise")
slope_summary_nbioe  <- summary(slope_results_nbioe , infer = c(TRUE, TRUE))
slope_summary_nbioe 
cld(emtrends(nbioe_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))
################ N bioe plot #####################################################################
Centered_Nbioe_plot <- ggplot(data = stat_sunflower_nbioe_clean, aes(x = Sdt_DD, y = centered_nbioe)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_nbioe_350, aes(x = Sdt_DD  , y = centered_nbioe), col = 'cyan2', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_nbioe_560, aes(x = Sdt_DD  , y = centered_nbioe), col = 'cyan3', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_nbioe_1050, aes(x = Sdt_DD  , y = centered_nbioe), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_nbioe_1680, aes(x = Sdt_DD  , y = centered_nbioe), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype ="dashed") +
 
  theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  
  ylab(expression('∆'* italic('ρ')['bioenergetics'] * ' (%)')) +
  xlab(expression (italic('Drought Severity Index')))
Centered_Nbioe_plot

###################################################################################################################
##################################################################################################################

###########################################################################################################################
########################## Analysis for Nstrucutre #####################################################
baseline_Nstrucutre_values <- stat_sunflower[stat_sunflower$time == "t1", c("n_trt", "ID_trt", "plant", "Nstrucutre")]
colnames(baseline_Nstrucutre_values)[4] <- "baseline_Nstrucutre"
stat_sunflower <- merge(stat_sunflower, baseline_Nstrucutre_values, 
                        by = c("n_trt", "ID_trt", "plant"), suffixes = c("", "_baseline"))
names(stat_sunflower)
stat_sunflower$centered_Nstrucutre <- ((stat_sunflower$Nstrucutre - stat_sunflower$baseline_Nstrucutre)/stat_sunflower$baseline_Nstrucutre)*100
hist(stat_sunflower$centered_Nstrucutre)

####################### Lmer for  centered Nstructure #################
stat_sunflower_t234 = subset(stat_sunflower, time!="t1")
hist(stat_sunflower_t234$centered_Nstrucutre) ## Normal dist

Nstrucutre_lmer <- lmer(centered_Nstrucutre~ Cumul_N*Sdt_DD +
                        (1|plant),
                      data = stat_sunflower_t234)
outlierTest(Nstrucutre_lmer)  ########## No outliers 
plot(Nstrucutre_lmer)
Anova(Nstrucutre_lmer) ## Only strong drought effects
AIC(Nstrucutre_lmer) # 1665.755
vif(Nstrucutre_lmer) # no colinerarity
plot(resid(Nstrucutre_lmer) ~ fitted(Nstrucutre_lmer)) ## good
r.squaredGLMM(Nstrucutre_lmer) ## marginal : 0.13, conditional: 0.31

summary(Nstrucutre_lmer) # minimal variability in the intercepts and slopes across plants: variability in chi is likely being captured by the fixed effects 
qqnorm(residuals(Nstrucutre_lmer))
qqline(residuals(Nstrucutre_lmer))
densityPlot(residuals(Nstrucutre_lmer))
shapiro.test(residuals(Nstrucutre_lmer))

cumul_N_levels_nstr <- c(350, 560, 1050, 1680)
Sdt_DD_range = seq(min(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   max(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   0.001)

slope_350_nstr <- test(emtrends(Nstrucutre_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_nstr <- test(emtrends(Nstrucutre_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_nstr <- test(emtrends(Nstrucutre_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_nstr <- test(emtrends(Nstrucutre_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_nstr <- summary(emmeans(Nstrucutre_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_nstr <- summary(emmeans(Nstrucutre_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_nstr <- summary(emmeans(Nstrucutre_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_nstr <- summary(emmeans(Nstrucutre_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_nstr_350 <- intercept_350_nstr + slope_350_nstr * Sdt_DD_range
reg_nstr_560 <- intercept_560_nstr + slope_560_nstr* Sdt_DD_range
reg_nstr_1050 <- intercept_1050_nstr + slope_1050_nstr * Sdt_DD_range
reg_nstr_1680 <- intercept_1680_nstr + slope_1680_nstr * Sdt_DD_range

regline_nstr_350 <- data.frame(Sdt_DD = Sdt_DD_range, centered_nstr= reg_nstr_350, Cumul_N = 350)
regline_nstr_560 <- data.frame(Sdt_DD = Sdt_DD_range, centered_nstr = reg_nstr_560, Cumul_N = 525)
regline_nstr_1050 <- data.frame(Sdt_DD = Sdt_DD_range, centered_nstr = reg_nstr_1050, Cumul_N = 1050)
regline_nstr_1680 <- data.frame(Sdt_DD = Sdt_DD_range, centered_nstr = reg_nstr_1680, Cumul_N = 1680)

slope_results_nstr <- emtrends(Nstrucutre_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_nstr <- contrast(slope_results_nstr, method = "pairwise")
slope_summary_nstr <- summary(slope_results_nstr, infer = c(TRUE, TRUE))
slope_summary_nstr
cld(emtrends(Nstrucutre_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))
################ N structure plot #####################################################################
Centered_Nstructure_plot <- ggplot(data = stat_sunflower_t234, aes(x = Sdt_DD, y = centered_Nstrucutre)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_nstr_350, aes(x = Sdt_DD  , y = centered_nstr), col = 'cyan2', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_nstr_560, aes(x = Sdt_DD  , y = centered_nstr), col = 'cyan3', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_nstr_1050, aes(x = Sdt_DD  , y = centered_nstr), col = 'darkcyan', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_nstr_1680, aes(x = Sdt_DD  , y = centered_nstr), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype="dashed") +
  
  
  theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  
  ylab(expression('∆'* italic('ρ')['structure'] * ' (%)')) +
  xlab(expression (italic('Drought Severity Index')))
Centered_Nstructure_plot

###########################################################################################################################
########################## Analysis for leaf PNUE #####################################################
names(stat_sunflower)
stat_sunflower$PNUE = stat_sunflower$anet_420/stat_sunflower$Narea
hist(stat_sunflower$PNUE)

baseline_PNUE_values <- stat_sunflower[stat_sunflower$time == "t1", c("n_trt", "ID_trt", "plant", "PNUE")]
colnames(baseline_PNUE_values)[4] <- "baseline_PNUE"
stat_sunflower <- merge(stat_sunflower, baseline_PNUE_values, 
                        by = c("n_trt", "ID_trt", "plant"), suffixes = c("", "_baseline"))
names(stat_sunflower)
stat_sunflower$centered_PNUE <- ((stat_sunflower$PNUE - stat_sunflower$baseline_PNUE)/stat_sunflower$baseline_PNUE)*100
hist(stat_sunflower$centered_PNUE)

####################### Lmer for  centerd PNUE #################
stat_sunflower_t234 = subset(stat_sunflower, time!="t1")
hist(stat_sunflower_t234$centered_PNUE) ## Normal dist

PNUE_lmer <- lmer(centered_PNUE~ Cumul_N*Sdt_DD +
                          (1|plant),
                        data = stat_sunflower_t234)
outlierTest(PNUE_lmer)  ########## No outliers 
plot(PNUE_lmer)
Anova(PNUE_lmer) 
AIC(PNUE_lmer) 
vif(PNUE_lmer) 
plot(resid(PNUE_lmer) ~ fitted(PNUE_lmer))
r.squaredGLMM(PNUE_lmer) 

summary(PNUE_lmer) 
qqnorm(residuals(PNUE_lmer))
qqline(residuals(PNUE_lmer))
densityPlot(residuals(PNUE_lmer))
shapiro.test(residuals(PNUE_lmer))

cumul_N_levels_PNUE<- c(350, 560, 1050, 1680)

Sdt_DD_range = seq(min(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   max(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   0.001)

slope_350_PNUE <- test(emtrends(PNUE_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_PNUE <- test(emtrends(PNUE_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_PNUE <- test(emtrends(PNUE_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_PNUE <- test(emtrends(PNUE_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_PNUE <- summary(emmeans(PNUE_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_PNUE <- summary(emmeans(PNUE_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_PNUE <- summary(emmeans(PNUE_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_PNUE <- summary(emmeans(PNUE_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_PNUE_350 <- intercept_350_PNUE + slope_350_PNUE * Sdt_DD_range
reg_PNUE_560 <- intercept_560_PNUE + slope_560_PNUE* Sdt_DD_range
reg_PNUE_1050 <- intercept_1050_PNUE + slope_1050_PNUE * Sdt_DD_range
reg_PNUE_1680 <- intercept_1680_PNUE + slope_1680_PNUE * Sdt_DD_range

regline_PNUE_350 <- data.frame(Sdt_DD = Sdt_DD_range, centered_PNUE= reg_PNUE_350, Cumul_N = 350)
regline_PNUE_560 <- data.frame(Sdt_DD = Sdt_DD_range, centered_PNUE = reg_PNUE_560, Cumul_N = 525)
regline_PNUE_1050 <- data.frame(Sdt_DD = Sdt_DD_range, centered_PNUE = reg_PNUE_1050, Cumul_N = 1050)
regline_PNUE_1680 <- data.frame(Sdt_DD = Sdt_DD_range, centered_PNUE = reg_PNUE_1680, Cumul_N = 1680)

slope_results_PNUE <- emtrends(PNUE_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_PNUE <- contrast(slope_results_PNUE, method = "pairwise")
slope_summary_PNUE <- summary(slope_results_PNUE, infer = c(TRUE, TRUE))
slope_summary_PNUE
cld(emtrends(PNUE_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))
################ N structure plot #####################################################################
Centered_PNUE_plot <- ggplot(data = stat_sunflower_t234, aes(x = Sdt_DD, y = centered_PNUE)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_PNUE_350, aes(x = Sdt_DD  , y = centered_PNUE), col = 'cyan2', lwd = 2, alpha = 0.8, , linetype="dashed") +
  geom_line(data = regline_PNUE_560, aes(x = Sdt_DD  , y = centered_PNUE), col = 'cyan3', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_PNUE_1050, aes(x = Sdt_DD  , y = centered_PNUE), col = 'darkcyan', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_PNUE_1680, aes(x = Sdt_DD  , y = centered_PNUE), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype="dashed") +
  
  
  theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  
  ylab(expression('∆'* italic('NUE')['leaf'] * ' (%)')) +
  xlab(expression (italic('Drought Severity Index')))
Centered_PNUE_plot

########################## Analysis for leaf WUE #####################################################
names(stat_sunflower)
stat_sunflower$LWUE = stat_sunflower$anet_420/stat_sunflower$gs_420
hist(stat_sunflower$LWUE)

baseline_LWUE_values <- stat_sunflower[stat_sunflower$time == "t1", c("n_trt", "ID_trt", "plant", "LWUE")]
colnames(baseline_LWUE_values)[4] <- "baseline_LWUE"
stat_sunflower <- merge(stat_sunflower, baseline_LWUE_values, 
                        by = c("n_trt", "ID_trt", "plant"), suffixes = c("", "_baseline"))
names(stat_sunflower)
stat_sunflower$centered_LWUE <- ((stat_sunflower$LWUE - stat_sunflower$baseline_LWUE)/stat_sunflower$baseline_LWUE)*100
hist(stat_sunflower$centered_LWUE)

####################### Lmer for  centerd PNUE #################
stat_sunflower_t234 = subset(stat_sunflower, time!="t1")
hist(stat_sunflower_t234$centered_LWUE) ## Normal dist
stat_sunflower_t234_LWUE = subset(stat_sunflower_t234, centered_LWUE<200)
nrow(stat_sunflower_t234_LWUE)
LWUE_lmer <- lmer(centered_LWUE~ Cumul_N*Sdt_DD +
                    (1|plant),
                  data = stat_sunflower_t234_LWUE)
outlierTest(LWUE_lmer)  
plot(LWUE_lmer)
Anova(LWUE_lmer) 
AIC(LWUE_lmer)
vif(LWUE_lmer) 
plot(resid(LWUE_lmer) ~ fitted(LWUE_lmer)) 
r.squaredGLMM(LWUE_lmer) 

summary(LWUE_lmer) 
qqnorm(residuals(LWUE_lmer))
qqline(residuals(LWUE_lmer))
densityPlot(residuals(LWUE_lmer))
shapiro.test(residuals(LWUE_lmer))

cumul_N_levels_LWUE<- c(350, 560, 1050, 1680)
Sdt_DD_range = seq(min(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   max(stat_sunflower_t234$Sdt_DD, na.rm = TRUE), 
                   0.001)

slope_350_LWUE <- test(emtrends(LWUE_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_LWUE <- test(emtrends(LWUE_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_LWUE <- test(emtrends(LWUE_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_LWUE <- test(emtrends(LWUE_lmer, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_LWUE <- summary(emmeans(LWUE_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_LWUE <- summary(emmeans(LWUE_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_LWUE <- summary(emmeans(LWUE_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_LWUE <- summary(emmeans(LWUE_lmer, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_LWUE_350 <- intercept_350_LWUE + slope_350_LWUE * Sdt_DD_range
reg_LWUE_560 <- intercept_560_LWUE + slope_560_LWUE* Sdt_DD_range
reg_LWUE_1050 <- intercept_1050_LWUE + slope_1050_LWUE * Sdt_DD_range
reg_LWUE_1680 <- intercept_1680_LWUE + slope_1680_LWUE * Sdt_DD_range

regline_LWUE_350 <- data.frame(Sdt_DD = Sdt_DD_range, centered_LWUE= reg_LWUE_350, Cumul_N = 350)
regline_LWUE_560 <- data.frame(Sdt_DD = Sdt_DD_range, centered_LWUE = reg_LWUE_560, Cumul_N = 525)
regline_LWUE_1050 <- data.frame(Sdt_DD = Sdt_DD_range, centered_LWUE = reg_LWUE_1050, Cumul_N = 1050)
regline_LWUE_1680 <- data.frame(Sdt_DD = Sdt_DD_range, centered_LWUE = reg_LWUE_1680, Cumul_N = 1680)

slope_results_LWUE <- emtrends(LWUE_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_LWUE <- contrast(slope_results_LWUE, method = "pairwise")
slope_summary_LWUE <- summary(slope_results_LWUE, infer = c(TRUE, TRUE))
slope_summary_LWUE
cld(emtrends(LWUE_lmer, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))
################ N structure plot #####################################################################
Centered_LWUE_plot <- ggplot(data = stat_sunflower_t234_LWUE, aes(x = Sdt_DD, y = centered_LWUE)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_LWUE_350, aes(x = Sdt_DD  , y = centered_LWUE), col = 'cyan2', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_LWUE_560, aes(x = Sdt_DD  , y = centered_LWUE), col = 'cyan3', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_LWUE_1050, aes(x = Sdt_DD  , y = centered_LWUE), col = 'darkcyan', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_LWUE_1680, aes(x = Sdt_DD  , y = centered_LWUE), col = 'lightcyan4', lwd = 2, alpha = 0.8) +
  
  
  theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  
  ylab(expression('∆'* italic('WUE')['leaf'] * ' (%)')) +
  xlab(expression (italic('Drought Severity Index')))
Centered_LWUE_plot
###########################################################################################################################

############################### Export Table Satistics ###########################################

################ gs ##########################
anova_results_gs <- as.data.frame(Anova(gs_lmer))
slope_summary_gs <- as.data.frame(slope_summary_gs)
anova_results_gs$model <- "gs"
slope_summary_gs$model <- "gs"
################ Chi ##########################
anova_results_chi <- as.data.frame(Anova(chi_lmer))
slope_summary_chi <- as.data.frame(slope_summary_chi)
anova_results_chi$model <- "Chi"
slope_summary_chi$model <- "Chi"
################ betaleaf ##########################
anova_results_betaleaf <- as.data.frame(Anova(Betaleaf_lmer))
slope_summary_betaleaf <- as.data.frame(slope_summary_Betaleaf)
anova_results_betaleaf$model <- "betaleaf"
slope_summary_betaleaf$model <- "betaleaf"
################ Nmass  ##########################
anova_results_Nmass <- as.data.frame(Anova(Nmass_lmer))
slope_summary_Nmass <- as.data.frame(slope_summary_Nmass)
anova_results_Nmass$model <- "Nmass"
slope_summary_Nmass$model <- "Nmass"
################ Narea ##########################slope_summary_Narea_df <- as.data.frame(slope_summary_Narea)
anova_results_Narea <- as.data.frame(Anova(Narea_lmer))
slope_summary_Narea <- as.data.frame(slope_summary_Narea)
anova_results_Narea$model <- "Narea"
slope_summary_Narea$model <- "Narea"
################ LMA ##########################
anova_results_LMA <- as.data.frame(Anova(LMA_lmer))
slope_summary_LMA <- as.data.frame(slope_summary_LMA)
anova_results_LMA$model <- "LMA"
slope_summary_LMA$model <- "LMA"
################ Vcmax ##########################
anova_results_Vcmax <- as.data.frame(Anova(Vcmax_lmer))
slope_summary_Vcmax <- as.data.frame(slope_summary_Vcmax)
anova_results_Vcmax$model <- "Vcmax"
slope_summary_Vcmax$model <- "Vcmax"
################ Jmax ##########################
anova_results_Jmax <- as.data.frame(Anova(Jmax_lmer))
slope_summary_jmax <- as.data.frame(slope_summary_jmax)
anova_results_Jmax$model <- "Jmax"
slope_summary_jmax$model <- "Jmax"
################ Nrubisco ##########################
anova_results_nrubisco <- as.data.frame(Anova(nrubisco_lmer))
slope_summary_nrubisco <- as.data.frame(slope_summary_nrubisco)
anova_results_nrubisco$model <- "nrubisco"
slope_summary_nrubisco$model <- "nrubisco"
################# nbioe ################
anova_results_nbioe <- as.data.frame(Anova(nbioe_lmer))
slope_summary_nbioe <- as.data.frame(slope_summary_nbioe)
anova_results_nbioe$model <- "nbioe"
slope_summary_nbioe$model <- "nbioe"
############# Nstructure ###################
anova_results_nstr <- as.data.frame(Anova(Nstrucutre_lmer))
slope_summary_nstr <- as.data.frame(slope_summary_nstr)
anova_results_nstr$model <- "nstr"
slope_summary_nstr$model <- "nstr"

##############Anet###############
anova_results_Anet <- as.data.frame(Anova(Anet_lmer))
anova_results_Anet$model <- "Anet"
slope_summary_anet <- as.data.frame(slope_summary_anet)
slope_summary_anet$model <- "Anet"
##############PNUE###############
anova_results_PNUE <- as.data.frame(Anova(PNUE_lmer))
anova_results_PNUE$model <- "PNUE"
slope_summary_PNUE <- as.data.frame(slope_summary_PNUE)
slope_summary_PNUE$model <- "PNUE"
##############LWUE###############
anova_results_LWUE <- as.data.frame(Anova(LWUE_lmer))
anova_results_LWUE$model <- "WUE"
slope_summary_LWUE <- as.data.frame(slope_summary_LWUE)
slope_summary_LWUE$model <- "WUE"

################# Combine all anova tables ############################
combined_anova_df <- rbind(anova_results_gs, anova_results_chi,anova_results_betaleaf, anova_results_Nmass, anova_results_Narea,
                           anova_results_LMA,anova_results_Vcmax,anova_results_Jmax,anova_results_nrubisco,
                           anova_results_nbioe, anova_results_nstr,
                           anova_results_LWUE,anova_results_PNUE,
                           anova_results_Anet)  

combined_anova_df   

write.csv(combined_anova_df, "../output/combined_anova_df.csv", row.names = TRUE)



combined_slopes <- rbind(slope_summary_gs, slope_summary_chi,slope_summary_betaleaf, slope_summary_Nmass, slope_summary_Narea,
                         slope_summary_LMA,slope_summary_Vcmax,slope_summary_jmax,slope_summary_nrubisco,
                         slope_summary_nbioe, slope_summary_nstr,
                         slope_summary_LWUE,slope_summary_PNUE,
                          slope_summary_anet)  

combined_slopes 

write.csv(combined_slopes, "../output/combined_slopes.csv", row.names = TRUE)
#####################################################################################




