#########################################
# packages necessary
#########################################
library(plotly)
library(lme4)
library(car)
library(emmeans)
library(ggeffects)
library(RColorBrewer)
library(multcompView)
library(nlme)
library(marginaleffects)
library(piecewiseSEM)
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

############## load data #################################################################### 
stat_all <- read.csv("../output/sunflower.csv")
names(stat_all)

stat_t4 = subset(stat_all, time =="t4")
table(stat_t4$time)
table(stat_t4$Cumul_N)

#############################  Analysis  #############################################
names(stat_t4)
######################## chi ######################################################################################### 

names(stat_t4)
gs_lmer_t4 <- lm(gs_420~ Cumul_N*Sdt_DD,
                    data = stat_t4)
Anova(gs_lmer_t4, type = "II")
shapiro.test(residuals(gs_lmer_t4)) 
outlierTest(gs_lmer_t4)  

Anova(gs_lmer_t4) # DD effects 
AIC(gs_lmer_t4) # 19.29
vif(gs_lmer_t4) # no colinerarity
plot(resid(gs_lmer_t4) ~ fitted(gs_lmer_t4)) ## good
r.squaredGLMM(gs_lmer_t4) ## marginal : 0.17, conditional: 0.17
summary(gs_lmer_t4) 
qqnorm(residuals(gs_lmer_t4))
qqline(residuals(gs_lmer_t4))
densityPlot(residuals(gs_lmer_t4))

residuals <- resid(gs_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
Sdt_DD_range_t4 = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                      max(stat_t4$Sdt_DD, na.rm = TRUE), 
                      0.001)
slope_350_gs_t4 <- test(emtrends(gs_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_gs_t4  <- test(emtrends(gs_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_gs_t4  <- test(emtrends(gs_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_gs_t4  <- test(emtrends(gs_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_gs_t4  <- summary(emmeans(gs_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_gs_t4  <- summary(emmeans(gs_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_gs_t4  <- summary(emmeans(gs_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_gs_t4 <- summary(emmeans(gs_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_gs_350_t4 <- intercept_350_gs_t4  + slope_350_gs_t4  * Sdt_DD_range_t4
reg_gs_560_t4 <- intercept_560_gs_t4  + slope_560_gs_t4  * Sdt_DD_range_t4
reg_gs_1050_t4 <- intercept_1050_gs_t4  + slope_1050_gs_t4  * Sdt_DD_range_t4
reg_gs_1680_t4 <- intercept_1680_gs_t4  + slope_1680_gs_t4  * Sdt_DD_range_t4

regline_gs_350_t4 <- data.frame(Sdt_DD = Sdt_DD_range_t4, gs = reg_gs_350_t4, Cumul_N = 350)
regline_gs_560_t4 <- data.frame(Sdt_DD = Sdt_DD_range_t4, gs = reg_gs_560_t4, Cumul_N = 560)
regline_gs_1050_t4 <- data.frame(Sdt_DD = Sdt_DD_range_t4, gs = reg_gs_1050_t4, Cumul_N = 1050)
regline_gs_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range_t4, gs = reg_gs_1680_t4, Cumul_N = 1680)

slope_results_gs_t4 <- emtrends(gs_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_gs_t4 <- contrast(slope_results_gs_t4, method = "pairwise")
slope_summary_gs_t4 <- summary(slope_results_gs_t4, infer = c(TRUE, TRUE))
slope_summary_gs_t4

slope_test_gs_t4 <- cld(emtrends(gs_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

###################################### Plot #######################
gs_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = gs_420)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_gs_350_t4, aes(x = Sdt_DD_range_t4, y = gs), col = 'cyan2', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_gs_560_t4, aes(x = Sdt_DD_range_t4, y = gs), col = 'cyan3', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_gs_1050_t4, aes(x = Sdt_DD_range_t4 , y = gs), col = 'darkcyan', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_gs_1680_t4, aes(x = Sdt_DD_range_t4 , y = gs), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype ="dashed") +
  
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
  
  ylab(expression(g[s] ~ (µmol ~ CO[2] ~ m^{-2} ~ s^{-1}))) +
  xlab(expression (italic('Drought Severity Index')))
gs_plot_t4 

######################## chi ######################################################################################### 
hist(stat_t4$chi)
nrow(stat_t4)
chi_lmer_t4 <- lm(chi~ Cumul_N*Sdt_DD, 
                  data = stat_t4)

shapiro.test(residuals(chi_lmer_t4)) # 0.44
outlierTest(chi_lmer_t4)  ### No outliers

Anova(chi_lmer_t4, type = "II") ## N and DD effects 
AIC(chi_lmer_t4) # - 226.248
vif(chi_lmer_t4) # no colinerarity
plot(resid(chi_lmer_t4) ~ fitted(chi_lmer_t4)) ## good
r.squaredGLMM(chi_lmer_t4) ## marginal : 0.28, conditional: 0.28
summary(chi_lmer_t4) # minimal variability in the intercepts and slopes across plants: variability in chi is likely being captured by the fixed effects 
qqnorm(residuals(chi_lmer_t4))
qqline(residuals(chi_lmer_t4))
densityPlot(residuals(chi_lmer_t4))

residuals <- resid(chi_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)

Sdt_DD_range_t4 = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_chi_t4 <- test(emtrends(chi_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_chi_t4  <- test(emtrends(chi_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_chi_t4  <- test(emtrends(chi_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_chi_t4  <- test(emtrends(chi_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_chi_t4  <- summary(emmeans(chi_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_chi_t4  <- summary(emmeans(chi_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_chi_t4  <- summary(emmeans(chi_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_chi_t4 <- summary(emmeans(chi_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_chi_350_t4 <- intercept_350_chi_t4  + slope_350_chi_t4  * Sdt_DD_range_t4
reg_chi_560_t4 <- intercept_560_chi_t4  + slope_560_chi_t4  * Sdt_DD_range_t4
reg_chi_1050_t4 <- intercept_1050_chi_t4  + slope_1050_chi_t4  * Sdt_DD_range_t4
reg_chi_1680_t4 <- intercept_1680_chi_t4  + slope_1680_chi_t4  * Sdt_DD_range_t4

regline_chi_350_t4 <- data.frame(Sdt_DD = Sdt_DD_range_t4, chi = reg_chi_350_t4, Cumul_N = 350)
regline_chi_560_t4 <- data.frame(Sdt_DD = Sdt_DD_range_t4, chi = reg_chi_560_t4, Cumul_N = 560)
regline_chi_1050_t4 <- data.frame(Sdt_DD = Sdt_DD_range_t4, chi = reg_chi_1050_t4, Cumul_N = 1050)
regline_chi_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range_t4, chi = reg_chi_1680_t4, Cumul_N = 1680)

slope_results_chi_t4 <- emtrends(chi_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_chi_t4 <- contrast(slope_results_chi_t4, method = "pairwise")
slope_summary_chi_t4 <- summary(slope_results_chi_t4, infer = c(TRUE, TRUE))
slope_summary_chi_t4

slope_test_chi_t4 <- cld(emtrends(chi_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

###################################### Plot #######################
Chi_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = chi)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_chi_350_t4, aes(x = Sdt_DD_range_t4, y = chi), col = 'cyan2', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_chi_560_t4, aes(x = Sdt_DD_range_t4, y = chi), col = 'cyan3', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_chi_1050_t4, aes(x = Sdt_DD_range_t4 , y = chi), col = 'darkcyan', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_chi_1680_t4, aes(x = Sdt_DD_range_t4 , y = chi), col = 'lightcyan4', lwd = 2, alpha = 0.8) +
  
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
  
  ylab(expression(chi * " (Pa Pa"^{-1} * ")"))+
  xlab(expression (italic('Drought Severity Index')))
Chi_plot_t4 

######################## Beta leaf ######################################################################################### 
hist(stat_t4$betaleaf)
hist(log(stat_t4$betaleaf))
nrow(stat_t4)
betaleaf_lmer_t4 <- lm(log(betaleaf)~ Cumul_N*Sdt_DD,
                    data = stat_t4)
Anova(betaleaf_lmer_t4, type = "II") ## drought effects 

plot(betaleaf_lmer_t4)
shapiro.test(residuals(betaleaf_lmer_t4)) # 0.64
outlierTest(betaleaf_lmer_t4)  ### No outliers
Anova(betaleaf_lmer_t4) ## N and DD effects 
AIC(betaleaf_lmer_t4) #14.99
vif(betaleaf_lmer_t4) # no colinerarity
plot(resid(betaleaf_lmer_t4) ~ fitted(betaleaf_lmer_t4)) ## good
r.squaredGLMM(betaleaf_lmer_t4) ## marginal : 0.27, conditional: 0.27
summary(betaleaf_lmer_t4) 
qqline(residuals(betaleaf_lmer_t4))
densityPlot(residuals(betaleaf_lmer_t4))

residuals <- resid(betaleaf_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)

Sdt_DD_range_t4 = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                      max(stat_t4$Sdt_DD, na.rm = TRUE), 
                      0.001)
slope_350_betaleaf_t4 <- test(emtrends(betaleaf_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_betaleaf_t4  <- test(emtrends(betaleaf_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_betaleaf_t4  <- test(emtrends(betaleaf_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_betaleaf_t4  <- test(emtrends(betaleaf_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_betaleaf_t4  <- summary(emmeans(betaleaf_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_betaleaf_t4  <- summary(emmeans(betaleaf_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_betaleaf_t4  <- summary(emmeans(betaleaf_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_betaleaf_t4 <- summary(emmeans(betaleaf_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_betaleaf_350_t4 <- exp(intercept_350_betaleaf_t4  + slope_350_betaleaf_t4  * Sdt_DD_range_t4)
reg_betaleaf_560_t4 <- exp(intercept_560_betaleaf_t4  + slope_560_betaleaf_t4  * Sdt_DD_range_t4)
reg_betaleaf_1050_t4 <- exp(intercept_1050_betaleaf_t4  + slope_1050_betaleaf_t4  * Sdt_DD_range_t4)
reg_betaleaf_1680_t4 <- exp(intercept_1680_betaleaf_t4  + slope_1680_betaleaf_t4  * Sdt_DD_range_t4)

regline_betaleaf_350_t4 <- data.frame(Sdt_DD = Sdt_DD_range_t4, betaleaf = reg_betaleaf_350_t4, Cumul_N = 350)
regline_betaleaf_560_t4 <- data.frame(Sdt_DD = Sdt_DD_range_t4, betaleaf = reg_betaleaf_560_t4, Cumul_N = 560)
regline_betaleaf_1050_t4 <- data.frame(Sdt_DD = Sdt_DD_range_t4, betaleaf = reg_betaleaf_1050_t4, Cumul_N = 1050)
regline_betaleaf_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range_t4, betaleaf = reg_betaleaf_1680_t4, Cumul_N = 1680)

slope_results_betaleaf_t4 <- emtrends(betaleaf_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_betaleaf_t4 <- contrast(slope_results_betaleaf_t4, method = "pairwise")
slope_summary_betaleaf_t4 <- summary(slope_results_betaleaf_t4, infer = c(TRUE, TRUE))
slope_summary_betaleaf_t4

slope_test_betaleaf_t4 <- cld(emtrends(betaleaf_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

###################################### Plot #######################
betaleaf_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = betaleaf)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_betaleaf_350_t4, aes(x = Sdt_DD_range_t4, y = betaleaf), col = 'cyan2', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_betaleaf_560_t4, aes(x = Sdt_DD_range_t4, y = betaleaf), col = 'cyan3', lwd = 2, alpha = 0.8, linetype ="dashed") +
  geom_line(data = regline_betaleaf_1050_t4, aes(x = Sdt_DD_range_t4 , y = betaleaf), col = 'darkcyan', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_betaleaf_1680_t4, aes(x = Sdt_DD_range_t4 , y = betaleaf), col = 'lightcyan4', lwd = 2, alpha = 0.8) +
  
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
  
  ylab(expression(italic('β') ['leaf'])) + 
  xlab(expression (italic('Drought Severity Index')))
betaleaf_plot_t4 

####################### Nmass #####################################################################################
hist(stat_t4$Nmass)
hist(log(stat_t4$Nmass))
Nmass_lmer_t4 <- lm(log(Nmass)~ Cumul_N*Sdt_DD,
                    data = stat_t4)
plot(Nmass_lmer_t4)
shapiro.test(residuals(Nmass_lmer_t4)) # 0.529
outlierTest(Nmass_lmer_t4)  ### No outliers
Anova(Nmass_lmer_t4, type="II") # N and DD effects
AIC(Nmass_lmer_t4) # 6.67
vif(Nmass_lmer_t4) # no colinerarity
plot(resid(Nmass_lmer_t4) ~ fitted(Nmass_lmer_t4)) ## good
r.squaredGLMM(Nmass_lmer_t4) ## marginal : 0.48, conditional: 0.48
summary(Nmass_lmer_t4) # minimal variability in the intercepts and slopes across plants: variability in chi is likely being captured by the fixed effects 
qqnorm(residuals(Nmass_lmer_t4))
qqline(residuals(Nmass_lmer_t4))
densityPlot(residuals(Nmass_lmer_t4))

residuals <- resid(Nmass_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
Sdt_DD_range_t4 = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_Nmass_t4 <- test(emtrends(Nmass_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_Nmass_t4 <- test(emtrends(Nmass_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_Nmass_t4 <- test(emtrends(Nmass_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_Nmass_t4 <- test(emtrends(Nmass_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_Nmass_t4 <- summary(emmeans(Nmass_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_Nmass_t4 <- summary(emmeans(Nmass_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_Nmass_t4 <- summary(emmeans(Nmass_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_Nmass_t4<- summary(emmeans(Nmass_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_Nmass_350_t4 <- exp(intercept_350_Nmass_t4 + slope_350_Nmass_t4 * Sdt_DD_range_t4)
reg_Nmass_560_t4  <- exp(intercept_560_Nmass_t4 + slope_560_Nmass_t4 * Sdt_DD_range_t4)
reg_Nmass_1050_t4  <- exp(intercept_1050_Nmass_t4 + slope_1050_Nmass_t4 * Sdt_DD_range_t4)
reg_Nmass_1680_t4  <- exp(intercept_1680_Nmass_t4 + slope_1680_Nmass_t4 * Sdt_DD_range_t4)

regline_Nmass_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range_t4, Nmass = reg_Nmass_350_t4 , Cumul_N = 350)
regline_Nmass_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range_t4, Nmass = reg_Nmass_560_t4 , Cumul_N = 560)
regline_Nmass_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range_t4, Nmass = reg_Nmass_1050_t4 , Cumul_N = 1050)
regline_Nmass_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range_t4, Nmass = reg_Nmass_1680_t4 , Cumul_N = 1680)

slope_results_Nmass_t4 <- emtrends(Nmass_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_Nmass_t4 <- contrast(slope_results_Nmass_t4, method = "pairwise")
slope_summary_Nmass_t4 <- summary(slope_results_Nmass_t4, infer = c(TRUE, TRUE))
slope_summary_Nmass_t4

slope_test_Nmass_t4 <- cld(emtrends(Nmass_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

############################ plot ##############################################################
Nmass_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = Nmass)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_Nmass_350_t4 , aes(x = Sdt_DD  , y = Nmass), col = 'cyan2', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Nmass_560_t4 , aes(x = Sdt_DD  , y = Nmass), col = 'cyan3', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Nmass_1050_t4 , aes(x = Sdt_DD  , y = Nmass), col = 'darkcyan', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Nmass_1680_t4 , aes(x = Sdt_DD  , y = Nmass), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype ="dashed") +
  
  
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
  
  ylab(expression(N[mass] ~ (gN ~ g^{-1})))+
  xlab(expression (italic('Drought Severity Index')))
Nmass_plot_t4
#######################################################################################################################

####################### Narea #####################################################################################
hist(stat_t4$Narea)
hist(log(stat_t4$Narea))
Narea_lmer_t4 <- lm(log(Narea)~ Cumul_N*Sdt_DD,
                      data = stat_t4)
plot(Narea_lmer_t4)
shapiro.test(residuals(Narea_lmer_t4)) # 0.29
outlierTest(Narea_lmer_t4)  ### No outliers
Anova(Narea_lmer_t4) ## N effects only
AIC(Narea_lmer_t4) # -6.61
vif(Narea_lmer_t4) # no colinerarity
plot(resid(Narea_lmer_t4) ~ fitted(Narea_lmer_t4)) ## good
r.squaredGLMM(Narea_lmer_t4) ## marginal : 0.41, conditional: 0.41
summary(Narea_lmer_t4) # minimal variability in the intercepts and slopes across plants: variability in chi is likely being captured by the fixed effects 
qqnorm(residuals(Narea_lmer_t4))
qqline(residuals(Narea_lmer_t4))
densityPlot(residuals(Narea_lmer_t4))

residuals <- resid(Narea_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_Narea_t4 <- test(emtrends(Narea_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_Narea_t4 <- test(emtrends(Narea_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_Narea_t4 <- test(emtrends(Narea_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_Narea_t4 <- test(emtrends(Narea_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_Narea_t4 <- summary(emmeans(Narea_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_Narea_t4 <- summary(emmeans(Narea_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_Narea_t4 <- summary(emmeans(Narea_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_Narea_t4<- summary(emmeans(Narea_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_Narea_350_t4 <- exp(intercept_350_Narea_t4 + slope_350_Narea_t4 * Sdt_DD_range)
reg_Narea_560_t4  <- exp(intercept_560_Narea_t4 + slope_560_Narea_t4 * Sdt_DD_range)
reg_Narea_1050_t4  <- exp(intercept_1050_Narea_t4 + slope_1050_Narea_t4 * Sdt_DD_range)
reg_Narea_1680_t4  <- exp(intercept_1680_Narea_t4 + slope_1680_Narea_t4 * Sdt_DD_range)

regline_Narea_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Narea = reg_Narea_350_t4 , Cumul_N = 350)
regline_Narea_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Narea = reg_Narea_560_t4 , Cumul_N = 560)
regline_Narea_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Narea = reg_Narea_1050_t4 , Cumul_N = 1050)
regline_Narea_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, Narea = reg_Narea_1680_t4 , Cumul_N = 1680)

slope_results_Narea_t4 <- emtrends(Narea_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_Narea_t4  <- contrast(slope_results_Narea_t4 , method = "pairwise")
slope_summary_Narea_t4  <- summary(slope_results_Narea_t4, infer = c(TRUE, TRUE))
slope_summary_Narea_t4

slope_test_Narea_t4 <- cld(emtrends(Narea_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

############################ plot ################################
Narea_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = Narea)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_Narea_350_t4 , aes(x = Sdt_DD  , y = Narea), col = 'cyan2', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Narea_560_t4 , aes(x = Sdt_DD  , y = Narea), col = 'cyan3', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Narea_1050_t4 , aes(x = Sdt_DD  , y = Narea), col = 'darkcyan', lwd = 2, alpha = 0.8,linetype ="dashed" ) +
  geom_line(data = regline_Narea_1680_t4 , aes(x = Sdt_DD  , y = Narea), col = 'lightcyan4', lwd = 2, alpha = 0.8,linetype ="dashed") +
  
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
  
  ylab(expression(N[area] ~ (gN ~ m^{-2}))) +
  xlab(expression (italic('Drought Severity Index')))
Narea_plot_t4
#######################################################################################################################

####################### LMA #####################################################################################
hist(stat_t4$LMA)
hist(log(stat_t4$LMA))
LMA_lmer_t4 <- lm(log(LMA)~ Cumul_N*Sdt_DD,
                      data = stat_t4)
plot(LMA_lmer_t4)
shapiro.test(residuals(LMA_lmer_t4)) # 0.65
outlierTest(LMA_lmer_t4)  ### No outliers
Anova(LMA_lmer_t4, type = "II") # Only drought effects
AIC(LMA_lmer_t4) # - 64.59
vif(LMA_lmer_t4) # no colinerarity
plot(resid(LMA_lmer_t4) ~ fitted(LMA_lmer_t4)) ## good
r.squaredGLMM(LMA_lmer_t4) ## marginal : 0.28, conditional: 0.28
summary(LMA_lmer_t4) 
qqnorm(residuals(LMA_lmer_t4))
qqline(residuals(LMA_lmer_t4))
densityPlot(residuals(LMA_lmer_t4))

residuals <- resid(LMA_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_LMA_t4 <- test(emtrends(LMA_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_LMA_t4 <- test(emtrends(LMA_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_LMA_t4 <- test(emtrends(LMA_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_LMA_t4 <- test(emtrends(LMA_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_LMA_t4 <- summary(emmeans(LMA_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_LMA_t4 <- summary(emmeans(LMA_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_LMA_t4 <- summary(emmeans(LMA_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_LMA_t4<- summary(emmeans(LMA_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_LMA_350_t4 <- exp(intercept_350_LMA_t4 + slope_350_LMA_t4 * Sdt_DD_range)
reg_LMA_560_t4  <- exp(intercept_560_LMA_t4 + slope_560_LMA_t4 * Sdt_DD_range)
reg_LMA_1050_t4  <- exp(intercept_1050_LMA_t4 + slope_1050_LMA_t4 * Sdt_DD_range)
reg_LMA_1680_t4  <- exp(intercept_1680_LMA_t4 + slope_1680_LMA_t4 * Sdt_DD_range)

regline_LMA_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, LMA = reg_LMA_350_t4 , Cumul_N = 350)
regline_LMA_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, LMA = reg_LMA_560_t4 , Cumul_N = 560)
regline_LMA_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, LMA = reg_LMA_1050_t4 , Cumul_N = 1050)
regline_LMA_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, LMA = reg_LMA_1680_t4 , Cumul_N = 1680)

slope_results_LMA_t4 <- emtrends(LMA_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_LMA_t4  <- contrast(slope_results_LMA_t4 , method = "pairwise")
slope_summary_LMA_t4  <- summary(slope_results_LMA_t4 , infer = c(TRUE, TRUE))
slope_summary_LMA_t4 

slope_test_LMA_t4 <- cld(emtrends(LMA_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

######################### Plot ################################
LMA_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = LMA)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_LMA_350_t4 , aes(x = Sdt_DD  , y = LMA), col = 'cyan2', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_LMA_560_t4 , aes(x = Sdt_DD  , y = LMA), col = 'cyan3', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_LMA_1050_t4 , aes(x = Sdt_DD  , y = LMA), col = 'darkcyan', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_LMA_1680_t4 , aes(x = Sdt_DD  , y = LMA), col = 'lightcyan4', lwd = 2, alpha = 0.8) +
  
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
  
  ylab(expression(LMA ~ (g ~ m^{-2}))) +
  xlab(expression (italic('Drought Severity Index')))
LMA_plot_t4
#######################################################################################################################
####################### Vcmax #####################################################################################
hist(log(stat_t4$vcmax_tleaf))
Vcmax_lmer_t4 <- lm(log(vcmax_tleaf)~ Cumul_N*Sdt_DD,
                    data = stat_t4)
shapiro.test(residuals(Vcmax_lmer_t4)) # 0.92
outlierTest(Vcmax_lmer_t4)  ### No outliers

Anova(Vcmax_lmer_t4, type ="II") # N effect only
AIC(Vcmax_lmer_t4) #4.24
vif(Vcmax_lmer_t4) # no colinerarity
plot(resid(Vcmax_lmer_t4) ~ fitted(Vcmax_lmer_t4)) ## good
r.squaredGLMM(Vcmax_lmer_t4) ## marginal : 0.39, conditional: 0.39
summary(Vcmax_lmer_t4) 
qqnorm(residuals(Vcmax_lmer_t4))
qqline(residuals(Vcmax_lmer_t4))
densityPlot(residuals(Vcmax_lmer_t4))

residuals <- resid(Vcmax_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_Vcmax_t4 <- test(emtrends(Vcmax_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_Vcmax_t4 <- test(emtrends(Vcmax_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_Vcmax_t4 <- test(emtrends(Vcmax_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_Vcmax_t4 <- test(emtrends(Vcmax_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_Vcmax_t4 <- summary(emmeans(Vcmax_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_Vcmax_t4 <- summary(emmeans(Vcmax_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_Vcmax_t4 <- summary(emmeans(Vcmax_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_Vcmax_t4<- summary(emmeans(Vcmax_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_Vcmax_350_t4 <- exp(intercept_350_Vcmax_t4 + slope_350_Vcmax_t4 * Sdt_DD_range)
reg_Vcmax_560_t4  <- exp(intercept_560_Vcmax_t4 + slope_560_Vcmax_t4 * Sdt_DD_range)
reg_Vcmax_1050_t4  <- exp(intercept_1050_Vcmax_t4 + slope_1050_Vcmax_t4 * Sdt_DD_range)
reg_Vcmax_1680_t4  <- exp(intercept_1680_Vcmax_t4 + slope_1680_Vcmax_t4 * Sdt_DD_range)

regline_Vcmax_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Vcmax = reg_Vcmax_350_t4 , Cumul_N = 350)
regline_Vcmax_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Vcmax = reg_Vcmax_560_t4 , Cumul_N = 560)
regline_Vcmax_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Vcmax = reg_Vcmax_1050_t4 , Cumul_N = 1050)
regline_Vcmax_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, Vcmax = reg_Vcmax_1680_t4 , Cumul_N = 1680)

slope_results_Vcmax_t4 <- emtrends(Vcmax_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_Vcmax_t4 <- contrast(slope_results_Vcmax_t4, method = "pairwise")
slope_summary_Vcmax_t4 <- summary(slope_results_Vcmax_t4, infer = c(TRUE, TRUE))
slope_summary_Vcmax_t4

slope_test_Vcmax_t4 <- cld(emtrends(Vcmax_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

###################### Plot##############################
Vcmax_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = vcmax_tleaf)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_Vcmax_350_t4 , aes(x = Sdt_DD  , y = Vcmax), col = 'cyan2', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_Vcmax_560_t4 , aes(x = Sdt_DD  , y = Vcmax), col = 'cyan3', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_Vcmax_1050_t4 , aes(x = Sdt_DD  , y = Vcmax), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_Vcmax_1680_t4 , aes(x = Sdt_DD  , y = Vcmax), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype="dashed") +
  
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
  
  ylab(expression(V[cmax] ~ (µmol ~ CO[2] ~ m^{-2} ~ s^{-1}))) +
  xlab(expression (italic('Drought Severity Index')))
Vcmax_plot_t4
#######################################################################################################################

####################### Jmax #####################################################################################
names(stat_t4)
hist(log(stat_t4$jmax_tleaf))
Jmax_lmer_t4 <- lm(log(jmax_tleaf)~ Cumul_N*Sdt_DD,
                      data = stat_t4)
plot(Jmax_lmer_t4) # no outliers
shapiro.test(residuals(Jmax_lmer_t4)) 
outlierTest(Jmax_lmer_t4)  ### No outliers

Anova(Jmax_lmer_t4, type = "II") # N effects only
AIC(Jmax_lmer_t4) #9.85
vif(Jmax_lmer_t4) # no colinerarity
plot(resid(Jmax_lmer_t4) ~ fitted(Jmax_lmer_t4)) ## good
r.squaredGLMM(Jmax_lmer_t4) ## marginal : 0.40, conditional: 0.40
summary(Jmax_lmer_t4) # minimal variability in the intercepts and slopes across plants: variability in chi is likely being captured by the fixed effects 
qqnorm(residuals(Jmax_lmer_t4))
qqline(residuals(Jmax_lmer_t4))
densityPlot(residuals(Jmax_lmer_t4))

residuals <- resid(Jmax_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)

Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_Jmax_t4 <- test(emtrends(Jmax_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_Jmax_t4 <- test(emtrends(Jmax_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_Jmax_t4 <- test(emtrends(Jmax_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_Jmax_t4 <- test(emtrends(Jmax_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_Jmax_t4 <- summary(emmeans(Jmax_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_Jmax_t4 <- summary(emmeans(Jmax_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_Jmax_t4 <- summary(emmeans(Jmax_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_Jmax_t4<- summary(emmeans(Jmax_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_Jmax_350_t4 <- exp(intercept_350_Jmax_t4 + slope_350_Jmax_t4 * Sdt_DD_range)
reg_Jmax_560_t4  <- exp(intercept_560_Jmax_t4 + slope_560_Jmax_t4 * Sdt_DD_range)
reg_Jmax_1050_t4  <- exp(intercept_1050_Jmax_t4 + slope_1050_Jmax_t4 * Sdt_DD_range)
reg_Jmax_1680_t4  <- exp(intercept_1680_Jmax_t4 + slope_1680_Jmax_t4 * Sdt_DD_range)

regline_Jmax_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Jmax = reg_Jmax_350_t4 , Cumul_N = 350)
regline_Jmax_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Jmax = reg_Jmax_560_t4 , Cumul_N = 560)
regline_Jmax_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Jmax = reg_Jmax_1050_t4 , Cumul_N = 1050)
regline_Jmax_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, Jmax = reg_Jmax_1680_t4 , Cumul_N = 1680)

slope_results_Jmax_t4 <- emtrends(Jmax_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_Jmax_t4 <- contrast(slope_results_Jmax_t4, method = "pairwise")
slope_summary_Jmax_t4 <- summary(slope_results_Jmax_t4, infer = c(TRUE, TRUE))
slope_summary_Jmax_t4

slope_test_Jmax_t4 <- cld(emtrends(Jmax_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

Jmax_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = jmax_tleaf)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_Jmax_350_t4 , aes(x = Sdt_DD  , y = Jmax), col = 'cyan2', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_Jmax_560_t4 , aes(x = Sdt_DD  , y = Jmax), col = 'cyan3', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_Jmax_1050_t4 , aes(x = Sdt_DD  , y = Jmax), col = 'darkcyan', lwd = 2, alpha = 0.8,linetype="dashed") +
  geom_line(data = regline_Jmax_1680_t4 , aes(x = Sdt_DD  , y = Jmax), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype="dashed") +
  
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
  
  ylab(expression('J'[max] * ' (µmol m' ^ '-2 ' * 's' ^ '-1' * ')')) +
  xlab(expression (italic('Drought Severity Index')))
Jmax_plot_t4
#######################################################################################################################

###################### Anet #####################################################################################
names(stat_t4)
hist(log(stat_t4$anet_420))
hist(stat_t4$anet_420)
Anet_lmer_t4 <- lm(log(anet_420)~ Cumul_N*Sdt_DD,
                     data = stat_t4)
plot(Anet_lmer_t4)
densityPlot(residuals(Anet_lmer_t4))
shapiro.test(residuals(Anet_lmer_t4)) # 0.814
outlierTest(Anet_lmer_t4)  ### two outliers to remove 

residuals_data <- data.frame(
  row_index = 1:nrow(stat_t4),
  residuals = rstudent(Anet_lmer_t4)
)
hist(residuals_data$residuals)
threshold <- 3
outliers <- abs(residuals_data$residuals) > threshold
# Filter the dataset to remove outliers
stat_t4_anet_clean <- stat_t4[!outliers, ]
nrow(stat_t4_anet_clean) # 46 , three removed
hist(stat_t4_anet_clean$anet_420) # better dist

############### Refit ##########################################
Anet_lmer_t4 <- lm(log(anet_420)~ Cumul_N*Sdt_DD,
                     data = stat_t4_anet_clean)
plot(Anet_lmer_t4)
densityPlot(residuals(Anet_lmer_t4))
shapiro.test(residuals(Anet_lmer_t4)) # 0.755
outlierTest(Anet_lmer_t4)  ### two outliers to remove 
Anova(Anet_lmer_t4) # N  effects only and slight interactions 
AIC(Anet_lmer_t4) # -11.44
vif(Anet_lmer_t4) # no colinerarity
plot(resid(Anet_lmer_t4) ~ fitted(Anet_lmer_t4)) ## good
r.squaredGLMM(Anet_lmer_t4) ## marginal : 0.39, conditional: 0.39
summary(Anet_lmer_t4) 
qqnorm(residuals(Anet_lmer_t4))
qqline(residuals(Anet_lmer_t4))
residuals <- resid(Anet_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_Anet_t4<- test(emtrends(Anet_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_Anet_t4 <- test(emtrends(Anet_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_Anet_t4 <- test(emtrends(Anet_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_Anet_t4 <- test(emtrends(Anet_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_Anet_t4<- summary(emmeans(Anet_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_Anet_t4 <- summary(emmeans(Anet_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_Anet_t4 <- summary(emmeans(Anet_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_Anet_t4 <- summary(emmeans(Anet_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_Anet_350_t4 <- exp(intercept_350_Anet_t4 + slope_350_Anet_t4* Sdt_DD_range)
reg_Anet_560_t4  <- exp(intercept_560_Anet_t4 + slope_560_Anet_t4 * Sdt_DD_range)
reg_Anet_1050_t4  <- exp(intercept_1050_Anet_t4 + slope_1050_Anet_t4 * Sdt_DD_range)
reg_Anet_1680_t4  <- exp(intercept_1680_Anet_t4 + slope_1680_Anet_t4 * Sdt_DD_range)

regline_Anet_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Anet = reg_Anet_350_t4 , Cumul_N = 350)
regline_Anet_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Anet = reg_Anet_560_t4 , Cumul_N = 560)
regline_Anet_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Anet = reg_Anet_1050_t4 , Cumul_N = 1050)
regline_Anet_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, Anet = reg_Anet_1680_t4 , Cumul_N = 1680)

slope_results_Anet_t4 <- emtrends(Anet_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_Anet_t4 <- contrast(slope_results_Anet_t4, method = "pairwise")
slope_summary_Anet_t4 <- summary(slope_results_Anet_t4, infer = c(TRUE, TRUE))
slope_summary_Anet_t4

slope_test_anet_t4 <- cld(emtrends(Anet_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

############################################ plot ############################
Anet_plot_t4 <- ggplot(data = stat_t4_anet_clean, aes(x = Sdt_DD, y = anet_420)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_Anet_350_t4 , aes(x = Sdt_DD  , y = Anet), col = 'cyan2', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Anet_560_t4 , aes(x = Sdt_DD  , y = Anet), col = 'cyan3', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Anet_1050_t4 , aes(x = Sdt_DD  , y = Anet), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_Anet_1680_t4 , aes(x = Sdt_DD  , y = Anet), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype="dashed") +
  
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
  
  ylab(expression(A[net] ~ (µmol ~ CO[2] ~ m^{-2} ~ s^{-1}))) +
  xlab(expression (italic('Drought Severity Index')))
Anet_plot_t4
######################################################################################################################

###################### Chlorophyll #####################################################################################
names(stat_t4)
hist(log(stat_t4$chl_mmolm2))
hist(stat_t4$chl_mmolm2)
Chl_lmer_t4 <- lm(log(chl_mmolm2)~ Cumul_N*Sdt_DD,
                    data = stat_t4)
plot(Chl_lmer_t4)
shapiro.test(residuals(Chl_lmer_t4)) # 0.83
outlierTest(Chl_lmer_t4)  ### No outliers
Anova(Chl_lmer_t4) # N effect 
AIC(Chl_lmer_t4) # -9.47
vif(Chl_lmer_t4) # no colinerarity
plot(resid(Chl_lmer_t4) ~ fitted(Chl_lmer_t4)) ## good
r.squaredGLMM(Chl_lmer_t4) ## marginal : 0.10, conditional: 0.10
summary(Chl_lmer_t4) 
qqnorm(residuals(Chl_lmer_t4))
qqline(residuals(Chl_lmer_t4))
densityPlot(residuals(Chl_lmer_t4))

residuals <- resid(Chl_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_Chl_t4<- test(emtrends(Chl_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_Chl_t4 <- test(emtrends(Chl_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_Chl_t4 <- test(emtrends(Chl_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_Chl_t4 <- test(emtrends(Chl_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_Chl_t4 <- summary(emmeans(Chl_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_Chl_t4 <- summary(emmeans(Chl_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_Chl_t4 <- summary(emmeans(Chl_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_Chl_t4<- summary(emmeans(Chl_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_Chl_350_t4 <- exp(intercept_350_Chl_t4 + slope_350_Chl_t4 * Sdt_DD_range)
reg_Chl_560_t4  <- exp(intercept_560_Chl_t4 + slope_560_Chl_t4 * Sdt_DD_range)
reg_Chl_1050_t4  <- exp(intercept_1050_Chl_t4 + slope_1050_Chl_t4 * Sdt_DD_range)
reg_Chl_1680_t4  <- exp(intercept_1680_Chl_t4 + slope_1680_Chl_t4 * Sdt_DD_range)

regline_Chl_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Chl = reg_Chl_350_t4 , Cumul_N = 350)
regline_Chl_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Chl = reg_Chl_560_t4 , Cumul_N = 560)
regline_Chl_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Chl = reg_Chl_1050_t4 , Cumul_N = 1050)
regline_Chl_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, Chl = reg_Chl_1680_t4 , Cumul_N = 1680)

slope_results_Chl_t4 <- emtrends(Chl_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_Chl_t4 <- contrast(slope_results_Chl_t4, method = "pairwise")
slope_summary_Chl_t4 <- summary(slope_results_Chl_t4, infer = c(TRUE, TRUE))
slope_summary_Chl_t4

slope_test_chl_t4 <- cld(emtrends(Chl_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

Chl_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = chl_mmolm2)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan3", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_Chl_350_t4 , aes(x = Sdt_DD  , y = Chl), col = 'cyan2', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_Chl_560_t4 , aes(x = Sdt_DD  , y = Chl), col = 'cyan3', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_Chl_1050_t4 , aes(x = Sdt_DD  , y = Chl), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_Chl_1680_t4 , aes(x = Sdt_DD  , y = Chl), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype="dashed") +
  
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
  
  ylab(expression('Chl'[area] * ' (mmol m' ^ '-2 ' * ')')) +
  xlab(expression (italic('Drought Severity Index')))
Chl_plot_t4
######################################################################################################################

####################### propo nrubisco #####################################################################################
names(stat_t4)
hist(log(stat_t4$nrubisco))
hist(stat_t4$nrubisco)
Nrubisco_lmer_t4 <- lm(nrubisco~ Cumul_N*Sdt_DD,
                     data = stat_t4)
plot(Nrubisco_lmer_t4)
shapiro.test(residuals(Nrubisco_lmer_t4)) # 0.42
outlierTest(Nrubisco_lmer_t4)  ### No outliers

Anova(Nrubisco_lmer_t4, type="II") # nothing 
AIC(Nrubisco_lmer_t4) #-54.40
vif(Nrubisco_lmer_t4) # no colinerarity
plot(resid(Nrubisco_lmer_t4) ~ fitted(Nrubisco_lmer_t4)) ## good
r.squaredGLMM(Nrubisco_lmer_t4) ## marginal : 0.017, conditional: 0.01
summary(Nrubisco_lmer_t4) 
qqnorm(residuals(Nrubisco_lmer_t4))
qqline(residuals(Nrubisco_lmer_t4))
densityPlot(residuals(Nrubisco_lmer_t4))

residuals <- resid(Nrubisco_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
trends_results_nrubisco_t4 <- emtrends(Nrubisco_lmer_t4, ~ Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
print(trends_results_nrubisco_t4)
trends_results_df_nrubisco_t4 <- as.data.frame(trends_results_nrubisco_t4)

Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_nrubisco_t4 <- test(emtrends(Nrubisco_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_nrubisco_t4 <- test(emtrends(Nrubisco_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_nrubisco_t4 <- test(emtrends(Nrubisco_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_nrubisco_t4 <- test(emtrends(Nrubisco_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_nrubisco_t4 <- summary(emmeans(Nrubisco_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_nrubisco_t4 <- summary(emmeans(Nrubisco_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_nrubisco_t4 <- summary(emmeans(Nrubisco_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_nrubisco_t4<- summary(emmeans(Nrubisco_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_nrubisco_350_t4 <- intercept_350_nrubisco_t4 + slope_350_nrubisco_t4 * Sdt_DD_range
reg_nrubisco_560_t4  <- intercept_560_nrubisco_t4 + slope_560_nrubisco_t4 * Sdt_DD_range
reg_nrubisco_1050_t4  <- intercept_1050_nrubisco_t4 + slope_1050_nrubisco_t4 * Sdt_DD_range
reg_nrubisco_1680_t4  <- intercept_1680_nrubisco_t4 + slope_1680_nrubisco_t4 * Sdt_DD_range

regline_nrubisco_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, nrubisco = reg_nrubisco_350_t4 , Cumul_N = 350)
regline_nrubisco_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, nrubisco = reg_nrubisco_560_t4 , Cumul_N = 560)
regline_nrubisco_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, nrubisco = reg_nrubisco_1050_t4 , Cumul_N = 1050)
regline_nrubisco_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, nrubisco = reg_nrubisco_1680_t4 , Cumul_N = 1680)

slope_results_nrubisco_t4 <- emtrends(Nrubisco_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_nrubisco_t4 <- contrast(slope_results_nrubisco_t4, method = "pairwise")
slope_summary_nrubisco_t4 <- summary(slope_results_nrubisco_t4, infer = c(TRUE, TRUE))
slope_summary_nrubisco_t4

cld(emtrends(Nrubisco_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

nrubisco_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = nrubisco)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_nrubisco_350_t4 , aes(x = Sdt_DD  , y = nrubisco), col = 'cyan2', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_nrubisco_560_t4 , aes(x = Sdt_DD  , y = nrubisco), col = 'cyan3', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_nrubisco_1050_t4 , aes(x = Sdt_DD  , y = nrubisco), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_nrubisco_1680_t4 , aes(x = Sdt_DD  , y = nrubisco), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype="dashed") +
  
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
  
  ylab(expression('ρ'[rubisco] * ' (gN gN' ^ '-1' * ')')) +
  xlab(expression (italic('Drought Severity Index')))
nrubisco_plot_t4
#######################################################################################################################

####################### Prop N bioenergetics #####################################################################################
names(stat_t4)
hist(log(stat_t4$nbioe))
hist(stat_t4$nbioe)
Nbioe_lmer_t4 <- lm(nbioe~ Cumul_N*Sdt_DD,
                         data = stat_t4)
plot(Nbioe_lmer_t4)
shapiro.test(residuals(Nbioe_lmer_t4)) # 0.702
outlierTest(Nbioe_lmer_t4)  ### No outliers
Anova(Nbioe_lmer_t4, type ="II") # nothing 
AIC(Nbioe_lmer_t4) #-219.728
vif(Nbioe_lmer_t4) # no colinerarity
plot(resid(Nbioe_lmer_t4) ~ fitted(Nbioe_lmer_t4)) ## good
r.squaredGLMM(Nbioe_lmer_t4) ## marginal : 0.087, conditional: 0.087
summary(Nbioe_lmer_t4) 
qqnorm(residuals(Nbioe_lmer_t4))
qqline(residuals(Nbioe_lmer_t4))
densityPlot(residuals(Nbioe_lmer_t4))

residuals <- resid(Nbioe_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
trends_results_nbioe_t4 <- emtrends(Nbioe_lmer_t4, ~ Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
print(trends_results_nbioe_t4)
trends_results_df_nbioe_t4 <- as.data.frame(trends_results_nbioe_t4)

Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_nbioe_t4 <- test(emtrends(Nbioe_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_nbioe_t4 <- test(emtrends(Nbioe_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_nbioe_t4 <- test(emtrends(Nbioe_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_nbioe_t4 <- test(emtrends(Nbioe_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_nbioe_t4 <- summary(emmeans(Nbioe_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_nbioe_t4 <- summary(emmeans(Nbioe_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_nbioe_t4 <- summary(emmeans(Nbioe_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_nbioe_t4<- summary(emmeans(Nbioe_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_nbioe_350_t4 <- intercept_350_nbioe_t4 + slope_350_nbioe_t4 * Sdt_DD_range
reg_nbioe_560_t4  <- intercept_560_nbioe_t4 + slope_560_nbioe_t4 * Sdt_DD_range
reg_nbioe_1050_t4  <- intercept_1050_nbioe_t4 + slope_1050_nbioe_t4 * Sdt_DD_range
reg_nbioe_1680_t4  <- intercept_1680_nbioe_t4 + slope_1680_nbioe_t4 * Sdt_DD_range

regline_nbioe_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, nbioe = reg_nbioe_350_t4 , Cumul_N = 350)
regline_nbioe_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, nbioe = reg_nbioe_560_t4 , Cumul_N = 560)
regline_nbioe_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, nbioe = reg_nbioe_1050_t4 , Cumul_N = 1050)
regline_nbioe_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, nbioe = reg_nbioe_1680_t4 , Cumul_N = 1680)

slope_results_nbioe_t4 <- emtrends(Nbioe_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_nbioe_t4 <- contrast(slope_results_nbioe_t4, method = "pairwise")
slope_summary_nbioe_t4 <- summary(slope_results_nbioe_t4, infer = c(TRUE, TRUE))
slope_summary_nbioe_t4

cld(emtrends(Nbioe_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

nbioe_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = nbioe)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_nbioe_350_t4 , aes(x = Sdt_DD  , y = nbioe), col = 'cyan2', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_nbioe_560_t4 , aes(x = Sdt_DD  , y = nbioe), col = 'cyan3', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_nbioe_1050_t4 , aes(x = Sdt_DD  , y = nbioe), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_nbioe_1680_t4 , aes(x = Sdt_DD  , y = nbioe), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype="dashed") +
  
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
  
  ylab(expression('ρ'[bioenergetics] * ' (gN gN' ^ '-1' * ')')) +
  xlab(expression (italic('Drought Severity Index')))
nbioe_plot_t4
#######################################################################################################################

###################### N light harvesting #####################################################################################
names(stat_t4)
hist(log(stat_t4$nlightharvesting))
hist(stat_t4$nlightharvesting)
Nlh_lmer_t4 <- lm(nlightharvesting~ Cumul_N*Sdt_DD,
                    data = stat_t4)
shapiro.test(residuals(Nlh_lmer_t4)) # 0.69
outlierTest(Nlh_lmer_t4)  ### No outliers
Anova(Nlh_lmer_t4, type ="II") # N effect only
AIC(Nlh_lmer_t4) # -275.0029
vif(Nlh_lmer_t4) # no colinerarity
plot(resid(Nlh_lmer_t4) ~ fitted(Nlh_lmer_t4)) ## good
r.squaredGLMM(Nlh_lmer_t4) ## marginal : 0.23, conditional: 0.23
summary(Nlh_lmer_t4) 
qqnorm(residuals(Nlh_lmer_t4))
qqline(residuals(Nlh_lmer_t4))
densityPlot(residuals(Nlh_lmer_t4))

residuals <- resid(Nlh_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
trends_results_Nlh_t4 <- emtrends(Nlh_lmer_t4, ~ Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
print(trends_results_Nlh_t4)
trends_results_df_Nlh_t4 <- as.data.frame(trends_results_Nlh_t4)

Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_Nlh_t4<- test(emtrends(Nlh_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_Nlh_t4 <- test(emtrends(Nlh_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_Nlh_t4 <- test(emtrends(Nlh_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_Nlh_t4 <- test(emtrends(Nlh_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_Nlh_t4<- summary(emmeans(Nlh_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_Nlh_t4 <- summary(emmeans(Nlh_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_Nlh_t4 <- summary(emmeans(Nlh_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_Nlh_t4 <- summary(emmeans(Nlh_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_Nlh_350_t4 <- intercept_350_Nlh_t4 + slope_350_Nlh_t4* Sdt_DD_range
reg_Nlh_560_t4  <- intercept_560_Nlh_t4 + slope_560_Nlh_t4 * Sdt_DD_range
reg_Nlh_1050_t4  <- intercept_1050_Nlh_t4 + slope_1050_Nlh_t4 * Sdt_DD_range
reg_Nlh_1680_t4  <- intercept_1680_Nlh_t4 + slope_1680_Nlh_t4 * Sdt_DD_range

regline_Nlh_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Nlh = reg_Nlh_350_t4 , Cumul_N = 350)
regline_Nlh_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Nlh = reg_Nlh_560_t4 , Cumul_N = 560)
regline_Nlh_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Nlh = reg_Nlh_1050_t4 , Cumul_N = 1050)
regline_Nlh_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, Nlh = reg_Nlh_1680_t4 , Cumul_N = 1680)

slope_results_Nlh_t4 <- emtrends(Nlh_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_Nlh_t4 <- contrast(slope_results_Nlh_t4, method = "pairwise")
slope_summary_Nlh_t4 <- summary(slope_results_Nlh_t4, infer = c(TRUE, TRUE))
slope_summary_Nlh_t4

cld(emtrends(Nlh_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

Nlh_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = nlightharvesting)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_Nlh_350_t4 , aes(x = Sdt_DD  , y = Nlh), col = 'cyan2', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_Nlh_560_t4 , aes(x = Sdt_DD  , y = Nlh), col = 'cyan3', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_Nlh_1050_t4 , aes(x = Sdt_DD  , y = Nlh), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_Nlh_1680_t4 , aes(x = Sdt_DD  , y = Nlh), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype="dashed") +
  
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
  
  ylab(expression('ρ'[lightharvesting] * ' (gN gN' ^ '-1' * ')')) +
  xlab(expression (italic('Drought Severity Index')))
Nlh_plot_t4

####################### Nstructure #####################################################################################
names(stat_t4)
hist(log(stat_t4$Nstrucutre))
hist(stat_t4$Nstrucutre)
Nstr_lmer_t4 <- lm(log(Nstrucutre)~ Cumul_N*Sdt_DD,
                      data = stat_t4)
shapiro.test(residuals(Nstr_lmer_t4)) # 0.8036
outlierTest(Nstr_lmer_t4)  ### No outliers

Anova(Nstr_lmer_t4, type ="II") # N and drought 
AIC(Nstr_lmer_t4) # 16.10
vif(Nstr_lmer_t4) # no colinerarity
plot(resid(Nstr_lmer_t4) ~ fitted(Nstr_lmer_t4)) ## good
r.squaredGLMM(Nstr_lmer_t4) ## marginal : 0.48, conditional: 0.48
summary(Nstr_lmer_t4) 
qqnorm(residuals(Nstr_lmer_t4))
qqline(residuals(Nstr_lmer_t4))
densityPlot(residuals(Nstr_lmer_t4))

residuals <- resid(Nstr_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
trends_results_Nstr_t4 <- emtrends(Nstr_lmer_t4, ~ Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
print(trends_results_Nstr_t4)
trends_results_df_Nstr_t4 <- as.data.frame(trends_results_Nstr_t4)

Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_Nstr_t4<- test(emtrends(Nstr_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_Nstr_t4 <- test(emtrends(Nstr_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_Nstr_t4 <- test(emtrends(Nstr_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_Nstr_t4 <- test(emtrends(Nstr_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_Nstr_t4 <- summary(emmeans(Nstr_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_Nstr_t4 <- summary(emmeans(Nstr_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_Nstr_t4 <- summary(emmeans(Nstr_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_Nstr_t4<- summary(emmeans(Nstr_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_Nstr_350_t4 <- exp(intercept_350_Nstr_t4 + slope_350_Nstr_t4 * Sdt_DD_range)
reg_Nstr_560_t4  <- exp(intercept_560_Nstr_t4 + slope_560_Nstr_t4 * Sdt_DD_range)
reg_Nstr_1050_t4  <- exp(intercept_1050_Nstr_t4 + slope_1050_Nstr_t4 * Sdt_DD_range)
reg_Nstr_1680_t4  <- exp(intercept_1680_Nstr_t4 + slope_1680_Nstr_t4 * Sdt_DD_range)

regline_Nstr_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Nstr = reg_Nstr_350_t4 , Cumul_N = 350)
regline_Nstr_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Nstr = reg_Nstr_560_t4 , Cumul_N = 560)
regline_Nstr_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, Nstr = reg_Nstr_1050_t4 , Cumul_N = 1050)
regline_Nstr_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, Nstr = reg_Nstr_1680_t4 , Cumul_N = 1680)

slope_results_Nstr_t4 <- emtrends(Nstr_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_Nstr_t4 <- contrast(slope_results_Nstr_t4, method = "pairwise")
slope_summary_Nstr_t4 <- summary(slope_results_Nstr_t4, infer = c(TRUE, TRUE))
slope_summary_Nstr_t4

cld(emtrends(Nstr_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

Nstr_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = Nstrucutre)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_Nstr_350_t4 , aes(x = Sdt_DD  , y = Nstr), col = 'cyan2', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Nstr_560_t4 , aes(x = Sdt_DD  , y = Nstr), col = 'cyan3', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Nstr_1050_t4 , aes(x = Sdt_DD  , y = Nstr), col = 'darkcyan', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_Nstr_1680_t4 , aes(x = Sdt_DD  , y = Nstr), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype="dashed") +
  
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
  
  ylab(expression('ρ'[structure] * ' (gN gN' ^ '-1' * ')')) +
  xlab(expression (italic('Drought Severity Index')))
Nstr_plot_t4
#######################################################################################################################
names(stat_t4)
stat_t4$PNUE = stat_t4$anet_420/stat_t4$Narea ###### leaf N use efficiency
stat_t4$LUE = stat_t4$anet_420/stat_t4$gs_420 ################# leaf water use efficiency
stat_t4_PNUE = subset(stat_t4, PNUE<40)
hist(stat_t4_PNUE$PNUE)
hist(log(stat_t4$PNUE))

hist(stat_t4$LUE)
stat_t4_LUE = subset(stat_t4, LUE<150)
hist(stat_t4_LUE$LUE)

hist(log(stat_t4_LUE$LUE))

LUE_lmer_t4 <- lm(log(LUE)~ Cumul_N*Sdt_DD,
                     data = stat_t4)
Anova(LUE_lmer_t4, type ="II") # N and drought effects, no interaction  
shapiro.test(residuals(LUE_lmer_t4)) 
outlierTest(LUE_lmer_t4)  ### No outliers
AIC(LUE_lmer_t4) # 58.97
vif(LUE_lmer_t4) # no colinerarity
plot(resid(LUE_lmer_t4) ~ fitted(LUE_lmer_t4)) ## good
r.squaredGLMM(LUE_lmer_t4) ## marginal : 0.48, conditional: 0.48
summary(LUE_lmer_t4) 
qqnorm(residuals(LUE_lmer_t4))
qqline(residuals(LUE_lmer_t4))
densityPlot(residuals(LUE_lmer_t4))

residuals <- resid(LUE_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_LUE_t4<- test(emtrends(LUE_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_LUE_t4 <- test(emtrends(LUE_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_LUE_t4 <- test(emtrends(LUE_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_LUE_t4 <- test(emtrends(LUE_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_LUE_t4 <- summary(emmeans(LUE_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_LUE_t4 <- summary(emmeans(LUE_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_LUE_t4 <- summary(emmeans(LUE_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_LUE_t4<- summary(emmeans(LUE_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_LUE_350_t4 <- exp(intercept_350_LUE_t4 + slope_350_LUE_t4 * Sdt_DD_range)
reg_LUE_560_t4  <- exp(intercept_560_LUE_t4 + slope_560_LUE_t4 * Sdt_DD_range)
reg_LUE_1050_t4  <- exp(intercept_1050_LUE_t4 + slope_1050_LUE_t4 * Sdt_DD_range)
reg_LUE_1680_t4  <- exp(intercept_1680_LUE_t4 + slope_1680_LUE_t4 * Sdt_DD_range)

regline_LUE_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, LUE = reg_LUE_350_t4 , Cumul_N = 350)
regline_LUE_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, LUE = reg_LUE_560_t4 , Cumul_N = 560)
regline_LUE_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, LUE = reg_LUE_1050_t4 , Cumul_N = 1050)
regline_LUE_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, LUE = reg_LUE_1680_t4 , Cumul_N = 1680)

slope_results_LUE_t4 <- emtrends(LUE_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_LUE_t4 <- contrast(slope_results_LUE_t4, method = "pairwise")
slope_summary_LUE_t4 <- summary(slope_results_LUE_t4, infer = c(TRUE, TRUE))
slope_summary_LUE_t4

cld(emtrends(LUE_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

WUE_plot_t4 <- ggplot(data = stat_t4_LUE, aes(x = Sdt_DD, y = LUE)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_LUE_350_t4 , aes(x = Sdt_DD  , y = LUE), col = 'cyan2', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_LUE_560_t4 , aes(x = Sdt_DD  , y = LUE), col = 'cyan3', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_LUE_1050_t4 , aes(x = Sdt_DD  , y = LUE), col = 'darkcyan', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_LUE_1680_t4 , aes(x = Sdt_DD  , y = LUE), col = 'lightcyan4', lwd = 2, alpha = 0.8) +
  
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
  
  ylab(expression(WUE[leaf] ~ (µmol ~ CO[2] ~ µmol^{-1} ~ H[2]*O ~ s^{-1}))) +
  xlab(expression (italic('Drought Severity Index')))
WUE_plot_t4

############################## PNUE ##################################################

PNUE_lmer_t4 <- lm(log(PNUE)~ Cumul_N*Sdt_DD,
                     data = stat_t4)
Anova(PNUE_lmer_t4, type = "II") # Nothing 

shapiro.test(residuals(PNUE_lmer_t4)) 
outlierTest(PNUE_lmer_t4)  ### one outlier, it is ok

AIC(PNUE_lmer_t4) # 29.06
vif(PNUE_lmer_t4) # no colinerarity
plot(resid(PNUE_lmer_t4) ~ fitted(PNUE_lmer_t4)) ## good
r.squaredGLMM(PNUE_lmer_t4) ## marginal : 0.11, conditional: 0.11
summary(PNUE_lmer_t4) 
qqnorm(residuals(PNUE_lmer_t4))
qqline(residuals(PNUE_lmer_t4))
densityPlot(residuals(PNUE_lmer_t4))

residuals <- resid(PNUE_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_PNUE_t4<- test(emtrends(PNUE_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_PNUE_t4 <- test(emtrends(PNUE_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_PNUE_t4 <- test(emtrends(PNUE_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_PNUE_t4 <- test(emtrends(PNUE_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_PNUE_t4 <- summary(emmeans(PNUE_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_PNUE_t4 <- summary(emmeans(PNUE_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_PNUE_t4 <- summary(emmeans(PNUE_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_PNUE_t4<- summary(emmeans(PNUE_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_PNUE_350_t4 <- exp(intercept_350_PNUE_t4 + slope_350_PNUE_t4 * Sdt_DD_range)
reg_PNUE_560_t4  <- exp(intercept_560_PNUE_t4 + slope_560_PNUE_t4 * Sdt_DD_range)
reg_PNUE_1050_t4  <- exp(intercept_1050_PNUE_t4 + slope_1050_PNUE_t4 * Sdt_DD_range)
reg_PNUE_1680_t4  <- exp(intercept_1680_PNUE_t4 + slope_1680_PNUE_t4 * Sdt_DD_range)

regline_PNUE_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, PNUE = reg_PNUE_350_t4 , Cumul_N = 350)
regline_PNUE_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, PNUE = reg_PNUE_560_t4 , Cumul_N = 560)
regline_PNUE_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, PNUE = reg_PNUE_1050_t4 , Cumul_N = 1050)
regline_PNUE_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, PNUE = reg_PNUE_1680_t4 , Cumul_N = 1680)

slope_results_PNUE_t4 <- emtrends(PNUE_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_PNUE_t4 <- contrast(slope_results_PNUE_t4, method = "pairwise")
slope_summary_PNUE_t4 <- summary(slope_results_PNUE_t4, infer = c(TRUE, TRUE))
slope_summary_PNUE_t4

PNUE_plot_t4 <- ggplot(data = stat_t4_PNUE, aes(x = Sdt_DD, y = PNUE)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_PNUE_350_t4 , aes(x = Sdt_DD  , y = PNUE), col = 'cyan2', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_PNUE_560_t4 , aes(x = Sdt_DD  , y = PNUE), col = 'cyan3', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_PNUE_1050_t4 , aes(x = Sdt_DD  , y = PNUE), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_PNUE_1680_t4 , aes(x = Sdt_DD  , y = PNUE), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype="dashed") +
  
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
  
  ylab(expression(NUE[leaf] ~ (µmol ~ CO[2] ~ g ~ N^{-1}))) +
  xlab(expression (italic('Drought Severity Index')))
PNUE_plot_t4

###################### plant_surface_area #####################################################################################
names(stat_t4)
hist(log(stat_t4$plant_surface_area))
PSA_lmer_t4 <- lm(log(plant_surface_area)~ Cumul_N*Sdt_DD,
                     data = stat_t4)
shapiro.test(residuals(PSA_lmer_t4)) # 0.067
outlierTest(PSA_lmer_t4)  ### 
Anova(PSA_lmer_t4, type = "II") # N and drought*N effects
AIC(PSA_lmer_t4) # 14.49
vif(PSA_lmer_t4) # no colinerarity
plot(resid(PSA_lmer_t4) ~ fitted(PSA_lmer_t4)) ## good
r.squaredGLMM(PSA_lmer_t4) ## marginal : 0.62, conditional: 0.62
summary(PSA_lmer_t4) 
qqnorm(residuals(PSA_lmer_t4))
qqline(residuals(PSA_lmer_t4))
densityPlot(residuals(PSA_lmer_t4))

residuals <- resid(PSA_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_PSA_t4<- test(emtrends(PSA_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_PSA_t4 <- test(emtrends(PSA_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_PSA_t4 <- test(emtrends(PSA_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_PSA_t4 <- test(emtrends(PSA_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_PSA_t4<- summary(emmeans(PSA_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_PSA_t4 <- summary(emmeans(PSA_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_PSA_t4 <- summary(emmeans(PSA_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_PSA_t4 <- summary(emmeans(PSA_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_PSA_350_t4 <- exp(intercept_350_PSA_t4 + slope_350_PSA_t4* Sdt_DD_range)
reg_PSA_560_t4  <- exp(intercept_560_PSA_t4 + slope_560_PSA_t4 * Sdt_DD_range)
reg_PSA_1050_t4  <- exp(intercept_1050_PSA_t4 + slope_1050_PSA_t4 * Sdt_DD_range)
reg_PSA_1680_t4  <- exp(intercept_1680_PSA_t4 + slope_1680_PSA_t4 * Sdt_DD_range)

regline_PSA_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, PSA = reg_PSA_350_t4 , Cumul_N = 350)
regline_PSA_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, PSA = reg_PSA_560_t4 , Cumul_N = 560)
regline_PSA_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, PSA = reg_PSA_1050_t4 , Cumul_N = 1050)
regline_PSA_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, PSA = reg_PSA_1680_t4 , Cumul_N = 1680)

slope_results_PSA_t4 <- emtrends(PSA_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_PSA_t4 <- contrast(slope_results_PSA_t4, method = "pairwise")
slope_summary_PSA_t4 <- summary(slope_results_PSA_t4, infer = c(TRUE, TRUE))
slope_summary_PSA_t4

cld(emtrends(PSA_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

PSA_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = plant_surface_area)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_PSA_350_t4 , aes(x = Sdt_DD  , y = PSA), col = 'cyan2', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_PSA_560_t4 , aes(x = Sdt_DD  , y = PSA), col = 'cyan3', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_PSA_1050_t4 , aes(x = Sdt_DD  , y = PSA), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_PSA_1680_t4 , aes(x = Sdt_DD  , y = PSA), col = 'lightcyan4', lwd = 2, alpha = 0.8) +
  
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
  
  ylab(expression('Total Leaf Area (cm' ^ '2' * ')')) +
  xlab(expression (italic('Drought Severity Index')))
PSA_plot_t4
######################################################################################################################

###################### Biomass #####################################################################################
names(stat_t4)
hist(log(stat_t4$biomass_all))
Biomass_lmer_t4 <- lm(log(biomass_all)~ Cumul_N*Sdt_DD,
                    data = stat_t4)
shapiro.test(residuals(Biomass_lmer_t4)) # 0.147
outlierTest(Biomass_lmer_t4)  ### p = 0.33

Anova(Biomass_lmer_t4, type = "II") # N and drought*N effects
AIC(Biomass_lmer_t4) # 46.72
vif(Biomass_lmer_t4) # no colinerarity
plot(resid(Biomass_lmer_t4) ~ fitted(Biomass_lmer_t4)) ## good
r.squaredGLMM(Biomass_lmer_t4) ## marginal : 0.53, conditional: 0.53
summary(Biomass_lmer_t4) 
qqnorm(residuals(Biomass_lmer_t4))
qqline(residuals(Biomass_lmer_t4))
densityPlot(residuals(Biomass_lmer_t4))

residuals <- resid(Biomass_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
trends_results_biomass_t4 <- emtrends(PSA_lmer_t4, ~ Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
print(trends_results_biomass_t4)
trends_results_df_biomass_t4 <- as.data.frame(trends_results_biomass_t4)

Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_biomass_t4 <- test(emtrends(Biomass_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_biomass_t4 <- test(emtrends(Biomass_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_biomass_t4 <- test(emtrends(Biomass_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_biomass_t4 <- test(emtrends(Biomass_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_biomass_t4 <- summary(emmeans(Biomass_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_biomass_t4 <- summary(emmeans(Biomass_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_biomass_t4 <- summary(emmeans(Biomass_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_biomass_t4 <- summary(emmeans(Biomass_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_biomass_350_t4 <- exp(intercept_350_biomass_t4 + slope_350_biomass_t4* Sdt_DD_range)
reg_biomass_560_t4  <- exp(intercept_560_biomass_t4 + slope_560_biomass_t4 * Sdt_DD_range)
reg_biomass_1050_t4  <- exp(intercept_1050_biomass_t4 + slope_1050_biomass_t4 * Sdt_DD_range)
reg_biomass_1680_t4  <- exp(intercept_1680_biomass_t4 + slope_1680_biomass_t4 * Sdt_DD_range)

regline_biomass_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, biomass = reg_biomass_350_t4 , Cumul_N = 350)
regline_biomass_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, biomass = reg_biomass_560_t4 , Cumul_N = 560)
regline_biomass_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, biomass = reg_biomass_1050_t4 , Cumul_N = 1050)
regline_biomass_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, biomass = reg_biomass_1680_t4 , Cumul_N = 1680)

slope_results_Biomass_t4 <- emtrends(Biomass_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_Biomass_t4 <- contrast(slope_results_Biomass_t4, method = "pairwise")
slope_summary_Biomass_t4 <- summary(slope_results_Biomass_t4, infer = c(TRUE, TRUE))
slope_summary_Biomass_t4

cld(emtrends(Biomass_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

biomass_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = biomass_t4)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_biomass_350_t4 , aes(x = Sdt_DD  , y = biomass), col = 'cyan2', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_biomass_560_t4 , aes(x = Sdt_DD  , y = biomass), col = 'cyan3', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_biomass_1050_t4 , aes(x = Sdt_DD  , y = biomass), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_biomass_1680_t4 , aes(x = Sdt_DD  , y = biomass), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype="dashed") +
  
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
  
  ylab(expression('Biomass (g)')) +
  xlab(expression (italic('Drought Severity Index')))
biomass_plot_t4
######################################################################################################################

###################### Ratios #####################################################################################
names(stat_t4)
plot(stat_t4$N_total_leaves,stat_t4$Abg_t4)
hist(stat_t4$RS_all)
RS_lmer_t4 <- lm(RS_all~ Cumul_N*Sdt_DD,
                        data = stat_t4)
outlierTest(RS_lmer_t4)
Anova(RS_lmer_t4, type = "II") # Nothing 


######################## Cost transpiration ############################################

Cost_trans_roots_lmer_t4 <- lm(Ecost~ Cumul_N*Sdt_DD,
                           data = stat_t4)

Anova(Cost_trans_roots_lmer_t4, type ="II") # Nothing
summary(Cost_trans_roots_lmer_t4)
shapiro.test(residuals(Cost_trans_roots_lmer_t4)) 
outlierTest(Cost_trans_roots_lmer_t4)  ### p = 0.11
AIC(Cost_trans_roots_lmer_t4) # 94.13
vif(Cost_trans_roots_lmer_t4) # no colinerarity
plot(resid(Cost_trans_roots_lmer_t4) ~ fitted(Cost_trans_roots_lmer_t4)) ## good
r.squaredGLMM(Cost_trans_roots_lmer_t4) ## marginal : 0.05, conditional: 0.05
summary(Cost_trans_roots_lmer_t4) 
qqnorm(residuals(Cost_trans_roots_lmer_t4))
qqline(residuals(Cost_trans_roots_lmer_t4))
densityPlot(residuals(Cost_trans_roots_lmer_t4))

residuals <- resid(Cost_trans_roots_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_cost_trans_t4 <- test(emtrends(Cost_trans_roots_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_cost_trans_t4 <- test(emtrends(Cost_trans_roots_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_cost_trans_t4 <- test(emtrends(Cost_trans_roots_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_cost_trans_t4 <- test(emtrends(Cost_trans_roots_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_cost_trans_t4 <- summary(emmeans(Cost_trans_roots_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_cost_trans_t4 <- summary(emmeans(Cost_trans_roots_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_cost_trans_t4 <- summary(emmeans(Cost_trans_roots_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_cost_trans_t4 <- summary(emmeans(Cost_trans_roots_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_cost_trans_350_t4 <- intercept_350_cost_trans_t4 + slope_350_cost_trans_t4* Sdt_DD_range
reg_cost_trans_560_t4  <- intercept_560_cost_trans_t4 + slope_560_cost_trans_t4 * Sdt_DD_range
reg_cost_trans_1050_t4  <- intercept_1050_cost_trans_t4 + slope_1050_cost_trans_t4 * Sdt_DD_range
reg_cost_trans_1680_t4  <- intercept_1680_cost_trans_t4 + slope_1680_cost_trans_t4 * Sdt_DD_range

regline_cost_trans_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, costrans = reg_cost_trans_350_t4 , Cumul_N = 350)
regline_cost_trans_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, costrans = reg_cost_trans_560_t4 , Cumul_N = 560)
regline_cost_trans_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, costrans = reg_cost_trans_1050_t4 , Cumul_N = 1050)
regline_cost_trans_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, costrans = reg_cost_trans_1680_t4 , Cumul_N = 1680)

slope_results_cost_trans_t4 <- emtrends(Cost_trans_roots_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_cost_trans_t4 <- contrast(slope_results_cost_trans_t4, method = "pairwise")
slope_summary_cost_trans_t4 <- summary(slope_results_cost_trans_t4, infer = c(TRUE, TRUE))
slope_summary_cost_trans_t4

cld(emtrends(Cost_trans_roots_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

transp_cost_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = Ecost)) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_cost_trans_350_t4 , aes(x = Sdt_DD  , y = costrans), col = 'cyan2', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_cost_trans_560_t4 , aes(x = Sdt_DD  , y = costrans), col = 'cyan3', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_cost_trans_1050_t4 , aes(x = Sdt_DD  , y = costrans), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_cost_trans_1680_t4 , aes(x = Sdt_DD  , y = costrans), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype="dashed") +
  
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
  
  ylab(expression(italic(E[Cost]) * ' ' * (g * ' roots' * ' ' * g * ' water' ^ '-1')))+
  xlab(expression (italic('Drought Severity Index')))
  transp_cost_plot_t4

############################################ Cost Nitrogen ##################################
hist(stat_t4$cost_nitrogen)
Cost_nitrogen_lmer_t4 <- lm(Ncost ~ Cumul_N*Sdt_DD,
                           data = stat_t4)
Anova(Cost_nitrogen_lmer_t4, type ="II") # N and interaction effects  
shapiro.test(residuals(Cost_nitrogen_lmer_t4)) # 0.22
outlierTest(Cost_nitrogen_lmer_t4)  ### 
AIC(Cost_nitrogen_lmer_t4) # 517.788
vif(Cost_nitrogen_lmer_t4) # no colinerarity
plot(resid(Cost_nitrogen_lmer_t4) ~ fitted(Biomass_lmer_t4)) ## good
r.squaredGLMM(Cost_nitrogen_lmer_t4) ## marginal : 0.238, conditional: 0.24
summary(Cost_nitrogen_lmer_t4) 
qqnorm(residuals(Cost_nitrogen_lmer_t4))
qqline(residuals(Cost_nitrogen_lmer_t4))
densityPlot(residuals(Cost_nitrogen_lmer_t4))

residuals <- resid(Cost_nitrogen_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)

Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_costN_t4 <- test(emtrends(Cost_nitrogen_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_costN_t4 <- test(emtrends(Cost_nitrogen_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_costN_t4 <- test(emtrends(Cost_nitrogen_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_costN_t4 <- test(emtrends(Cost_nitrogen_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_costN_t4 <- summary(emmeans(Cost_nitrogen_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_costN_t4 <- summary(emmeans(Cost_nitrogen_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_costN_t4 <- summary(emmeans(Cost_nitrogen_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_costN_t4 <- summary(emmeans(Cost_nitrogen_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_costN_350_t4 <- intercept_350_costN_t4 + slope_350_costN_t4* Sdt_DD_range
reg_costN_560_t4  <- intercept_560_costN_t4 + slope_560_costN_t4 * Sdt_DD_range
reg_costN_1050_t4  <- intercept_1050_costN_t4 + slope_1050_costN_t4 * Sdt_DD_range
reg_costN_1680_t4  <- intercept_1680_costN_t4 + slope_1680_costN_t4 * Sdt_DD_range

regline_costN_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, costN = reg_costN_350_t4 , Cumul_N = 350)
regline_costN_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, costN = reg_costN_560_t4 , Cumul_N = 560)
regline_costN_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, costN = reg_costN_1050_t4 , Cumul_N = 1050)
regline_costN_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, costN = reg_costN_1680_t4 , Cumul_N = 1680)

slope_results_costN_t4 <- emtrends(Cost_nitrogen_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_costN_t4 <- contrast(slope_results_costN_t4, method = "pairwise")
slope_summary_costN_t4 <- summary(slope_results_costN_t4, infer = c(TRUE, TRUE))
slope_summary_costN_t4

cld(emtrends(Cost_nitrogen_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

N_cost_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = Ncost )) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_costN_350_t4 , aes(x = Sdt_DD  , y = costN), col = 'cyan2', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_costN_560_t4 , aes(x = Sdt_DD  , y = costN), col = 'cyan3', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_costN_1050_t4 , aes(x = Sdt_DD  , y = costN), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_costN_1680_t4 , aes(x = Sdt_DD  , y = costN), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype="dashed") +
  
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
  
  ylab(expression(italic(N[Cost]) * ' ' * (g * ' roots' * ' ' * g * N * ' leaves'^{-1}))) +
  xlab(expression (italic('Drought Severity Index')))
N_cost_plot_t4
######################################################################################################################
############################## betaplant ####################################################
names(stat_t4)

beta_lmer_t4 <- lm(log(betaplant)~ Cumul_N*Sdt_DD,
                    data = stat_t4)
Anova(beta_lmer_t4, type ="II") # N effects and drought
summary(beta_lmer_t4) 

summary(beta_lmer_t4)
shapiro.test(residuals(beta_lmer_t4)) 
outlierTest(beta_lmer_t4)  
AIC(beta_lmer_t4) 
vif(beta_lmer_t4) # no colinerarity
plot(resid(beta_lmer_t4) ~ fitted(beta_lmer_t4)) ## good
r.squaredGLMM(beta_lmer_t4) 
summary(beta_lmer_t4) 
qqnorm(residuals(beta_lmer_t4))
qqline(residuals(beta_lmer_t4))
densityPlot(residuals(beta_lmer_t4))

residuals <- resid(beta_lmer_t4)
hist(residuals, breaks = 20, main = "Histogram of Residuals") ## not bad

cumul_N_levels <- c(350,560, 1050, 1680)
Sdt_DD_range = seq(min(stat_t4$Sdt_DD, na.rm = TRUE), 
                   max(stat_t4$Sdt_DD, na.rm = TRUE), 
                   0.001)
slope_350_beta_t4 <- test(emtrends(beta_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 350)))[1,2]
slope_560_beta_t4 <- test(emtrends(beta_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 560)))[1,2]
slope_1050_beta_t4 <- test(emtrends(beta_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1050)))[1,2]
slope_1680_beta_t4 <- test(emtrends(beta_lmer_t4, ~1, var = "Sdt_DD", at = list(Cumul_N = 1680)))[1,2]

intercept_350_beta_t4 <- summary(emmeans(beta_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 350)))[1,2]
intercept_560_beta_t4 <- summary(emmeans(beta_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 560)))[1,2]
intercept_1050_beta_t4 <- summary(emmeans(beta_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1050)))[1,2]
intercept_1680_beta_t4 <- summary(emmeans(beta_lmer_t4, ~1, at = list(Sdt_DD = 0,Cumul_N = 1680)))[1,2]

reg_beta_350_t4 <- exp(intercept_350_beta_t4 + slope_350_beta_t4* Sdt_DD_range)
reg_beta_560_t4  <- exp(intercept_560_beta_t4 + slope_560_beta_t4 * Sdt_DD_range)
reg_beta_1050_t4  <- exp(intercept_1050_beta_t4 + slope_1050_beta_t4 * Sdt_DD_range)
reg_beta_1680_t4  <- exp(intercept_1680_beta_t4 + slope_1680_beta_t4 * Sdt_DD_range)

regline_beta_350_t4  <- data.frame(Sdt_DD = Sdt_DD_range, beta = reg_beta_350_t4 , Cumul_N = 350)
regline_beta_560_t4  <- data.frame(Sdt_DD = Sdt_DD_range, beta = reg_beta_560_t4 , Cumul_N = 560)
regline_beta_1050_t4  <- data.frame(Sdt_DD = Sdt_DD_range, beta = reg_beta_1050_t4 , Cumul_N = 1050)
regline_beta_1680_t4 <- data.frame(Sdt_DD = Sdt_DD_range, beta = reg_beta_1680_t4 , Cumul_N = 1680)

slope_results_beta_t4 <- emtrends(beta_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels))
slope_comparison_beta_t4 <- contrast(slope_results_beta_t4, method = "pairwise")
slope_summary_beta_t4 <- summary(slope_results_beta_t4, infer = c(TRUE, TRUE))
slope_summary_beta_t4

cld(emtrends(beta_lmer_t4, ~Cumul_N, var = "Sdt_DD", at = list(Cumul_N = cumul_N_levels)))

############################ Beta plant plot #################################################################
Beta_plot_t4 <- ggplot(data = stat_t4, aes(x = Sdt_DD, y = betaplant )) + 
  
  geom_point(aes(shape = factor(drought), color = Cumul_N), size = 3, alpha=0.9, width=20) +
  scale_shape_manual(values = c(17, 15, 16), name = "Drought Scenario") +
  scale_color_gradient(low = "cyan2", high = "lightcyan4", name = "N level") +
  
  geom_line(data = regline_beta_350_t4 , aes(x = Sdt_DD  , y = beta), col = 'cyan2', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_beta_560_t4 , aes(x = Sdt_DD  , y = beta), col = 'cyan3', lwd = 2, alpha = 0.8) +
  geom_line(data = regline_beta_1050_t4 , aes(x = Sdt_DD  , y = beta), col = 'darkcyan', lwd = 2, alpha = 0.8, linetype="dashed") +
  geom_line(data = regline_beta_1680_t4 , aes(x = Sdt_DD  , y = beta), col = 'lightcyan4', lwd = 2, alpha = 0.8, linetype="dashed") +
  
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
        panel.grid.major = element_line(colour = "white")) +
  
  ylab(expression(italic('β') ['plant'])) + 
  xlab(expression (italic('Drought Severity Index')))
Beta_plot_t4

########################################################################################
names(stat_t4)
hist(stat_t4$ratio_J_V)
ratio_J_V_lmer_t4 <- lm(ratio_J_V~ Cumul_N*Sdt_DD,
                     data = stat_t4)
Anova(ratio_J_V_lmer_t4, type="II") # Nothing 

############################### Export Table Satistics ###########################################

################ gs ##########################
anova_results_gs_t4 <- as.data.frame(Anova(gs_lmer_t4))
slope_summary_gs_t4 <- as.data.frame(slope_summary_gs_t4)
anova_results_gs_t4$model <- "gs"
slope_summary_gs_t4$model <- "gs"
################ Chi ##########################
anova_results_chi_t4 <- as.data.frame(Anova(chi_lmer_t4))
slope_summary_chi_t4 <- as.data.frame(slope_summary_chi_t4)
anova_results_chi_t4$model <- "Chi"
slope_summary_chi_t4$model <- "Chi"
################ betaleaf ##########################
anova_results_betaleaf_t4 <- as.data.frame(Anova(betaleaf_lmer_t4))
slope_summary_betaleaf_t4 <- as.data.frame(slope_summary_betaleaf_t4)
anova_results_betaleaf_t4$model <- "betaleaf"
slope_summary_betaleaf_t4$model <- "betaleaf"
################ Nmass  ##########################
anova_results_Nmass_t4 <- as.data.frame(Anova(Nmass_lmer_t4))
slope_summary_Nmass_t4 <- as.data.frame(slope_summary_Nmass_t4)
anova_results_Nmass_t4$model <- "Nmass"
slope_summary_Nmass_t4$model <- "Nmass"
################ Narea ##########################slope_summary_Narea_df <- as.data.frame(slope_summary_Narea)
anova_results_Narea_t4 <- as.data.frame(Anova(Narea_lmer_t4))
slope_summary_Narea_t4 <- as.data.frame(slope_summary_Narea_t4)
anova_results_Narea_t4$model <- "Narea"
slope_summary_Narea_t4$model <- "Narea"
################ LMA ##########################
anova_results_LMA_t4 <- as.data.frame(Anova(LMA_lmer_t4))
slope_summary_LMA_t4 <- as.data.frame(slope_summary_LMA_t4)
anova_results_LMA_t4$model <- "LMA"
slope_summary_LMA_t4$model <- "LMA"
################ Vcmax ##########################
anova_results_Vcmax_t4 <- as.data.frame(Anova(Vcmax_lmer_t4))
slope_summary_Vcmax_t4 <- as.data.frame(slope_summary_Vcmax_t4)
anova_results_Vcmax_t4$model <- "Vcmax"
slope_summary_Vcmax_t4$model <- "Vcmax"
################ Jmax ##########################
anova_results_Jmax_t4 <- as.data.frame(Anova(Jmax_lmer_t4))
slope_summary_Jmax_t4 <- as.data.frame(slope_summary_Jmax_t4)
anova_results_Jmax_t4$model <- "Jmax"
slope_summary_Jmax_t4$model <- "Jmax"
################ Nrubisco ##########################
anova_results_nrubisco_t4 <- as.data.frame(Anova(Nrubisco_lmer_t4))
slope_summary_nrubisco_t4 <- as.data.frame(slope_summary_nrubisco_t4)
anova_results_nrubisco_t4$model <- "nrubisco"
slope_summary_nrubisco_t4$model <- "nrubisco"
################# nbioe ################
anova_results_nbioe_t4 <- as.data.frame(Anova(Nbioe_lmer_t4))
slope_summary_nbioe_t4 <- as.data.frame(slope_summary_nbioe_t4)
anova_results_nbioe_t4$model <- "nbioe"
slope_summary_nbioe_t4$model <- "nbioe"
############# Nstructure ###################
anova_results_nstr_t4 <- as.data.frame(Anova(Nstr_lmer_t4))
slope_summary_Nstr_t4 <- as.data.frame(slope_summary_Nstr_t4)
anova_results_nstr_t4$model <- "nstr"
slope_summary_Nstr_t4$model <- "nstr"
############# Chlorophyll ###################
anova_results_Chl_t4 <- as.data.frame(Anova(Chl_lmer_t4))
slope_summary_Chl_t4 <- as.data.frame(slope_summary_Chl_t4)
anova_results_Chl_t4$model <- "Chl"
slope_summary_Chl_t4$model <- "Chl"
###############N light harvesting ####################
anova_results_Nlh_t4 <- as.data.frame(Anova(Nlh_lmer_t4))
slope_summary_Nlh_t4 <- as.data.frame(slope_summary_Nlh_t4)
anova_results_Nlh_t4$model <- "Nlh"
slope_summary_Nlh_t4$model <- "Nlh"
##############Anet###############
anova_results_Anet_t4 <- as.data.frame(Anova(Anet_lmer_t4))
anova_results_Anet_t4$model <- "Anet"
slope_summary_Anet_t4 <- as.data.frame(slope_summary_Anet_t4)
slope_summary_Anet_t4$model <- "Anet"
##############WUE###############
anova_results_WUE_t4 <- as.data.frame(Anova(LUE_lmer_t4))
anova_results_WUE_t4$model <- "WUE"
slope_summary_WUE_t4 <- as.data.frame(slope_summary_LUE_t4)
slope_summary_WUE_t4$model <- "WUE"
##############NUE###############
anova_results_PNUE_t4 <- as.data.frame(Anova(PNUE_lmer_t4))
anova_results_PNUE_t4$model <- "PNUE"
slope_summary_PNUE_t4 <- as.data.frame(slope_summary_PNUE_t4)
slope_summary_PNUE_t4$model <- "PNUE"
############# PSA ###################
anova_results_PSA_t4 <- as.data.frame(Anova(PSA_lmer_t4))
slope_summary_PSA_t4 <- as.data.frame(slope_summary_PSA_t4)
anova_results_PSA_t4$model <- "PSA"
slope_summary_PSA_t4$model <- "PSA"
############# Biomass ###################
anova_results_Biomass_t4 <- as.data.frame(Anova(Biomass_lmer_t4))
slope_summary_Biomass_t4 <- as.data.frame(slope_summary_Biomass_t4)
anova_results_Biomass_t4$model <- "Biomass"
slope_summary_Biomass_t4$model <- "Biomass"
############# Cost Nitrogen  ###################
anova_results_costN_t4 <- as.data.frame(Anova(Cost_nitrogen_lmer_t4))
slope_summary_costN_t4 <- as.data.frame(slope_summary_costN_t4)
anova_results_costN_t4$model <- "Ncost"
slope_summary_costN_t4$model <- "Ncost"
############# Cost transpiration  ###################
anova_results_costtrans_t4 <- as.data.frame(Anova(Cost_trans_roots_lmer_t4))
slope_summary_costtrans_t4 <- as.data.frame(slope_summary_cost_trans_t4)
anova_results_costtrans_t4$model <- "Ecost"
slope_summary_costtrans_t4$model <- "Ecost"
############# beta plant  ###################
anova_results_beta_t4 <- as.data.frame(Anova(beta_lmer_t4))
slope_summary_beta_t4 <- as.data.frame(slope_summary_beta_t4)
anova_results_beta_t4$model <- "beta plant"
slope_summary_beta_t4$model <- "beta plant"
##################### combine ###############################################################
combined_anova_df_t4 <- rbind(anova_results_gs_t4,anova_results_chi_t4, anova_results_betaleaf_t4, 
                          anova_results_Nmass_t4, anova_results_Narea_t4,
                          anova_results_LMA_t4,anova_results_Vcmax_t4,anova_results_Jmax_t4,anova_results_nrubisco_t4,
                          anova_results_nbioe_t4, anova_results_nstr_t4,anova_results_Nlh_t4,anova_results_Chl_t4,
                          anova_results_WUE_t4, anova_results_PNUE_t4,
                          anova_results_Anet_t4,anova_results_PSA_t4, anova_results_Biomass_t4, 
                           anova_results_costN_t4, anova_results_costtrans_t4,anova_results_beta_t4)  

combined_anova_df_t4   

write.csv(combined_anova_df_t4, "../output/combined_anova_df_t4.csv", row.names = TRUE)


combined_slopes_t4 <- rbind(slope_summary_gs_t4,slope_summary_chi_t4, slope_summary_betaleaf_t4, slope_summary_Nmass_t4, 
                      slope_summary_Narea_t4,slope_summary_LMA_t4,slope_summary_Vcmax_t4,slope_summary_Jmax_t4,slope_summary_nrubisco_t4,
                      slope_summary_nbioe_t4, slope_summary_Nstr_t4,slope_summary_Nlh_t4,slope_summary_Chl_t4,
                      slope_summary_WUE_t4,slope_summary_PNUE_t4,
                      slope_summary_Anet_t4,slope_summary_PSA_t4,slope_summary_Biomass_t4, 
                      slope_summary_costN_t4, slope_summary_costtrans_t4, slope_summary_beta_t4)  
combined_slopes_t4 

write.csv(combined_slopes_t4, "../output/combined_slopes_t4.csv", row.names = TRUE)

####################### Structural equation model ################################

sunflower_beta <- psem(
  lme2 <- lm(Ncost~ Sdt_DD + Cumul_N , data = stat_t4),
  lme3 <- lm(Ecost~Sdt_DD, data = stat_t4),
  lme8 <- lm(betaplant~Ncost + Ecost, data = stat_t4),
  lme8 <- lm(betaleaf~betaplant + Sdt_DD + Cumul_N , data = stat_t4),
  lme7 <- lm(chi~betaleaf, data = stat_t4),
  lme8 <- lm(Narea~chi, data = stat_t4), 
  lme8 <- lm(vcmax_tleaf~Narea, data = stat_t4)
)
summary(sunflower_beta)
plot(sunflower_beta)



