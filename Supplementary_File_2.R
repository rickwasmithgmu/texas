#### Demographic Data Analysis ####
# Rick WA Smith
# Austin W Reynolds
# 14 August 2023

### Table of Contents ###
#1. Dependencies 
#2. User-defined functions
#3. Data wrangling
#4. Rural/Urban 
#5. Farmer status
#6. Race
#7. Sex
#8. Farmer + Race
#9. Farmer + Sex
#10. Marital status
#11. Immigrant status
#12. Multiple Testing Corrections
#13. Survival Analysis
#14. Figures



################################################################################
#1. Dependencies
################################################################################
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(survival)
library(survminer)
library(ggfortify)
library(timereg)



################################################################################
#2. User-defined functions
################################################################################
#Create custom split violin plotting function 
# from jan-glx Stack Overflow https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}



################################################################################
#3. Data Wrangling
################################################################################

#Read and Attach Data
data_path <- "/Users/austin_reynolds/Downloads/TX_Vital_Records_FINAL.csv"
data <- read.csv(data_path, fileEncoding = "UTF-8-BOM")
##Assign new variable for transforming some variables
data_transformed <- data
attach(data)

# Root Transform Age Data to Approximate Normal Dist.
## Replace Age Values in 'data_transformed' with Cube Rooted Values
data_transformed$age_yrs <- (max(data_transformed$age_yrs + 1) - data_transformed$age_yrs)^(1 / 2)

## Visualize Original vs Transformed Distributions
hist(data$age_yrs) # Original Data Dist
hist(data_transformed$age_yrs) # Cube Rooted Dist - LEFT BOUNDED



################################################################################
#4. Rural / Urban
################################################################################
#Life Span by Rural/Urban divide
## Calculate Group Means
lifespan_rural <- aggregate(data$age_yrs, list(data$rur_urb), FUN = mean)

## Mean Reduction in Life Span
lifespan_rural_reduction <- abs(lifespan_rural[1, ][2] - lifespan_rural[2, ][2])

## Non-Parametric Test for Life Span by Farmer Status
kruskal_rural <- kruskal.test(age_yrs ~ rur_urb, data = data)
kruskal_rural_pval <- kruskal_rural$p.value

## ANOVA Using Root Transformed Data
anova_rural <- aov(age_yrs ~ rur_urb, data = data_transformed)
summary(anova_rural)
#plot(anova_rural)
anova_rural <- unlist(summary(anova_rural)) # get actual p value
anova_rural_pval <- anova_rural["Pr(>F)1"] # print actual p value

## Visualize Variation in Life Span by Rural/Urban Divide
data$rur_urb <- factor(data$rur_urb, levels = c("U", "R")) #assign rural and urban classes to factors

violin_plot_rural <- ggplot(data = subset(data, !is.na(rur_urb)), aes(x = rur_urb, y = age_yrs, fill = rur_urb)) +
  geom_violin(alpha = 0.5) +
  scale_fill_manual(values = c("cadetblue", "burlywood4", "black")) +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, colour = "black") +
  labs(title = "", x = "", y = "Age in Years") + # change plot and axis titles
  scale_x_discrete(labels = c("R" = "Rural", "U" = "Urban")) +
  theme_light() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 16))



################################################################################
#5. Farmer Status
################################################################################
# Life Span by Farmer Status
## Calculate Group Means
lifespan_farmer <- aggregate(data$age_yrs, list(data$farmer_at_death), FUN = mean)

## Mean Reduction in Life Span
lifespan_farmer_reduction <- abs(lifespan_farmer[1, ][2] - lifespan_farmer[3, ][2])

## Non-Parametric Test for Life Span by Farmer Status
kruskal_farmer <- kruskal.test(age_yrs ~ farmer_at_death, data = data)
kruskal_farmer_pval <- kruskal_farmer$p.value

## ANOVA Using Root Transformed Data
anova_farmer <- aov(age_yrs ~ farmer_at_death, data = data_transformed)
summary(anova_farmer)
#plot(anova_farmer)
anova_farmer <- unlist(summary(anova_farmer)) # get actual p value
anova_farmer_pval <- anova_farmer["Pr(>F)1"] # print actual p value

##Visualize Variation in Life Span by Farmer Status
violin_plot_farmer<-ggplot(data = subset(data, !is.na(farmer_at_death)), aes(x = farmer_at_death, y = age_yrs, fill=farmer_at_death)) +
  geom_violin(alpha = 0.5) +
  scale_fill_manual(values=c("cadetblue","burlywood4", "black")) +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, colour = "black") +
  theme_light()+
  labs(title="", x="", y="")+
  scale_x_discrete(labels=c("D" = "Urban", "N" = "Rural Non-Farmer", "Y"="Farmer"))+
  theme(legend.position = "none")+
  theme(text = element_text(size = 16))



################################################################################
#6. Race
################################################################################
# Life Span by Race
## Calculate Group Means
lifespan_race <- aggregate(data$age_yrs, list(data$race), FUN = mean)

## Mean Reduction in Life Span
lifespan_race_reduction <- abs(lifespan_race[2, ][2] - lifespan_race[1, ][2])

## Non-parametric Test for Life Span by Race
kruskal_race <- kruskal.test(age_yrs ~ race, data = data)
kruskal_race_pval <- kruskal_race$p.value

## ANOVA Using Cube Root Transformed Data
anova_race <- aov(age_yrs ~ race, data = data_transformed)
summary(anova_race)
#plot(anova_race)
anova_race <- unlist(summary(anova_race)) # get actual p value
anova_race_pval <- anova_race["Pr(>F)1"] # print actual p value

## Visualize Variation in Life Span by Race
data$race <- factor(data$race, levels = c("W", "B"))
violin_plot_race <- ggplot(data = subset(data, !is.na(race)), aes(x = race, y = age_yrs, fill = race)) +
  geom_violin(alpha = 0.5) +
  scale_fill_manual(values = c("cadetblue", "burlywood4")) +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, colour = "black") +
  labs(title = "", x = "", y = "Age in Years") + # change plot and axis titles
  scale_x_discrete(labels = c("B" = "Black", "W" = "White")) +
  guides(fill = "none") +
  theme_light()



################################################################################
#7. Sex
################################################################################
# Life Span by Sex
## Calculate Group Means
lifespan_sex <- aggregate(data$age_yrs, list(data$sex), FUN = mean)

## Mean Reduction in Life Span
lifespan_sex_reduction <- abs(lifespan_sex[2, ][2] - lifespan_sex[1, ][2])

## Non-parametric Test for Life Span by Sex
kruskal_sex <- kruskal.test(age_yrs ~ sex, data = data)
kruskal_sex_pval <- kruskal_sex$p.value

## Anova Using Root Transformed Data
anova_sex <- aov(age_yrs ~ sex, data = data_transformed)
summary(anova_sex)
#plot(anova_sex)
anova_sex <- unlist(summary(anova_sex)) # get actual p values
anova_sex_pval <- anova_sex["Pr(>F)1"] # print actual p value

## Visualize Variation in Life Span by Sex
violin_plot_sex <- ggplot(data = subset(data, !is.na(sex)), aes(x = sex, y = age_yrs, fill = sex)) +
  geom_violin(alpha = 0.5) +
  scale_fill_manual(values = c("cadetblue", "burlywood4")) +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, colour = "black") +
  labs(title = "", x = "", y = "Age in Years") +
  scale_x_discrete(labels = c("F" = "Female", "M" = "Male")) +
  guides(fill = "none") +
  theme_light()



################################################################################
#8. Farmer + Race
################################################################################
# Life Span by Farmer Status and Race
## Anova using cuberoot transformed data
anova_farmer_race <- aov(age_yrs ~ farmer_at_death * race, data = data_transformed)
summary(anova_farmer_race)
#plot(anova_farmer_race)
anova_farmer_race <- unlist(summary(anova_farmer_race)) # get actual p values
anova_farmer_race_pval1 <- anova_farmer_race["Pr(>F)1"] # print actual p value for farming status
anova_farmer_race_pval2 <- anova_farmer_race["Pr(>F)2"] # print actual p value for race
anova_farmer_race_pval3 <- anova_farmer_race["Pr(>F)3"] # print actual p value for farming*race

## Visualize Life Span by Farmer Status and Race
data$farmer_at_death <- factor(data$farmer_at_death, levels = c("D", "N", "Y"))

violin_plot_farmer_race <- ggplot(data = subset(data, !is.na(farmer_at_death)), aes(x = farmer_at_death, y = age_yrs, fill = race)) +
  geom_split_violin(alpha = 0.5) +
  scale_fill_manual(values = c("cadetblue", "burlywood4")) +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, colour = "black", position = position_dodge(0.5)) +
  labs(title = "", x = "", y = "") +
  scale_x_discrete(labels = c("D" = "Urban", "N" = "Rural Non-Farmer", "Y" = "Farmer")) +
  theme_light()



################################################################################
#9. Farmer + Sex
################################################################################
#Vizualize Lifespan by Farmer Status and Sex
violin_plot_farmer_sex <- ggplot(data = subset(data, !is.na(farmer_at_death)), aes(x = farmer_at_death, y = age_yrs, fill = sex)) +
  geom_split_violin(alpha = 0.5) +
  scale_fill_manual(values = c("cadetblue", "burlywood4")) +
  stat_summary(fun = "median", geom = "crossbar", width = 0.5, colour = "black", position = position_dodge(0.5)) +
  labs(title = "", x = "", y = "") +
  scale_x_discrete(labels = c("D" = "Dallas", "N" = "Rural Non-Farmer", "Y" = "Farmer")) +
  theme_light()



################################################################################
#10. Marital Status
################################################################################
# Lifespan by Marital Status
## Group Means
lifespan_marital <- aggregate(data$age_yrs, list(data$marital_status), FUN = mean)

## Reduction in Life Span
lifespan_marital_reduction <- abs(lifespan_marital[1, ][2] - lifespan_marital[3, ][2])

## Non-parametric Test for Life Span by Marital Status
kruskal_marital <- kruskal.test(age_yrs ~ marital_status, data = data)
kruskal_marital_pval <- kruskal_marital$p.value

## Anova Using Root Transformed Data
anova_marital <- aov(age_yrs ~ marital_status, data = data_transformed)
summary(anova_marital)
#plot(anova_marital)
anova_marital <- unlist(summary(anova_marital)) # get actual p value
anova_marital_pval <- anova_marital["Pr(>F)1"] # print actual p value

## Visualize Variation in Life Span by Marital Status
data$marital_status <- factor(data$marital_status, levels = c("S", "M", "D", "W"))

violin_plot_marital <- ggplot(data = subset(data, !is.na(marital_status)), aes(x = marital_status, y = age_yrs, fill = sex)) +
  geom_split_violin(alpha = 0.5) +
  scale_fill_manual(values = c("cadetblue", "burlywood4")) +
  stat_summary(fun = "median", geom = "crossbar", width = 0.25, colour = "black", position = position_dodge(0.25)) +
  labs(title = "", x = "", y = "Age in Years") + # change plot and axis titles
  scale_x_discrete(labels = c("S" = "Single", "M" = "Married", "D" = "Divorced", "W" = "Widowed")) +
  theme_light()



################################################################################
#11. Immigrant Status
################################################################################
# Lifespan by Immigrant Status
## Calculate Group Means
lifespan_immigrant <- aggregate(data$age_yrs, list(data$imm_parents), FUN = mean)

## Mean Reduction in Life Span
lifespan_immigrant_reduction <- abs(lifespan_immigrant[2, ][2] - lifespan_immigrant[1, ][2])

## Non-parametric Test for Life Span by Immigrant Status
kruskal_immigrant <- kruskal.test(age_yrs ~ imm_parents, data = data)
kruskal_immigrant_pval <- kruskal_immigrant$p.value

## Anova Using Root Transformed Data
anova_immigrant <- aov(age_yrs ~ imm_parents, data = data_transformed)
summary(anova_immigrant)
#plot(anova_immigrant)
anova_immigrant <- unlist(summary(anova_immigrant)) # get actual p value
anova_immigrant_pval <- anova_immigrant["Pr(>F)1"] # print actual p value



################################################################################
#12. Multiple Testing Corrections
################################################################################
#Bonferroni corrections for multiple testing
kruskal_pvals <- as.vector(c(kruskal_rural_pval, kruskal_farmer_pval, kruskal_race_pval, 
                             kruskal_sex_pval, kruskal_marital_pval, kruskal_immigrant_pval)) # create a vector of pvalues for kts
anova_pvals <- as.vector(c(anova_rural_pval, anova_farmer_pval, anova_race_pval, anova_farmer_race_pval1, 
                           anova_farmer_race_pval2, anova_farmer_race_pval3, anova_sex_pval, anova_marital_pval, anova_immigrant_pval)) # create a vector of pvalues for aovs
kruskal_pvals.adj <- p.adjust(kruskal_pvals, method = "bonferroni", n = length(kruskal_pvals))
anova_pvals.adj <- p.adjust(anova_pvals, method = "bonferroni", n = length(anova_pvals))



################################################################################
#13. Survival Analysis
################################################################################
# we start by asking if being a tenant farmer reduces one's probability of survival
# relative to rural non-farmers and urban residents
# to do this we plot a Kaplan-Meier survival curve
survivor.fit.farming <- survfit(Surv(age_yrs) ~ farmer_at_death, data = data)
survivor.fit.farmingplus <- survfit(Surv(age_yrs) ~ farmer_at_death + sex + race, data = data)
print(survivor.fit.farming)
survival_plot <- ggsurvplot(survivor.fit.farming, 
                  data = data,
                  conf.int = TRUE,
                  pval = TRUE,
                  font.x = c(),
                  font.y = c(),
                  font.legend = c(),
                  xlab = "Age in Years",
                  ylab = "Survival Probability",
                  legend.labs = c("Urban", "Rural Non Farmer", "Farmer"),
                  legend.position = "right",
                  palette = c("cadetblue", "burlywood4", "black")
)
survival_plot <- survival_plot$plot # extract pathworkable plot from ggsurvplot object

# We find that tenant farmers have significantly shorter survival times compared 
# to urban residents and rural non-farmers (median survival time = 64 years, 
# versus 70 and 69 years respectively; p<0.0001)


# We next used a Cox proportional hazards model to assess the association between 
# age at death and farmer status, adjusted for both sex (Male, Female) and race (White, Black)
farmer_cox <- coxph(Surv(age_yrs) ~ farmer_at_death + sex + race, data = data)
summary(farmer_cox)

# We find that being a tenant farmers has a significant positive association with 
# earlier age at death (AHR = 1.29; 95% CI 1.16-1.43 ; p = 2.8e-6). 
# We also find that being a man is positively associated with earlier age at death
# (AHR = 1.41; 95% CI 1.29-1.54; p = 6.41e-15), while being white is negatively 
# associated with earlier age at death (AHR = 0.51; 95% CI 0.45-0.57; p < 2e-16)


# We assess the proportional hazards assumption using Schoenfeld residuals
residuals.zph <- cox.zph(farmer_cox)
print(residuals.zph)
par(mfrow=c(1,3))
plot(residuals.zph)
# and find that farming status and race violate this assumption, suggesting that 
# their effects on age at death do vary over time


# To evaluate the time-dependent effects of each variable on survival we fit a 
# Aalen's additive regression
aa.fit <- aareg(Surv(age_yrs) ~ farmer_at_death + sex + race, data = data)
summary(aa.fit)
#The plot shows that over time, the effect of being a farmer on age at death increases, 
#as does being a man. The effect of being white on age at death decreases through time.

#Assess goodness of fit using predicted survival curves
#start by bootstrapping the aalen model to get confidence intervals
aalen.fit <- aalen(Surv(age_yrs) ~ farmer_at_death + sex + race, data = data,
                   n.sim=500, resample.iid=1)
summary(aalen.fit)

#next create a model matrix containing each possible combination of farmer, sex, race variables
newdata <- as.data.frame(model.matrix(age_yrs ~ farmer_at_death + sex + race, data = data)) %>% distinct()
newdata$farmer_at_death<-c("Y","N","N","Y","N","Y","Y","N","D","D","D","D")
newdata$sex<-c("M","M","F","F","F","F","M","M","M","M","F","F")
newdata$race<-c("W","W","B","W","W","B","B","B","B","W","W","B")

#finally predict new data from bootstrapped model
aalen_fit_predictions<-predict.aalen(aalen.fit,newdata)
#plot predicted survival curve within actual curve + CI for each possible combination


################################################################################
#14. FIGURES
################################################################################

#Figure 1: Map of Study Area in Texas
texas_counties <- map_data("county", region = "texas")

vec1 <- c("ellis","navarro","hill")
vec2 <- c("dallas")

dt <- texas_counties %>% mutate(selected = ifelse(subregion %in% vec1, "1", "0"))
dt$selected[which(dt$subregion %in% vec2)] <- "2"

scale_fill_manual(breaks = c("2", "1", "0"),)

p <- ggplot(data = dt,
            mapping = aes(x = long, y = lat,
                          group = group, fill = selected))
p + geom_polygon(color = "gray90", size = 0.1) +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
  guides(fill = FALSE)+
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank()) +
  scale_fill_manual(breaks = c("2", "1", "0"),values=c("black", "cadetblue", "darkgray"))



#Figure 2: Age at Death across Urban and Rural Residence.
violin_plot_rural+violin_plot_farmer+survival_plot+plot_annotation(tag_levels = 'I')



#Figure3: Age at Death across I) Sex, II) Sex and Proximity to Farm Labor, and III) Sex and Marital Status
((violin_plot_sex+violin_plot_farmer_sex)/violin_plot_marital)+plot_annotation(tag_levels = 'I')



#Figure 4: Age at Death across I) Racial categories, and II) Race and Proximity to Farm Labor
violin_plot_race+violin_plot_farmer_race+plot_annotation(tag_levels = 'I')



#Figure S1: Cumulative hazards with their 95% confidence intervals estimated with Aalen’s additive hazards model for farmer status, sex, and race.
autoplot(aa.fit) + geom_hline(yintercept = 0)


#Figure S2: Comparison of empirical (black solid lines) K-M survival curves with those estimated from the Aalen’s model (red solid lines) and 95% confidence interval (red dash lines).
par(mfrow=c(4,3))
#DFB - Dallas, Female, Black
plot(aalen_fit_predictions,col=1, se = 0, uniform = 0, specific.comps =12)
lines(survivor.fit.farmingplus[1], col = 2)
title(main = "Urban")
#NFB - rural non-farmer, female, black and so on...
plot(aalen_fit_predictions,col=1, se = 0, uniform = 0, specific.comps =3)
lines(survivor.fit.farmingplus[5], col = 2)
title(main = "Rural Non-Farmer")
#YFB
plot(aalen_fit_predictions,col=1, se = 0, uniform = 0, specific.comps =6)
lines(survivor.fit.farmingplus[9], col = 2)
title(main = "Rural Farmer")
#DMB
plot(aalen_fit_predictions,col=1, se = 0, uniform = 0, specific.comps =9)
lines(survivor.fit.farmingplus[3], col = 2)
#NMB
plot(aalen_fit_predictions,col=1, se = 0, uniform = 0, specific.comps =8)
lines(survivor.fit.farmingplus[7], col = 2)
#YMB
plot(aalen_fit_predictions,col=1, se = 0, uniform = 0, specific.comps =7)
lines(survivor.fit.farmingplus[11], col = 2)
#DFW
plot(aalen_fit_predictions,col=1, se = 0, uniform = 0, specific.comps =11)
lines(survivor.fit.farmingplus[2], col = 2)
#NFW
plot(aalen_fit_predictions,col=1, se = 0, uniform = 0, specific.comps =5)
lines(survivor.fit.farmingplus[6], col = 2)
#YFW
plot(aalen_fit_predictions,col=1, se = 0, uniform = 0, specific.comps =4)
lines(survivor.fit.farmingplus[10], col = 2)
#DMW
plot(aalen_fit_predictions,col=1, se = 0, uniform = 0, specific.comps =10)
lines(survivor.fit.farmingplus[4], col = 2)
#NMW
plot(aalen_fit_predictions,col=1, se = 0, uniform = 0, specific.comps =2)
lines(survivor.fit.farmingplus[8], col = 2)
#YMW
plot(aalen_fit_predictions,col=1, se = 0, uniform = 0, specific.comps =1)
lines(survivor.fit.farmingplus[12], col = 2)

