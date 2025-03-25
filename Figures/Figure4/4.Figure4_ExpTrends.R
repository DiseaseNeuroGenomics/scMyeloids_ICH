# Load necessary libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(emmeans)
library(ggpubr)

library(ggplot2)
library(dplyr)
# Load dataset
data <- read.csv("/Users/dimitrioskyriakis/Documents/Projects/scMyeloids_ICH/Paper_Tables/7.Mice_Experiments_Results.csv")
data$OFVelocity_Day3 <- as.numeric(data$OFVelocity_Day3)
data$OFVelocity_Day7 <- as.numeric(data$OFVelocity_Day7)
data$OFDistance_Day3 <- as.numeric(data$OFDistance_Day3)
data$OFDistance_Day7 <- as.numeric(data$OFDistance_Day7)
data$OFInZoneOZ_Day3 <- as.numeric(data$OFInZoneOZ_Day3)
data$OFInZoneOZ_Day7 <- as.numeric(data$OFInZoneOZ_Day7)
data$OFInZoneCorners_Day3 <- as.numeric(data$OFInZoneCorners_Day3)
data$OFInZoneCorners_Day7 <- as.numeric(data$OFInZoneCorners_Day7)
data$OFInZoneBoarders_Day3 <- as.numeric(data$OFInZoneBoarders_Day3)
data$OFInZoneBoarders_Day7 <- as.numeric(data$OFInZoneBoarders_Day7)
data$RotarodLatencyMeans_Day3 <- as.numeric(data$RotarodLatencyMeans_Day3)
data$RotarodLatencyMeans_Day7 <- as.numeric(data$RotarodLatencyMeans_Day7)
data$RotarodLatencyMedians_Day3 <- as.numeric(data$RotarodLatencyMedians_Day3)
data$RotarodLatencyMedians_Day7 <- as.numeric(data$RotarodLatencyMedians_Day7)
data$RotarodRPMMeans_Day3 <- as.numeric(data$RotarodRPMMeans_Day3)
data$RotarodRPMMeans_Day7 <- as.numeric(data$RotarodRPMMeans_Day7)
data$RotarodRPMMedians_Day3 <- as.numeric(data$RotarodRPMMedians_Day3)
data$RotarodRPMMedians_Day7 <- as.numeric(data$RotarodRPMMedians_Day7)
data$X <- NULL
data$X.1 <- NULL
data$X.2 <- NULL
data$X.3 <- NULL
##########
### ALL
##########


colnames(data)



# ======================================= RotarodRPMMedians ====================================================
# ======================================= RotarodRPMMedians ====================================================
# ======================================= RotarodRPMMedians ====================================================
# ======================================= RotarodRPMMedians ====================================================
# Reshape data into long format for analysis
df <- data %>% filter(Condition %in% c('Control','RA','Everolimus'))%>%
  select(Animal.Code, Condition, Gender,Dose,
         RotarodRPMMedians_BSL, RotarodRPMMedians_Day3, RotarodRPMMedians_Day7) %>%
  pivot_longer(cols = starts_with("RotarodRPMMedians"),
               names_to = "Timepoint",
               values_to = "RotarodRPMMedians_Score")
df$Timepoint[df$Timepoint=='RotarodRPMMedians_BSL'] <- 'BSL'
df$Timepoint[df$Timepoint=='RotarodRPMMedians_Day3'] <- 'Day3'
df$Timepoint[df$Timepoint=='RotarodRPMMedians_Day7'] <- 'Day7'
df$Timepoint <- factor(df$Timepoint, levels = c("BSL", "Day3", "Day7"),ordered = TRUE)
df$Timepoint_Numeric <- as.numeric(df$Timepoint)
df$Gender <- as.factor(df$Gender)
df$Condition <- as.factor(df$Condition)
df$Dose_scaled <- scale(df$Dose)


hist(df$RotarodRPMMedians_Score[df$Condition =='Control'],breaks=20)

# ====================== CHECK DISTRIBUTION MODEL ==============================
gaus_model <- lmer(RotarodRPMMedians_Score ~ Timepoint   + Dose_scaled + (1 | Animal.Code),
                   data = df%>% filter(Condition %in% c('Control')))
gamma_model <- glmer(RotarodRPMMedians_Score ~ Timepoint   + Dose_scaled + (1 | Animal.Code),
                     data = df%>% filter(Condition %in% c('Control')),family = Gamma(link = 'log'))
poisson_model <- glmer(RotarodRPMMedians_Score ~ Timepoint   + Dose_scaled + (1 | Animal.Code),
                       data = df%>% filter(Condition %in% c('Control')),family = poisson(link = 'log'))
BIC(poisson_model,gaus_model,gamma_model)
# ------------------------------------------------------------------------------
# ANOVA table
# Perform ANOVA
model <- glmer(RotarodRPMMedians_Score ~ Timepoint * Condition  + Dose_scaled + (1 | Animal.Code),
                       data = df%>% filter(Condition %in% c('Control','RA')),family = Gamma(link = 'log'))
anova_results <- anova(model)
summary(model)

emmeans_model <- emmeans(model, pairwise ~ Condition | Timepoint, adjust = "tuckey")
summary(emmeans_model)


# ANOVA table
# Perform ANOVA
model <- glmer(RotarodRPMMedians_Score ~ Timepoint * Condition  + Dose_scaled + (1 | Animal.Code),
               data = df%>% filter(Condition %in% c('Control','Everolimus')),family = Gamma(link = 'log'))
anova_results <- anova(model)
summary(model)
emmeans_model <- emmeans(model, pairwise ~ Condition | Timepoint, adjust = "tuckey")
summary(emmeans_model)


model <- glmer(RotarodRPMMedians_Score ~ Timepoint_Numeric * Condition  + Dose_scaled + (1 | Animal.Code),
               data = df%>% filter(Condition %in% c('Control','RA','Everolimus')),family = Gamma(link = 'log'))
emtrends_model <- emtrends(model, pairwise ~ Condition, var = "Timepoint_Numeric", adjust = "bonferroni")
summary(emtrends_model)


p <- ggplot(df, aes(x = Timepoint, y = RotarodRPMMedians_Score, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(title = "Rotarod Over Time by Condition", x = "Timepoint", y = "Rotarod") +
  scale_color_manual(values = c("Control" = "black", "Everolimus" = "#CC79A7", "RA" = "#E69F00")) +
  theme_bw() +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 6)
  )

p

ggsave("scMyeloids_ICH/Paper_Figures/Figure4/RotarodRPMMedians_All_Lines.pdf", plot = p, width = 6, height = 4, device = cairo_pdf)

# ======================================= RotarodLatencyMedians ====================================================
# ======================================= RotarodLatencyMedians ====================================================
# ======================================= RotarodLatencyMedians ====================================================
# ======================================= RotarodLatencyMedians ====================================================
# Reshape data into long format for analysis
data$RotarodLatencyMedians
df <- data %>% filter(Condition %in% c('Control','RA','Everolimus'))%>%
  select(Animal.Code, Condition, Gender,Dose,contains('RotarodLatencyMedians')) %>%
  pivot_longer(cols = starts_with("RotarodLatencyMedians"),
               names_to = "Timepoint",
               values_to = "RotarodLatencyMedians_Score")

df$Timepoint[df$Timepoint=='RotarodLatencyMedians_BSL'] <- 'BSL'
df$Timepoint[df$Timepoint=='RotarodLatencyMedians_Day3'] <- 'Day3'
df$Timepoint[df$Timepoint=='RotarodLatencyMedians_Day7'] <- 'Day7'
df$Timepoint <- factor(df$Timepoint, levels = c("BSL", "Day3", "Day7"),ordered = TRUE)
df$Gender <- as.factor(df$Gender)
df$Condition <- as.factor(df$Condition)
df$Dose_scaled <- scale(df$Dose)
df$Timepoint_Numeric <- as.numeric(df$Timepoint)

hist(df$RotarodLatencyMedians_Score[df$Condition =='Control'])

# ====================== CHECK DISTRIBUTION MODEL ==============================
gaus_model <- lmer(RotarodLatencyMedians_Score ~ Timepoint   + Dose_scaled + (1 | Animal.Code),
                   data = df%>% filter(Condition %in% c('Control')))
gamma_model <- glmer(RotarodLatencyMedians_Score ~ Timepoint   + Dose_scaled + (1 | Animal.Code),
                     data = df%>% filter(Condition %in% c('Control')),family = Gamma(link = 'log'))
poisson_model <- glmer(RotarodLatencyMedians_Score ~ Timepoint   + Dose_scaled + (1 | Animal.Code),
                       data = df%>% filter(Condition %in% c('Control')),family = poisson(link = 'log'))
BIC(poisson_model,gaus_model,gamma_model)
# ------------------------------------------------------------------------------

# ANOVA table
# Perform ANOVA
model <- lmer(RotarodLatencyMedians_Score ~ Timepoint * Condition  + Dose + (1 | Animal.Code),
              data = df%>% filter(Condition %in% c('Control','RA')))
anova_results <- anova(model)
summary(model)

emmeans_model <- emmeans(model, pairwise ~ Condition | Timepoint, adjust = "tuckey")
summary(emmeans_model)

# ANOVA table
# Perform ANOVA
model <- lmer(RotarodLatencyMedians_Score ~ Timepoint * Condition  + Dose + (1 | Animal.Code),
              data = df%>% filter(Condition %in% c('Control','Everolimus')))
anova_results <- anova(model)
summary(model)


emmeans_model <- emmeans(model, pairwise ~ Condition | Timepoint, adjust = "tuckey")
summary(emmeans_model)



model <- lmer(RotarodLatencyMedians_Score ~ Timepoint_Numeric * Condition  + Dose_scaled + (1 | Animal.Code),
               data = df%>% filter(Condition %in% c('Control','RA','Everolimus')))
emtrends_model <- emtrends(model, pairwise ~ Condition, var = "Timepoint_Numeric", adjust = "bonferroni")
summary(emtrends_model)



p <- ggplot(df, aes(x = Timepoint, y = RotarodLatencyMedians_Score, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(title = "RotarodLatencyMedians_Score", x = "Timepoint", y = "RotarodLatencyMedians") +
  scale_color_manual(values = c("Control" = "black", "Everolimus" = "#CC79A7", "RA" = "#E69F00")) +
  theme_bw() +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 6)
  )


ggsave("scMyeloids_ICH/Paper_Figures/Figure4/RotarodLatencyMedians_All_Lines.pdf", plot = p, width = 6, height = 4, device = cairo_pdf)


# ======================================= MNSS ====================================================
# ======================================= MNSS ====================================================
# ======================================= MNSS ====================================================
# ======================================= MNSS ====================================================
df <- data %>% filter(Condition %in% c('Control','RA','Everolimus'))%>%
  select(Animal.Code, Condition, Gender,Dose,contains('MNSS')) %>%
  pivot_longer(cols = starts_with("MNSS"),
               names_to = "Timepoint",
               values_to = "MNSS_Score")
df$Timepoint[df$Timepoint=='MNSS_BSL'] <- 'BSL'
df$Timepoint[df$Timepoint=='MNSS_2h'] <- '2h'
df$Timepoint[df$Timepoint=='MNSS_Day1'] <- 'Day1'
df$Timepoint[df$Timepoint=='MNSS_Day3'] <- 'Day3'
df$Timepoint[df$Timepoint=='MNSS_Day7'] <- 'Day7'
df$Timepoint <- factor(df$Timepoint,levels = c("BSL",'2h',"Day1", "Day3", "Day7"),ordered = TRUE)
df$Gender <- as.factor(df$Gender)
df$Condition <- as.factor(df$Condition)
df$Dose_scaled <- scale(df$Dose)
df$Timepoint_Numeric <- as.numeric(df$Timepoint)

hist(df$MNSS_Score[df$Condition =='Control'])

df$MNSS_Score
# ====================== CHECK DISTRIBUTION MODEL ==============================
gaus_model <- lmer(MNSS_Score+1 ~ Timepoint   + Dose_scaled + (1 | Animal.Code),
                   data = df%>% filter(Condition %in% c('Control')))
gamma_model <- glmer(MNSS_Score+1 ~ Timepoint   + Dose_scaled + (1 | Animal.Code),
                     data = df%>% filter(Condition %in% c('Control')),family = Gamma(link = 'log'))
poisson_model <- glmer(MNSS_Score+1 ~ Timepoint   + Dose_scaled + (1 | Animal.Code),
                       data = df%>% filter(Condition %in% c('Control')),family = poisson(link = 'log'))
BIC(poisson_model,gaus_model,gamma_model)
# ------------------------------------------------------------------------------


# ANOVA table
model <- lmer(MNSS_Score ~ Timepoint * Condition  + Dose_scaled + (1 | Animal.Code),
              data = df%>% filter(Condition %in% c('Control','RA')))
anova_results <- anova(model)
summary(model)
emmeans_model <- emmeans(model, pairwise ~ Condition | Timepoint, adjust = "tuckey")
summary(emmeans_model)
# ANOVA table
# Perform ANOVA
model <- lmer(MNSS_Score ~ Timepoint * Condition  + Dose_scaled + (1 | Animal.Code),
              data = df%>% filter(Condition %in% c('Control','Everolimus')))
anova_results <- anova(model)
summary(model)
emmeans_model <- emmeans(model, pairwise ~ Condition | Timepoint, adjust = "tuckey")
summary(emmeans_model)


model <- glmer(MNSS_Score+1 ~ Timepoint_Numeric * Condition  + Dose_scaled + (1 | Animal.Code),
              data = df%>% filter(Condition %in% c('Control','RA','Everolimus')),family = Gamma(link = 'log'))
emtrends_model <- emtrends(model, pairwise ~ Condition, var = "Timepoint_Numeric", adjust = "bonferroni")
summary(emtrends_model)


p <- ggplot(df, aes(x = Timepoint, y = MNSS_Score, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(title = "MNSS_Score", x = "Timepoint", y = "MNSS") +
  scale_color_manual(values = c("Control" = "black", "Everolimus" = "#CC79A7", "RA" = "#E69F00")) +
  theme_bw() +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 6)
  )


ggsave("scMyeloids_ICH/Paper_Figures/Figure4/MNSS_Score_All_Lines.pdf", plot = p, width = 6, height = 4, device = cairo_pdf)


