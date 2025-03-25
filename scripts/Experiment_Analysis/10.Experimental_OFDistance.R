setwd("/Users/dimitrioskyriakis/Downloads/")


# ============== RUN GAUSSIAN ===============
# ============== RUN GAUSSIAN ===============
# ============== RUN GAUSSIAN ===============

# Load necessary libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(emmeans)

# Load dataset
data <- read.csv("Modified_ICH_Microglia_EXP.csv")
data$OFDistance_Day3 <- as.numeric(data$OFDistance_Day3)
data$OFDistance_Day7 <- as.numeric(data$OFDistance_Day7)
data$X <- NULL
data$X.1 <- NULL
data$X.2 <- NULL
data$X.3 <- NULL
##########
### ALL
##########

colnames(data)
# Reshape data into long format for analysis
long_data <- data %>% filter(Condition %in% c('Control','RA'))%>%
  select(Animal.Code, Condition, Gender,Dose,contains('OFDistance_')) %>%
  pivot_longer(cols = starts_with("OFDistance"),
               names_to = "Timepoint",
               values_to = "OFDistance_Score")

# Convert variables to factors
long_data$Timepoint <- factor(long_data$Timepoint,
                              levels = c("OFDistance_BSL", "OFDistance_Day3", "OFDistance_Day7"),
                              ordered = TRUE)
long_data$Gender <- as.factor(long_data$Gender)
long_data$Condition <- as.factor(long_data$Condition)

hist(long_data$OFDistance_Score)

# Fit a linear mixed-effects model including Gender as a covariate
gaus_model <- lmer(OFDistance_Score ~ Timepoint + Condition + Dose + (1 | Animal.Code), data = long_data)
summary(gaus_model)

gaus_model <- lmer(OFDistance_Score ~ Timepoint * Condition + Dose + (1 | Animal.Code), data = long_data)
summary(gaus_model)

# ANOVA table
anova(gaus_model)

# Post-hoc comparisons
posthoc <- emmeans(gaus_model, pairwise ~  Condition | Timepoint, adjust = "bonferroni")
summary(posthoc)

# Visualization of MNSS score trends by Condition and Gender
ggplot(long_data, aes(x = Timepoint, y = OFDistance_Score, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  facet_wrap(~ Gender) +
  labs(title = "OFDistance_ Score Over Time by Condition and Gender",
       x = "Timepoint",
       y = "OFDistance_Score") +
  theme_minimal()



# ========================= Everolimus vs Controls ================================
# ========================= Everolimus vs Controls ================================
# ========================= Everolimus vs Controls ================================
### Everolimus vs Controls

# Filter only Everolimus and Control groups
data_Everolimus <- long_data %>% filter(Condition %in% c("Everolimus", "Control"))
data_Everolimus$Gender <- as.factor(data_Everolimus$Gender)
data_Everolimus$Condition <- as.factor(data_Everolimus$Condition)
data_Everolimus$Timepoint <- factor(data_Everolimus$Timepoint,
                              levels = c("OFDistance_BSL", "OFDistance_Day3", "OFDistance_Day7"),
                              ordered = TRUE)

data_Everolimus$Condition
# Fit a linear mixed-effects model for Everolimus vs Control
model <- lmer(OFDistance_Score ~ Timepoint * Condition + Gender + (1 | Animal.Code), data = data_Everolimus)

# ANOVA table
anova(model)

# Post-hoc comparisons (Everolimus vs Control across timepoints)
posthoc <- emmeans(model, pairwise ~ Timepoint | Condition, adjust = "bonferroni")
summary(posthoc)

# Visualization of MNSS score trends for Everolimus vs Control
ggplot(data_Everolimus, aes(x = Timepoint, y = MNSS_Score, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  facet_wrap(~ Gender) +
  labs(title = "MNSS Scores Over Time: Everolimus vs Control",
       x = "Timepoint",
       y = "MNSS Score") +
  theme_minimal()
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# ============================= RA vs Controls ==================================
# ============================= RA vs Controls ==================================
# ============================= RA vs Controls ==================================
# Filter only RA and Control groups
data_filtered <- data %>% filter(Condition %in% c("RA", "Control"))

# Reshape data into long format for analysis
long_data <- data_filtered %>%
  select(Animal.Code, Condition, Gender, MNSS_BSL, MNSS_2h, MNSS_Day1, MNSS_Day3, MNSS_Day7) %>%
  pivot_longer(cols = starts_with("MNSS"),
               names_to = "Timepoint",
               values_to = "MNSS_Score")

# Convert variables to factors
long_data$Timepoint <- factor(long_data$Timepoint,
                              levels = c("MNSS_BSL", "MNSS_2h", "MNSS_Day1", "MNSS_Day3", "MNSS_Day7"),
                              ordered = TRUE)
long_data$Gender <- as.factor(long_data$Gender)
long_data$Condition <- as.factor(long_data$Condition)

# Fit a linear mixed-effects model for RA vs Control
model <- lmer(MNSS_Score ~ Timepoint * Condition + Gender + (1 | Animal.Code), data = long_data)

# ANOVA table
anova(model)

# Post-hoc comparisons (RA vs Control across timepoints)
posthoc <- emmeans(model, pairwise ~ Timepoint | Condition, adjust = "bonferroni")
summary(posthoc)

# Visualization of MNSS score trends for RA vs Control
ggplot(long_data, aes(x = Timepoint, y = MNSS_Score, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  facet_wrap(~ Gender) +
  labs(title = "MNSS Scores Over Time: RA vs Control",
       x = "Timepoint",
       y = "MNSS Score") +
  theme_minimal()
