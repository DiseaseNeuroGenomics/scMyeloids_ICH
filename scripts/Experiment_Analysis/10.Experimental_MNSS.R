setwd("/Users/dimitrioskyriakis/Downloads/")

# Load necessary libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(emmeans)

# Load dataset
data <- read.csv("Modified_ICH_Microglia_EXP.csv")

##########
### ALL
##########


# Reshape data into long format for analysis
long_data <- data %>%filter(Condition %in% c('Control','RA'))%>%
  select(Animal.Code, Condition,Dose, Gender, MNSS_BSL, MNSS_2h, MNSS_Day1, MNSS_Day3, MNSS_Day7) %>%
  pivot_longer(cols = starts_with("MNSS"),
               names_to = "Timepoint",
               values_to = "MNSS_Score")

# Convert variables to factors
long_data$Timepoint <- factor(long_data$Timepoint,
                              levels = c("MNSS_BSL", "MNSS_2h", "MNSS_Day1", "MNSS_Day3", "MNSS_Day7"),
                              ordered = TRUE)
long_data$Gender <- as.factor(long_data$Gender)
long_data$Condition <- as.factor(long_data$Condition)


# Fit a linear mixed-effects model including Gender as a covariate
gaus_model <- lmer(MNSS_Score ~ Timepoint * Condition + Dose + (1 | Animal.Code), data = long_data)
gamma_model <- glmer(MNSS_Score ~ Timepoint * Condition + Dose + (1|Animal.Code),
                     family = poisson(link = "log"), data = long_data)
BIC(gaus_model,gamma_model)
summary(gaus_model)
summary(gamma_model)

# ANOVA table
anova(gamma_model)

# Post-hoc comparisons
posthoc <- emmeans(gamma_model, pairwise ~  Condition | Timepoint, adjust = "bonferroni")
summary(posthoc)
posthoc <- emmeans(gaus_model, pairwise ~  Condition | Timepoint, adjust = "bonferroni")
summary(posthoc)






# Reshape data into long format for analysis
long_data <- data %>%filter(Condition %in% c('Control','Everolimus'))%>%
  select(Animal.Code, Condition,Dose, Gender, MNSS_BSL, MNSS_2h, MNSS_Day1, MNSS_Day3, MNSS_Day7) %>%
  pivot_longer(cols = starts_with("MNSS"),
               names_to = "Timepoint",
               values_to = "MNSS_Score")

# Convert variables to factors
long_data$Timepoint <- factor(long_data$Timepoint,
                              levels = c("MNSS_BSL", "MNSS_2h", "MNSS_Day1", "MNSS_Day3", "MNSS_Day7"),
                              ordered = TRUE)
long_data$Gender <- as.factor(long_data$Gender)
long_data$Condition <- as.factor(long_data$Condition)


# Fit a linear mixed-effects model including Gender as a covariate
gaus_model <- lmer(MNSS_Score ~ Timepoint * Condition + Dose + (1 | Animal.Code), data = long_data)
gamma_model <- glmer(MNSS_Score ~ Timepoint * Condition + Dose + (1|Animal.Code),
                     family = poisson(link = "log"), data = long_data)
BIC(gaus_model,gamma_model)
summary(gaus_model)
summary(gamma_model)

# ANOVA table
anova(gamma_model)

# Post-hoc comparisons
posthoc <- emmeans(gamma_model, pairwise ~  Condition | Timepoint, adjust = "bonferroni")
summary(posthoc)
posthoc <- emmeans(gaus_model, pairwise ~  Condition | Timepoint, adjust = "bonferroni")
summary(posthoc)


























# Fit a linear mixed-effects model including Gender as a covariate
model <- lmer(MNSS_Score ~ Timepoint * Condition + Gender + (1 | Animal.Code), data = long_data)

# ANOVA table
anova(model)

# Post-hoc comparisons
posthoc <- emmeans(model, pairwise ~ Timepoint | Condition, adjust = "bonferroni")
summary(posthoc)

# Visualization of MNSS score trends by Condition and Gender
ggplot(long_data, aes(x = Timepoint, y = MNSS_Score, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  facet_wrap(~ Gender) +
  labs(title = "MNSS Scores Over Time by Condition and Gender",
       x = "Timepoint",
       y = "MNSS Score") +
  theme_minimal()


### Everolimus vs Controls

# Filter only Everolimus and Control groups
data_filtered <- data %>% filter(Condition %in% c("Everolimus", "Control"))

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

long_data
# Fit a linear mixed-effects model for Everolimus vs Control
model <- lmer(MNSS_Score ~ Timepoint * Condition + Gender + (1 | Animal.Code), data = long_data)

# ANOVA table
anova(model)

# Post-hoc comparisons (Everolimus vs Control across timepoints)
posthoc <- emmeans(model, pairwise ~ Timepoint | Condition, adjust = "bonferroni")
summary(posthoc)

# Visualization of MNSS score trends for Everolimus vs Control
ggplot(long_data, aes(x = Timepoint, y = MNSS_Score, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  facet_wrap(~ Gender) +
  labs(title = "MNSS Scores Over Time: Everolimus vs Control",
       x = "Timepoint",
       y = "MNSS Score") +
  theme_minimal()

### RA
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


#### RA but dose as covariate
# Filter only RA and Control groups
data_filtered <- data %>% filter(Condition %in% c("RA", "Control"))

# Reshape data into long format for analysis
long_data <- data_filtered %>%
  select(Animal.Code, Condition, Dose, MNSS_BSL, MNSS_2h, MNSS_Day1, MNSS_Day3, MNSS_Day7) %>%
  pivot_longer(cols = starts_with("MNSS"),
               names_to = "Timepoint",
               values_to = "MNSS_Score")

# Convert variables to factors
long_data$Timepoint <- factor(long_data$Timepoint,
                              levels = c("MNSS_BSL", "MNSS_2h", "MNSS_Day1", "MNSS_Day3", "MNSS_Day7"),
                              ordered = TRUE)
long_data$Condition <- as.factor(long_data$Condition)

# Fit a linear mixed-effects model for RA vs Control with Dose as a covariate
model <- lmer(MNSS_Score ~ Timepoint * Condition + Dose + (1 | Animal.Code), data = long_data)

# ANOVA table
anova(model)

# Post-hoc comparisons (RA vs Control across timepoints)
posthoc <- emmeans(model, pairwise ~ Timepoint | Condition, adjust = "bonferroni")
summary(posthoc)

# Visualization of MNSS score trends for RA vs Control with Dose as a covariate
ggplot(long_data, aes(x = Timepoint, y = MNSS_Score, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  labs(title = "MNSS Scores Over Time: RA vs Control (Adjusted for Dose)",
       x = "Timepoint",
       y = "MNSS Score") +
  theme_minimal()

### Everolimus adjusted for dose
# Filter only Everolimus and Control groups
data_filtered <- data %>% filter(Condition %in% c("Everolimus", "Control"))

# Reshape data into long format for analysis
long_data <- data_filtered %>%
  select(Animal.Code, Condition, Dose, MNSS_BSL, MNSS_2h, MNSS_Day1, MNSS_Day3, MNSS_Day7) %>%
  pivot_longer(cols = starts_with("MNSS"),
               names_to = "Timepoint",
               values_to = "MNSS_Score")

# Convert variables to factors
long_data$Timepoint <- factor(long_data$Timepoint,
                              levels = c("MNSS_BSL", "MNSS_2h", "MNSS_Day1", "MNSS_Day3", "MNSS_Day7"),
                              ordered = TRUE)
long_data$Condition <- as.factor(long_data$Condition)

# Fit a linear mixed-effects model for Everolimus vs Control with Dose as a covariate
model <- lmer(MNSS_Score ~ Timepoint * Condition + Dose + (1 | Animal.Code), data = long_data)

# ANOVA table
anova(model)

# Post-hoc comparisons (Everolimus vs Control across timepoints)
posthoc <- emmeans(model, pairwise ~ Timepoint | Condition, adjust = "bonferroni")
summary(posthoc)

# Visualization of MNSS score trends for Everolimus vs Control with Dose as a covariate
ggplot(long_data, aes(x = Timepoint, y = MNSS_Score, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  labs(title = "MNSS Scores Over Time: Everolimus vs Control (Adjusted for Dose)",
       x = "Timepoint",
       y = "MNSS Score") +
  theme_minimal()

### Everolimus, RA, vs control
# Filter only Everolimus, RA, and Control groups
data_filtered <- data %>% filter(Condition %in% c("Everolimus", "RA", "Control"))

# Reshape data into long format for analysis
long_data <- data_filtered %>%
  select(Animal.Code, Condition, Dose, MNSS_BSL, MNSS_2h, MNSS_Day1, MNSS_Day3, MNSS_Day7) %>%
  pivot_longer(cols = starts_with("MNSS"),
               names_to = "Timepoint",
               values_to = "MNSS_Score")

# Convert variables to factors
long_data$Timepoint <- factor(long_data$Timepoint,
                              levels = c("MNSS_BSL", "MNSS_2h", "MNSS_Day1", "MNSS_Day3", "MNSS_Day7"),
                              ordered = TRUE)
long_data$Condition <- as.factor(long_data$Condition)

# Fit a linear mixed-effects model for Everolimus & RA vs Control with Dose as a covariate
model <- lmer(MNSS_Score ~ Timepoint * Condition + Dose + (1 | Animal.Code), data = long_data)

# ANOVA table
anova(model)

# Post-hoc comparisons (Everolimus & RA vs Control across timepoints)
posthoc <- emmeans(model, pairwise ~ Timepoint | Condition, adjust = "bonferroni")
summary(posthoc)

# Visualization of MNSS score trends for Everolimus, RA vs Control (adjusted for Dose)
ggplot(long_data, aes(x = Timepoint, y = MNSS_Score, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  labs(title = "MNSS Scores Over Time: Everolimus & RA vs Control (Adjusted for Dose)",
       x = "Timepoint",
       y = "MNSS Score") +
  theme_minimal()

### combine treatment everolimus vs control
# Recode Everolimus & RA as "Treated" and keep Control as is
data <- data %>%
  mutate(Treatment_Group = case_when(
    Condition %in% c("Everolimus", "RA") ~ "Treated",
    Condition == "Control" ~ "Control"
  ))

# Filter only Treated and Control groups
data_filtered <- data %>% filter(Treatment_Group %in% c("Treated", "Control"))

# Reshape data into long format for analysis
long_data <- data_filtered %>%
  select(Animal.Code, Treatment_Group, Dose, MNSS_BSL, MNSS_2h, MNSS_Day1, MNSS_Day3, MNSS_Day7) %>%
  pivot_longer(cols = starts_with("MNSS"),
               names_to = "Timepoint",
               values_to = "MNSS_Score")

# Convert variables to factors
long_data$Timepoint <- factor(long_data$Timepoint,
                              levels = c("MNSS_BSL", "MNSS_2h", "MNSS_Day1", "MNSS_Day3", "MNSS_Day7"),
                              ordered = TRUE)
long_data$Treatment_Group <- as.factor(long_data$Treatment_Group)

# Fit a linear mixed-effects model for Treated vs Control with Dose as a covariate
model <- lmer(MNSS_Score ~ Timepoint * Treatment_Group + Dose + (1 | Animal.Code), data = long_data)

# ANOVA table
anova(model)

# Post-hoc comparisons (Treated vs Control across timepoints)
posthoc <- emmeans(model, pairwise ~ Timepoint | Treatment_Group, adjust = "bonferroni")
summary(posthoc)

# Visualization of MNSS score trends for Treated vs Control (adjusted for Dose)
ggplot(long_data, aes(x = Timepoint, y = MNSS_Score, color = Treatment_Group, group = Treatment_Group)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  labs(title = "MNSS Scores Over Time: Treated vs Control (Adjusted for Dose)",
       x = "Timepoint",
       y = "MNSS Score") +
  theme_minimal()
