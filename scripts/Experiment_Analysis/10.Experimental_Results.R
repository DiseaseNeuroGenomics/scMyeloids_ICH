setwd("/Users/dimitrioskyriakis/Downloads/")

# Load necessary libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(emmeans)

# Load dataset
data <- read.csv("../../../../../Downloads/Modified_ICH_Microglia_Experimental_Data.csv")
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

Died_earlier <- data %>%filter(is.na(MNSS_Day7) )
table(Died_earlier$Condition,Died_earlier$Gender)

Without_exclude <- data %>%filter(Animal.Code %!in% false_list)

Died_earlier_plus_exclude <- data %>%filter(!is.na(MNSS_Day7) & Animal.Code %!in% false_list)
table(Died_earlier_plus_exclude$Condition,Died_earlier_plus_exclude$Gender)


table(data$Condition,data$Gender)
false_list <- c("BNM", "CNM",'D1+3R','A1+3R','NNM','KNM death during operation','N3R','B1R','E1+3R')
length(false_list)
data_sure <- data %>% filter(Animal.Code %!in% false_list)
table(data_sure$Condition,data_sure$Gender)
data$Animal.Code

length(unique(data_sure$Animal.Code))
length(unique(data$Animal.Code))



#data_f <- data %>% drop_na(MNSS_Day7)
data_f<- data
data_f

df_mnss <- data_f %>%
  select(Animal.Code, Gender, Condition, Date.of.Birth, Gender.1, Batch,
         Surgery.Day, Death.at.Day, Drug.Tx.Tube, Dose,
         contains("MNSS"),contains("BW"))

df_rotarod <- data_f %>%
  select(Animal.Code, Gender, Condition, Date.of.Birth, Gender.1, Batch,
         Surgery.Day, Death.at.Day, Drug.Tx.Tube, Dose,
         contains("Rotarod"),BW_BSL, ,BW_Day3, BW_Day7)

df_of <- data_f %>%
  select(Animal.Code, Gender, Condition, Date.of.Birth, Gender.1, Batch,
         Surgery.Day, Death.at.Day, Drug.Tx.Tube, Dose,
         contains("OF"),contains("BW"))

colnames(df_mnss)
colnames(df_rotarod)

# Define timepoints present in the column names
of_timepoints <- c("BSL", "Day3","Day7")
mnss_timepoints <- c("BSL", "2h", "Day1", "Day3","Day7")


df_mnss_long <- df_mnss %>%
  pivot_longer(
    cols = matches(paste0(".*_(", paste(mnss_timepoints, collapse = "|"), ")$")),  # Selects all time-dependent columns
    names_to = c(".value", "Timepoint"),
    names_pattern = "(.+?)_(BSL|2h|Day[0-9]+)"  # Captures the variable name and timepoint
  )
df_mnss_long



df_rotarod_long <- df_rotarod %>%
  pivot_longer(
    cols = matches(paste0(".*_(", paste(of_timepoints, collapse = "|"), ")$")),  # Selects all time-dependent columns
    names_to = c(".value", "Timepoint"),
    names_pattern = "(.+?)_(BSL|2h|Day[0-9]+)"  # Captures the variable name and timepoint
  )



df_mnss_long$Timepoint <- factor(df_mnss_long$Timepoint, levels = c("BSL", "2h", "Day1", "Day3", "Day7"),ordered = TRUE)
df_mnss_long$Timepoint_Numeric <- as.numeric(df_mnss_long$Timepoint)
df_mnss_long$Gender <- as.factor(df_mnss_long$Gender)
df_mnss_long$Condition <- as.factor(df_mnss_long$Condition)

df_rotarod_long$Timepoint <- factor(df_rotarod_long$Timepoint, levels = c("BSL", "Day3", "Day7"),ordered = TRUE)
df_rotarod_long$Timepoint_Numeric <- as.numeric(df_rotarod_long$Timepoint)
df_rotarod_long$Gender <- as.factor(df_rotarod_long$Gender)
df_rotarod_long$Condition <- as.factor(df_rotarod_long$Condition)

df_rotarod_long


##########
### ALL
##########
library(tidyverse)
# Define the "not in" operator
`%!in%` <- function(x, table) {
  !(x %in% table)
}



df_mnss_long

#%%  ======================= ALL SAMPLES ======================================
#%%  ======================= ALL SAMPLES ======================================
#%%  ======================= ALL SAMPLES ======================================

hist(df_rotarod_long$RotarodRPMMedians,breaks=30)

model <- lmer(MNSS ~ Timepoint + Condition  + Dose + Gender+ (1 | Animal.Code),
                data = df_mnss_long)
# ANOVA table
# Perform ANOVA
anova_results <- anova(model)
summary(model)
# Compute p-values
anova_results$"Pr(>F)" <- pf(anova_results$"F value", anova_results$npar, df.residual(model), lower.tail = FALSE)
# Print the updated ANOVA table with p-values
anova_results



# =================== RotarodRPMMedians Pairwise Comparisons ==================
model <- lmer(RotarodRPMMedians ~ Timepoint * Condition + Dose  + (1|Gender) + (1|Animal.Code),
              data = df_rotarod_long)
summary(model)

emmeans_model <- emmeans(model, pairwise ~ Condition | Timepoint, adjust = "bonferroni")
summary(emmeans_model)

model <- glmer(RotarodRPMMedians ~ Timepoint * Condition + Dose + (1|Gender) + (1|Animal.Code),
               family = Gamma(link = "log"), data = df_rotarod_long)

summary(model)
emmeans_model <- emmeans(model, pairwise ~ Condition | Timepoint, adjust = "bonferroni")
summary(emmeans_model)


# ======================== RotarodRPMMedians Trends ============================
poisson_model <- glmer(RotarodRPMMedians ~ Timepoint_Numeric * Condition + Dose + BW + (1|Gender) + (1|Animal.Code),
               family =  poisson(link = "log"), data = df_rotarod_long)
emtrends_model <- emtrends(poisson_model, pairwise ~ Condition, var = "Timepoint_Numeric", adjust = "bonferroni")
summary(emtrends_model)

Gamma_model <- glmer(RotarodRPMMedians ~ Timepoint_Numeric * Condition + Dose + BW + (1|Gender) + (1|Animal.Code),
               family =  Gamma(link = "log"), data = df_rotarod_long)
emtrends_model <- emtrends(Gamma_model, pairwise ~ Condition, var = "Timepoint_Numeric", adjust = "bonferroni")
summary(emtrends_model)

gaussian_model <- lmer(RotarodRPMMedians ~ Timepoint_Numeric * Condition + Dose + BW + (1|Gender) + (1|Animal.Code),
              data = df_rotarod_long)
emtrends_model <- emtrends(gaussian_model, pairwise ~ Condition, var = "Timepoint_Numeric", adjust = "bonferroni")
summary(emtrends_model)

# Compare AIC
AIC(Gamma_model, gaussian_model,poisson_model)
# Compare BIC
BIC(Gamma_model, gaussian_model,poisson_model)

# =============== RotarodRPMMedians Pairwise Comparisons  SCALED ================
df_rotarod_long$RotarodRPMMedians_scaled <- scale(df_rotarod_long$RotarodRPMMedians)
df_rotarod_long$Timepoint_Numeric_scaled <- scale(df_rotarod_long$Timepoint_Numeric)
df_rotarod_long$Dose_scaled <- scale(df_rotarod_long$Dose)
df_rotarod_long$BW_scaled <- scale(df_rotarod_long$BW)
# ======================== RotarodRPMMedians Trends ============================
model <- lmer(RotarodRPMMedians_scaled ~ Timepoint_Numeric_scaled * Condition +
                Dose_scaled + BW_scaled +
                (1|Gender) + (1|Animal.Code),
              data = df_rotarod_long)
emtrends_model <- emtrends(model, pairwise ~ Condition, var = "Timepoint_Numeric_scaled", adjust = "bonferroni")
summary(emtrends_model)

# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------





hist(df_mnss_long_clean$MNSS,breaks=10)
table(df_mnss_long_clean$MNSS)
# =================== MNSS Pairwise Comparisons ==================
df_mnss_long_clean <- df_mnss_long%>%na.omit()
poisson_model <- glmer(MNSS ~ Timepoint * Condition + Dose + BW + (1|Gender) + (1|Animal.Code),
               family = poisson(link = "log"), data = df_mnss_long_clean)
emmeans_model <- emmeans(poisson_model, pairwise ~ Condition | Timepoint, adjust = "tukey")
summary(emmeans_model)

gamma_model <- glmer(MNSS+1 ~ Timepoint * Condition + Dose + BW + (1|Gender) + (1|Animal.Code),
               family = Gamma(link = "log"), data = df_mnss_long_clean)
emmeans_model <- emmeans(gamma_model, pairwise ~ Condition | Timepoint, adjust = "tukey")
summary(emmeans_model)

gauss_model <- lmer(MNSS ~ Timepoint * Condition + Dose + BW + (1|Gender) + (1|Animal.Code),
              data = df_mnss_long_clean)
emmeans_model <- emmeans(gauss_model, pairwise ~ Condition | Timepoint, adjust = "bonferroni")
summary(emmeans_model)

BIC(gamma_model, gauss_model,poisson_model)


# ======================== MNSS Trends ============================
poisson_model <- glmer(MNSS ~ Timepoint_Numeric * Condition + Dose + BW + (1|Gender) + (1|Animal.Code),
               family = poisson(link = "log"), data = df_mnss_long_clean)
emtrends_model <- emtrends(poisson_model, pairwise ~ Condition, var = "Timepoint_Numeric", adjust = "tukey")
summary(emtrends_model)

gamma_model <- glmer(MNSS+1 ~ Timepoint_Numeric * Condition + Dose + BW + (1|Gender) + (1|Animal.Code),
               family = Gamma(link = "log"), data = df_mnss_long_clean)
emtrends_model <- emtrends(gamma_model, pairwise ~ Condition, var = "Timepoint_Numeric", adjust = "tukey")
summary(emtrends_model)

gauss_model <- lmer(MNSS ~ Timepoint_Numeric * Condition + Dose + BW + (1|Gender) + (1|Animal.Code),
                data = df_mnss_long_clean)
emtrends_model <- emtrends(gauss_model, pairwise ~ Condition, var = "Timepoint_Numeric", adjust = "tukey")
summary(emtrends_model)

BIC(gamma_model, gauss_model,poisson_model)


# =============== MNSS Pairwise Comparisons  SCALED ================
df_mnss_long_clean$MNSS_scaled <- scale(df_mnss_long_clean$MNSS)
df_mnss_long_clean$Timepoint_Numeric_scaled <- scale(df_mnss_long_clean$Timepoint_Numeric)
df_mnss_long_clean$Dose_scaled <- scale(df_mnss_long_clean$Dose)
df_mnss_long_clean$BW_scaled <- scale(df_mnss_long_clean$BW)
# ======================== MNSS Trends ============================
model <- lmer(MNSS_scaled ~ Timepoint_Numeric_scaled * Condition +
                Dose_scaled + BW_scaled +
                (1|Gender) + (1|Animal.Code),
              data = df_mnss_long_clean)
emtrends_model <- emtrends(model, pairwise ~ Condition, var = "Timepoint_Numeric_scaled", adjust = "tukey")
summary(emtrends_model)
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------








#%%  ======================= Exclude Animal IDS ======================================
#%%  ======================= Exclude Animal IDS ======================================
#%%  ======================= Exclude Animal IDS ======================================
# List of Animal.Code values that should be FALSE
false_list <- c("BNM", "CNM",'D1+3R','A1+3R','NNM','KNM','N3R','B1R','E1+3R')
df_mnss_long_sure <- df_mnss_long %>% filter(Animal.Code %!in% false_list)
df_rotarod_long_sure <- df_rotarod_long   %>% filter(Animal.Code %!in% false_list)
df_rotarod_long_sure

df_rotarod_long_sure_clean <- df_rotarod_long_sure#%>%na.omit()

# =================== RotarodRPMMedians Pairwise Comparisons ==================
model <- glmer(RotarodRPMMedians ~ Timepoint * Condition + Dose + BW + (1|Gender) + (1|Animal.Code),
               family = Gamma(link = "log"), data = df_rotarod_long_sure_clean)
emmeans_model <- emmeans(model, pairwise ~ Condition | Timepoint, adjust = "bonferroni")
summary(emmeans_model)

# ======================== RotarodRPMMedians Trends ============================
model <- glmer(RotarodRPMMedians ~ Timepoint_Numeric * Condition + Dose + BW + (1|Gender) + (1|Animal.Code),
               family = Gamma(link = "log"), data = df_rotarod_long_sure_clean)
emtrends_model <- emtrends(model, pairwise ~ Condition, var = "Timepoint_Numeric", adjust = "bonferroni")
summary(emtrends_model)


# =============== RotarodRPMMedians Pairwise Comparisons  SCALED ================
df_rotarod_long_sure$RotarodRPMMedians_scaled <- scale(df_rotarod_long_sure$RotarodRPMMedians)
df_rotarod_long_sure$Timepoint_Numeric_scaled <- scale(df_rotarod_long_sure$Timepoint_Numeric)
df_rotarod_long_sure$Dose_scaled <- scale(df_rotarod_long_sure$Dose)
df_rotarod_long_sure$BW_scaled <- scale(df_rotarod_long_sure$BW)
# ======================== RotarodRPMMedians Trends ============================
model <- lmer(RotarodRPMMedians_scaled ~ Timepoint_Numeric_scaled * Condition +
                Dose_scaled + BW_scaled +
                (1|Gender) + (1|Animal.Code),
              data = df_rotarod_long_sure_clean)
emtrends_model <- emtrends(model, pairwise ~ Condition, var = "Timepoint_Numeric_scaled", adjust = "tukey")
summary(emtrends_model)
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------


# =================== MNSS Pairwise Comparisons ==================
df_mnss_long_clean <- df_mnss_long_sure%>%na.omit()
model <- glmer(MNSS ~ Timepoint * Condition + Dose + BW + (1|Gender) + (1|Animal.Code),
               family = poisson(link = "log"), data = df_mnss_long_clean)
emmeans_model <- emmeans(model, pairwise ~ Condition | Timepoint, adjust = "tukey")
summary(emmeans_model)

# ======================== MNSS Trends ============================
model <- glmer(MNSS ~ Timepoint_Numeric * Condition + Dose + BW + (1|Gender) + (1|Animal.Code),
               family = poisson(link = "log"), data = df_mnss_long_clean)
emtrends_model <- emtrends(model, pairwise ~ Condition, var = "Timepoint_Numeric", adjust = "tukey")
summary(emtrends_model)

# =============== MNSS Pairwise Comparisons  SCALED ================
df_mnss_long_clean$MNSS_scaled <- scale(df_mnss_long_clean$MNSS)
df_mnss_long_clean$Timepoint_Numeric_scaled <- scale(df_mnss_long_clean$Timepoint_Numeric)
df_mnss_long_clean$Dose_scaled <- scale(df_mnss_long_clean$Dose)
df_mnss_long_clean$BW_scaled <- scale(df_mnss_long_clean$BW)
# ======================== MNSS Trends ============================
model <- lmer(MNSS_scaled ~ Timepoint_Numeric_scaled * Condition +
                Dose_scaled + BW_scaled +
                (1|Gender) + (1|Animal.Code),
              data = df_mnss_long_clean)
emtrends_model <- emtrends(model, pairwise ~ Condition, var = "Timepoint_Numeric_scaled", adjust = "tukey")
summary(emtrends_model)
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------






# ====================== SEP MALE FEMALE ===========================
# ====================== SEP MALE FEMALE ===========================
# ====================== SEP MALE FEMALE ===========================
males_df_mnss_long_clean <- df_mnss_long_clean%>%filter(Gender=="M")
females_df_mnss_long_clean <- df_mnss_long_clean%>%filter(Gender=="F")

model <- lmer(MNSS_scaled ~ Timepoint_Numeric_scaled * Condition +
                Dose_scaled + BW_scaled + (1|Animal.Code),
              data = males_df_mnss_long_clean)
emtrends_model <- emtrends(model, pairwise ~ Condition, var = "Timepoint_Numeric_scaled", adjust = "tukey")
summary(emtrends_model)

model <- lmer(MNSS_scaled ~ Timepoint_Numeric_scaled * Condition +
                Dose_scaled + BW_scaled + (1|Animal.Code),
              data = females_df_mnss_long_clean)
emtrends_model <- emtrends(model, pairwise ~ Condition, var = "Timepoint_Numeric_scaled", adjust = "tukey")
summary(emtrends_model)
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------


ggplot(df_mnss_long, aes(x = Timepoint, y = MNSS, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(title = "MNSS Scores Over Time",
       x = "Timepoint",
       y = "MNSS Score") +
  theme_minimal()

ggplot(df_mnss_long%>%na.omit(), aes(x = Timepoint, y = MNSS, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  facet_wrap(~ Gender) +
  labs(title = "MNSS Scores Over Time",
       x = "Timepoint",
       y = "MNSS Score") +
  theme_minimal()



# Visualization of MNSS score trends for RA vs Control
ggplot(df_mnss_long_clean, aes(x = Timepoint, y = MNSS, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(title = "Subset MNSS Scores Over Time",
       x = "Timepoint",
       y = "MNSS Score") +
  theme_minimal()

ggplot(df_mnss_long_clean, aes(x = Timepoint, y = MNSS, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  facet_wrap(~ Gender) +
  labs(title = "MNSS Scores Over Time: RA vs Control",
       x = "Timepoint",
       y = "MNSS Score") +
  theme_minimal()

# Visualization of RotarodRPMMedians score trends for RA vs Control
ggplot(df_rotarod_long_sure_clean, aes(x = Timepoint, y = RotarodRPMMedians, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  facet_wrap(~ Gender) +
  labs(title = "RotarodRPMMedians Scores Over Time: RA vs Control",
       x = "Timepoint",
       y = "RotarodRPMMedians Score") +
  theme_minimal()
ggplot(df_rotarod_long_sure_clean, aes(x = Timepoint, y = RotarodRPMMedians, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(title = "RotarodRPMMedians Scores Over Time: RA vs Control",
       x = "Timepoint",
       y = "RotarodRPMMedians Score") +
  theme_minimal()
df_rotarod_long_sure








# =========================== PCA ============================
# =========================== PCA ============================
# =========================== PCA ============================
# =========================== PCA ============================
# Select only the columns for PCA
pca_data <- data_f %>%
  select(contains("MNSS_D"), contains("OFVe"), contains("Rotarod")) %>%
  drop_na() %>%
  mutate(across(everything(), as.numeric))
# Standardize the data (scale)
pca_data_scaled <- scale(pca_data)
# Step 3: Perform PCA
pca_result <- prcomp(pca_data_scaled, center = TRUE, scale. = TRUE)
# Step 4: Extract PCA scores
pca_scores <- as.data.frame(pca_result$x)
# Step 5: Plot the first two principal components (PC1 vs PC2)
library(ggplot2)
cond_df <- data_f %>%
  select(Condition,contains("MNSS_D"), contains("OFVe"), contains("Rotarod")) %>%
  drop_na()

ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = cond_df$Condition)) +  # Color by Condition (or any other grouping variable)
  theme_bw() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2")




Died_earlier_plus_exclude


df_mnss_final <- Died_earlier_plus_exclude %>%
  select(Animal.Code, Gender, Condition, Date.of.Birth, Gender.1, Batch,
         Surgery.Day, Death.at.Day, Drug.Tx.Tube, Dose,
         contains("MNSS"),contains("BW"))
mnss_timepoints <- c("BSL", "2h", "Day1", "Day3","Day7")


df_mnss_final_long <- df_mnss_final %>%
  pivot_longer(
    cols = matches(paste0(".*_(", paste(mnss_timepoints, collapse = "|"), ")$")),  # Selects all time-dependent columns
    names_to = c(".value", "Timepoint"),
    names_pattern = "(.+?)_(BSL|2h|Day[0-9]+)" )
df_mnss_final_long$Timepoint <- factor(df_mnss_final_long$Timepoint, levels = c("BSL", "2h", "Day1", "Day3", "Day7"),ordered = TRUE)
df_mnss_final_long$Timepoint_Numeric <- as.numeric(df_mnss_final_long$Timepoint)
df_mnss_final_long$Gender <- as.factor(df_mnss_final_long$Gender)
df_mnss_final_long$Condition <- as.factor(df_mnss_final_long$Condition)



# =================== MNSS Pairwise Comparisons ==================
model <- glmer(MNSS ~ Timepoint * Condition + Dose + BW + (1|Gender) + (1|Animal.Code),
               family = poisson(link = "log"), data = df_mnss_final_long)
emmeans_model <- emmeans(model, pairwise ~ Condition | Timepoint, adjust = "bonferroni")
summary(emmeans_model)

# ======================== MNSS Trends ============================
model <- glmer(MNSS ~ Timepoint_Numeric * Condition + Dose + BW + (1|Gender) + (1|Animal.Code),
               family = poisson(link = "log"), data = df_mnss_final_long)
emtrends_model <- emtrends(model, pairwise ~ Condition, var = "Timepoint_Numeric", adjust = "bonferroni")
summary(emtrends_model)

# =============== MNSS Pairwise Comparisons  SCALED ================
df_mnss_final_long$MNSS_scaled <- scale(df_mnss_final_long$MNSS)
df_mnss_final_long$Timepoint_Numeric_scaled <- scale(df_mnss_final_long$Timepoint_Numeric)
df_mnss_final_long$Dose_scaled <- scale(df_mnss_final_long$Dose)
df_mnss_final_long$BW_scaled <- scale(df_mnss_final_long$BW)
# ======================== MNSS Trends ============================
model <- lmer(MNSS_scaled ~ Timepoint_Numeric_scaled * Condition +
                Dose_scaled + BW_scaled +
                (1|Gender) + (1|Animal.Code),
              data = df_mnss_final_long)
emtrends_model <- emtrends(model, pairwise ~ Condition, var = "Timepoint_Numeric_scaled", adjust = "tukey")
summary(emtrends_model)
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------

ggplot(df_mnss_final_long, aes(x = Timepoint, y = MNSS, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(title = "MNSS Scores Over Time: RA vs Control",
       x = "Timepoint",
       y = "MNSS Score") +
  theme_minimal()

ggplot(df_mnss_final_long, aes(x = Timepoint, y = MNSS, color = Condition, group = Condition)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  facet_wrap(~ Gender) +
  labs(title = "MNSS Scores Over Time",
       x = "Timepoint",
       y = "MNSS Score") +
  theme_minimal()

males_df_mnss_long_clean <- df_mnss_final_long%>%filter(Gender=="M")
females_df_mnss_long_clean <- df_mnss_final_long%>%filter(Gender=="F")

# Fit a linear mixed-effects model including Gender as a covariate
model <- lmer(MNSS ~ Timepoint * Condition +Dose+BW+ (1 | Animal.Code), data = males_df_mnss_long_clean)
# ANOVA table
anova(model)
# Post-hoc comparisons
posthoc <- emmeans(model, pairwise ~ Condition |Timepoint, adjust = "bonferroni")
summary(posthoc)

# Fit a linear mixed-effects model including Gender as a covariate
model <- lmer(MNSS ~ Timepoint * Condition +Dose+BW+ (1 | Animal.Code), data = females_df_mnss_long_clean)
plot(model)
# ANOVA table
anova(model)
# Post-hoc comparisons
posthoc <- emmeans(model, pairwise ~ Condition |Timepoint, adjust = "bonferroni")
summary(posthoc)




# ================================= ROC CURVE AND PREDICTIONS ==================================
# ================================= ROC CURVE AND PREDICTIONS ==================================
# ================================= ROC CURVE AND PREDICTIONS ==================================
# Load required libraries
library(pROC)
library(caret)

cond_compare <- 'RA'
cond_compare <- 'Everolimus'

df <- Died_earlier_plus_exclude %>% filter(Condition %in% c('Control',cond_compare))
# Assuming your dataframe is named df
# Convert categorical variables to factors
df$Condition <- as.factor(df$Condition)  # RA (1) vs. Control (0)
df$Gender <- as.factor(df$Gender)
df$Animal.Code <- as.factor(df$Animal.Code)

# Fit the logistic regression model
logit_model <- glm(Condition ~ Gender + Dose  + MNSS_2h + MNSS_Day1 + MNSS_Day3 + MNSS_Day7 +
                     RotarodRPMMedians_BSL + RotarodRPMMedians_Day3 + RotarodRPMMedians_Day7,
                   data = df, family = "binomial")
# Summarize the model
summary(logit_model)
plot(logit_model)

# Predict probabilities
probabilities <- predict(logit_model, type = "response")
# Predict classes
predicted_classes <- ifelse(probabilities > 0.5, cond_compare, "Control")
predicted_classes <- factor(predicted_classes, levels = levels(df$Condition))

length(predicted_classes)
length(df$Condition)
# Create confusion matrix
confusion_matrix <- table(Actual = df$Condition, Predicted = predicted_classes)
print(confusion_matrix)

# Convert confusion matrix to a data frame for plotting
confusion_df <- as.data.frame(confusion_matrix)

# Plot the confusion matrix
ggplot(confusion_df, aes(x = Predicted, y = Actual, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), color = "black", size = 6) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Confusion Matrix",
       x = "Predicted Condition",
       y = "Actual Condition",
       fill = "Count") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"))

# Extract values from the confusion matrix
TN <- confusion_matrix[1, 1] # True Negatives
FP <- confusion_matrix[1, 2] # False Positives
FN <- confusion_matrix[2, 1] # False Negatives
TP <- confusion_matrix[2, 2] # True Positives

# Calculate metrics
accuracy <- (TN + TP) / (TN + FP + FN + TP)
precision <- TP / (TP + FP)
recall <- TP / (TP + FN)
f1_score <- 2 * (precision * recall) / (precision + recall)

# Print metrics
cat("Accuracy:", accuracy, "\n")
cat("Precision:", precision, "\n")
cat("Recall:", recall, "\n")
cat("F1-Score:", f1_score, "\n")

# Create ROC curve
roc_curve <- roc(df$Condition, probabilities)
# Plot ROC curve
plot(roc_curve, main = "ROC Curve for Control vs RA (then Everolimus)", col = "blue", lwd = 2)
auc_value <- auc(roc_curve)
legend("bottomright", legend = paste("AUC =", round(auc_value, 2)), bty = "n")

