# --------------------------------------------------------------------------------
# Cytokine Activity Analysis and Differential Enrichment Score Comparison in Clusters
# --------------------------------------------------------------------------------
# Description:
# This script performs an analysis on cytokine activity across two clusters.
# The analysis includes loading data uploaded to the IREA platform, calculating
# cytokine activation scores for each cluster, and visualizing differences in enrichment
# scores between the clusters to assess activation and significance levels.

# Steps in the Analysis:
# 1. Data Import:
#    - Load cytokine activity data, including differential expression values
#      for two clusters.
#    - This data allows for calculating and comparing cytokine activation scores
#      across clusters.

# 2. Enrichment Score Comparison:
#    - Import enrichment scores (or activation scores) for each cluster.
#    - These scores provide a quantitative measure of cytokine activity for each
#      cluster, enabling comparison between them.
#    - Statistical significance of each cytokine’s differential expression is also
#      examined to highlight key cytokines with meaningful activation differences.

# 3. Visualization:
#    - Use ggplot2 to create a comparative plot of enrichment scores for cytokines
#      across clusters.
#    - Each cytokine is plotted with its enrichment score, with points sized by the
#      degree of activation and labeled to indicate statistical significance.
#    - This visual representation provides insights into which cytokines show
#      significant activity differences between the clusters.

# Goal:
# This script aims to identify cytokines with significant activation differences
# between clusters, providing insights into the molecular drivers of cluster-specific
# phenotypes based on cytokine activity.

# --------------------------------------------------------------------------------
# Load necessary libraries
# --------------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
setwd('~/')
save_Fig_dir <- 'Figure2/'

#---------------------------------------------------------------------
# HELPER  Variables
#---------------------------------------------------------------------
TICK_Colors = list(

  ### Interferons
  'IFNa1'='#20aa8b',
  'IFNb'='#20aa8b',
  'IFNe'='#20aa8b',
  'IFNk'='#20aa8b',
  'IFNg'='#20aa8b',
  'IFNl'='#20aa8b',
  'IFNl1'='#20aa8b',
  'IFNl2'='#20aa8b',

  ### IL1 group
  'IL1a'= '#922a6e',
  'IL1b'= '#922a6e',
  'IL1ra'= '#922a6e',
  'IL18'= '#922a6e',
  'IL33'= '#922a6e',
  'IL36a'= '#922a6e',
  'IL36RA'= '#922a6e',

  ### Common g chain / IL13
  'IL2'='#F6C76D',
  'IL4'='#F6C76D',
  'IL13'='#F6C76D',
  'IL15'='#F6C76D',
  'IL7'='#F6C76D',
  'TSLP'='#F6C76D',
  'IL9'='#F6C76D',
  'IL21'='#F6C76D',

  ### Common b chain
  'IL3'='#8668c0',
  'IL5'='#8668c0',
  'GM-CSF'='#8668c0',

  ### IL6 / IL12
  'IL6'='#198591',
  'IL11'='#198591',
  'IL27'='#198591',
  'IL30'='#198591',
  'IL31'='#198591',
  'LIF'='#198591',
  'OSM'='#198591',
  'Cardiotrophin-1'='#198591',
  'Neuropoietin'='#198591',
  'IL12'='#198591',
  'IL23'='#198591',
  'IL-Y'='#198591',

  ### IL10
  'IL10'='#df70b0',
  'IL19'='#df70b0',
  'IL20'='#df70b0',
  'IL22'='#df70b0',
  'IL24'='#df70b0',

  ### IL17
  'IL17A'='#096954',
  'IL17B'='#096954',
  'IL17C'='#096954',
  'IL17D'='#096954',
  'IL17E'='#096954',
  'IL17F'='#096954',

  ### Growth Factor
  'Flt3l'='#472d8c',
  'IL34'='#472d8c',
  'M-CSF'='#472d8c',
  'G-CSF'='#472d8c',
  'SCF'='#472d8c',
  'EGF'='#472d8c',
  'VEGF'='#472d8c',
  'FGF-basic'='#472d8c',
  'HGF'='#472d8c',
  'IGF-I'='#472d8c',

  ### TNF
  'LTA1-B2'='#53b8c4',
  'LTA2-B1'='#53b8c4',
  'TNFa'='#53b8c4',
  'OX40L'='#53b8c4',
  'CD40L'='#53b8c4',
  'FasL'='#53b8c4',
  'CD27L'='#53b8c4',
  'CD30L'='#53b8c4',
  '41BBL'='#53b8c4',
  'TRAIL'='#53b8c4',
  'RANKL'='#53b8c4',
  'TWEAK'='#53b8c4',
  'APRIL'='#53b8c4',
  'BAFF'='#53b8c4',
  'LIGHT'='#53b8c4',
  'TL1A'='#53b8c4',
  'GITRL'='#53b8c4',

  ### Complement
  'C3a'='#b7572e',
  'C5a'='#b7572e',

  ### Other
  'Prolactin'='#91a160',
  'Leptin'='#91a160',
  'Adiponectin'='#91a160',
  'Resistin'='#91a160',
  'TGF-beta-1'='#91a160',
  'GDNF'='#91a160',
  'Persephin'='#91a160',
  'Noggin'='#91a160',
  'Decorin'='#91a160',
  'TPO'='#91a160'

)




# Read a CSV file named 'MTC_123_df_irea.csv' from the 'Downloads' folder
# into a data frame called 'Control'.
Control <- read.csv('MTC_123_df_irea.csv')

# Read a CSV file named 'MTC_456_df_irea.csv' from the 'Downloads' folder
# into a data frame called 'CASE'.
CASE <- read.csv('MTC_456_df_irea.csv')

# Add a new column 'Cond' to the 'Control' data frame, setting its value to "MTC_123"
# to label it as data from the control condition.
Control$Cond <- "MTC_123"

# Add a new column 'Cond' to the 'CASE' data frame, setting its value to "MTC_456"
# to label it as data from the case condition.
CASE$Cond <- "MTC_456"

# Combine (row-bind) the 'Control' and 'CASE' data frames into a single data frame 'df',
# stacking the rows from both datasets to create a unified dataset.
df <- rbind(Control, CASE)

# Calculate the negative log10 of the adjusted p-value (padj) + 1,
# storing the result in a new column 'Log10Qval' in 'df'.
df$Log10Qval <- -log10(df$padj + 1)

# Copy values from 'Log10Qval' to initialize a new column 'directed' in 'df'.
df$directed <- df$Log10Qval

# Create a 'multiplier' variable that is -1 if 'Cond' is 'MTC_123' and 1 otherwise.
# This allows for direction-based modification of 'Log10Qval'.
multiplier <- ifelse(df$Cond == 'MTC_123', -1, 1)

# Apply the 'multiplier' to 'Log10Qval' and store the result in 'directed',
# creating a signed version of 'Log10Qval' based on the condition in 'Cond'.
df$directed <- df$Log10Qval * multiplier

# Create a 'Significant' column, marking rows with adjusted p-value (padj) < 0.05
# as "+" (significant), and "" (not significant) otherwise.
df$Significant <- ifelse(df$padj < 0.05, '+', "")

# Create a data frame named 'my_dataframe' with two columns:
# 'Cytokine' contains the names from the TICK_Colors object,
# 'Groups' contains the values from the TICK_Colors object as a vector.
my_dataframe <- data.frame('Cytokine' = names(TICK_Colors),
                           'Groups' = as.vector(unlist(c(TICK_Colors))))

# Filter 'my_dataframe' to include only rows where the 'Cytokine' column
# matches any value in the 'Cytokine' column of the 'df' data frame.
# Store this filtered data frame in 'my_dataframe_sub'.
my_dataframe_sub <- my_dataframe %>% filter(Cytokine %in% df$Cytokine)

# Set the row names of 'my_dataframe_sub' to the values in the 'Cytokine' column.
rownames(my_dataframe_sub) <- my_dataframe_sub$Cytokine

# Check which cytokines in 'my_dataframe_sub' are present in the 'Cytokine' column of 'df'.
my_dataframe_sub$Cytokine %in% df$Cytokine

# Check which cytokines in the 'df' data frame are present in the 'my_dataframe_sub' data frame.
as.vector(df$Cytokine) %in% my_dataframe_sub$Cytokine

# Create a unique ordering of cytokines from the 'Cytokine' column of 'my_dataframe_sub'.
order_ticks <- unique(my_dataframe_sub$Cytokine)

# Convert 'Cytokine' column in 'df' to a factor with levels defined by 'order_ticks'.
# This ensures 'df$Cytokine' has a consistent ordering based on 'order_ticks'.
df$Cytokine <- factor(df$Cytokine, levels = order_ticks)

# Order rows in 'df' based on the levels of the 'Cytokine' factor, saving the result as 'df_ord'.
df_ord <- df[order(df$Cytokine), ]

# Add a 'Groups' column to 'df_ord' by mapping the 'Cytokine' values in 'df_ord'
# to the 'Groups' column in 'my_dataframe_sub'.
df_ord$Groups <- my_dataframe_sub[df_ord$Cytokine, ]

# Add a new column 'ES' (Enrichment Score) to 'df_ord' by copying values from the 'Enrichment.Score' column.
df_ord$ES <- df_ord$Enrichment.Score

# Display the last few rows of the updated 'df_ord' to check final data structure.
tail(df_ord)


# Create a ggplot object named 'Fig2f' based on the 'df_ord' data frame.
# Map 'Cytokine' to the x-axis, 'Cond' to the y-axis,
# 'ES' to both color and fill aesthetics, the absolute value of 'ES' to point size,
# and 'Significant' to the text label for points.
Fig2f <- ggplot(df_ord, aes(x = Cytokine,
                            y = Cond,
                            color = ES,
                            fill = ES,
                            size = abs(ES),
                            label = Significant)) +

  # Add a scatter plot layer with points for each data entry.
  geom_point() +

  # Overlay text labels on points with alpha transparency and white color for contrast.
  geom_text(alpha = 1, size = 4, colour = "white") +

  # Rotate the axis labels to enhance readability of x-axis labels.
  RotatedAxis() + theme_cowplot() +

  # Adjust text styles: set global text size to 12,
  # rotate x-axis labels 90 degrees with adjustments for alignment,
  # and set y-axis text size to 12.
  theme(text = element_text(size = 12),
        axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12)) +

  # Center the plot title and position the legend at the bottom.
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'bottom') +

  # Define a gradient color scale for 'ES' with 'navy' for low values,
  # 'white' for mid values, and 'firebrick3' for high values.
  scale_color_gradient2(low = 'navy', mid = 'white', high = 'firebrick3') +

  # Remove labels for both x and y axes.
  xlab("") + ylab("") +

  # Set x-axis text color based on the 'Groups' column in 'my_dataframe_sub'
  # for visual distinction of each cytokine group.
  theme(axis.text.x = element_text(colour = my_dataframe_sub$Groups))



# Display the resulting plot 'Fig2f'.
pdf(paste0(save_Fig_dir,'2024_12_Figure_2F_GL.pdf'),width = 16,height = 3.5)
Fig2f
dev.off()

png(paste0(save_Fig_dir,'2024_12_Figure_2F_GL.png'),width = 16,height = 3.5)
Fig2f
dev.off()