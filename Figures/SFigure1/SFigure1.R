library(ggplot2)
library(Seurat)
library(cowplot)
library(patchwork)
mydata <- readRDS('/Users/dimitrioskyriakis/Documents/Projects/2025_ICH_Project/Data_backup/2.2024_11_Final_Immune.rds')

Real_donor_ids <- sort(unique(mydata@meta.data$donor))
Crypto_ids <- sprintf("Donor_%02d", 1:length(Real_donor_ids))
names(Crypto_ids) <- Real_donor_ids
Crypto_ids

df <- mydata@meta.data
# Use lapply to replace real IDs with encrypted IDs
df$donor <- unlist(lapply(df$donor, function(x) Crypto_ids[x]))
df[,c('orig.ident','donor')]
mydata$donor <-  df$donor
g1 <- VlnPlot(mydata,features=c("nGenes","nCounts"),stack=T,group.by="donor")&theme_bw()
g2 <- VlnPlot(mydata,features=c("log.genes","log.umis"),stack=T,group.by="donor")&theme_bw()
g3 <- ggplot(mydata@meta.data,aes(x=donor,y=0.2)) +
  geom_tile(aes(fill = Kit), colour = "black")+ theme_bw()+coord_flip()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),
        axis.text.x=element_blank(),axis.ticks.x=element_blank())+ylab("Kit")+
  theme(axis.title.x= element_text(angle=90, vjust=.5, hjust=1),
        plot.margin = unit(c(0,0,0,0), "cm"))+ 
  scale_fill_grey(start=0.8, end=0.2, na.value = NA) 
g4 <- ggplot(mydata@meta.data,aes(x=donor,y=0.2)) +
  geom_tile(aes(fill = sex), colour = "black")+ theme_bw()+coord_flip()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),axis.ticks.x=element_blank())+ylab("sex")+
  theme(axis.title.x= element_text(angle=90, vjust=.5, hjust=1),
        plot.margin = unit(c(0,0,0,0), "cm"))+ 
  scale_fill_manual(values=c("#E69F00", "#56B4E9"), na.value = NA)
g5 <- ggplot(mydata@meta.data,aes(x=donor,y=0.2)) +
  geom_tile(aes(fill = race), colour = "black")+ theme_bw()+coord_flip()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),axis.ticks.x=element_blank())+ylab("race")+
  scale_fill_manual(values=c("#FFFE99","#F5C085","#BDADD4","#7FC97F"), na.value = NA)+
  theme(axis.title.x= element_text(angle=90, vjust=.5, hjust=1),
        plot.margin = unit(c(0,0,0,0), "cm"))
g6 <- ggplot(mydata@meta.data,aes(x=donor,y=0.2)) +
  geom_tile(aes(fill = as.factor(mRS)), colour = "black")+ theme_bw()+coord_flip()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),axis.ticks.x=element_blank())+ylab("mRS")+
  scale_fill_manual(values=c("#4D4D4D","#999999","#E0E0E0","#F9DAC7","#EF8962","#B3242B"), na.value = NA)+
  theme(axis.title.x= element_text(angle=90, vjust=.5, hjust=1),
        plot.margin = unit(c(0,0,0,0), "cm"))
g7 <- ggplot(mydata@meta.data,aes(x=donor))+geom_bar()+theme_bw()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank())+
  coord_flip()

sfig1A <- g1+g3+g4+g5 +g6+g7+ 
  plot_layout(guides = "collect",widths = c(2,0.2,0.2,0.2,0.2,1.5))& theme(legend.position = 'right')




Categories <- c('sex','race','BL','CL','dx','mRS')

library(dplyr)
df_unique <- df[,c('orig.ident','donor',Categories)] %>%
  distinct(donor, .keep_all = TRUE)

df_unique_ord <- df_unique[order(df_unique$donor),]
rownames(df_unique_ord) <- df_unique_ord$donor
df_unique_ord


# Create a function to generate a pie chart for each category
create_pie_chart <- function(category) {
  # Subset data for the category
  category_data <- df_unique_ord[[category]]
  
  # Summarize data
  summary_data <- as.data.frame(table(category_data))
  colnames(summary_data) <- c(category, "Count")
  
  # Calculate percentages
  summary_data$Percentage <- round(100 * summary_data$Count / sum(summary_data$Count))
  
  # Create pie chart
  colors_list <- 

  p <- ggplot(summary_data, aes(x = "", y = Count, fill = category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    geom_text(aes(label = paste0(Percentage, "%")), position = position_stack(vjust = 0.5),size = 3) +
    labs(title = paste("Pie Chart:", category), fill = category) +
    theme_void() +
    theme(legend.position = "right")+
    theme(
      legend.position = "right",
      legend.text = element_text(size = 6)  # Adjust the size of legend text here
    )
  if(category =='sex'){p <- p+scale_fill_manual(values=c("#E69F00", "#56B4E9"), na.value = NA)}
  if(category =='race'){p <- p+scale_fill_manual(values=c("#FFFE99","#F5C085","#BDADD4","#7FC97F"), na.value = NA)}
  if(category =='dx'){p <- p+scale_fill_manual(values=c("#E78AC3","#8DA1CB", "#EF8C61"), na.value = NA)}
  if(category =='mRS'){p <- p+scale_fill_manual(values=c("#4D4D4D","#999999","#E0E0E0","#F9DAC7","#EF8962","#B3242B"), na.value = NA)}
  p 
}

# Generate pie charts for all categories
plots <- lapply(Categories, create_pie_chart)

p12 <- (plots[[1]] + plots[[2]])
p34 <- (plots[[3]] + plots[[4]])
p56 <-  (plots[[5]] + plots[[6]])


options(repr.plot.width=8, repr.plot.height=6)
sfig1A
sfig1B <- p12/p34/p56
sfig1B

pdf('SFigures/SFigure1/SupFigure1_A.pdf', width = 8, height = 6)
options(repr.plot.width=8, repr.plot.height=6)
sfig1A +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 5)  # Adjust the size of legend text here
  )
dev.off()

pdf('SFigures/SFigure1/SupFigure1_B.pdf', width = 8, height = 6)
options(repr.plot.width=8, repr.plot.height=6)
sfig1B
dev.off()


pdf('SFigures/SFigure1/SupFigure1.pdf', width = 12, height = 6)
options(repr.plot.width=12, repr.plot.height=6)
sfig1A + sfig1B + plot_layout(widths = c(2,0.2,0.2,0.2,0.2,1.5,5))
dev.off()
