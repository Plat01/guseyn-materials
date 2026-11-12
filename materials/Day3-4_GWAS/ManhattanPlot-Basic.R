
library(tidyverse)


# To erase objects from the working space - Clean up of the memory
rm(list = ls())
graphics.off()

#gwasResultsKr3 <- read.table("GLM_Aphano-KR3_Stats.txt", sep = "\t", dec = ".", header = TRUE)
gwasResultsKr3 <- read.table("CMLM_Aphano-KR3_240K_chr3.hmp.TASSEL_Filter_+_Filter_+_AphanoRRI_prTassel_+_MtrPopStruct_Admixture_K8_stats.txt", sep = "\t", dec = ".", header = TRUE)
## gwasResultsKr7 <- read.table("PodCoiling_Chrom7_Stats-MLM-TASSEL.csv.txt", sep = "\t", dec = ".", header = TRUE)

## gwasResults <- rbind(gwasResultsKr3, gwasResultsKr7)
gwasResults <- gwasResultsKr3


# We need to prepare the X axis. Indeed we do not want to display the position of SNP in bp, but just show the chromosome name instead.

axisdf = gwasResults %>% group_by(Chr) %>% summarize(center=( max(Pos) + min(Pos) ) / 2 )

#Ready to make the plot using ggplot2:

ggplot(gwasResults, aes(x = Pos, y = -log10(p))) +
    
    # Show all points
    geom_point( aes(color=as.factor(Chr)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("firebrick", "darkgreen" ), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = paste("Chrom.", axisdf$Chr), breaks= axisdf$center ) +
#    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )


