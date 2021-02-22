# Loading packages-----------------------------------------------------------------
library(adegraphics)
library(adegenet)
library(factoextra)
library(ggplot2)
library(pegas)

# Set up, data load, formatting-------------------------------------------
setwd ("C:/Users/parentg/Documents/R/Projects/2020-Limnology-and-Oceanography_Calanus-comment/Script/Results")

geno<- read.loci("C:/Users/parentg/Documents/R/Projects/2020-Limnology-and-Oceanography_Calanus-comment/Data/indel_all.txt", allele.sep = "-", loci.sep = "\t", col.loci = 2:7, col.pop = 8, na.strings = "NA")
info <- read.delim("C:/Users/parentg/Documents/R/Projects/2020-Limnology-and-Oceanography_Calanus-comment/Data/indel_all.txt", header=TRUE)
geno.gi <- loci2genind(geno, ploidy=2)

# PCA---------------------------------------------------------------------
geno.giNA<-tab(geno.gi, NA.method="mean")

pca <- dudi.pca(geno.giNA, center=TRUE, scale=FALSE, scannf = FALSE, nf = 2)
(k <- 100 * pca$eig/sum(pca$eig))

eig.val <- get_eigenvalue(pca)
info$Spp<-as.factor(info$Spp)

fig.PCA<-ggplot(pca$li, aes(x=Axis1,y=Axis2))+
               geom_point(size = 2, shape = 19, alpha=0.6, aes(colour=info$Spp))+
               xlab("PC 1 (52.0%)") + ylab("PC 2 (22.1%)")+
               theme_bw() + theme(panel.grid.major = element_blank())

fig.PCA + theme_classic() + scale_color_manual(values = c("yellow", "blue","green3"))+
  scale_x_continuous(expand=c(0,0),limits=c(-4, 6), breaks=seq(-4,6, by=2))+scale_y_continuous(expand=c(0,0),limits=c(-4, 3), breaks=seq(-4, 3, by=1))+
  theme(axis.text.x = element_text(color = "gray20", size = 12),
        axis.text.y = element_text(color = "gray20", size = 12),  
        axis.title.x = element_text(color = "gray20", size = 12),
        axis.title.y = element_text(color = "gray20", size = 12),legend.position = "none")
ggsave("Fig1_final.tiff", width = 10, height = 10, units = "cm", dpi=300)
