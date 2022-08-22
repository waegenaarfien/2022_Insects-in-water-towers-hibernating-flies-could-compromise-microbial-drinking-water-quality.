#### 16S DATA OTU TABLE ####

# Set working directory 
setwd("/media/projects1/Fien/Vliegenexperimenten/December")

# load  packages 
library(VennDiagram)
library(scales)   #for percent formatting
library(ggplot2)  #for adequate plotting
library(plyr)     #data wrangling (mapvalues)
library(tidyr)    #tidy data
library(phyloseq) #for microbiome census data processing
library(ade4)     #for ecological calculations
library(splitstackshape) #for csplit
library(knitr)    #for simple markdown tables
library(xtable)   #for more advanced tables
library(ape)      #dealing with phylogenetic trees
library("openxlsx") #handling Excel files without Java deps
library(readxl)   #faster handling of excel files
library(data.table) #data wrangling
library(SPECIES)  #alpha diversity estimators
library(parallel) #parallel computation in R
library(RCM) # advanced ordination
library(RColorBrewer)
library(car) #for Anova data analysis
library(stringr)
library(cowplot)
library(lattice)
library(gridExtra)
library(grid)
library("tibble")
library("edgeR")
library("magrittr")
library(gganimate)
library(ggrepel)
library(dplyr)
library(vegan)
library(microbiome)
#install.packages("writexl")
library(writexl)

#### DATA TRANSFORMATION ####

# Read OTU table (samples are columns, OTUs are rows)
######################## DWG ##############################"
otu.tbl <- read_xlsx(paste0("/media/projects1/Fien/Vliegenexperimenten/December/Genus_import.xlsx"),sheet = 1,col_types=NULL)
x <- otu.tbl$Genus
otu.tbl <- otu.tbl[,-1]
rownames(otu.tbl) <- x

# Make metadata file from sample names
metadata.s <- data.frame(do.call(rbind, lapply(strsplit(colnames(otu.tbl),"_"), rbind)))
metadata <- cbind(colnames(otu.tbl),metadata.s)
colnames(metadata) <- c("Full_name", "Sample","Replicate")

# Make phyloseq object 
otumat.ns <- as.matrix(otu.tbl)
OTU       <- otu_table(otumat.ns,taxa_are_rows = TRUE)
physeqobj <- phyloseq(OTU)
physeqobj.meta <- merge_phyloseq(physeqobj,metadata)


#### PLOT THEME ####

# Standardised function for all figures so they are publishable
papertheme <- function(ggplotobject){
  library(ggplot2)  
  ggplotnewobject <- ggplotobject + 
    theme(axis.ticks.x = element_line(colour = "#696969",size=0.5), 
          axis.ticks.y =element_line(colour = "#696969",size=0.5),
          axis.ticks.length=unit(0.2,"cm"),
          panel.background = element_rect(fill = "white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.line.x  = element_blank(),
          axis.line.y  = element_blank(),
          plot.title = element_text(size = 16, family="sans", face="bold"),
          axis.text.x = element_text(size = 14,colour="#696969",family="sans"),
          axis.text.y = element_text(size = 14,colour="#696969",family="sans"),
          axis.title.x = element_text(size =16,color="#696969",family="sans"),
          axis.title.y = element_text(size=16,color="#696969",family="sans"),
          strip.text.x = element_text(size = 14),
          legend.key = element_blank(),
          legend.text = element_text(size=12,family="sans", color="#696969"),
          legend.title=element_text(size=13,family="sans", color="#696969", face = "bold"),
          strip.background = element_blank(),
          axis.line=element_blank(),
          panel.border = element_rect(colour = "#696969",size=0.5, fill = "transparent"),
          plot.margin = margin(0.5, 0.5, 0.5, 0, "cm"))+
    guides(shape = guide_legend(override.aes = list(size = 1.5)),
           color = guide_legend(override.aes = list(size = 1.5))) +
    theme(legend.title = element_text(size = 10), 
          legend.text  = element_text(size = 10),
          legend.key.size = unit(0.1, "lines"))
  return(ggplotnewobject)}


# Load colour palette for all the graphs except bar data
values <- brewer.pal(8, "RdYlBu")

# Define colours
n <- 19 #n <- 19
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
speciesPalette <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


#### BARPLOT ALL DATA ####
# Make relative abundances as percentages
otus_perc <- otu.tbl / rep(colSums(otu.tbl), each = nrow(otu.tbl))*100

# Make dataframe for plotting
otus <- data.frame(otus_perc)
otus$Genus <- rownames(otu.tbl)
otus_m <- data.frame(na.omit(reshape2::melt(otus)))
otus_m$variable <- as.character(otus_m$variable)
otus_m$value <- as.numeric(otus_m$value)
metadata.s <- data.frame(do.call(rbind, lapply(strsplit((otus_m$variable),"_"), rbind)))
otus_m <- cbind(otus_m,metadata.s)
colnames(otus_m) <- c("Genus","Full_name","Abundance","Sample","Replicate")

# Add absolute abundances 
cell_dens_LE <- 247262.22
cell_dens_DE <- 710500.00
cell_dens_BL <- 46251.11
cell_dens_ST <- 421.1214

for (i in 1:length(otus_m$Sample)){
  if (otus_m$Sample[i] == "LE"){
    otus_m$AbsAb[i] <- cell_dens_LE*otus_m$Abundance[i]/100}
  if (otus_m$Sample[i] == "DE"){
    otus_m$AbsAb[i] <- cell_dens_DE*otus_m$Abundance[i]/100}
  if (otus_m$Sample[i] == "ST"){
    otus_m$AbsAb[i] <- cell_dens_ST*otus_m$Abundance[i]/100}
  if (otus_m$Sample[i] == "BL"){
    otus_m$AbsAb[i] <- cell_dens_BL*otus_m$Abundance[i]/100}
}

write_xlsx(x = otus_m, path = "otus_m.xlsx", col_names = TRUE)

statistic_dataset <- otus_m

# Ordenen obv abundanties 
otus_m <- otus_m[order(-otus_m$Abundance),]

# Genus based
TopGenera <- otus_m$Genus[1:25]

# Add function defining "not in"
`%!in%` <- purrr::compose(`!`, `%in%`)

#Change to others (all genera less than 25 most abundant)
otus_m[otus_m$Genus %!in% TopGenera,]$Genus <- "Others"

# Reorder levels 
old.lvl <- levels(as.factor(otus_m$Genus))
otus_m$Genus <- factor(otus_m$Genus, levels=c("Others",sort(old.lvl[old.lvl!="Others"], decreasing=F)))

# Change names in metadata
for (i in 1:length(otus_m$Sample)){
  if (otus_m$Sample[i] == "BL"){
    otus_m$Sample[i] <- "Blank"
  }
  if (otus_m$Sample[i] == "DE"){
    otus_m$Sample[i] <- "Dead flies"
  }  
  if (otus_m$Sample[i] == "LE"){
    otus_m$Sample[i] <- "Living flies"
  }
  if (otus_m$Sample[i] == "ST"){
    otus_m$Sample[i] <- "Start"
  }
}

otus_m$Sample <- factor(otus_m$Sample,                                    # Change ordering manually
                  levels = c("Start", "Blank", "Dead flies", "Living flies"))
# Plot relative abundance
plot_all_1 <- ggplot(otus_m, aes(fill = Genus, y = Abundance)) +
  geom_bar(aes(x=Replicate), stat="identity", position="stack", colour="grey")+
  facet_grid(.~Sample,scales = "free_x")

plot_all_1 <- papertheme(plot_all_1) + 
  labs(y="Relative abundance (%)", x="Replicate") +
  theme(strip.text.y = element_text(face="bold", size=15),
        strip.text.x = element_text(face="bold", size=12),
        plot.title = element_text(hjust=0.5, size=22),
        panel.border = element_rect(fill=NA,color="#696969", size=0.5, linetype="solid"),
        strip.background = element_rect(fill=NA,color="#696969", size=0.5, linetype="solid"),
        legend.text = element_text(size=15,family="sans"), axis.text.x = element_text(angle = -45, hjust = 0),
        legend.key.size = unit(1,"line") ) +
  scale_fill_manual(values= c("grey", speciesPalette[1:25]), name= "Genus")+ 
  theme_bw(base_size = 20, base_family = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background = element_rect(colour = "black", fill = "#EFEFEE"))

plot_all_1

dpi = 300
png("figures/Barplot_OTUS_RunA_all_top25.png",width=15*dpi,height=8*dpi,res=dpi)
plot_all_1
dev.off()

# Plot absolute abundance
plot_all_2 <- ggplot(otus_m, aes(fill = Genus, y=AbsAb)) +
  geom_bar(aes(x=Replicate), stat="identity", position="stack", colour="grey")+
  facet_grid(.~Sample,scales = "free_x")

plot_all_2 <- papertheme(plot_all_2) + 
  labs(y="Absolute abundance (cells/mL)", x="Replicate") +
  theme(strip.text.y = element_text(face="bold", size=15),
        strip.text.x = element_text(face="bold", size=12),
        plot.title = element_text(hjust=0.5, size=22),
        panel.border = element_rect(fill=NA,color="#696969", size=0.5, linetype="solid"),
        strip.background = element_rect(fill=NA,color="#696969", size=0.5, linetype="solid"),
        legend.text = element_text(size=15,family="sans"), axis.text.x = element_text(angle = -45, hjust = 0),
        legend.key.size = unit(1,"line") ) +
  scale_fill_manual(values= c("grey", speciesPalette[1:25]), name= "Genus")+ 
  theme_bw(base_size = 20, base_family = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background = element_rect(colour = "black", fill = "#EFEFEE"))

plot_all_2

dpi = 300
png("figures/Barplot_OTUS_RunA_all_top25_AbsAb.png",width=15*dpi,height=8*dpi,res=dpi)
plot_all_2
dev.off()

# Merge rel and abs abundances
plot_rel_abs <- ggarrange(plot_all_1, plot_all_2, labels = c("A", "B"),font.label = list(size = 20, face = "bold", color ="black"),
                          common.legend = TRUE, legend = "right")
plot_rel_abs
# Export
dpi=300
png("figures/plot_rel_abs.tif",width=16*dpi,height=7*dpi,res=dpi)
plot_rel_abs
dev.off()


#### BETA DIVERSITY ####

# PcoA
set.seed(4798)
phylo.pcoa <- ordinate(physeqobj.meta, "PCoA", "bray")

# Calculate coordinates
beta.div.co <- data.frame(phylo.pcoa[["vectors"]][,c(1,2)])
colnames(beta.div.co) <- c("Axis1","Axis2")
var <-  round(100*phylo.pcoa[["values"]][["Eigenvalues"]]/(sum(phylo.pcoa[["values"]][["Eigenvalues"]])),1)

x <- c(brewer.pal(12,"Paired"),brewer.pal(12, "Set3"))

# Change names in metadata
for (i in 1:length(metadata$Sample)){
  if (metadata$Sample[i] == "BL"){
    metadata$Sample[i] <- "Blank"
  }
  if (metadata$Sample[i] == "DE"){
    metadata$Sample[i] <- "Dead flies"
  }  
  if (metadata$Sample[i] == "LE"){
    metadata$Sample[i] <- "Living flies"
  }
  if (metadata$Sample[i] == "ST"){
    metadata$Sample[i] <- "Start"
  }
}

#x <- c("#D57E7E","#A2CDCD","#C6D57E","#FFD24C")
#x <- c("#E7E0C9", "#C1CFC0","#6B7AA1", "#11324D")
x <- c("#6B7AA1","#DAB88B", "#C1CFC0", "#11324D")


# Plot
pcoa.plot <- ggplot() +
  geom_point(aes(x=beta.div.co$Axis1, y=beta.div.co$Axis2,fill=metadata$Sample,color=metadata$Sample), size = 4) +
  labs(x= paste0("Axis1 (",var[1], "%)"), y=paste0("Axis2 (",var[2], "%)"), color = "Condition", fill = "Condition", size = 6)+
  #scale_fill_manual(values=x)+
  scale_color_manual(values =x)+
  scale_shape()+
  coord_fixed(ratio = 1)+
  theme_bw(base_size = 30, base_family = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background = element_rect(colour = "black", fill = "#EFEFEE"))
#pcoa.plot <- papertheme(pcoa.plot)
pcoa.plot
# Export plot 
dpi <- 300
png("figures/diversity.png",width=12*dpi,height=8*dpi,res=dpi)
pcoa.plot
dev.off()

####### STATISTICS 

# Check normality: NOT NORMAL --> KRUSKALL WALLIS TEST
library(ggpubr)
ggqqplot(statistic_dataset$Abundance)
ggqqplot(log10(statistic_dataset$Abundance))

# Kruskal test --> significant difference between all the groups
library(dplyr)
kruskal.test(Abundance ~ Sample, data = otus_m)

# Which groups are different: Pairwise Wilcoxon Rank Sum Tests; correction with bonferoni
pairwise.wilcox.test(otus_m$Abundance, otus_m$Sample,p.adjust.method = "bonf", paired = FALSE)

##### alpha diversity #####
alpha <- estimate_richness(physeqobj.meta)

# Make one dataframe
alpha <- cbind(alpha, metadata)

# Average and sd for biological triplates
alpha_averages<- aggregate(x=alpha$InvSimpson,
                          by=list(alpha$Sample),
                          FUN= function(x) c(avg = mean(x), sdev = sd(x)))
alpha_averages<- data.frame(cbind(alpha_averages[,1:1],alpha_averages$x[,1:2]))
colnames(alpha_averages) <- c("Sample", "avg", "sd")

alpha_averages$Sample <- factor(alpha_averages$Sample,                                    # Change ordering manually
                        levels = c("Start", "Blank", "Dead flies", "Living flies"))

for (i in 1:length(alpha_averages$Sample)){
  if (is.na(alpha_averages$sd[i]) == TRUE){
    alpha_averages$sd[i] <- 0
  }
}

alpha_averages$avg <- as.numeric(alpha_averages$avg)
alpha_averages$sd <- as.numeric(alpha_averages$sd)

colours <- c("#6B7AA1","#11324D","#DAB88B", "#C1CFC0")

alpha_plot <- ggplot(data = alpha_averages, aes(x = Sample, y = avg, color = Sample, fill = Sample, legend = FALSE))+
  geom_point(size = 5, alpha = 1)+
  #geom_line(aes(x= Time,y = avg, group = Condition, colour = Condition), size = 2, show.legend = FALSE)+
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd, x = Sample), width=0.02, color="black")+
  #geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=Condition), colour=NA, alpha=0.3 )+
  scale_fill_manual(values = colours)+
  scale_colour_manual(values = colours)+
  labs(x="Condition", y = expression("Diversity index D"[2]*"(A.U.)"), legend = FALSE)+
  theme_bw(base_size = 30, base_family = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background = element_rect(colour = "black", fill = "#EFEFEE"),legend.position = "none")
alpha_plot
# Export plot 
dpi <- 300
png("figures/diversity_alpha_otus.png",width=12*dpi,height=8*dpi,res=dpi)
alpha_plot
dev.off()

# Merge alpha and beta figure
plot_A_B <- ggarrange(alpha_plot, pcoa.plot, labels = c("A", "B"),font.label = list(size = 20, face = "bold", color ="black"),
                          common.legend = TRUE, legend = "right")
plot_A_B
# Export
dpi=300
png("figures/plot_Seq_A_B.png",width=20*dpi,height=7*dpi,res=dpi)
plot_A_B
dev.off()
