setwd("/media/projects1/Fien/Vliegenexperimenten/December")

# Packages
library("Phenoflow") 
library("flowViz")
library("ggplot2")
library("flowAI")
library("cowplot")
library("flowAI")
library("scales") # for plotting time series with ggplot
library("xlsx") # for working with excel files

# Seed
set.seed(777)

# Loading FCM data and metadata: SG
fcsfiles <- list.files(path = "/media/projects1/Fien/Vliegenexperimenten/December/SG", recursive = TRUE, pattern = ".fcs", full.names = TRUE)
flowData <- read.flowSet(files = fcsfiles, transformation = FALSE, emptyValue = FALSE)

# Asinh transformation 
flowData_transformed <- transform(flowData,`BL1-H`= asinh(`BL1-H`), 
                                  `SSC-H`= asinh(`SSC-H`), 
                                  `BL3-H`= asinh(`BL3-H`), 
                                  `FSC-H`= asinh(`FSC-H`),
                                  `SSC-W`= asinh(`SSC-W`),
                                  `FSC-A`= asinh(`FSC-A`))
param <- c("BL1-H","BL3-H","SSC-H","FSC-H", "SSC-W","FSC-A")
remove(flowData)

# Extract metadata from sample names (remove first and last characters)
flowCore::sampleNames(flowData_transformed) <- substring(flowCore::sampleNames(flowData_transformed), 1, nchar(flowCore::sampleNames(flowData_transformed))-4)

# Extract metadata
metadata <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_transformed),"_"), rbind)))
colnames(metadata) <- c("Time", "Condition", "Bio_Replicate","Tech_Replicate")

# Evaluation of the gate T1 
sqrcut <- matrix(c(7.7,7.7,8.6,15,15,4,5.5,7.3,13.6,4), ncol = 2, nrow = 5)
colnames(sqrcut) <- c("BL1-H","BL3-H")
polyGate <- polygonGate(.gate = sqrcut, filterId = "Cells")
xyplot(`BL3-H` ~ `BL1-H`, 
       data = flowData_transformed[1:10], #sample nummer
       filter = polyGate,
       scales = list(y = list(limits = c(3,15)), x = list(limits = c(5,15))),
       axis = axis.default, 
       nbin = 125, 
       par.strip.text = list(col = "white", font = 2, cex = 0.5), 
       smooth = FALSE)

# Extract cell counts --> all cells SG
counts <- flowCore::filter(flowData_transformed, polyGate) 
counts <- toTable(summary(counts))

# Extract volumes (as uL) --> calculate densities all cells 
counts$vol <- as.numeric(flowCore::fsApply(flowData_transformed, FUN = function(x) x@description$`$VOL`))/1000
metadata$density <- (counts$true/counts$vol)*1000


# Adapt levels
metadata$Time <- as.factor(metadata$Time)
levels(metadata$Time) <- c("0","1","2","5","7")
metadata$Time <- as.numeric(as.character(metadata$Time))

# Average cell density and sd for biological triplates
metadata.sum <- aggregate(x=metadata$density,
                          by=list(metadata$Time,metadata$Condition, metadata$Bio_Replicate),
                          FUN= function(x) c(avg = mean(x), sdev = sd(x)))
metadata.sum <- data.frame(cbind(metadata.sum[,1:3],metadata.sum$x[,1:2]))
colnames(metadata.sum) <- c("Time","Condition", "Bio_Replicate","Avg.dens","Sd.dens")

metadata.sum.sum <- aggregate(x=metadata.sum$Avg.dens,
                              by=list(metadata.sum$Time,metadata.sum$Condition),
                              FUN= function(x) c(avg = mean(x), sdev = sd(x)))
metadata.sum.sum <- data.frame(cbind(metadata.sum.sum[,1:2],metadata.sum.sum$x[,1:2]))
colnames(metadata.sum.sum) <- c("Time","Condition","Avg.dens","Sd.dens")

metadata.sum.sum$Condition <- as.factor(metadata.sum.sum$Condition)
levels(metadata.sum.sum$Condition) <- c("Blank","Dead flies","Living flies")

write_xlsx(x = metadata.sum.sum, path = "metadata.sum.sum.xlsx", col_names = TRUE)

# Plot cell densities 
#colours <- c("#D57E7E","#A2CDCD","#C6D57E","#FFD24C")
#colours <- c("#E7E0C9", "#C1CFC0","#6B7AA1", "#11324D")
colours <- c("#6B7AA1","#DAB88B", "#C1CFC0", "#11324D")


plot_1 <- ggplot(data=metadata.sum.sum)+
  geom_point(aes(x=Time,y=Avg.dens, colour = Condition, fill = Condition), shape = 21, size=5, alpha = 1)+
  geom_line(aes(x=Time,y=Avg.dens, group = Condition, colour = Condition), size = 2, show.legend = FALSE)+
  #facet_grid(rows = vars(Place), scales = "free_y")+
  geom_errorbar(aes(ymin=Avg.dens-Sd.dens,ymax=Avg.dens+Sd.dens, x = Time),width=0.2)+
  scale_y_continuous(trans = log10_trans(),labels=scientific,limits = c(10000,1000000)) +
  scale_fill_manual(values = colours)+
  scale_colour_manual(values = colours)+
  labs(fill = "Condition", x="Time (days)", y = "Average density (cells/mL)")+
  theme_bw(base_size = 30, base_family = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background = element_rect(colour = "black", fill = "#EFEFEE"))
plot_1
# Export
dpi=300
png("figures/plot_TotalCellCounts.png",width=12*dpi,height=8*dpi,res=dpi)
plot_1
dev.off()


################################# Fingerprinting ###########################
### Calculate fingerprint with bw = 0.01 
flowData_transformed <- Subset(flowData_transformed, polyGate)
summary <- fsApply(x = flowData_transformed, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
maxval <- max(summary[,"BL1-H"]) #Replace with the column representing the green fluorescence channel (e.g. "FITC-H")
mytrans <- function(x) x/maxval
flowData_transformed <- transform(flowData_transformed,`BL1-H`= mytrans(`BL1-H`), 
                                  `SSC-H`= mytrans(`SSC-H`), 
                                  `BL3-H`= mytrans(`BL3-H`), 
                                  `FSC-H`= mytrans(`FSC-H`),
                                  `SSC-W`= mytrans(`SSC-W`),
                                  `FSC-A`= mytrans(`FSC-A`))

### Randomly resample to the lowest sample size
flowData_transformed <- FCS_resample(flowData_transformed, replace=TRUE)

# Rename
metadata$Condition <- as.factor(metadata$Condition)
levels(metadata$Condition) <- c("Blank","Dead flies","Living flies")
metadata$Time <- as.factor(metadata$Time)
levels(metadata$Time) <- c("D0","D1","D2","D5","D7")

#flowData_transformed_0_7 <- flowData_transformed[metadata$Time == "0" | metadata$Time == "7",]  
#metadata_0_7 <- metadata[metadata$Time == "0" | metadata$Time == "7",]  


fbasis <- flowBasis(flowData_transformed, param, nbin=128, 
                    bw=0.01,normalize=function(x) x)

# Beta-diversity assessment of fingerprint
beta.div <- beta_div_fcm(fbasis, ord.type="PCoA")

# Plot ordination 
plot_beta <- plot_beta_fcm(beta.div, color = metadata$Condition) + 
  scale_fill_manual(values = colours)+
  scale_colour_manual(values = colours)+
  labs(color = "Condition", title = "") +
  geom_text(aes(label=metadata$Time),hjust=1.5, vjust=0)+
  geom_point(size = 5, alpha = 1)+
  scale_shape()+
  coord_fixed(ratio = 1)+
  theme_bw(base_size = 30, base_family = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background = element_rect(colour = "black", fill = "#EFEFEE"))
plot_beta

# Export
dpi=300
png("figures/diversity_RES.png",width=12*dpi,height=8*dpi,res=dpi)
plot_beta
dev.off()

### A AND B FIGURE 
#install.packages("ggpubr")
library(ggpubr)

plot_FCM <- ggarrange(plot_1, plot_beta, labels = c("A", "B"),font.label = list(size = 20, face = "bold", color ="black"),
                   common.legend = TRUE, legend = "right")
plot_FCM
# Export
dpi=300
png("figures/plot_FCM.tif",width=16*dpi,height=7*dpi,res=dpi)
plot_FCM
dev.off()


### STATISTIC on cell counts #### (ppt Inez)
D7 <- metadata[metadata$Time == "D7",]

# Check normality: actually not needed because n > 30
library(ggpubr)
ggqqplot(metadata$density)
chechnormality <- ggqqplot(D7$density)

# Check if the variances are equal --> YES; ONE WAY ANOVA
bartlett.test(density ~ Condition, data = D7)

# One Way Anova
library(dplyr)
density_aov <- aov(density ~ Condition, data = metadata)
summary(density_aov)


### STATISTIC on dissimilarity between groups (distances) and within groups  #######

# Bray-Curtis dissimilarity matrix
BrayCurtis <- vegan::vegdist(fbasis@basis, method = "bray")

# Add total sample name to metadata
metadata$Fullname <- flowCore::sampleNames(flowData_transformed) 

# Checking assumptions PERMANOVA: Check homogeneity of group dispersions
# Between all of the seperate groups (time and condition) --> YES
betadisper <- vegan::betadisper(BrayCurtis, group = metadata$Fullname)
vegan::permutest(betadisper)
plot(betadisper)

# Over time
betadisper <- vegan::betadisper(BrayCurtis, group = metadata$Time)
vegan::permutest(betadisper)
plot(betadisper)

# Between Conditions
betadisper <- vegan::betadisper(BrayCurtis, group = metadata$Condition)
vegan::permutest(betadisper)
plot(betadisper)

### Permanova can be applied
permanova <- vegan::adonis2(BrayCurtis ~ Condition, data = metadata, permutations = permute::how(nperm = 1000))
permanova
str(permanova)

## ANOSIM
#anosim <- vegan::anosim(BrayCurtis, metadata$Condition, distance= "bray", permutations= 1000)
#anosim
#plot(anosim)

# Between Conditions (Blanc, Living flies) --> p = 0.002; R2 = 0.06316
BrayCurtisSubsetA <- vegan::vegdist(fbasis@basis[(metadata$Condition == "Blanc") | (metadata$Condition == "Living flies"),],  method = "bray")
GroupsA <- metadata[metadata$Condition == "Blanc" | metadata$Condition == "Living flies",]
permanovaA <- vegan::adonis2(BrayCurtisSubsetA ~ Condition, data = GroupsA, permutations = permute::how(nperm = 1000))
permanovaA

# Between Conditions (Blanc, Dead flies) --> p = 9.999 e-05; R2 = 0.22737
BrayCurtisSubsetB <- vegan::vegdist(fbasis@basis[(metadata$Condition == "Blanc") | (metadata$Condition == "Dead flies"),],  method = "bray")
GroupsB <- metadata[metadata$Condition == "Blanc" | metadata$Condition == "Dead flies",]
permanovaB <- vegan::adonis2(BrayCurtisSubsetB ~ Condition, data = GroupsB, permutations = permute::how(nperm = 1000))
permanovaB

# Between Conditions (Dead flies, Living flies) --> p = 2e-04; R2 = 0.12734
BrayCurtisSubsetC <- vegan::vegdist(fbasis@basis[(metadata$Condition == "Dead flies") | (metadata$Condition == "Living flies"),],  method = "bray")
GroupsC <- metadata[metadata$Condition == "Dead flies" | metadata$Condition == "Living flies",]
permanovaC <- vegan::adonis2(BrayCurtisSubsetC ~ Condition, data = GroupsC, permutations = permute::how(nperm = 1000))
permanovaC

# OR pairwise test
pairwise_permanova <- pairwiseAdonis::pairwise.adonis(BrayCurtis, metadata$Condition, sim.method= "bray", perm = 1000) 
pairwise_permanova

# Alpha diversity
### Calculate Diversity from normalized fingerprint 
Diversity.fbasis <- Diversity(fbasis,d=3,plot=FALSE, R=999)

# Adapt levels
metadata$Time <- as.factor(metadata$Time)
levels(metadata$Time) <- c("0","1","2","5","7")
metadata$Time <- as.numeric(as.character(metadata$Time))

# Make one dataframe
alpha <- cbind(Diversity.fbasis, metadata)

library(dplyr)
alpha_averages <- alpha %>% group_by(Time, Condition) %>% summarize(avg=mean(D2, na.rm = TRUE),sd=sd(D2, na.rm = TRUE))

alpha <- ggplot(data = alpha_averages, aes(x = Time, y = avg, color = Condition))+
  geom_point(size = 5, alpha = 1)+
  geom_line(aes(x= Time,y = avg, group = Condition, colour = Condition), size = 2, show.legend = FALSE)+
  #geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2, x = metadata$Time), width=0.02, color="black")+
  geom_ribbon(aes(ymin=avg-sd, ymax=avg+sd, fill=Condition), colour=NA, alpha=0.3 )+
  scale_fill_manual(values = colours)+
  scale_colour_manual(values = colours)+
  labs(color = "Condition", x="Time (days)", y = expression("Diversity index D"[2]*"(A.U.)"))+
  theme_bw(base_size = 30, base_family = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background = element_rect(colour = "black", fill = "#EFEFEE"))
alpha

# Export
dpi=300
png("figures/alpha_diversity.png",width=12*dpi,height=8*dpi,res=dpi)
alpha
dev.off()



##################### TOC ################"

setwd("/media/projects1/Fien/Vliegenexperimenten/December")

# Importing and reading the data
library(openxlsx)
library("xlsx")
library(ggplot2)
TOC <- read.xlsx("/media/projects1/Fien/Vliegenexperimenten/December/TOC_graph.xlsx",sheet = 1,colNames=TRUE)

# Average ATP and sd for biological triplates
metadata.sum <- aggregate(x=TOC$TOC,
                          by=list(TOC$Time, TOC$Condition),
                          FUN= function(x) c(avg = mean(x), sdev = sd(x)))
metadata.sum <- data.frame(cbind(metadata.sum[,1:2],metadata.sum$x[,1:2]))
colnames(metadata.sum) <- c("Time","Condition","Avg.TOC","Sd.TOC")

# Adapt levels
metadata.sum$Time <- as.factor(metadata.sum$Time)
levels(metadata.sum$Time) <- c("0","1","2","5","7")
metadata.sum$Time <- as.numeric(as.character(metadata.sum$Time))
metadata.sum$Condition <- as.factor(metadata.sum$Condition)
levels(metadata.sum$Condition) <- c("Blank","Dead flies","Living flies")

plot_TOC <- ggplot(data=metadata.sum)+
  geom_point(aes(x=Time,y=Avg.TOC, colour = Condition, fill = Condition), shape = 21, size=5, alpha=1)+
  geom_line(aes(x=Time,y=Avg.TOC, group = Condition, colour = Condition), size = 2, show.legend = FALSE)+
  #facet_grid(rows = vars(Place), cols = vars(BiolRepl), scales = "fixed")+
  geom_errorbar(aes(ymin=Avg.TOC-Sd.TOC,ymax=Avg.TOC+Sd.TOC, x = Time),width=0.2)+
  #ylim(0,0.12)+
  #scale_y_continuous(trans = log10_trans(),labels=scientific) +
  scale_fill_manual(values = colours)+
  scale_colour_manual(values = colours)+
  labs(colour = "Condition", x="Time (days)", y = "TOC (pbb)")+
  theme_bw(base_size = 30, base_family = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background = element_rect(colour = "black", fill = "#EFEFEE"))
plot_TOC
# Export
dpi=300
png("figures/TOC.png",width=12*dpi,height=8*dpi,res=dpi)
plot_TOC
dev.off()

write.xlsx(metadata.sum, "metadata.sum.sum_TOC.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

#write.xlsx(x = metadata.sum.sum, path = "metadata.sum.sum_TOC.xlsx", col_names = TRUE)

############################### ATP ANALYSIS EXPERIMENT DECEMBER ###########################
setwd("/media/projects1/Fien/Vliegenexperimenten/December")

# Importing and reading the data
library(openxlsx)
library(ggplot2)
ATP <- read.xlsx("/media/projects1/Fien/Vliegenexperimenten/December/ATP/ATP_flyexp3.xlsx",sheet = 1,colNames=TRUE)

# Average ATP and sd for biological triplates
metadata.sum <- aggregate(x=ATP$ATP,
                          by=list(ATP$Time, ATP$Condition, ATP$BiolRepl),
                          FUN= function(x) c(avg = mean(x), sdev = sd(x)))
metadata.sum <- data.frame(cbind(metadata.sum[,1:3],metadata.sum$x[,1:2]))
colnames(metadata.sum) <- c("Time","Condition","BiolRepl","Avg.ATP","Sd.ATP")

# Average ATP and sd for biological triplates
metadata.sum.sum <- aggregate(x=metadata.sum$Avg.ATP,
                              by=list(metadata.sum$Time, metadata.sum$Condition),
                              FUN= function(x) c(avg = mean(x), sdev = sd(x)))
metadata.sum.sum <- data.frame(cbind(metadata.sum.sum[,1:2],metadata.sum.sum$x[,1:2]))
colnames(metadata.sum.sum) <- c("Time","Condition","Avg.ATP","Sd.ATP")

metadata.sum.sum$Condition <- as.factor(metadata.sum.sum$Condition)
levels(metadata.sum.sum$Condition) <- c("Blank","Dead flies","Living flies")

plot_ATP <- ggplot(data=metadata.sum.sum)+
  geom_point(aes(x=Time,y=Avg.ATP, colour = Condition, fill = Condition), shape = 21, size=5, alpha=1)+
  geom_line(aes(x=Time,y=Avg.ATP, group = Condition, colour = Condition), size = 2, show.legend = FALSE)+
  geom_errorbar(aes(ymin=Avg.ATP-Sd.ATP,ymax=Avg.ATP+Sd.ATP, x = Time),width=0.2)+
  #ylim(0,0.12)+
  #scale_y_continuous(trans = log10_trans(),labels=scientific) +
  scale_fill_manual(values = colours)+
  scale_colour_manual(values = colours)+
  labs(colour = "Condition", x="Time (days)", y = "ATP (ng/L)")+
  theme_bw(base_size = 30, base_family = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background = element_rect(colour = "black", fill = "#EFEFEE"))
plot_ATP
# Export
dpi=300
png("figures/ATP_sum.png",width=12*dpi,height=8*dpi,res=dpi)
plot_ATP
dev.off()

write_xlsx(x = metadata.sum.sum, path = "metadata.sum.sum_ATP.xlsx", col_names = TRUE)

# Merge ATP and TOC figure
plot_TOC_ATP <- ggarrange(plot_TOC, plot_ATP, labels = c("A", "B"),font.label = list(size = 20, face = "bold", color ="black"),
                      common.legend = TRUE, legend = "right")
plot_TOC_ATP
# Export
dpi=300
png("figures/plot_TOC_ATP.png",width=16*dpi,height=7*dpi,res=dpi)
plot_TOC_ATP
dev.off()
