rm(list=ls())
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4955027/pdf/f1000research-5-10545.pdf

install.packages("pacman")
library("pacman")
pacman::p_load(ggplot2,dada2,phyloseq,gridExtra,grid,lattice,
               ggrepel,DECIPHER,tidyverse,ggforce,phangorn,ggpubr)

# Use the color-blind-friendly palette from http://www.cookbook-r.com
theme_set(theme_bw()) 
color_pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

setwd("D:/1research/5 IZZY RYGB/Reanalysis_200916/Phyloseq")
getwd()

set.seed(100)

###Finish dada2###

### Step 1
#Construct phylogenetic tree

# Open the RData image from Dada2 to get sequences counts "seqtab.nochim", or
# you can also import from excel use these two lines
seqtab.nochim <- read.csv(file.choose(), header=T, row.names = 1)
seqtab.nochim <- as.matrix(seqtab.nochim)

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

save.image("16S_phyloseq_1.Rdata")

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order 
fit = pml(treeNJ, data=phang.align) 

# negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2) 
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
# This step takes very long time (several hours!)
detach("package:phangorn", unload=TRUE)

save.image("16S_phyloseq_2.Rdata")

### Step 2
# Combine data into a phyloseq object

# Import metadata
samdf <- read.csv(file.choose(), header=T) # The first Column is named SampleID
rownames(samdf)=samdf$SampleID #Use the SampleID column as the row names

all(rownames(seqtab.nochim) %in% rownames(samdf)) # check if data matches metadata

### Step 3
#Silva taxonomy file

taxa <- read.csv(file.choose(), header=T, row.names=1) # use corrected taxa (combine, change NA names)
taxa <- as.matrix(taxa)

ps_silva_raw <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                         sample_data(samdf), tax_table(taxa),phy_tree(fitGTR$tree))
ps_silva_raw

ps_silva = ps_silva_raw

#ps_silva <- prune_samples(sample_names(ps_silva) != "SHAMMIX", ps_silva)
#ps_silva

#correctedsamplenames = read.csv(file.choose(), header=F)
#correctedsamplenames = correctedsamplenames[,1]
#sample_names(ps_silva) <- correctedsamplenames
#sample_names(ps_silva)

#ps_silva@sam_data[["SampleID"]] = correctedsamplenames
#ps_silva@sam_data[["SampleID"]]

# Remove some samples with extremely few reads (eg. <1000)
# Method 1
# ps_silva <- prune_samples(sample_sums(ps_silva) >= 1000, ps_silva)

# Method 2
# ps_silva <- prune_samples(sample_names(ps_silva) != "MNS5", ps_silva)
# ps_silva <- subset_samples(ps_silva, Age != "2")

# Method 3
# outliers <- c("DUC2", "F6D165", "M3D175", "M4D175", "M5D175", "M6D175")
# outliers <- read.csv(file.choose(), header=F)
# outliers <- as.vector(as.matrix(outliers))
# ps_silva <- prune_samples(!(sample_names(ps_silva) %in% outliers), ps_silva)
# ps_silva

min(sample_sums(ps_silva)) # what is the minimum sequencing depth?

n_distinct(sample_data(ps_silva)$SampleID) # how many samples are there?

# remove control group which are not used in this paper
ps_silva <- ps_silva %>%
  subset_samples(Class %in% c(1:14))

### Step 4
# Taxonomic Filtering

# Show available ranks in the dataset
rank_names(ps_silva)

# Create table, number of features for each phyla
table(tax_table(ps_silva)[, "Phylum"], exclude = NULL)

# delete ambiguous phylum annotation
ps_silva <- subset_taxa(ps_silva, !is.na(Phylum) & !Phylum %in% c("Bacteria_unclassified", "Eukaryota_unclassified"))
ps_silva

# Prevalence control
ps_silva <- filter_taxa(ps_silva, function(x) sum(x>3)>0.065*nsamples(ps_silva), TRUE) 
ps_silva

# Method 1: Julie's criteria, a taxa should have at least 10 reads in one sample
#ps_silva <- filter_taxa(ps_silva, function(x) sum(x>9)>0, TRUE) 
#ps_silva

# Method 2: any taxa at least with non zero reads in 5% of samples
# ps_silva <- filter_taxa(ps_silva, function(x) sum(x>0)>0.05*nsamples(ps_silva), TRUE)
# ps_silva

# Method 3: any taxa with at least 4 counts in at lease 10% of samples
# ps_silva <- filter_taxa(ps_silva, function(x) sum(x>3)>0.1*nsamples(ps_silva), TRUE)
# ps_silva

# After excluding samples in phyloseq, you need to remove taxa that are
# no longer present in any of remaining samples.
ps_silva <- ps_silva %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE) # keep only taxa with more than 0 total counts
ps_silva

save.image("16S_phyloseq_5.Rdata")

# Step 4
# Transform data to proportions
ps_silva.prop <- transform_sample_counts(ps_nocontrol, function(otu) otu/sum(otu))

## Visualize alpha-diversity. 
# It is suggested to use untrimmed taxa (with singletons and rare taxa)
# Use raw reads to plot
# alpha_methods = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")

plot_richness(ps_raw_nocontrol, x="Group_Day", measures=c("Observed","Shannon","Simpson"), color="Class")+
  geom_boxplot()
# scale_x_discrete(limits=c("F-1", "F0", "F1C", "F1E", "F1R", "F1S", "F5C", "F5E", "F5R", "F5S", "F10C", "F10E", "F10R", "F10S", "DUE", "DUR", "DUS", "JEC", "JEE", "JER", "JES", "ILC", "ILE", "ILR", "ILS", "CEE", "CER", "CES", "COC", "COE", "COR", "COS", "MNC", "MNE", "MNR", "MNS", "HR", "HS1", "HS2"))
# scale_x_discrete is used to re-order the x axis

# Use rarefied reads to plot
ps_rar <- rarefy_even_depth(ps_raw_nocontrol, 
                            sample.size = min(sample_sums(ps_raw_nocontrol)), 
                            rngseed = 935, # set this so that subsampling is reproducible
                            replace = F, 
                            trimOTUs = T)

plot_richness(ps_rar, x="Group_Day", measures=c("Observed","Shannon"), color="Class")+
  geom_boxplot()
#scale_x_discrete(limits=c("F-1", "F0", "F1C", "F1E", "F1R", "F1S", "F5C", "F5E", "F5R", "F5S", "F10C", "F10E", "F10R", "F10S", "DUE", "DUR", "DUS", "JEC", "JEE", "JER", "JES", "ILC", "ILE", "ILR", "ILS", "CEE", "CER", "CES", "COC", "COE", "COR", "COS", "MNC", "MNE", "MNR", "MNS", "HR", "HS1", "HS2"))

alpha_dataframe = estimate_richness(ps_rar)
write.csv(alpha_dataframe, "alpha.csv")

# line plot: not ready to use
p<- ggplot(df2, aes(x=dose, y=len, group=supp, color=supp)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                position=position_dodge(0.05))+
  labs(title="Tooth length per dose", x="Dose (mg)", y = "Length")+
  scale_color_manual(values=c('#999999','#E69F00'))






#1 PCoA with weighted Unifrac - All samples
ordu = ordinate(ps_nocontrol, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_nocontrol, ordu, color="Day", shape="Day") + 
  scale_shape_manual(values=c(0,1,9,2,5:8,0,1,9,2,5:8))+
  scale_color_manual(values = c(rep("blue",8),rep("red",8)))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  # ggtitle("MDS/PCoA on weighted-UniFrac distance")+
  # guides(fill = guide_legend(override.aes = aes(label = "")))+
  # guides(fill = guide_legend(override.aes = list(shape = 21)))+
  stat_ellipse(geom = "polygon", type="norm", alpha=0.4,aes(fill = Day))
  scale_fill_manual(values = c("blue","red","yellow","orange","gray","cyan","green","purple","black"))+
  ggsave("wunifrac_all.png", height=5, width=5, units='in', dpi=600)



  #1 PCoA with weighted Unifrac - All samples
  ordu = ordinate(ps_nocontrol, "PCoA", "unifrac", weighted=TRUE)
  plot_ordination(ps_nocontrol, ordu, color="Day", label = "NewID") +
    scale_color_manual(values = c(rep(c("blue","red","yellow","orange","gray","cyan","green"),1),"blue","red"))+
    #scale_shape_manual(values=c(0,1,9,2,5:8,0,1,9,2,5:8))+
    geom_point(size=4)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title = element_blank()) +
    #scale_fill_discrete(breaks=c("1", "3", "6", "9", "16", "pre-atb", "post-atb", "RYGB_donor","Sham_donor"))+
    stat_ellipse(geom = "polygon", type="norm", alpha=0.2, aes(fill = Day),linetype=2)+
    scale_fill_manual("",breaks=c("1", "3", "6", "9", "16", "pre-atb", "post-atb", "RYGB_donor","Sham_donor"),
                      values = c(rep(c("blue","red","yellow","orange","gray","cyan","green"),2)),
                      labels= c("D1", "D3", "D6", "D9", "D16", "pre-atb", "post-atb", "RYGB_donor","Sham_donor"))
                      
  
  
#2 PCoA with weighted Unifrac - All samples excluding post ATB
ps_expostATB <- subset_samples(ps_nocontrol, Class != 2)
ps_expostATB <- ps_expostATB %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE) # keep only taxa with more than 0 total counts
ps_expostATB

ordu = ordinate(ps_expostATB, "PCoA", "wunifrac")
plot_ordination(ps_expostATB, ordu, color="Day", label = "NewID") + 
  scale_color_manual(values = c("orange", "yellow", "gray", "green","cyan","black","pink", "purple"))+
  geom_point(size=2) + guides(color=FALSE) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.position = "none") +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.2, aes(fill = Day), linetype=2)+
  scale_fill_manual(breaks=c("pre-atb","1", "3", "6", "9", "16", "RYGB_donor","Sham_donor"),
                    values = c("black", "orange", "yellow", "gray", "green", "cyan", "pink", "purple"),
                    labels= c("pre-atb","D1", "D3", "D6", "D9", "D16", "RYGB_donor","Sham_donor"))+
  ggsave("PCoA_wunifrac_all_expostATB_label.png", height=5, width=5, units='in', dpi=600)

plot_ordination(ps_expostATB, ordu, color="Day") + 
  scale_color_manual(values = c("orange", "yellow", "gray", "green","cyan","black","pink", "purple"))+
  geom_point(size=2)+ guides(color=FALSE, shape=FALSE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.position = "none") +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.2, aes(fill = Day),linetype=2)+
  scale_fill_manual(breaks=c("pre-atb","1", "3", "6", "9", "16", "RYGB_donor","Sham_donor"),
                    values = c("black", "orange", "yellow", "gray", "green", "cyan", "pink", "purple"),
                    labels= c("pre-atb","D1", "D3", "D6", "D9", "D16", "RYGB_donor","Sham_donor"))+
  ggsave("PCoA_wunifrac_all_expostATB.png", height=5, width=5, units='in', dpi=600)


#2 PCoA with weighted Unifrac - pre-ATB
ps_preATB <- subset_samples(ps_nocontrol,  Class %in% c(1,13,14))
ps_preATB <- ps_preATB %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE) # keep only taxa with more than 0 total counts
ps_preATB

ordu = ordinate(ps_preATB, "PCoA", "wunifrac")
plot_ordination(ps_preATB, ordu, color="Group", label = "NewID") + 
  scale_color_manual(values = c("pink","black", "purple","black"))+
  geom_point(size=3) + guides(color=FALSE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.position = "none") +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.2, aes(fill = Day), linetype=2)+
  scale_fill_manual(breaks= c("RYGB_donor","RYGB_recipient","Sham_donor","Sham_recipient"),
                    values= c("pink","black", "purple","black"),
                    labels= c("RYGB_donor","RYGB_recipient","Sham_donor","Sham_recipient"))+
  ggsave("PCoA_wunifrac_preATB_label.png", height=5, width=5, units='in', dpi=600)

plot_ordination(ps_preATB, ordu, color="Day") + 
  scale_color_manual(values = c("orange", "yellow", "gray", "green","cyan","black","pink", "purple"))+
  geom_point(size=2)+ guides(color=FALSE, shape=FALSE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.position = "none") +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.2, aes(fill = Day),linetype=2)+
  scale_fill_manual(breaks=c("pre-atb","1", "3", "6", "9", "16", "RYGB_donor","Sham_donor"),
                    values = c("black", "orange", "yellow", "gray", "green", "cyan", "pink", "purple"),
                    labels= c("pre-atb","D1", "D3", "D6", "D9", "D16", "RYGB_donor","Sham_donor"))+
  ggsave("PCoA_wunifrac_all_expostATB.png", height=5, width=5, units='in', dpi=600)












































#3 PCoA with weighted Unifrac - RYGB samples
ps_RYGB <- subset_samples(ps_nocontrol, Group %in% c("RYGB_donor", "RYGB_recipient"))
ps_RYGB <- ps_RYGB %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE) # keep only taxa with more than 0 total counts
ps_RYGB

ordu = ordinate(ps_RYGB, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_RYGB, ordu, color="Group_Day", shape="Group_Day") + 
  scale_shape_manual(values=c(0,1,9,2,5:8))+
  scale_color_manual(values = c(rep("blue",8)))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  ggsave("wunifrac_RYGB.png", height=5, width=5, units='in', dpi=600)

plot_ordination(ps_RYGB, ordu, color="Group_Day", label= "NewID", shape="Group_Day") + 
  scale_shape_manual(values=c(0,1,9,2,5:8))+
  scale_color_manual(values = c(rep("blue",8)))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  ggsave("wunifrac_RYGB_label.png", height=5, width=5, units='in', dpi=600)

#4 PCoA with weighted Unifrac - RYGB samples exclude post ATB
ps_RYGB_expostATB <- subset_samples(ps_RYGB,  Class != 2)
ps_RYGB_expostATB <- ps_RYGB_expostATB %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE) # keep only taxa with more than 0 total counts
ps_RYGB_expostATB

ordu = ordinate(ps_RYGB_expostATB, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_RYGB_expostATB, ordu, color="Group_Day", shape="Group_Day") + 
  scale_shape_manual(values=c(0,1,2,5:8))+
  scale_color_manual(values = c(rep("blue",7)))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  ggsave("wunifrac_RYGB_expostATB.png", height=5, width=5, units='in', dpi=600)

plot_ordination(ps_RYGB_expostATB, ordu, color="Group_Day", label= "NewID", shape="Group_Day") + 
  scale_shape_manual(values=c(0,1,2,5:8))+
  scale_color_manual(values = c(rep("blue",7)))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  ggsave("wunifrac_RYGB_expostATB_label.png", height=5, width=5, units='in', dpi=600)

#5 PCoA with weighted Unifrac - SHAM samples
ps_SHAM <- subset_samples(ps_nocontrol, Group %in% c("Sham_donor", "Sham_recipient"))
ps_SHAM <- ps_SHAM %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE) # keep only taxa with more than 0 total counts
ps_SHAM

ordu = ordinate(ps_SHAM, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_SHAM, ordu, color="Group_Day", shape="Group_Day") + 
  scale_shape_manual(values=c(0,1,9,2,5:8))+
  scale_color_manual(values = c(rep("red",8)))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_SHAM.png", height=5, width=5, units='in', dpi=600)

plot_ordination(ps_SHAM, ordu, color="Group_Day", label= "NewID", shape="Group_Day") + 
  scale_shape_manual(values=c(0,1,9,2,5:8))+
  scale_color_manual(values = c(rep("red",8)))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_SHAM_label.png", height=5, width=5, units='in', dpi=600)

#6 PCoA with weighted Unifrac - SHAM samples exclude post ATB
ps_SHAM_expostATB <- subset_samples(ps_SHAM,  Class != 2)
ps_SHAM_expostATB <- ps_SHAM_expostATB %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE) # keep only taxa with more than 0 total counts
ps_SHAM_expostATB

ordu = ordinate(ps_SHAM_expostATB, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_SHAM_expostATB, ordu, color="Group_Day", shape="Group_Day") + 
  scale_shape_manual(values=c(0,1,2,5:8))+
  scale_color_manual(values = c(rep("red",7)))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_SHAM_expostATB.png", height=5, width=5, units='in', dpi=600)

plot_ordination(ps_SHAM_expostATB, ordu, color="Group_Day", label= "NewID", shape="Group_Day") + 
  scale_shape_manual(values=c(0,1,2,5:8))+
  scale_color_manual(values = c(rep("red",7)))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_SHAM_expostATB_label.png", height=5, width=5, units='in', dpi=600)

#7 PCoA with weighted Unifrac - preATB (Donor, RYGBr SHAMr)
ps_preATB <- subset_samples(ps_nocontrol,  Class %in% c(1,13,14))
ps_preATB <- ps_preATB %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE) # keep only taxa with more than 0 total counts
ps_preATB

ordu = ordinate(ps_preATB, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_preATB, ordu, color="Group_Day", shape="Group_Day") + 
  scale_shape_manual(values=c(15,16,15,16))+
  scale_color_manual(values = c("blue","blue","red","red")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_preATB.png", height=5, width=5, units='in', dpi=600)

plot_ordination(ps_preATB, ordu, color="Group_Day", label= "NewID", shape="Group_Day") + 
  scale_shape_manual(values=c(15,16,15,16))+
  scale_color_manual(values = c("blue","blue","red","red")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_preATB_label.png", height=5, width=5, units='in', dpi=600)

#8 PCoA with weighted Unifrac - postATB (Donor, RYGBr SHAMr)
ps_postATB <- subset_samples(ps_nocontrol,  Class %in% c(2,13,14))
ps_postATB <- ps_postATB %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE) # keep only taxa with more than 0 total counts
ps_postATB

ordu = ordinate(ps_postATB, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_postATB, ordu, color="Group_Day", shape="Group_Day") + 
  scale_shape_manual(values=c(15,16,15,16))+
  scale_color_manual(values = c("blue","blue","red","red")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_postATB.png", height=5, width=5, units='in', dpi=600)

plot_ordination(ps_postATB, ordu, color="Group_Day", label= "NewID", shape="Group_Day") + 
  scale_shape_manual(values=c(15,16,15,16))+
  scale_color_manual(values = c("blue","blue","red","red")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_postATB_label.png", height=5, width=5, units='in', dpi=600)

#9 PCoA with weighted Unifrac - D1 (Donor, RYGBr SHAMr)
ps_D1 <- subset_samples(ps_nocontrol,  Class %in% c(3,4,13,14))
ps_D1 <- ps_D1 %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE) # keep only taxa with more than 0 total counts
ps_D1

ordu = ordinate(ps_D1, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_D1, ordu, color="Group_Day", shape="Group_Day") + 
  scale_shape_manual(values=c(15,16,15,16))+
  scale_color_manual(values = c("blue","blue","red","red")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_D1.png", height=5, width=5, units='in', dpi=600)

plot_ordination(ps_D1, ordu, color="Group_Day", label= "NewID", shape="Group_Day") + 
  scale_shape_manual(values=c(15,16,15,16))+
  scale_color_manual(values = c("blue","blue","red","red")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_D1_label.png", height=5, width=5, units='in', dpi=600)

#10 PCoA with weighted Unifrac - D3 (Donor, RYGBr SHAMr)
ps_D3 <- subset_samples(ps_nocontrol,  Class %in% c(5,6,13,14))
ps_D3 <- ps_D3 %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE) # keep only taxa with more than 0 total counts
ps_D3

ordu = ordinate(ps_D3, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_D3, ordu, color="Group_Day", shape="Group_Day") + 
  scale_shape_manual(values=c(15,16,15,16))+
  scale_color_manual(values = c("blue","blue","red","red")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_D3.png", height=5, width=5, units='in', dpi=600)

plot_ordination(ps_D3, ordu, color="Group_Day", label= "NewID", shape="Group_Day") + 
  scale_shape_manual(values=c(15,16,15,16))+
  scale_color_manual(values = c("blue","blue","red","red")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_D3_label.png", height=5, width=5, units='in', dpi=600)

#11 PCoA with weighted Unifrac - D6 (Donor, RYGBr SHAMr)
ps_D6 <- subset_samples(ps_nocontrol,  Class %in% c(7,8,13,14))
ps_D6 <- ps_D6 %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE) # keep only taxa with more than 0 total counts
ps_D6

ordu = ordinate(ps_D6, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_D6, ordu, color="Group_Day", shape="Group_Day") + 
  scale_shape_manual(values=c(15,16,15,16))+
  scale_color_manual(values = c("blue","blue","red","red")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_D6.png", height=5, width=5, units='in', dpi=600)

plot_ordination(ps_D6, ordu, color="Group_Day", label= "NewID", shape="Group_Day") + 
  scale_shape_manual(values=c(15,16,15,16))+
  scale_color_manual(values = c("blue","blue","red","red")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_D6_label.png", height=5, width=5, units='in', dpi=600)

#12 PCoA with weighted Unifrac - D9 (Donor, RYGBr SHAMr)
ps_D9 <- subset_samples(ps_nocontrol,  Class %in% c(9,10,13,14))
ps_D9 <- ps_D9 %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE) # keep only taxa with more than 0 total counts
ps_D9

ordu = ordinate(ps_D9, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_D9, ordu, color="Group_Day", shape="Group_Day") + 
  scale_shape_manual(values=c(15,16,15,16))+
  scale_color_manual(values = c("blue","blue","red","red")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_D9.png", height=5, width=5, units='in', dpi=600)

plot_ordination(ps_D9, ordu, color="Group_Day", label= "NewID", shape="Group_Day") + 
  scale_shape_manual(values=c(15,16,15,16))+
  scale_color_manual(values = c("blue","blue","red","red")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_D9_label.png", height=5, width=5, units='in', dpi=600)

#12 PCoA with weighted Unifrac - D16 (Donor, RYGBr SHAMr)
ps_D16 <- subset_samples(ps_nocontrol,  Class %in% c(11,12,13,14))
ps_D16 <- ps_D16 %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE) # keep only taxa with more than 0 total counts
ps_D16

ordu = ordinate(ps_D16, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_D16, ordu, color="Group_Day", shape="Group_Day") + 
  scale_shape_manual(values=c(15,16,15,16))+
  scale_color_manual(values = c("blue","blue","red","red")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_D16.png", height=5, width=5, units='in', dpi=600)

plot_ordination(ps_D16, ordu, color="Group_Day", label= "NewID", shape="Group_Day") + 
  scale_shape_manual(values=c(15,16,15,16))+
  scale_color_manual(values = c("blue","blue","red","red")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()) + 
  ggsave("wunifrac_D16_label.png", height=5, width=5, units='in', dpi=600)

save.image("16S_phyloseq_5.Rdata")



library(readr)
library(tidyverse)
library(dplyr)
source("D:/2Learning/R/1MyScripts/ancom_v2.1.R")

# ANCOM_ASV pre vs post
ps_pre_post_ANCOM_ASV <- subset_samples(ps_nocontrol,  Class %in% c(1,2))
ps_pre_post_ANCOM_ASV <- ps_pre_post_ANCOM_ASV %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE) # keep only taxa with more than 0 total counts
ps_pre_post_ANCOM_ASV

otu_data=t(as.data.frame(otu_table(ps_pre_post_ANCOM_ASV)))
taxa_data=as.data.frame(tax_table(ps_nocontrol))
taxa_data$taxa_id=rownames(taxa_data)
meta_data=as.matrix(sample_data(ps_pre_post_ANCOM_ASV))

# Step 1: Data preprocessing
# Change SampleID column name in "sample_var"
feature_table = otu_data; sample_var = "SampleID"; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data2 = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM
# Change group column name in "main_var"
main_var = "Group"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL
res = ANCOM(feature_table, meta_data2, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
mean_difference = res$fig$data$x
data_out = left_join(cbind(res$out, mean_difference), taxa_data, by = "taxa_id")
data_out = data_out[order(data_out$W, decreasing = TRUE),]
write_csv(data_out, "ANCOM_ASV pre vs post.csv")

# Step 3: Volcano Plot
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig = res$fig +  
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)+
  ggsave("ANCOM_ASV pre vs post.png", height=5, width=5, units='in', dpi=600)
fig


# ANCOM_Family pre vs post
ps_pre_post_ANCOM_Family <- subset_samples(ps_nocontrol,  Class %in% c(1,2))
ps_pre_post_ANCOM_Family <- ps_pre_post_ANCOM_Family %>%
  tax_glom("Family")
ps_pre_post_ANCOM_Family <- ps_pre_post_ANCOM_Family %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE) # keep only taxa with more than 0 total counts
ps_pre_post_ANCOM_Family

otu_data=t(as.data.frame(otu_table(ps_pre_post_ANCOM_Family)))
taxa_data=as.data.frame(tax_table(ps_nocontrol))
taxa_data$taxa_id=rownames(taxa_data)
meta_data=as.matrix(sample_data(ps_pre_post_ANCOM_Family))

# Step 1: Data preprocessing
feature_table = otu_data; sample_var = "SampleID"; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data2 = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM
# Change group column name in "main_var"
main_var = "Group"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL
res = ANCOM(feature_table, meta_data2, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
mean_difference = res$fig$data$x
data_out = left_join(cbind(res$out, mean_difference), taxa_data, by = "taxa_id")
data_out = data_out[order(data_out$W, decreasing = TRUE),]
write_csv(data_out, "ANCOM_Family pre vs post.csv")

# Step 3: Volcano Plot
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig = res$fig +  
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)+
  ggsave("ANCOM_Family pre vs post.png", height=5, width=5, units='in', dpi=600)
fig













# PCoA of Faecal samples
ps_faeces<- subset_samples(ps_silva, Type != "T")
out.pcoa.logfaeces <- ordinate(ps_faeces,  method = "MDS", distance = "bray")
evalsfaeces <- out.pcoa.logfaeces$values[,1]
plot_ordination(ps_faeces, out.pcoa.logfaeces, shape = "ExpGroup") +
  labs(col = "Group", shape = "ExpGroup") +
  scale_shape_manual(values=1:9) +
  coord_fixed(sqrt(evalsfaeces[2] / evalsfaeces[1])) + 
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = factor(Group), label = factor(SampleID)))

# Remove F2 and human samples
ps_faecesR_ex_F2_H<- subset_samples(ps_faeces, Group != "F0")
ps_faecesR_ex_F2_H<- subset_samples(ps_faecesR_ex_F2_H, Group != "HR")
ps_faecesR_ex_F2_H<- subset_samples(ps_faecesR_ex_F2_H, Group != "HS1")
ps_faecesR_ex_F2_H<- subset_samples(ps_faecesR_ex_F2_H, Group != "HS2")

out.pcoa.logfaeces_ex_F2_H <- ordinate(ps_faecesR_ex_F2_H,  method = "MDS", distance = "bray")
evalsfaecesR_F2_H <- out.pcoa.logfaeces_ex_F2_H$values[,1]
plot_ordination(ps_faecesR_ex_F2_H, out.pcoa.logfaeces_ex_F2_H, shape = "ExpGroup") +
  labs(col = "Group", shape = "ExpGroup") +
  scale_shape_manual(values=1:9) +
  coord_fixed(sqrt(evalsfaecesR_F2_H[2] / evalsfaecesR_F2_H[1])) + 
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = Group, label = SampleID))

# weighted Unifrac PCoA of faecal samples after removing F2 and human samples
ordu = ordinate(ps_faecesR_ex_F2_H, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_faecesR_ex_F2_H, ordu, shape="ExpGroup") + 
  scale_shape_manual(values=1:9) +
  labs(col = "Group", shape = "ExpGroup") +
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = Group, label = SampleID))

# PCoA of F1
ps_F1<- ps_silva %>%
  subset_samples(Group %in% c("F1C", "F1R", "F1S", "F1E")) 
ps_F1

out.pcoa.logF1 <- ordinate(ps_F1,  method = "MDS", distance = "bray")
evalsF1 <- out.pcoa.logF1$values[,1]
plot_ordination(ps_F1, out.pcoa.logF1, shape = "ExpGroup") +
  labs(col = "Group", shape = "ExpGroup") +
  scale_shape_manual(values=1:9) +
  coord_fixed(sqrt(evalsF1[2] / evalsF1[1])) + 
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = Group, label = SampleID))

# weighted Unifrac PCoA F1
ordu = ordinate(ps_F1, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_F1, ordu, shape="ExpGroup") + 
  scale_shape_manual(values=1:9) +
  labs(col = "Group", shape = "ExpGroup") +
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = Group, label = SampleID))

# PCoA of F5
ps_F5 <- ps_silva %>%
  subset_samples(Group %in% c("F5C", "F5R", "F5S", "F5E")) 

out.pcoa.logF5 <- ordinate(ps_F5,  method = "MDS", distance = "bray")
evalsF5 <- out.pcoa.logF5$values[,1]
plot_ordination(ps_F5, out.pcoa.logF5, shape = "ExpGroup") +
  labs(col = "Group", shape = "ExpGroup") +
  scale_shape_manual(values=1:9) +
  coord_fixed(sqrt(evalsF5[2] / evalsF5[1])) + 
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = Group, label = SampleID))

# weighted Unifrac PCoA F5
ordu = ordinate(ps_F5, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_F5, ordu, shape="ExpGroup") + 
  scale_shape_manual(values=1:9) +
  labs(col = "Group", shape = "ExpGroup") +
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = Group, label = SampleID))

# PCoA of F10
ps_F10 <- ps_silva %>%
  subset_samples(Group %in% c("F10C", "F10R", "F10S", "F10E")) 

out.pcoa.logF10 <- ordinate(ps_F10,  method = "MDS", distance = "bray")
evalsF10 <- out.pcoa.logF10$values[,1]
plot_ordination(ps_F10, out.pcoa.logF10, shape = "ExpGroup") +
  labs(col = "Group", shape = "ExpGroup") +
  scale_shape_manual(values=1:9) +
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = Group, label = SampleID))

# weighted Unifrac PCoA F10
ordu = ordinate(ps_F10, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_F10, ordu, shape="ExpGroup") + 
  scale_shape_manual(values=1:9) +
  labs(col = "Group", shape = "ExpGroup") +
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = Group, label = SampleID))

# PCoA of tissue samples
ps_tissue<- subset_samples(ps_silva, Type != "F")
out.pcoa.logtissue <- ordinate(ps_tissue,  method = "MDS", distance = "bray")
evalstissue <- out.pcoa.logtissue$values[,1]
plot_ordination(ps_tissue, out.pcoa.logtissue, shape = "ExpGroup") +
  labs(col = "Group", shape = "ExpGroup") +
  scale_shape_manual(values=1:9) +
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = Group, label = SampleID))

# weighted Unifrac PCoA of tissue samples
ordu = ordinate(ps_tissue, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_tissue, ordu, shape="ExpGroup") + 
  scale_shape_manual(values=1:9) +
  labs(col = "Group", shape = "ExpGroup") +
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = Group, label = SampleID))

# PCoA of colon
ps_colon <- ps_silva %>%
  subset_samples(Group %in% c("COC", "COR", "COS", "COE")) 

out.pcoa.logcolon <- ordinate(ps_colon,  method = "MDS", distance = "bray")
evalscolon <- out.pcoa.logcolon$values[,1]
plot_ordination(ps_colon, out.pcoa.logcolon, shape = "ExpGroup") +
  labs(col = "Group", shape = "ExpGroup") +
  scale_shape_manual(values=1:4) +
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = Group, label = SampleID))

# weighted Unifrac PCoA colon
ordu = ordinate(ps_colon, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_colon, ordu, shape="ExpGroup") + 
  scale_shape_manual(values=1:4) +
  labs(col = "Group", shape = "ExpGroup") +
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = Group, label = SampleID))

# PCoA of mucus
ps_mucus <- ps_silva %>%
  subset_samples(Group %in% c("MNC", "MNR", "MNS", "MNE")) 

out.pcoa.logmucus <- ordinate(ps_mucus,  method = "MDS", distance = "bray")
plot_ordination(ps_mucus, out.pcoa.logmucus, shape = "ExpGroup") +
  labs(col = "Group", shape = "ExpGroup") +
  scale_shape_manual(values=1:4) +
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = Group, label = SampleID))

# weighted Unifrac PCoA mucus
ordu = ordinate(ps_mucus, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_mucus, ordu, shape="ExpGroup") + 
  scale_shape_manual(values=1:4) +
  labs(col = "Group", shape = "ExpGroup") +
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = Group, label = SampleID))

# PCoA of ileum
ps_ileum <- ps_silva %>%
  subset_samples(Group %in% c("ILC", "ILR", "ILS", "ILE")) 

out.pcoa.logileum <- ordinate(ps_ileum,  method = "MDS", distance = "bray")
plot_ordination(ps_ileum, out.pcoa.logileum, shape = "ExpGroup") +
  labs(col = "Group", shape = "ExpGroup") +
  scale_shape_manual(values=1:4) +
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = Group, label = SampleID))

# weighted Unifrac PCoA ileum
ordu = ordinate(ps_ileum, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps_ileum, ordu, shape="ExpGroup") + 
  scale_shape_manual(values=1:4) +
  labs(col = "Group", shape = "ExpGroup") +
  geom_point(alpha=.7, size=2) + 
  geom_label_repel(aes(fill = Group, label = SampleID))

# DPCoA (double principal coordinates analysis), a phylogenetic ordination method which
# provides a biplot representation of both samples and taxonomic categories
out.dpcoa.log <- ordinate(ps_silva, method = "DPCoA")
evals <- out.dpcoa.log$eig
plot_ordination(ps_silva, out.dpcoa.log, color = "Group", label= "Group", shape = "ExpGroup") +
  labs(col = "Group", shape = "ExpGroup") +
  coord_fixed(sqrt(evals[2] / evals[1]))

plot_ordination(ps_silva, out.dpcoa.log, type = "Species", color = "Phylum") +
  coord_fixed(sqrt(evals[2] / evals[1]))

## DPCoA plot with proportion
ord.DPCoA.silva <- ordinate(ps_silva.prop, method="DPCoA")
plot_ordination(ps_silva.prop, ord.DPCoA.silva, color="Group", title="DPCoA")

save.image("16S_phyloseq_8.Rdata")

# Bar plot:
# make each group samples together
ps_grouped <- merge_samples(ps_silva, "Group" )
ps_grouped.prop <- transform_sample_counts(ps_grouped, function(OTU) OTU/sum(OTU))
theme_set(theme_bw())

P=plot_bar(ps_grouped.prop, fill = "Phylum") + geom_bar(stat = "identity") +
  facet_wrap(~Phylum, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) + xlab("Group") +
  scale_x_discrete(limits=c("F-1", "F0", "F1C", "F1E", "F1R", "F1S", "F5C", "F5E", "F5R", "F5S", "F10C", "F10E", "F10R", "F10S", "DUE", "DUR", "DUS", "JEC", "JEE", "JER", "JES", "ILC", "ILE", "ILR", "ILS", "CEE", "CER", "CES", "COC", "COE", "COR", "COS", "MNC", "MNE", "MNR", "MNS", "HR", "HS1", "HS2"))

# use below to show plots one by one
# facet_wrap_paginate(~Phylum, scales = "free_y", ncol=1, nrow=1, page=1)

plot_bar(ps_grouped.prop, fill = "Phylum") + geom_bar(stat = "identity") +
  facet_wrap_paginate(~Phylum, scales = "free_y", ncol=3, nrow=1, page=1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) + xlab("Group") +
  scale_x_discrete(limits=c("F-1", "F0", "F1C", "F1E", "F1R", "F1S", "F5C", "F5E", "F5R", "F5S", "F10C", "F10E", "F10R", "F10S", "DUE", "DUR", "DUS", "JEC", "JEE", "JER", "JES", "ILC", "ILE", "ILR", "ILS", "CEE", "CER", "CES", "COC", "COE", "COR", "COS", "MNC", "MNE", "MNR", "MNS"))

# delete ambiguous phylum annotation
ps_silva_8phylum <- subset_taxa(ps_silva, !Phylum %in% c("Deferribacteres","Fusobacteria", "Synergistetes"))
ps_silva_8phylum_prop <- transform_sample_counts(ps_silva_8phylum, function(OTU) OTU/sum(OTU))
plot_bar(ps_silva_8phylum_prop, x= "Group", fill = "Phylum") + geom_bar(stat = "identity") +
  facet_wrap_paginate(~Phylum, scales = "free_y", ncol=1, nrow=1, page=6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) + xlab("Group") 
  scale_x_discrete(limits=c("F-1", "F0", "F1C", "F1E", "F1R", "F1S", "F5C", "F5E", "F5R", "F5S", "F10C", "F10E", "F10R", "F10S", "DUE", "DUR", "DUS", "JEC", "JEE", "JER", "JES", "ILC", "ILE", "ILR", "ILS", "CEE", "CER", "CES", "COC", "COE", "COR", "COS", "MNC", "MNE", "MNR", "MNS"))
  

top20_silva <- names(sort(taxa_sums(ps_silva), decreasing=TRUE))[1:20]
ps_silva.top20 <- transform_sample_counts(ps_silva, function(OTU) OTU/sum(OTU))
ps_silva.top20 <- prune_taxa(top20_silva, ps_silva.top20)
plot_bar(ps_silva.top20, x="Group", fill="Phylum") + facet_wrap(~Phylum, scales="free_y")

# Heatmap
plot_heatmap(carbom_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Class", taxa.order = "Class", 
             trans=NULL, low="beige", high="red", na.value="beige")
# Plot tree
plot_tree(ps_silva, color = "Group", shape = "Phylum", label.tips = "Phylum", size = "abundance")

# Network analysis
plot_net(ps_silva, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="Class", point_label="Order") 

save.image("16S_phyloseq_5.Rdata")


# below are used to change fonts of elements in the figure
# p + theme_bw() + 
#  geom_text(mapping = aes(label = SampleID), size = 6, vjust = 1.5) + #label text size
#  theme(text = element_text(size = 16)) + #axis font
#  geom_point(size = 4) #label size

# merge_samples cannot work on Factor type, so you need to change Group type 
# merge_samples just calculate all group sum, not average. It doesn't matter if you use relative abundance
# If you use absolute abundance, you need to calculate average by the following
# aaaa <- aggregate(. ~ as.factor(sample_data(ps_silva)$Group), data = as.data.frame(otu_table(ps_silva)), mean)
# rownames(aaaa)=aaaa$`as.factor(sample_data(ps_silva)$Group)`
# aaaa <- aaaa[,-1]
# ps_silva2 <- phyloseq(otu_table(aaaa, taxa_are_rows = FALSE), tax_table(taxa), phy_tree(fitGTR$tree))

# Calculate weighted unifrac distances
# weighted_unifrac_silva <- UniFrac(ps_silva.prop, weighted = TRUE, normalized = TRUE, parallel = TRUE, fast = TRUE)
# write.table(as.matrix(weighted_unifrac_silva), "wunifrac_matrix_silva.csv")
# ord.nmds.wunifrac <- ordinate(ps_silva.prop, method="NMDS", distance="wunifrac")
# plot_ordination(ps_silva.prop, ord.nmds.wunifrac.silva, color="Group", label= "Group", 
#                 shape = "ExpGroup", title="Weighted unifrac NMDS")

# plot abundance
# here the plot Firmicutes 
plot_abundance = function(phyloseqproject, title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(phyloseqproject, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Group",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) 
}

plotBefore = plot_abundance(ps_silva,"")
plotAfter = plot_abundance(ps_silva.prop,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2,  plotBefore, plotAfter)