# -------------------------------------------------------------------------
#   Author: Priscilla A San Juan
#   Topic: Microbiome pipeline for possum microbiome
#   Experiment 1
# -------------------------------------------------------------------------

# LOAD PACKAGES  --------------------------------------------------------
req_pkg <- c("readr","microbiome","dplyr","phyloseq","vegan","ggplot2","utils", 
             "DESeq2","ape","forcats","DescTools","scales","tidyr","ResourceSelection", 
             "gridExtra","ggbeeswarm","mvabund","RColorBrewer","randomcoloR","MASS", 
             "viridis","ggpubr","decontam","ggbeeswarm","plotly", "stringr") 
# Load all required packages and show version
for(i in req_pkg){
  print(i)
  print(packageVersion(i))
  library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE)
}
library(agricolae)

# Set color palettes ------------------------------------------------------
library(randomcoloR)
colors_scheme_treatment <- c("pre"="#33A02C", "post"="#1F78B4")
palette_qual_all <- distinctColorPalette(25)
pie(rep(1, 25), col=palette_qual_all)

# 2 LOAD FILES ------------------------------------------------------------
# Load new otu and tax table
OTUexp1 <- as.matrix(fread("data_copied_from_possum_new_data_sept_2024/possum_new_otu_tab_15_August_2025_editednames.tsv"), rownames = "V1")
OTUexp1 <- otu_table(OTUexp1, taxa_are_rows = T)

TAXexp1 <- as.matrix(fread("data_copied_from_possum_new_data_sept_2024/possum_new_tax_tab_15_August_2025.tsv", header=T), rownames = "V1")
TAXexp1 <- tax_table(TAXexp1)

# Load two sample data sheets and full join
SAMexp1 <- read.delim("possum_metadata_Exp1.tsv", header = TRUE, sep = "\t")
row.names(SAMexp1) <- SAMexp1$submit_form_code

SAMexp2mod <- read.delim("possum_metadata_Exp2.tsv", header = TRUE, sep = "\t")
row.names(SAMexp2mod) <- SAMexp2mod$submit_form_code

# Joining data
sampledata_poss <- full_join(SAMexp1, SAMexp2mod, by = c("submit_form_code"))
row.names(sampledata_poss) <- sampledata_poss$submit_form_code
sampledata_poss <- sample_data(sampledata_poss)

# make phyloseq object
NewPossExp1 <- phyloseq(OTUexp1,TAXexp1,sampledata_poss)

# Statistics of read coverage across samples
max(sample_sums(NewPossExp1))
range(sample_sums(NewPossExp1))
hist(sample_sums(NewPossExp1), breaks = 20000, col = "orange", main = "Histogram of sample sequencing depth")

# Filter for nuisance taxa, such as chloroplast or other known contaminants
# Taxa cleaning to remove sequences that aligned to chloroplast and mitochondria
NewPossExp1 <- NewPossExp1 %>% 
  subset_taxa(
    Kingdom == "Bacteria" &
      Family != "Mitochondria" &
      Class != "Chloroplast" &
      Phylum != "Cyanobacteria" &
      Phylum != "Chloroplast" &
      Phylum != "Chloroflexi")

# 3 SUBSET DATA -----------------------------------------------------------
# Subset by new sample locations 
FranzJosefSamples <- prune_samples(sample_data(NewPossExp1)$Location.x == "Franz Josef", NewPossExp1)
FranzJosefSamples <- prune_taxa(taxa_sums(FranzJosefSamples)>1,FranzJosefSamples)

# Drop unknown
FranzJosefSamples <- subset_samples(FranzJosefSamples, treat.global.reformed != "unknown")

# Drop post.second (captive exposed)
FranzJosefSamples <-  subset_samples(FranzJosefSamples, treat.global.reformed != "post.second")

# Drop sample poss22_PosMic159 first
FranzJosefSamples <- prune_samples(sample_names(FranzJosefSamples) != "PosMic148", FranzJosefSamples)
FranzJosefSamples <- prune_taxa(taxa_sums(FranzJosefSamples)>1,FranzJosefSamples)


# 5 RICHNESS --------------------------------------------------------------
alpha.FJ <- plot_richness(FranzJosefSamples, 
                               x="treat.global.reformed", 
                               measures="Shannon",
                               color="treat.global.reformed") +
  geom_boxplot() + 
  scale_colour_manual(values = c("pre"="#33A02C",
                                 "post.toxin"="#1F78B4",
                                 "post.toxin.new"="#B2DF8A",
                                 "post.second"="#A6CEE3"), 
                      labels = c("pre" = "Wild possums", 
                                 "post.toxin" = "Survivors known", 
                                 "post.toxin.new"="Survivors unknown", 
                                 "post.second"="Captive exposed")) + 
  theme_light() + xlab("Treatment")

# Filter out category "post.second"
filtered_FJ <- alpha.FJ$data %>% filter(treat.global.reformed != "post.second")

filtered_FJ %>%
  mutate(treatOrder = fct_relevel(treat.global.reformed, "pre","post.toxin")) %>%
  ggplot(aes(x=treatOrder, y=value, col=treat.global.reformed)) +
  geom_boxplot() + 
  geom_quasirandom(size = 3.0) + 
  #scale_colour_manual(values = colors_scheme_treatment) +
  scale_x_discrete(labels = c("pre" = "Wild possums", "post.toxin" = "Survivors known")) +
  scale_colour_manual(values =  c("pre"="#33A02C", "post.toxin"="#1F78B4"),
                    labels = c("pre" = "Wild possums", "post.toxin" = "Survivors known")) +
  xlab("Treatment") + ylab("Shannon Diversity Index") +
  # labs(col = "Possum Group") +
  ggtitle("Franz Josef Possums") +
  annotate("text", x = 1, y = 5, label = "a", color = "black", size = 4) +
  annotate("text", x = 2, y = 5, label = "a", color = "black", size = 4) +
  annotate("text", x = 2.3, y = 2.6, label = "ANOVA, p=0.813", color = "black", size = 4) +
  theme_light() + theme(legend.position="none")

# ANOVA and Tukey test | anova better for categorical
aov.FJ <- aov(value~treat.global.reformed, data=alpha.FJ$data)
summary(aov.FJ)
HSD.FJ <- HSD.test(aov.FJ, "treat.global.reformed", group=T); HSD.FJ
tuk.FJ <- TukeyHSD(aov.FJ); tuk.FJ


# Bacterial Abundance -----------------------------------------------------
# Merge rare taxa to speed up examples

# Only pre and post samples
FranzJosefSamplesWildSurKnown <- subset_samples(FranzJosefSamples, treat.global.reformed != "post.second")
FranzJosefSamplesWildSurKnown <- prune_taxa(taxa_sums(FranzJosefSamplesWildSurKnown)>1,FranzJosefSamplesWildSurKnown)
FranzJosefSamplesWildSurKnown <- microbiome::transform(FranzJosefSamplesWildSurKnown, "compositional")

FJ.rel1 <- aggregate_rare(FranzJosefSamplesWildSurKnown, level = "Phylum", detection = 0/100, prevalence = 0/100)
FJ.rel2 <- aggregate_rare(FranzJosefSamplesWildSurKnown, level = "Family", detection = 0.01, prevalence = 10/100)
FJ.rel2b <- aggregate_rare(FranzJosefSamplesWildSurKnown, level = "Family", detection = 0.01, prevalence = 5/100)

FJ.rel3 <- aggregate_rare(FranzJosefSamplesWildSurKnown, level = "Genus", detection = 0.03, prevalence = 10/100)
physeq_Synergistes <- subset_taxa(FranzJosefSamplesWildSurKnown, Phylum == "Synergistetes")

physeq_Synergistes_melt <- physeq_Synergistes %>%
  psmelt
FJ.rel.melt <- FranzJosefSamplesWildSurKnown %>%
  psmelt
FJ_phy.bac <- FJ.rel1 %>%
  psmelt 
FJ_fam.bac <- FJ.rel2b %>%
  psmelt 
FJ_gen.bac <- FJ.rel3 %>%
  psmelt 

# Phyla
more_phyla_pal <- c("Actinobacteria"="#6EE6D7",
                    "Bacteroidetes"="#E0E9A3", 
                    "Elusimicrobia"="#E1C5DF", 
                    "Epsilonbacteraeota"="#DC8053",
                    "Firmicutes"="#C5EAD2", 
                    "Proteobacteria"="#5E6C92", 
                    "Fusobacteria"="#DB4DAE", 
                    "Lentisphaerae"="#9BD371",
                    "Synergistetes"="#7B7FD9",
                    "Spirochaetes"="#7BACDA",
                    "Myxococcota"="#7A3EDF",
                    "Tenericutes"="#75E7D0",
                    "Planctomycetes"="#B9828C",
                    "Verrucomicrobiota"="#E3DA4A",
                    "Acidobacteria"="grey", 
                    "Deferribacteres"="pink")
custom_labels <- c("pre" = "Wild possums", 
                    "post.toxin" = "Survivors known", 
                    "post.toxin.new"="Survivors unknown", 
                    "post.second"="Captive exposed")
FJphylaComp <- ggplot(FJ_phy.bac, aes(x = Eartag_or_name, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~treat.global.reformed, scales = "free_x", space = "free_x", labeller = as_labeller(custom_labels)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_light()+
  theme(axis.title.x = element_blank(), 
        strip.text = element_text(size = 10, face = "bold", color = "black",hjust = 0),
        axis.text.x = element_text(angle = 90)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance") +
  xlab("Sample name") + 
  ggtitle("Franz Josef Bacterial Phyla Composition"); FJphylaComp

# PhylaComp all samples combined by treatment to calculate F:B ratio
FJphylaCompByTreat <- ggplot(FJ_phy.bac, aes(x = treat.global.reformed, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~treat.global.reformed, scales = "free_x", space = "free_x", labeller = as_labeller(custom_labels)) +
  scale_fill_manual(values = more_phyla_pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_light()+
  theme(axis.title.x = element_blank(), 
        strip.text = element_text(size = 10, face = "bold", color = "black",hjust = 0),
        axis.text.x = element_blank()) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance") +
  xlab("Sample name") + 
  ggtitle("Bacterial Phyla Composition"); FJphylaCompByTreat

# Extract data from above figure to get relative abundance of F and B ratio
FB_data <- data.frame(Phyla=FJphylaCompByTreat$data$OTU, Sample=FJphylaCompByTreat$data$Sample, 
                      RelativeAbundance=FJphylaCompByTreat$data$Abundance, Treatment=FJphylaCompByTreat$data$treat.global.reformed)
# Select only Firmicutes and Bacteroidetes
FB_data_refined <- FB_data %>%
  filter(Phyla==c("Firmicutes", "Bacteroidetes"))

# Separate by treatment
Wild_only <- FB_data_refined %>%
  filter(Treatment=="pre")
Survivors_only <- FB_data_refined %>%
  filter(Treatment=="post.toxin")

# Calculate average relative abundance by treatment
WildFB <- Wild_only %>%
  group_by(Phyla) %>%
  summarise(avg=mean(RelativeAbundance))
SurvFB <- Survivors_only %>%
  group_by(Phyla) %>%
  summarise(avg=mean(RelativeAbundance))

more_fam_pal <- c("#6EE6D7","#E0E9A3", "#E1C5DF", "#DC8053",
                           "#C5EAD2", "#5E6C92", "#DB4DAE", "#9BD371",
                           "#7B7FD9","#7BACDA","#7A3EDF","#75E7D0",
                           "#B9828C","#E3DA4A","pink")
      
FJ_fam.bac$treat.global.reformed <- factor(FJ_fam.bac$treat.global.reformed, levels = c("pre","post.toxin","post.toxin.new","post.second"))                           
more_fam_pal <- c("Acidaminococcaceae"="chartreuse3",
                  "Akkermansiaceae"="#6EE6D7",
                  "Bacteroidaceae"="#E0E9A3", 
                  "Christensenellaceae"="#E1C5DF", 
                  "Clostridiaceae_1"="#DC8053", 
                  "Clostridiales_vadinBB60_group"="#C5EAD2",
                  "Eggerthellaceae"="#5E6C92", 
                  "Enterobacteriaceae"="pink",
                  "Erysipelotrichaceae"="#DB4DAE", 
                  "Lachnospiraceae"="#9BD371", 
                  "Muribaculaceae"="#7B7FD9",
                  "Other"="#7BACDA",
                  "Prevotellaceae"="#7A3EDF",
                  "Rikenellaceae"="#75E7D0",
                  "Ruminococcaceae"="#B9828C",
                  "Spirochaetaceae"="azure2",
                  "Synergistaceae"="#E3DA4A",
                  "Tannerellaceae"="coral")

# Custom strip labels
custom_labels <- c("pre" = "Wild possums", 
                   "post.toxin" = "Survivors known", 
                   "post.toxin.new"="Survivors unknown", 
                   "post.second"="Captive exposed")

FJfamComp <- ggplot(FJ_fam.bac, aes(x = Sample, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~treat.global.reformed, scales = "free_x", space = "free_x", labeller = as_labeller(custom_labels)) +
  scale_fill_manual(values = more_fam_pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
  theme_light()+ 
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.title.x = element_blank(), 
        strip.text = element_text(size = 10, face = "bold", color = "black",hjust = 0),
        axis.text.x = element_blank()) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance") +
  xlab("Sample name"); FJfamComp

more_gen_pal <- c("Alistipes"="#6EE6D7",
                  "Bacteroides"="#E0E9A3", 
                  "Blautia"="#E1C5DF",
                  "Christensenellaceae_R-7_group"="#DC8053",
                  "Candidatus_Stoquefichus"="#C5EAD2", 
                  "Faecalibacterium"="#5E6C92", 
                  "Lachnospiraceae_NK4A136_group"="#DB4DAE",
                  "Marvinbryantia"="#9BD371",
                  "Other"="#7B7FD9",
                  "Lachnospiraceae_UCG-008"="#7BACDA",
                  "Ruminiclostridium_9"="#7A3EDF",
                  "Ruminococcaceae_NK4A214_group"="#75E7D0",
                  "Ruminococcaceae_UCG-005"="#B9828C",
                  "Ruminococcaceae_UCG-013"="#E3DA4A",
                  "Ruminococcaceae_UCG-014"="pink",
                  "Ruminococcus_1"="#bf812d",
                  "Sarcina"="#35978f",
                  "Shuttleworthia"="palegreen3",
                  "Unknown"="grey")
FJ_gen.bac$treat.global <- factor(FJ_gen.bac$treat.global.reformed, levels = c("pre","post.toxin","post.toxin.new","post.second"))

FJgenComp <- ggplot(FJ_gen.bac, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~treat.global.reformed, scales = "free_x", space = "free_x", labeller = as_labeller(custom_labels)) +
  scale_fill_manual(values = more_gen_pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_light()+
  theme(axis.title.x = element_blank(), 
        strip.text = element_text(size = 10, face = "bold", color = "black",hjust = 0),
        axis.text.x = element_blank()) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance") +
  xlab("Sample name") + 
  ggtitle("Bacterial Genus Composition"); FJgenComp

# remove 0 reads
physeq_Synergistes_b <- prune_samples(sample_sums(physeq_Synergistes)>1,physeq_Synergistes)
physeq_Synergistes_b_melt <- physeq_Synergistes_b %>%
  psmelt
physeq_Synergistes_c <- prune_taxa(taxa_sums(physeq_Synergistes_b)>1,physeq_Synergistes_b)
physeq_Synergistes_c_melt <- physeq_Synergistes_c %>%
  psmelt

physeq_Synergistes_c_melt$treat.global <- factor(physeq_Synergistes_c_melt$treat.global.reformed, 
levels = c("pre","post.toxin","post.toxin.new","post.second"))

phy_Syn <- ggplot(physeq_Synergistes_c_melt, 
                  aes(x = treat.global, y = log(Abundance), col = treat.global)) + 
  geom_boxplot() + 
  geom_quasirandom(size = 3.0) + 
  scale_colour_manual(values = c("pre"="#33A02C", "post.toxin"="#1F78B4")) +
  theme_light()+
  theme(legend.position = "none") +
  scale_x_discrete(label = labels) +
  ylab("Number of Reads") +
  xlab("Treatment") + 
  ggtitle("Phylum Synergistetes Abundance"); phy_Syn

# Custom X-axis labels 
labels <- c("Wild possums", "Survivors known")
 
phyloseq::psmelt(NewPossExp1) %>%
  ggplot(data = ., aes(x = treat.global.reformed, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Phylum, scales = "free")

# 7 ORDINATION ------------------------------------------------------------
# A couple of standardization options before doing the analysis
data.otu.FJ <- otu_table(FranzJosefSamples)        

# A couple of pairwise distance metric 
distanceFJ <- phyloseq::distance(data.otu.FJ, method="bray") 

# Ordination
ordinationFJ <- ordinate(FranzJosefSamples, method="PCoA", distance="bray", k= 2, trymax=1000)

# Stress of distance matrix
stressplot(monoMDS(distanceFJ))

colors_scheme_treatment <- c("pre"="#33A02C",
                             "post.toxin"="#1F78B4",
                             "post.toxin.new"="#B2DF8A",
                             "post.second"="#A6CEE3")

# NMDS all samples colored by INSERT VARIABLE  
nmdsFJ <- plot_ordination((FranzJosefSamples), ordinationFJ, 
                                 # color = "Sample_ID_2024", 
                                  color = "treat.global.reformed", 
                                 # shape = "treat.global",
                                 ) + 
  geom_point(size=4, alpha=1) + 
  theme_classic() + 
  stat_ellipse(geom = "polygon", 
                aes(fill = treat.global.reformed), 
                alpha = 0.1,
                size=0.1) +
  scale_colour_manual(values = c("pre"="#33A02C",
                                 "post.toxin"="#1F78B4")) +
  scale_fill_manual(values = c("pre"="#33A02C",
                               "post.toxin"="#1F78B4")) +
  #geom_line(aes(group = Eartag), col="#7D7D7D7D") +
  annotate("text", 
           x = 0.38, 
           y = -0.65, 
           label = "PERMANOVA", 
           color = "black", size = 4) +
  annotate("text", 
           x = 0.42, 
           y = -0.7, 
           label = "p=0.001, r2=0.11", 
           color = "black", size = 4) +
  theme(axis.text = element_text(size=15), 
        axis.title = element_text(size=15),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black")); nmdsFJ + 
  theme(legend.position="bottom")


# 8  PermANOVA ------------------------------------------------------------
# Distance metric 
metaFJ <- as(sample_data(FranzJosefSamples),"data.frame")
distFJ <- phyloseq::distance(FranzJosefSamples,"bray")

# One-factor PERMANOVA
permFJtreatment <- adonis2(distFJ ~ metaFJ$treat.global.reformed, metaFJ)
print(permFJtreatment)

# 9 BETA DIVERSITY --------------------------------------------------------
data.otu.FJ <- otu_table(FranzJosefSamples)

# Convert to data frame and transpose and calculate using betadiver 
beta.paFJ <- betadiver(as.data.frame(t(data.otu.FJ)), "hk")          # presence/absence
betaFJ <- vegdist(as.data.frame(t(data.otu.FJ)), method = "bray")    # abundance

# Calculate beta diversity using distances from above and matching to sample data
dispersionFJ <- with(sample_data(FranzJosefSamples), betadisper(beta, treat.global))

# Permutation test
permFJ <- permutest(dispersionFJ, pairwise=T)

# Putting data in long form so ggplot can read
distanceFJ.treatment  <- data.frame(distance = dispersionFJ$distances, treat = sample_data(FranzJosefSamples)$treat.global)

# Plotting in ggplot
betaFJ.treatment <- ggplot(distanceFJ.treatment, aes(x = treat, y = distance, col = treat)) + xlab("") +
  scale_colour_manual(name="Key", values=colors_scheme_treatment) +
  ylab("Distance to Centroid") + 
  geom_boxplot() + 
  ggtitle("Franz Josef Possum Microbiome Betadiversity") +
  geom_quasirandom(size=6) + theme_classic() +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=20), 
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "none") 
betaFJ.treatment + aes(x = factor(treat, level = c('pre', 'post', 'post-tox'))) + xlab("Treatment")

# anova to see if statistically significantly different 
temp.aov <- aov(distance ~ treat, data = distanceFJ.treatment)
summary(temp.aov)
temp.HSD <- HSD.test(temp.aov, "treat", group = T)
temp.HSD
betaFJ.treatment$label <- temp.HSD$groups$M

# DESeq2 ------------------------------------------------------------------
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
FJ_DESeq <- phyloseq_to_deseq2(FranzJosefSamples, ~ treat.global)
geoMeans = apply(counts(FJ_DESeq), 1, gm_mean)
FJ_DESeq = estimateSizeFactors(FJ_DESeq, geoMeans = geoMeans)
FJ_DESeq = DESeq(FJ_DESeq, fitType="local")

res = results(FJ_DESeq, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(FranzJosefSamples)[rownames(sigtab), ], "matrix"))
head(sigtab)

# -------------------------------------------------------------------------
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
sigtabfam = subset(sigtab, !is.na(Family))

# Phylum colors
phyla_pal <- c("Actinobacteriota"="#6EE6D7","Bacteroidota"="#E0E9A3", "Campylobacterota"="#7A3EDF", 
               "Deinococcota"="#DC8053","Firmicutes"="#C5EAD2", "Proteobacteria"="#5E6C92", 
               "Deferribacterota"="#DB4DAE", "Patescibacteria"="#9BD371")

taxa_overrep_cap <- ggplot(sigtabfam, aes(y=Family, x=log2FoldChange, color=Phylum, alpha=0.5, stroke=1)) + 
  geom_vline(xintercept = 0, linetype = "dotted", size = 1) +
  geom_point(size=6) + ylab("Family") + xlab("log2FoldChange") +
  scale_color_manual(values = palette_qual_all) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust=0.5), 
        axis.text = element_text(size=14), axis.title = element_text(size=18),
        legend.title.align = 0.5, legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15)) + guides(color = guide_legend("Bacterial Phyla")) +
  ggtitle("Differential abundance of bacterial families")
taxa_overrep_cap
