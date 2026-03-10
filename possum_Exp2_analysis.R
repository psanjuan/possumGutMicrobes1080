# -------------------------------------------------------------------------
#   Author: Priscilla A San Juan
#   Topic: Microbiome pipeline for possum microbiome
#   Experiment 2
# -------------------------------------------------------------------------

# LOAD PACKAGES -----------------------------------------------------------
req_pkg <- c("readr","microbiome","dplyr","phyloseq","vegan","ggplot2","utils", 
             "DESeq2","ape","DescTools","scales","tidyr","ResourceSelection", 
             "gridExtra","ggbeeswarm","mvabund","RColorBrewer","randomcoloR","MASS", 
             "viridis","ggpubr","decontam","ggbeeswarm") 
# Load all required packages and show version
for(i in req_pkg){
  print(i)
  print(packageVersion(i))
  library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE)
}
library(data.table)


# Load new otu and tax table
OTUexp2 <- as.matrix(fread("data/possum_new_otu_tab_15_August_2025.tsv"), rownames = "V1")
OTUexp2 <- otu_table(OTUexp2, taxa_are_rows = T)

TAXexp2 <- as.matrix(fread("data/possum_new_tax_tab_15_August_2025.tsv", header=T), rownames = "V1")
TAXexp2 <- tax_table(TAXexp2)

SAMexp2 <- read.delim("data/possum_metadata_15_August_2025.tsv", header = TRUE, sep = "\t")
SAMexp2 <- sample_data(SAMexp2)
row.names(SAMexp2) <- SAMexp2$seq_name

# make phyloseq object
NewPossExp2 <- phyloseq(OTUexp2,TAXexp2,SAMexp2) 

# Filter for nuisance taxa, such as chloroplast or other known contaminants
# Taxa cleaning to remove sequences that aligned to chloroplast and mitochondria
NewPossExp2 <- NewPossExp2 %>% # 5342 taxa and 49 samples
  subset_taxa(
    Kingdom == "Bacteria" &
      Family != "Mitochondria" &
      Class != "Chloroplast" &
      Phylum != "Cyanobacteria" &
      Phylum != "Chloroplast" &
      Phylum != "Chloroflexi")

sample_sums(NewPossExp2)

# Look for skew in coverage across sample types
set.seed(711)
level_order <- c('sample', 'control') #set your own variable order

df = as.data.frame(sample_data(NewPossExp2))
df$LibrarySize = sample_sums(NewPossExp2)
df = df[order(df$LibrarySize),]
df$Index = seq(nrow(df))


# Plot ordered library size across 'batch' & coloured by 'replicate'
ggplot(data=df, aes(x=Index, y=LibrarySize, colour= sample_control))+
  geom_point()+
  facet_wrap(~ factor(sample_control, level = level_order)) +
  scale_y_continuous(trans='sqrt')

# DECONTAM: FILTER CONTAMS USING NEG CONTROL ------------------------------
library(decontam)
library(plotly)
ps <- NewPossExp2

# Inspect library sizes ---------------------------------------------------
df <- as.data.frame(sample_data(ps))
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
tem <- ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_control)) + geom_point()
ggplotly(tem)

# Use prevalence method to filter out suspected contaminants --------------
sample_data(ps)$is.neg <- sample_data(ps)$sample_control == "control"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold = 0.1)
hist(contamdf.prev$p) # looked at histogram to determine cut-off, range from 0.1-0.7 is the same
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$sample_control == "control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_control == "sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples ----------
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), 
                    pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Remove contams from phyloseq obj ----------------------------------------
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
ps.noncontam # 5323 taxa and 49 samples

# Remove the negative control from analysis -------------------------------
ps.noncontam <- prune_samples(sample_data(ps.noncontam)$sample_control != "control", ps.noncontam) 

# ready to use in phyloseq ------------------------------------------------
NewPossDecontam <- ps.noncontam # 5323 taxa and 47 samples

# every time you subset you must check for range in taxa sums and sample sums
range(taxa_sums(NewPossDecontam)) # taxa sums= 0 2600
# bc there is a zero we must set a min threshold (e.g. 10)
NewPossDecontam <- prune_taxa(taxa_sums(NewPossDecontam)>=1,NewPossDecontam) # 2073 taxa and 47 samples
range(sample_sums(NewPossDecontam))
range(taxa_sums(NewPossDecontam)) # 2073 taxa and 47 samples

# Need to drop sample poss22_PosMic159 - dropped in the sample sheet already!

# Alpha diversity ---------------------------------------------------------
library(agricolae)
library(forcats)

# drop samples with no eartag (possum individual id) and no treat.global.reformed
NewPossDecontamDropBlanksUnknown <- NewPossDecontam %>% # 2073 taxa and 25 samples
  subset_samples(sample_data(NewPossDecontam)$Eartag != "") %>%
  subset_samples(sample_data(NewPossDecontam)$treat.global.reformed != "")

NewPossDecontamDropBlanksUnknown <- NewPossDecontamDropBlanksUnknown %>% # 2073 taxa and 24 samples
  subset_samples(sample_data(NewPossDecontamDropBlanksUnknown)$treat.global.reformed != "unknown")

range(sample_sums(NewPossDecontamDropBlanksUnknown))
range(taxa_sums(NewPossDecontamDropBlanksUnknown))
NewPossDecontamDropBlanksUnknown <- prune_taxa(taxa_sums(NewPossDecontamDropBlanksUnknown)>=1,NewPossDecontamDropBlanksUnknown) # 1371 taxa and 24 samples

alphaNewTreatGlobal <- plot_richness(NewPossDecontamDropBlanksUnknown, 
                          x="treat.global.reformed", 
                          measures="Shannon",
                          color="treat.global.reformed") 
alphaNewTreatGlobalFig <- alphaNewTreatGlobal$data %>%
  mutate(treatOrder = fct_relevel(treat.global.reformed, 
                                  "pre",
                                  "post.toxin",
                                  "post.toxin.new",
                                  "post.second")) %>%
  ggplot(aes(x=treatOrder, y=value, col=treat.global.reformed)) +
  geom_boxplot() + 
  geom_quasirandom(size = 3.0) + 
  scale_x_discrete(labels = c("pre" = "Wild possums", 
                              "post.toxin" = "Survivors known", 
                              "post.toxin.new"="Survivors unknown", 
                              "post.second"="Captive exposed")) +
  scale_colour_manual(values = c("pre"="#33A02C",
                                 "post.toxin"="#1F78B4",
                                 "post.toxin.new"="#B2DF8A",
                                 "post.second"="#A6CEE3"), 
                      labels = c("pre" = "Wild possums", 
                                 "post.toxin" = "Survivors known", 
                                 "post.toxin.new"="Survivors unknown", 
                                 "post.second"="Captive exposed")) +
  xlab("Treatment") + 
  ylab("Shannon Diversity Index") + 
  ggtitle("Alpha Diversity of Possums") +
  theme_classic() + 
  theme(legend.position="none"); alphaNewTreatGlobalFig
alphaNewTreatGlobalFig$layers <- alphaNewTreatGlobalFig$layers[-1]

# ANOVA and Tukey test | anova better for categorical
aovNewTreatGlobal <- aov(value~treat.global.reformed, data=alphaNewTreatGlobal$data)
summary(aovNewTreatGlobal)
hsdNewTreatGlobal <- HSD.test(aovNewTreatGlobal, "treat.global.reformed", group=T); hsdNewTreatGlobal
tukeyTreatGlobal <- TukeyHSD(aovNewTreatGlobal); tukeyTreatGlobal

# Bacterial Abundance -----------------------------------------------------
# Merge rare taxa to speed up examples
NewPossDecontamRel <- microbiome::transform(NewPossDecontamDropBlanksUnknown, "compositional")
PS.rel1 <- aggregate_rare(NewPossDecontamRel, level = "Phylum", detection = 0/100, prevalence = 0/100)
PS.rel2 <- aggregate_rare(NewPossDecontamRel, level = "Family", detection = 0.01, prevalence = 10/100)
PS.rel3 <- aggregate_rare(NewPossDecontamRel, level = "Genus", detection = 0.03, prevalence = 10/100)

PS.rel.melt <- NewPossDecontamRel %>%
  psmelt
all_phylum.bac <- PS.rel1 %>%
  psmelt 
all_fam.bac <- PS.rel2 %>%
  psmelt 
all_gen.bac <- PS.rel3 %>%
  psmelt 
custom_labels2 <- c("pre" = "Wild possums", 
                    "post.toxin" = "Survivors known", 
                    "post.toxin.new"="Survivors unknown", 
                    "post.second"="Captive exposed")

# Phyla -------------------------------------------------------------------
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
                    "Verrucomicrobiota"="#E3DA4A")

all_phylum.bac$treat.global.reformed <- factor(all_phylum.bac$treat.global.reformed, 
                                               levels = c("pre",
                                                          "post.toxin",
                                                          "post.toxin.new",
                                                          "post.second"))
phylaComp <- ggplot(all_phylum.bac, 
                    aes(x = Sample, 
                        y = Abundance, 
                        fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~treat.global.reformed, 
             scales = "free_x", 
             space = "free_x", 
             labeller = as_labeller(custom_labels2)) +
  scale_fill_manual(values = more_phyla_pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_light()+
  theme(axis.title.x = element_blank(), 
        strip.text = element_text(size = 10, face = "bold", color = "black",hjust = 0),
        axis.text.x = element_blank()) +
        # axis.text.x = element_text(angle = 90)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance") +
  xlab("Sample name") + 
  ggtitle("Bacterial Phyla Composition"); phylaComp

# Family ------------------------------------------------------------------
more_fam_pal <- c("#6EE6D7",
                  "#E0E9A3", 
                  "#E1C5DF", 
                  "#DC8053",
                  "#C5EAD2", 
                  "#5E6C92", 
                  "#DB4DAE", 
                  "#9BD371",
                  "#7B7FD9",
                  "#7BACDA",
                  "#7A3EDF",
                  "#75E7D0",
                  "#B9828C",
                  "#E3DA4A",
                  "pink")

all_fam.bac$treat.global.reformed <- factor(all_fam.bac$treat.global.reformed, 
                                            levels = c("pre",
                                                       "post.toxin",
                                                       "post.toxin.new",
                                                       "post.second"))
famComp <- ggplot(all_fam.bac, 
                    aes(x = Sample, 
                        y = Abundance, 
                        fill = Family)) + 
  geom_bar(stat = "identity", 
           position = "fill") +
  facet_grid(~treat.global.reformed, 
             scales = "free_x", 
             space = "free_x", 
             labeller = as_labeller(custom_labels2)) +
  scale_fill_manual(values = more_fam_pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_light()+
  theme(axis.title.x = element_blank(), 
        strip.text = element_text(size = 10, face = "bold", color = "black",hjust = 0),
        axis.text.x = element_blank()) +
  # axis.text.x = element_text(angle = 90)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance") +
  xlab("Sample name") + 
  ggtitle("Bacterial Family Composition"); famComp

# Genus -------------------------------------------------------------------
more_gen_pal <- c("#6EE6D7",
                  "#E0E9A3",
                  "#E1C5DF",
                  "#DC8053",
                  "#C5EAD2", 
                  "#5E6C92", 
                  "#DB4DAE", 
                  "#9BD371",
                  "#7B7FD9",
                  "#7BACDA",
                  "#7A3EDF",
                  "coral",
                  "#B9828C",
                  "#E3DA4A",
                  "pink",
                  "#bf812d",
                  "#35978f",
                  "grey")
all_gen.bac$treat.global.reformed <- factor(all_gen.bac$treat.global.reformed, 
                                            levels = c("pre",
                                                       "post.toxin",
                                                       "post.toxin.new",
                                                       "post.second"))
genComp <- ggplot(all_gen.bac, 
                  aes(x = Sample, 
                      y = Abundance, 
                      fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~treat.global.reformed, 
             scales = "free_x", 
             space = "free_x",
             labeller = as_labeller(custom_labels2)) +
  scale_fill_manual(values = more_gen_pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.0))) +
  theme_light() + 
  # theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.title.x = element_blank(),
        strip.text = element_text(size = 10, face = "bold", color = "black",hjust = 0),
        axis.text.x = element_blank()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance") +
  xlab("Sample name"); genComp

# -------------------------------------------------------------------------
# BETA DIVERSITY: ORDINATION & PERMANOVA 
# Set seed for reproducibility
set.seed(2)

# Ordination --------------------------------------------------------------
data.ordinaPCoA.bray <- ordinate(NewPossDecontamRel, method="PCoA", distance="bray")
data.ordinaNMDS.bray <- ordinate(NewPossDecontamRel, method="NMDS", distance="bray")

PCoA.bray <- plot_ordination((NewPossDecontamRel),
                             data.ordinaPCoA.bray,
                             # color="Eartag",
                             # shape="treat.global.reformed",
                             color="treat.global.reformed") + 
  theme_classic() + 
  ggtitle("PCoA with Bray-Curtis Distance") +
  scale_colour_manual(values = c("post.second"="#A6CEE3",
                                 "post.toxin"="#1F78B4",
                                 "post.toxin.new"="#B2DF8A",
                                 "pre"="#33A02C"),
                      name = "Treatment",
                      labels = c("Wild possums",
                                 "Survivors known",
                                 "Survivors unknown",
                                 "Captive exposed")) +
  scale_fill_manual(values = c("post.second"="#A6CEE3",
                               "post.toxin"="#1F78B4",
                               "post.toxin.new"="#B2DF8A",
                               "pre"="#33A02C")) +
  # scale_shape_manual(values = c(16, 17, 15, 18, 8), 
  #                    breaks = c("pre",
  #                               "post.toxin",
  #                               "post.toxin.new",
  #                               "post.second")) +
  geom_line(aes(group = Eartag), col="#7D7D7D7D") +
  stat_ellipse(geom = "polygon", 
               aes(fill = treat.global.reformed), 
               alpha = 0.1,
               size=0.1) +
  geom_point(size=5, alpha=1) + 
  guides(colour = guide_legend(reverse=T)) +
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=16)); PCoA.bray 

colors14 <- c("#A6CEE3",
              "#1F78B4",
              "#B2DF8A",
              "#33A02C",
              "#FB9A99",
              "#E31A1C",  
              "#FDBF6F",
              "#FF7F00",
              "#CAB2D6",
              "#6A3D9A",
              "#F9E79F",
              "#B15928",
              "pink",
              "grey")

# PermANOVA ---------------------------------------------------------------
df.poss <- as(sample_data(NewPossDecontamRel),"data.frame")
dist.matrix.wha <- phyloseq::distance(NewPossDecontamRel,"bray")

# One-factor PERMANOVA
permanovaTreat <- adonis2(dist.matrix.wha~df.poss$treat.global.reformed,
                         df.poss,
                         permutations = 100000); permanovaTreat

dispersion <- betadisper(dist.matrix.wha, df.poss$treat.global.reformed)

# Run the permutational test with pairwise comparisons
pairwise_results <- permutest(dispersion, pairwise = TRUE, permutations = 999)

# Distance to centroid ----------------------------------------------------
df.poss$NewGroup <- forcats::fct_collapse(df.poss$treat.global.reformed, 
                      NoToxin=c("pre","post.toxin.new"),
                      Toxin=c("post.toxin","post.second"))
dispersion <- betadisper(dist.matrix.wha, df.poss$NewGroup)
distances_df <- data.frame(Distance = dispersion$distances, Group = dispersion$group)

library(dplyr)
betadiv <- ggplot(distances_df, 
                  aes(x = Group, 
                      y = Distance)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(y = "Distance to Centroid", x = "Group")
betadivfig <- betadiv$data %>%
  mutate(treatOrder = fct_relevel(Group, "NoToxin","Toxin")) %>%
  ggplot(aes(x=treatOrder, y=Distance, col=Group)) +
  geom_boxplot() + 
  geom_quasirandom(size = 3.0) + 
  scale_colour_manual(values = c("NoToxin"="#33A02C", "Toxin"="#1F78B4")) +
  xlab("Treatment") + 
  ylab("Distance to Centroid") + 
  annotate("text", x = 1, y = 0.69, label = "a", color = "black", size = 4) +
  annotate("text", x = 2, y = 0.69, label = "b", color = "black", size = 4) +
  annotate("text", x = 2, y = 0.62, label = "ANOVA, p=0.03", color = "black", size = 4) +
  theme_classic() + 
  theme(legend.position="none"); betadivfig

betadisper_result <- betadisper(bray_dist, groups)

# Step 4: Permutational test for homogeneity of dispersion
anova_result <- anova(dispersion); anova_result
pairwise_result <- permutest(dispersion, pairwise = TRUE); pairwise_result

# DESeq2 ------------------------------------------------------------------
# removing captive exposed
PrePostTox <-  NewPossDecontamDropBlanksUnknown %>%
  subset_samples(sample_data(NewPossDecontamDropBlanksUnknown)$treat.global.reformed != "post.second") 
# removing survivors unknown in order to compare wild and known survivors
PrePostTox <- PrePostTox %>%
  subset_samples(sample_data(PrePostTox)$treat.global.reformed != "post.toxin.new")
range(sample_sums(PrePostTox)) # 418 5259
range(taxa_sums(PrePostTox)) # 0 691
PrePostTox <- prune_taxa(taxa_sums(PrePostTox)>=1,PrePostTox) # 2 691

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
poss_DESeq <- phyloseq_to_deseq2(PrePostTox, ~ treat.global.reformed)
geoMeans = apply(counts(poss_DESeq), 1, gm_mean)
poss_DESeq = estimateSizeFactors(poss_DESeq, geoMeans = geoMeans)
poss_DESeq = DESeq(poss_DESeq, fitType="local")

res = results(poss_DESeq, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
# BUG!  object 'NewPossDecontamEartagTreat' not found
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(PrePostTox)[rownames(sigtab), ], "matrix"))
allsigtab = cbind(as(res, "data.frame"), as(tax_table(PrePostTox)[rownames(res), ], "matrix"))
head(sigtab)
head(allsigtab)

# -------------------------------------------------------------------------
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
sigtabfam = subset(sigtab, !is.na(Family))

taxa_overrep_fam <- ggplot(sigtabgen, aes(y=Family, x=log2FoldChange, color=Phylum, alpha=0.9, stroke=1)) + 
  geom_vline(xintercept = 0, linetype = "dotted", size = 1) +
  theme_light()+
  geom_point(size=6) + ylab("Bacterial Genus") + xlab("log2FoldChange") +
  scale_color_manual(values = more_phyla_pal) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust=0.5), 
        axis.text = element_text(size=14), axis.title = element_text(size=15),
        # legend.title.align = 0.5, 
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15)) + guides(color = guide_legend("Bacterial Family")) +
  ggtitle("Differential abundance of bacterial genera")
taxa_overrep_fam

# -------------------------------------------------------------------------

taxa_overrep_cap <- ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum, alpha=0.9, stroke=1)) + 
  geom_vline(xintercept = 0, linetype = "dotted", size = 1) +
  theme_light()+
  geom_point(size=6) + ylab("Bacterial Genus") + xlab("log2FoldChange") +
  scale_color_manual(values = more_phyla_pal) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust=0.5), 
        axis.text = element_text(size=14), axis.title = element_text(size=15),
        # legend.title.align = 0.5, 
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15)) + guides(color = guide_legend("Bacterial Family")) +
  ggtitle("Differential abundance of bacterial genera")
taxa_overrep_cap


