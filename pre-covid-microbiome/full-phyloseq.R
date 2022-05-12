
  #---- load packages and set seed ---- ####
```{r include=FALSE}
library(reshape2)
library(reshape)
library(stringi)
library(plyr)
library(dplyr)
library(tidyr)
library(DT)
library(data.table)
library(openxlsx)
library(rstatix)
library(gtools) 
library(gdata)
library(devtools)
library(readr)
library(patchwork)
library(plotly)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(ggsci)

library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(ANCOMBC)
library(qiime2R)
library(vegan)

set.seed(194175)
```



#---- functions ---- ####
```{r}
# function change N/A to zero
na.zero <- function (x) {x[is.na(x)] <- 0 ;return(x)  }
gm_mean = function(x, na.rm=TRUE){exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}

# faster alternative to the phyloseq psmelt() function
pseq.long <- function(x, transform.counts = "NULL") {
  asv_tab <- asv_tab_lng <- tax_tab <- NULL
  if (is.null(transform.counts)) {
    x <- x
  }
  else if (transform.counts == "log10") {
    x <- transform(x, "log10")
  }
  else if (transform.counts == "Z-OTU") {
    x <- transform(x, "Z", "OTU")
  }
  else if (transform.counts == "Z-Sample") {
    x <- transform(x, "Z", "Sample")
  }
  else if (transform.counts == "compositional") {
    x <- transform(x, "compositional", "OTU")
  }
  else {
    stop("Please provide appropriate transformation")
  }
  asv_tab <- as.data.frame(abundances(x)) # get asvs/otus
  asv_tab$asv_id <- rownames(asv_tab) # add a new column for ids
  #head(asv_tab)
  asv_tab_lng <- asv_tab %>%
    reshape2::melt()
  #head(asv_tab_lng)
  colnames(asv_tab_lng) <- c("asv_id", "sample", "abundance")
  # tax_tab <- as.data.frame(tax_table(x)) # get taxonomy note: can be slow
  tax_tab <- as.data.frame(as(x@tax_table, "matrix")) # get taxonomy note: can be slow
  # tax_tab <- as.data.frame(tax_tab)
  tax_tab$asv_id <- rownames(tax_tab) # add a new column for ids
  asv_tax_tab <- asv_tab_lng %>%
    right_join(tax_tab, by = "asv_id")
  
  meta_tg <- meta(x)
  meta_tg$sample <- rownames(meta_tg)
  asv_tax_meta_tab <- asv_tax_tab %>%
    right_join(meta_tg, by = "sample")
  
  return(asv_tax_meta_tab)
}
```


# ----  load all data ---- ####
```{r include=FALSE}
#according to https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121

OTU_table_L8 <- read_delim("/Users/mykajaapyoungapelian/covid-microbiome/metadataFiles/covid_microbiome_otu_table-dada2-L8.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # read otu table following dada2 processing

tax_table <- read_delim("/Users/mykajaapyoungapelian/covid-microbiome/metadataFiles/covid_microbiome_taxonomy.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # read taxonomy table text file

sampledata <- read_delim("/Users/mykajaapyoungapelian/covid-microbiome/metadataFiles/covid_microbiome_mapping_corrected.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # read mapping file with sequencing data

tax_table <- read_qza("/Users/mykajaapyoungapelian/covid-microbiome/metadataFiles/covid_microbiome_taxtab.qza") # read taxonomy qza file
tax_table <- parse_taxonomy(tax_table$data)
TAX <- data.frame(tax_table)

full_metadata <- read_delim("/Users/mykajaapyoungapelian/covid-microbiome/metadataFiles/covid_microbiome_full_metadata.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # full patient metadata


OTU <- column_to_rownames(OTU_table_L8, var = "OTUID")
MAP <- column_to_rownames(sampledata, var = "SampleID")
MAP2 <- column_to_rownames(full_metadata, var = "sampleID")
otuColNames <- colnames(OTU) # otuColNames
otuRowNames <- rownames(OTU) # otuRowNames
colnames(OTU) <- NULL # remove OTU row names and col names 
rownames(OTU) <- NULL # dim(OTU) 

mOTU <- as.matrix(OTU, rownames.force = TRUE) # force mOTU matrix to be numeric instead of character class
class(mOTU) <- 'numeric'
colnames(mOTU) <- otuColNames # add row names and col names back to matrix otuColNames
rownames(mOTU) <- otuRowNames
class(mOTU[1,]) 
mTAX <- as.matrix(TAX)
```


# ---- create phyloseq object & subset samples ---- ####
```{r include = FALSE}
p_OTU = otu_table(mOTU, taxa_are_rows = TRUE)
p_TAX = tax_table(mTAX)
p_MAP = sample_data(MAP)
p_MAP2 = sample_data(MAP2)
phyloseq <- merge_phyloseq(p_OTU, p_TAX, p_MAP, p_MAP2)


# filter samples to ensure quality
rank_names(phyloseq)
get_taxa_unique(phyloseq, "Kingdom")
get_taxa_unique(phyloseq)
phyloseq_filtered <- phyloseq
phyloseq_filtered <- subset_taxa(phyloseq_filtered, !Kingdom=="d__Eukaryota")
phyloseq_filtered <- subset_taxa(phyloseq_filtered, !Kingdom=="Unassigned")
phyloseq_filtered <- subset_taxa(phyloseq_filtered, !Genus=="d__Bacteria")
get_taxa_unique(phyloseq_filtered, "Kingdom")
get_taxa_unique(phyloseq_filtered, "Family")
length(get_taxa_unique(phyloseq_filtered, "Family"))# 69 families in total
phyloseq_filtered <- subset_taxa(phyloseq_filtered, !Family=="f__Mitochondria")
get_taxa_unique(phyloseq_filtered, "Family") # 100 families in total
length(get_taxa_unique(phyloseq_filtered, "Family")) # 69 total 
sort(sample_sums(phyloseq_filtered))

# cut-off of samples below 1000 reads 
phyloseq_filtered_1000 <- phyloseq_filtered
phyloseq_filtered_1000 <- prune_samples(sample_sums(phyloseq_filtered_1000)>=1000, phyloseq_filtered_1000)

# in case not removed in the mapping file, clean controls with the scripts below. 
phyloseq_filtered_1000 <- subset_samples(phyloseq_filtered_1000, !sample_ID_provider == "Mock_Standard_community")
phyloseq_filtered_1000 <- subset_samples(phyloseq_filtered_1000, !sample_ID_provider == "Mock_Log_community")
phyloseq_filtered_1000 <- subset_samples(phyloseq_filtered_1000, !sample_ID_provider == "Useq_Mock_Standard")
phyloseq_filtered_1000 <- subset_samples(phyloseq_filtered_1000, !sample_ID_provider == "Useq_Mock_Log")
phyloseq_filtered_1000 <- subset_samples(phyloseq_filtered_1000, !sample_ID_provider == "Zymobiomics_gDNA")
phyloseq_filtered_1000 <- subset_samples(phyloseq_filtered_1000, !sample_ID_provider == "Zymobiomics_gDNA_1")
phyloseq_filtered_1000 <- subset_samples(phyloseq_filtered_1000, !sample_ID_provider == "Zymobionics_gDNA_2")
phyloseq_filtered_1000 <- subset_samples(phyloseq_filtered_1000, !sample_ID_provider == "Swab_Fernanda")
phyloseq_filtered_1000 <- subset_samples(phyloseq_filtered_1000, !sample_ID_provider == "Swab_Shu")
phyloseq_filtered_1000 <- subset_samples(phyloseq_filtered_1000, !sample_ID_provider == "Swab_Paul")
phyloseq_filtered_1000 <- subset_samples(phyloseq_filtered_1000, !sample_ID_provider == "Neg")
phyloseq_filtered_1000 <- subset_samples(phyloseq_filtered_1000, !sample_ID_provider == "Useq_neg")
phyloseq_filtered_1000 <- subset_samples(phyloseq_filtered_1000, !sample_ID_provider == "Mock_gDNA")
phyloseq_filtered_1000 <- subset_samples(phyloseq_filtered_1000, !sample_ID_provider == "20-023100") # remove duplicate 
phyloseq_filtered_1000 <- prune_taxa(taxa_sums(phyloseq_filtered_1000)>0, phyloseq_filtered_1000)

#check final file , only samples with more than 1000 reads should be included
sort(sample_sums(phyloseq_filtered_1000))
phyloseq_filtered_1000
# View(phyloseq_filtered_1000$sample_data)

# save files
# write.table(otu_table(phyloseq_filtered_1000), file = "phyloseq_filtered_1000.txt", sep = "\t") 
# write.table(tax_table(phyloseq_filtered_1000), file = "phyloseq_filtered_1000_tax.txt", sep = "\t")

#change file name to be able to run scripts directly. decide which physeq cutoff file to use.
phyl_utr <- phyloseq_filtered_1000
sam_mat = as(sample_data(phyl_utr), "data.frame")
sort(sample_sums(phyl_utr))

# define colors
cols_group <- c('yes' = rgb(0, 115, 189, maxColorValue = 255),
                'no' = rgb(191, 0, 58, maxColorValue = 255))
```


# ---- recode relevant metadata ---- #### 
```{r include = FALSE}
# recode pos/neg metadata 
sample_data(phyl_utr)$covid = recode(sample_data(
  phyl_utr)$covid,
  `0` = 'negative',
  `1` = 'positive'
)

sam_mat = as(sample_data(phyl_utr), "data.frame")
colnames(sam_mat)

# nr of samples
length(unique(rownames(sam_mat)))

# nr groups in COVID
unique(sam_mat$covid)

# recode covid severity metadata
sample_data(phyl_utr)$severity = recode(sample_data(
  phyl_utr)$severity,
  `0` = 'negative',
  `1` = 'mild',
  `2` = 'moderate',
  `3` = 'severe',
  `4` = 'death',
)

# level severity classes
levels(sample_data(phyl_utr)$severity) <- c("negative", "mild", "moderate", "severe", "death")

sam_mat = as(sample_data(phyl_utr), "data.frame")
colnames(sam_mat)

# nr samples
length(unique(rownames(sam_mat)))

# nr groups in COVID severity class type
unique(sam_mat$severity)

# nr of groups in group negative class
length(rownames(unique(dplyr::filter(sam_mat, 
                                     grepl("negative", severity,
                                           ignore.case=TRUE))%>%
                         select(severity)))) 

# nr of groups in group mild class
length(rownames(unique(dplyr::filter(sam_mat,
                                     grepl("mild",severity,
                                           ignore.case=TRUE))%>%
                         select(severity)))) 

# nr of groups in group moderate class
length(rownames(unique(dplyr::filter(sam_mat,
                                     grepl("moderate",severity,
                                           ignore.case=TRUE))%>%
                         select(severity)))) 

# nr of groups in group severe class
length(rownames(unique(dplyr::filter(sam_mat,
                                     grepl("severe",severity,
                                           ignore.case=TRUE))%>%
                         select(severity)))) 

# nr of groups in group death class
length(rownames(unique(dplyr::filter(sam_mat, 
                                     grepl("death", severity,
                                           ignore.case=TRUE))%>%
                         select(severity)))) 

# recode hospital transfer metadata 
sample_data(phyl_utr)$hosp_transfer = recode(sample_data(
  phyl_utr)$hosp_transfer,
  `0` = 'no',
  `1` = 'yes',
  `#N/A` = NULL
)      

sam_mat = as(sample_data(phyl_utr), "data.frame")

# nr samples
length(unique(rownames(sam_mat)))

# nr groups in hosp_transfer
unique(sam_mat$hosp_transfer)

# nr groups in group 1 (hospital transferred samples)
length(rownames(unique(dplyr::filter(sam_mat, 
                                     grepl('yes', hosp_transfer,
                                           ignore.case=TRUE))%>%
                         select(hosp_transfer)))) 

# nr groups in group 2 (not hospital transferred samples)
length(rownames(unique(dplyr::filter(sam_mat,
                                     grepl('no',hosp_transfer,
                                           ignore.case=TRUE))%>%
                         select(hosp_transfer))))

# recode pre-antibiotic metadata 
sample_data(phyl_utr)$pre_antibiotic_use = recode(sample_data(
  phyl_utr)$pre_antibiotic_use,
  `0` = 'no',
  `1` = 'yes',
  `#N/A` = NULL
)


sam_mat = as(sample_data(phyl_utr), "data.frame")

# nr samples
length(unique(rownames(sam_mat)))

# nr groups in pre_antibiotic_use
unique(sam_mat$pre_antibiotic_use)

sam_mat = as(sample_data(phyl_utr), "data.frame")
```


# ---- write filtered otu table, taxonomy table, and metadata ---- ####
```{r}
# subset table and aggregate to genus level
physeq_gen <- phyl_utr 

physeq_gen <- subset_samples(physeq_gen, gender %in% c(0, 1))
physeq_gen <- subset_samples(physeq_gen, covid %in% c("negative", "positive"))
physeq_gen <- subset_samples(physeq_gen, hosp_transfer %in% c('no', NA))
physeq_gen <- subset_samples(physeq_gen, pre_antibiotic_use %in% c(NA))
# use this `physeq_gen` phyloseq object for the ANCOM-BC analysis following exclusion criteria

# create otu table from physeq_gen object 
physeq_gen_otu <- otu_table(physeq_gen)
physeq_gen_otu<- cbind(newColName = rownames(physeq_gen_otu), physeq_gen_otu)
rownames(physeq_gen_otu) <- 1:nrow(physeq_gen_otu)
physeq_gen_otu <- as.data.frame(physeq_gen_otu)
setnames(physeq_gen_otu, c("newColName"), c("featureID"))
#write.table(physeq_gen_otu, file = "physeq_gen_otu.txt", sep = "\t", row.names = FALSE, col.names = TRUE) # write otu table

# create taxonomy table from physeq_gen object 
physeq_gen_tax <- tax_table(physeq_gen)
physeq_gen_tax<- cbind(newColName = rownames(physeq_gen_tax), physeq_gen_tax)
rownames(physeq_gen_tax) <- 1:nrow(physeq_gen_tax)
physeq_gen_tax <- as.data.frame(physeq_gen_tax)
setnames(physeq_gen_tax, c("newColName"), c("featureID"))
#write.table(physeq_gen_tax, file = "physeq_gen_tax.txt", sep = "\t", row.names = FALSE, col.names = TRUE) # write taxonomy table 

# create metadata file from filtered physeq_gen object 
physeq_gen_metadata <- sample_data(physeq_gen)
physeq_gen_metadata <- cbind(newColName = rownames(physeq_gen_metadata), physeq_gen_metadata)
rownames(physeq_gen_metadata) <- 1:nrow(physeq_gen_metadata)
physeq_gen_metadata <- as.data.frame(physeq_gen_metadata)
setnames(physeq_gen_metadata, c("newColName"), c("sampleID"))
#write.table(physeq_gen_metadata, file = "physeq_gen_metadata.txt", sep = "\t", row.names = FALSE, col.names = TRUE) # write filtered metadata file

#physeq_gen <- phyloseq::tax_glom(physeq_gen, "Genus", NArm = FALSE)
```




# ---- write final mapping file ---- ####
```{r}
# read in taxonomy table and otu table made from phyloseq object 
otu_table_final <- read_delim("/Users/mykajaapyoungapelian/covid-microbiome/physeq_gen_otu.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
tax_table_final <- read_delim("/Users/mykajaapyoungapelian/covid-microbiome/physeq_gen_tax.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# create data.table objects and merge on OTU as key
otuDT <- data.table(otu_table_final, key = "featureID")
taxDT <- data.table(tax_table_final, key = "featureID")
otuTaxDT <- merge(otuDT, taxDT, by = "featureID")

# melt data.table to include sample and phylogenetic designation 
idVars <- names(otuTaxDT[, 2:(ncol(otuTaxDT)-7)]) # obtain all sample numbers by excluding taxonomy

otuTaxDTMelted <- melt(otuTaxDT, id.vars = c("featureID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), measure.vars = c(idVars))
setnames(otuTaxDTMelted, c("variable", "value"), c("sampleID", "reads"))

otuTaxDTMelted$taxa_merged <- paste(otuTaxDTMelted$Kingdom, otuTaxDTMelted$Phylum, otuTaxDTMelted$Class, otuTaxDTMelted$Order, otuTaxDTMelted$Family, otuTaxDTMelted$Genus,  sep = ";")

otuTaxDTMelted <- otuTaxDTMelted[, c("taxa_merged", "sampleID", "reads")]
otuTaxDTMelted <- cbind(newColName = rownames(otuTaxDTMelted), otuTaxDTMelted)
rownames(otuTaxDTMelted) <- 1:nrow(otuTaxDTMelted)
castedDT <- cast(otuTaxDTMelted, sampleID ~ taxa_merged, fun.aggregate = sum, value = "reads")

physeq_gen_final_mapping <- merge(physeq_gen_metadata, castedDT, by = "sampleID")
#write.table(physeq_gen_final_mapping, file = "physeq_gen_final_mapping.txt", sep = "\t", na = "", row.names = FALSE, col.names = TRUE) # write final mapping file
```



# ---- mean reads per sample and range ---- ####
```{r}
physeq_gen_otu_melted <- setDT(physeq_gen_otu)
m_vars <- colnames(physeq_gen_otu)[2:length(colnames(physeq_gen_otu))]

potu_melt <- melt(physeq_gen_otu, id.vars = "featureID", measure.vars = m_vars, value.name = "reads", "sample")
potu_melt$reads <- as.numeric(potu_melt$reads)

potu_melt <- potu_melt[reads > 0]

potu_melt[, c := round(mean(reads), 2), by = sample]

round(sum(potu_melt$reads)/145,2)

range(potu_melt$reads) # 1-9331

# calculate unique asv per sample, mean, and range
potu_melt[, .(count = uniqueN(featureID)), by = sample]
mean(potu_melt[, .(count = uniqueN(featureID)), by = sample]$count)
range(potu_melt[, .(count = uniqueN(featureID)), by = sample]$count) # 2 - 99
```



# ---- top 20 genus phyloseq objects ---- ####
```{r}
physeq_gen_1 <- microbiome::aggregate_taxa(physeq_gen, "Genus") # create phyloseq on genus 
physeq_gen_2 <- aggregate_top_taxa2(physeq_gen, top = 20,  "Genus") # top 20 genus
# changed to aggregate_top_taxa2 which is another version of aggregate_top_taxa in microbiome package

# create phyloseq object with relative abundance of phyloseq_gen_2
physeq_gen_RA <- microbiome::transform(physeq_gen_2, "compositional") #compositional
taxa_sums(physeq_gen_RA)
#write.table(taxa_sums(physeq_gen_RA), file = "taxa_sums.txt", sep = "\t") # write top 20 taxa relative abundance 

TopGenus1 <- names(sort(taxa_sums(physeq_gen_RA), TRUE)[1:21])
Genus1 <- prune_taxa(TopGenus1, physeq_gen_RA)
Genus1
taxa_sums(Genus1)
```


# ---- mean and sd of RA per taxa per severity ---- ####
```{r}
# Pick the core (>0.1% relative abundance in >10% of the samples)
physeq_gen_sum <- aggregate_rare(physeq_gen, level = "Genus", detection = 0.1/100, prevalence = 10/100) #reduce to genus

physeq_gen_neg <- subset_samples(physeq_gen_sum, severity %in% "negative")
physeq_gen_neg_sum <- microbiome::transform(physeq_gen_neg, "compositional")
taxa_summary(physeq_gen_neg_sum, level = "Genus")

physeq_gen_mild <- subset_samples(physeq_gen_sum, severity %in% "mild")
physeq_gen_mild_sum <- microbiome::transform(physeq_gen_mild, "compositional")
taxa_summary(physeq_gen_mild_sum, level = "Genus")

physeq_gen_moderate <- subset_samples(physeq_gen_sum, severity %in% "moderate")
physeq_gen_moderate_sum <- microbiome::transform(physeq_gen_moderate, "compositional")
taxa_summary(physeq_gen_moderate_sum, level = "Genus")

physeq_gen_severe <- subset_samples(physeq_gen_sum, severity %in% "severe")
physeq_gen_severe_sum <- microbiome::transform(physeq_gen_severe, "compositional")
taxa_summary(physeq_gen_severe_sum, level = "Genus")

physeq_gen_death <- subset_samples(physeq_gen_sum, severity %in% "death")
physeq_gen_death_sum <- microbiome::transform(physeq_gen_death, "compositional")
taxa_summary(physeq_gen_death_sum, level = "Genus")
```




# ---- palettes and colors ---- ####
```{r include = FALSE}
# palette choices need to reflect the length of top 20 phyloseq objects 
PaletteChoice = rev(rainbow(length(taxa_sums(physeq_gen_RA)))) #pick colors
PaletteChoice

# use ggsci library and add additional colors for visualizations. 
PaletteChoice4 = pal_npg('nrc')(10)
PaletteChoice4 = c(PaletteChoice4, "#8c6bb1", "#bebada", "#fb8072", "#1f78b4", "#fb9a99", "#01665e", "#8e0152", "#d1e5f0", "#a50f15", "#4d4d4d", "#fed976")
```


```{r}
# plot un-grouped taxa covid pos/neg
{plot.composition.relAbun <- plot_composition(physeq_gen_RA, 
                                              sample.sort = NULL,
                                              otu.sort = names(sort(
                                                taxa_sums(physeq_gen_RA), F)[1:21]),
                                              x.label = "sample",
                                              group_by = "covid")
  #add settings plot
  plot.composition.relAbun + 
    scale_fill_manual(values=PaletteChoice4) + 
    theme_bw() +
    ylab("Relative Abundance")+
    scale_x_discrete(labels = NULL) +
    theme(legend.justification  = "top",
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(size = 13, face = "bold"),
          axis.title.y = element_text(color="Grey1", size=11, face="bold"),
          axis.text.x =element_text(color="Grey1", size=11, angle = 90),
          axis.text=element_text(color="Grey1", size=11, face="bold"),
          panel.grid.major = element_blank()
    )
  
}
#ggsave(filename = "supplementaryFigure2.tiff", width = 20, height = 10, units = "cm")
```



# ---- sars-cov-2 positive and negative plots ungroups ---- ####
```{r}
#plot ungrouped taxa covid pos/neg - npj
{plot.composition.relAbun <- plot_composition(physeq_gen_RA, 
                                              sample.sort = NULL,
                                              otu.sort = names(sort(
                                                taxa_sums(physeq_gen_RA), F)[1:21]),
                                              x.label = "sample",
                                              group_by = "covid")
  #add settings plot
  plot.composition.relAbun + 
    scale_fill_manual(values=PaletteChoice4) + 
    theme_bw() +
    ggtitle("COVID - relative abundance on genus level") +
    ylab("Relative Abundance")+
    theme_bw() + 
    scale_x_discrete(labels = NULL) +
    theme(legend.justification  = "top",
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(size = 13, face = "bold"),
          axis.title.y = element_text(color="Grey1", size=11, face="bold"),
          axis.text.x =element_text(color="Grey1", size=11, angle = 90),
          axis.text=element_text(color="Grey1", size=11, face="bold"),
          panel.grid.major = element_blank(),
          axis.ticks = element_blank()
    )
  
}
#ggsave(filename = "plotTaxaCovidUngrouped2_npj.pdf", width = 36, height = 20, units = "cm")
````
# ---- top 20 bar plots ---- ####
```{r}
# stacked bar plot on relative abundance for SARS-CoV-2 positive and negative 
gg_barplot_RA <- plot_bar(physeq_gen_2, x="covid", fill = "Genus", 
                          facet_grid = ~Genus) + geom_bar(aes(color=Genus, fill=Genus), 
                                                          stat="identity", position="stack") + theme(legend.position = "none")

#create data table from ggplot data 
tabRA <- data.table(gg_barplot_RA$data) 
# add relative abundance column 
tabRA[, abundance_sum := sum(Abundance), by = Sample] # sum of abundance per sample
tabRA[, relative_abundance := Abundance/abundance_sum]
tabRA <- tabRA %>% relocate(relative_abundance, .before = Abundance)
tabRA$OTU <- factor(tabRA$OTU, levels = c("Actinomyces", "Alloprevotella", "Atopobium", "Corynebacterium", "Dolosigranulum", "Fusobacterium", "Gemella", "Granulicatella", "Haemophilus", "Lactobacillus", "Leptotrichia", "Megasphaera", "Moraxella", "Neisseria", "Porphyromonas", "Prevotella", "Rothia", "Staphylococcus", "Streptococcus", "Veillonella", "Other"))

covid_stacked2 <- ggplot(tabRA, aes(x = covid, y = relative_abundance, fill = OTU)) +
  geom_col(position = "fill") + 
  labs(fill = "Genus") +
  scale_fill_manual(values=PaletteChoice4)+
  xlab("SARS-CoV-2") + ylab("Relative Abundance") +
  theme(
    #legend.position = "None", 
    text = element_text(size = 13))
covid_stacked2
#ggsave(filename = "figure1.pdf", width = 20, height = 13, units = "cm")

```

```{r}
# stacked bar plot on relative abundance for SARS-CoV-2 severity classes 
gg_barplot_RA_severity <- plot_bar(physeq_gen_2, x="severity", fill = "Genus", 
                                   facet_grid = ~Genus) + geom_bar(aes(color=Genus, fill=Genus), 
                                                                   stat="identity", position="stack") + theme(legend.position = "none")

#create data table from ggplot data of class data.table
tabRASeverity <- data.table(gg_barplot_RA_severity$data) 
# add relative abundance column 
tabRASeverity[, abundance_sum := sum(Abundance), by = Sample] # sum of abundance per sample
tabRASeverity[, relative_abundance := Abundance/abundance_sum]
tabRASeverity <- tabRASeverity %>% relocate(relative_abundance, .before = Abundance)
# level severity 
tabRASeverity$severity <- factor(tabRASeverity$severity, levels = c("negative", "mild", "moderate", "severe", "death"))
# level unique taxa and add other to end 
tabRASeverity$OTU <- factor(tabRASeverity$OTU, levels = c("Actinomyces", "Alloprevotella", "Atopobium", "Corynebacterium", "Dolosigranulum", "Fusobacterium", "Gemella", "Granulicatella", "Haemophilus", "Lactobacillus", "Leptotrichia", "Megasphaera", "Moraxella", "Neisseria", "Porphyromonas", "Prevotella", "Rothia", "Staphylococcus", "Streptococcus", "Veillonella", "Other"))

severity_stacked <- ggplot(tabRASeverity, aes(x = severity, y = relative_abundance, fill = OTU)) +
  geom_col(position = "fill") + 
  labs(fill = "Genus") +
  scale_fill_manual(values=PaletteChoice4)+
  theme(text = element_text(size = 13)) +
  xlab("COVID-19 Severity") + ylab(NULL)
severity_stacked
#ggsave(filename = "figure2.tiff", width = 20, height = 16, units = "cm")
```

```{r}
# merge stacked bar charts for covid and severity 
covid_severity_stacked <- covid_stacked2 + severity_stacked +
  plot_annotation(theme = theme(plot.caption = element_text(hjust = 0)), tag_levels = c("a", "b"))
covid_severity_stacked
#ggsave(filename = "covid_severity_stacked.tiff", height = 6, width = 11)
```


# ---- severity class plots ungrouped ---- ####
```{r}
# subset table and aggregate to genus level
physeq_gen1 <- phyl_utr
physeq_gen1 <- subset_samples(physeq_gen1, severity %in% c("negative", "mild", "moderate", "severe", "death", NA))
physeq_gen1 <- phyloseq::tax_glom(physeq_gen1, "Genus", NArm = FALSE)
physeq_gen1_1 <- microbiome::aggregate_taxa(physeq_gen1, "Genus")
physeq_gen2_2 <- aggregate_top_taxa2(physeq_gen1, top = 20, "Genus") #reduce to genus

physeq_gen1_RA <- microbiome::transform(physeq_gen2_2, "compositional") #compositional
taxa_sums(physeq_gen1_RA)
sample_data(physeq_gen1_RA)$severity <- factor(sample_data(physeq_gen1_RA)$severity)
```

```{r}
# severity plot un-grouped
sample_data(physeq_gen1_RA)$severity <- factor(sample_data(physeq_gen1_RA)$severity, levels = c("negative", "mild", "moderate", "severe", "death"))

{plot.composition.relAbun <- plot_composition(physeq_gen1_RA,  
                                              sample.sort = NULL,
                                              otu.sort = names(sort(
                                                taxa_sums(physeq_gen1_RA), F)[1:21]),
                                              x.label = "severity",
                                              group_by = "severity")
  
  #add settings plot
  plot.composition.relAbun + 
    scale_fill_manual(values=c(PaletteChoice4)) + 
    theme_bw() +
    scale_x_discrete(labels = NULL) +
    ylab("Relative Abundance") +
    theme(legend.justification  = "top",
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(size = 13, face = "bold"),
          axis.title.y = element_text(color="Grey1", size=11, face="bold"),
          axis.text.x =element_text(color="Grey1", size=11, angle = 90),
          axis.text=element_text(color="Grey1", size=11, face="bold"),
          panel.grid.major = element_blank()
          
    )
  
}
#ggsave(filename = "supplementaryFigure3.tiff", width = 36, height = 18, units = "cm")
```



# ---- alpha diversity sars-cov-2 and severity plots ---- ####
```{r}
# chao diversity covid 
shan_cys <- physeq_gen

sample_data(shan_cys)$read_count <- readcount(shan_cys)
shan_cys <- cbind(meta(shan_cys),
                  microbiome::alpha(shan_cys, index = "chao1"))

# export chao index per sample
#write.table(shan_cys, file="chao_covid.txt", sep="\t", row.names = FALSE, col.names = TRUE)

#define colors
cols_group <- c('negative' = PaletteChoice4[20],
                'positive' = PaletteChoice4[8]) 


# plot chao diversity with color dots
alpha_chao1 = {ggplot(shan_cys, aes(covid, chao1)) +
    geom_violin(aes(color = covid)) +
    geom_jitter(aes(color = covid), alpha = 0.75) +
    scale_color_manual("SARS-CoV-2", values = cols_group) +
    theme_bw() +
    #stat_compare_means(comparisons = make_pairs(shan_cys$covid)) +
    #theme(legend.position = "None") +
    ylab("Chao1 Diversity") + xlab("SARS-CoV-2")
}
alpha_chao1

#alpha_chao1 + ylab("Chao1 Diversity") + xlab("SARS-CoV-2")
alpha_chao1

#calculate significance
print(wilcox_result <- wilcox.test(chao1~covid, data = shan_cys,  alternative = "two.sided", paired = F, correct = T, conf.int = T, conf.level = 0.95))
# data:  chao1 by covid
# W = 2125, p-value = 0.1751
# 95 percent confidence interval:
#  -10.000060   1.999996

# save plot as rds 
#write_rds(alpha_chao1, "chaoCovid.Rds")
# and read and plot afterwards:
#read_rds("chaoCovid.Rds")

#ggsave("alphaChao1Covid.tiff", height=4, width=5)
#rm(list = ls(all = T)) # remove and re-run analysis for each alpha diversity plot 
```

```{r}
# inverse simpson diversity (1/d) covid 
shan_cys1 <- physeq_gen

sample_data(shan_cys1)$read_count <- readcount(shan_cys1)
shan_cys1 <- cbind(meta(shan_cys1),
                   microbiome::alpha(shan_cys1, index = "diversity_inverse_simpson"))

# export inverse simpson index per sample
#write.table(shan_cys1, file="inverse_simpson_covid.txt", sep="\t", row.names = FALSE, col.names = TRUE)

# define colors
cols_group <- c('negative' = PaletteChoice4[20],
                'positive' = PaletteChoice4[8])

# plot 1/d  with color dots
alpha_inverse_simpson = {ggplot(shan_cys1, aes(covid, diversity_inverse_simpson)) +
    geom_violin(aes(color = covid)) +
    geom_jitter(aes(color = covid), alpha = 0.75) +
    
    scale_color_manual("SARS-CoV-2", values = cols_group) +
    theme_bw() +
    #theme(legend.position = "None") +
    ylab("Inverse Simpson Diversity (1/D)")+xlab("SARS-CoV-2")
  #stat_compare_means(comparisons = make_pairs(shan_cys1$covid))
}
#alpha_inverse_simpson+ylab("Inverse Simpson Diversity (1/D)")+xlab("SARS-CoV-2")
alpha_inverse_simpson

print(wilcox_result <- wilcox.test(diversity_inverse_simpson~covid, data = shan_cys1,  alternative = "two.sided", paired = F, correct = T, conf.int = T, conf.level = 0.95))

# data:  diversity_inverse_simpson by covid
# W = 2299, p-value = 0.5195
# 95 percent confidence interval:
#  -2.723478  1.338805

# save plot as rds 
#write_rds(alpha_inverse_simpson, "inverseSimpsonCovid.Rds")
# and read and plot afterwards:
#read_rds("inverseSimpsonCovid.Rds")

#ggsave("inverseSimpsonDiversity_2.tiff", height=4, width=5)
#rm(list = ls(all = T)) # remove and re-run analysis for each alpha diversity plot 
``` 

```{r}
# alpha diversity severity
# chao diversity severity
shan_cys <- physeq_gen
sample_data(shan_cys)$read_count <- readcount(shan_cys)
shan_cys <- cbind(meta(shan_cys),
                  microbiome::alpha(shan_cys, index = "chao1"))
# level severity 
shan_cys$severity <- factor(shan_cys$severity, levels = c("negative", "mild", "moderate", "severe", "death"))

# export chao index per sample
#write.table(shan_cys, file="chao_severity.txt", sep="\t", row.names = FALSE, col.names = TRUE)

# define colors
cols_group_npj <- c('negative' = PaletteChoice4[1],
                    'mild' = PaletteChoice4[2],
                    'moderate' = PaletteChoice4[3],
                    'severe' = PaletteChoice4[4],
                    'death' = PaletteChoice4[5])

# plot chao diversity with color dots
chaoSeverity = {ggplot(shan_cys, aes(severity, chao1)) +
    geom_violin(aes(color = severity)) +
    geom_jitter(aes(color = severity), alpha = 0.75) +
    scale_color_manual("Severity", values = cols_group_npj) +
    theme_bw()+
    #theme(legend.position = "None") +
    ylab( "Chao1 Diversity") + xlab("Severity")+
    stat_compare_means(comparisons = make_pairs(shan_cys$severity), method = "wilcox.test")
}
#shan_all+ylab( "Chao1 diversity")
chaoSeverity

# save plot as rds 
#write_rds(chaoSeverity, "chaoSeverity.Rds")
# and read and plot afterwards:
#read_rds("chaoSeverity.Rds")

#ggsave("alphaChao1Severity.tiff", height=4, width=5)
#rm(list = ls(all = T)) # remove and re-run analysis for each alpha diversity plot 
```

```{r}
# inverse simpson 1/d diversity severity
shan_cys1 <- physeq_gen
sample_data(shan_cys1)$read_count <- readcount(shan_cys1)
shan_cys1 <- cbind(meta(shan_cys1),
                   microbiome::alpha(shan_cys1, index = "diversity_inverse_simpson"))
# level severity 
shan_cys1$severity <- factor(shan_cys1$severity, levels = c("negative", "mild", "moderate", "severe", "death"))

# export inverse simpson index per sample
#write.table(shan_cys1, file="inverse_simpson_severity.txt", sep="\t", row.names = FALSE, col.names = TRUE)

# define colors
cols_group_npj <- c('negative' = PaletteChoice4[1],
                    'mild' = PaletteChoice4[2],
                    'moderate' = PaletteChoice4[3],
                    'severe' = PaletteChoice4[4],
                    'death' = PaletteChoice4[5])

# plot 1/d with color dots
inverseSimpsonSeverity = {ggplot(shan_cys1, aes(severity, diversity_inverse_simpson)) +
    geom_violin(aes(color = severity)) +
    geom_jitter(aes(color = severity), alpha = 0.75) +
    scale_color_manual("Severity", values = cols_group_npj) +
    theme_bw() +
    #theme(legend.position = "None") +
    ylab("Inverse Simpson Diversity (1/D)") + xlab("Severity")+
    stat_compare_means(comparisons = make_pairs(shan_cys1$severity), method = "wilcox.test")
}
#inverseSimpsonSeverity+ylab("Inverse Simpson Diversity (1/D)")
inverseSimpsonSeverity
# severe and death were significant p = 0.032 without wilcox.test 
# none significant with wilcox.test 

# save plot as rds 
#write_rds(inverseSimpsonSeverity, "inverseSimpsonSeverity.Rds")
# and read and plot afterwards:
#read_rds("inverseSimpsonSeverity.Rds")

#ggsave("inverseSimpsonSeverity.tiff", height=4, width=5)
#rm(list = ls(all = T)) # remove and re-run analysis for each alpha diversity plot 
```


```{r}
# plot all alpha diversity plots together 
chaoCovid <- readRDS("~/covid-microbiome/plots/alphaDiversity/alphaRds/chaoCovid.Rds")
inverseSimpsonCovid <- readRDS("~/covid-microbiome/plots/alphaDiversity/alphaRds/inverseSimpsonCovid.Rds")
chaoSeverity <- readRDS("~/covid-microbiome/plots/alphaDiversity/alphaRds/chaoSeverity.Rds")
inverseSimpsonSeverity <- readRDS("~/covid-microbiome/plots/alphaDiversity/alphaRds/inverseSimpsonSeverity.Rds")


covid_alpha_all <- chaoCovid + inverseSimpsonCovid + chaoSeverity + inverseSimpsonSeverity +
  plot_annotation(
    theme = theme(plot.caption = element_text(hjust = 0), plot.tag.position = "topleft"),
    tag_levels = c("a", "b", "c", "d"))

covid_alpha_all

#ggsave("~/covid-microbiome/plots/alphaDiversityFigure.tiff", height = 7, width = 11)
```



# ---- ANCOM-BC sars-cov-2 positive and negative script ---- ####
```{r}
# following the selection of covid pos/neg, no pre_antibiotic_use, and no transfers, use this phyloseq object for ANCOM-BC
ancom_motu <- physeq_gen

# optional check to see what is removed during ANCOM data processing
length(taxa_sums(ancom_motu))
ancom_motu <- prune_taxa(taxa_sums(ancom_motu) > 1, ancom_motu) #remove low read taxa
length(taxa_sums(ancom_motu))

# state what to test 
compare_group = "covid"     #group
tax_level = "Genus"         #level of taxonomy
used_formula = "covid" # for covariate correction add a second group 
# adding confounders (age, gender, sars-cov-2 infection status, etc.)

# re-level
sample_data(ancom_motu)$covid <- factor(sample_data(ancom_motu)$covid)
sample_data(ancom_motu)$covid <- relevel(sample_data(ancom_motu)$covid,"positive") 

#group sublevels (positive and negative, mild and death for severity, etc.)
comb1 = c("negative" , "positive")

# loop run ANCOMBC: I give featureID instead of OTUID to the column name so that we can combine it with the taxa_table
# make sure that the taxa table does not have a space between feature and ID

ancom_pval <- matrix(c("empty"),ncol=4,byrow=TRUE)
colnames(ancom_pval) <- c("featureID","grp1", "grp2", "pval") #
ancom_padj <- matrix(c("empty"),ncol=4,byrow=TRUE)
colnames(ancom_padj) <- c("featureID","grp1", "grp2", "padj")

group_list <- list(comb1)
for (comb in group_list){
  print(comb)
  
  #subset table and aggregate to level
  ancom_motu_temp <- subset_samples(ancom_motu, covid %in% c(comb[1], comb[2]))
  ancom_motu_temp = tax_glom(ancom_motu_temp, tax_level)
  
  #Run ancombc function
  #notice that we have to change the lib_cut to 0 read cutoff - change it according to the dataset
  out = ancombc(phyloseq = ancom_motu_temp, formula = used_formula,
                p_adj_method = "fdr", zero_cut = 0.90, lib_cut = 0,
                group = compare_group, struc_zero = TRUE, neg_lb = FALSE,
                tol = 1e-5, max_iter = 100, conserve = TRUE,
                alpha = 0.05, global = TRUE)
  
  #select data for table
  res = out$res$p_val %>% as.data.frame()
  res$OTUID = row.names(res)
  res = dplyr::filter(res, res[1] <0.05)
  res_pval = res
  colnames(res_pval)[colnames(res_pval) == paste(compare_group,comb[1], sep="")] <- "pval"
  ifelse(nrow(res_pval) == 0, print("none"), {
    res_pval$grp1 <- comb[1]
    res_pval$grp2 <- comb[2]
    res_pval<-dplyr::select(res_pval, OTUID, grp1, grp2, pval) %>% as.matrix()
    ancom_pval <-  rbind(ancom_pval, res_pval)})
  rm(res)
  
  res = out$res$q_val %>% as.data.frame()
  res$OTUID = row.names(res)
  res_zero = cbind(res,(out$zero_ind%>% as.data.frame()))
  res_padj = subset(res_zero, (out$res$diff_abn)=='TRUE') 
  res_padj=  dplyr::select (res_padj,
                            as.name(paste(compare_group,comb[1], sep="")), OTUID)
  #(dim(df)[1] == 0)
  colnames(res_padj)[colnames(res_padj) == paste(compare_group,comb[1], sep="")] <- "padj"
  ifelse(nrow(res_padj) == 0, print("none"), {
    res_padj$grp1 <- comb[1]
    res_padj$grp2 <- comb[2]
    res_padj<-dplyr::select(res_padj, OTUID, grp1, grp2, padj) %>% as.matrix()
    ancom_padj <- rbind(ancom_padj, res_padj)})
  rm(out)|rm(res)|rm(res_zero)|rm(res_pval)|rm(res_padj)|rm(comb)
}

{rm(compare_group)|rm(used_formula)|rm(group_list)|rm(comb1)} #clean
rm(ancom_motu_temp)
```

```{r}
# create objects which list the p-values and adjusted p-values of relevant taxa across pos/neg
ancom_pval_covid <- ancom_pval[-c(1),]
ancom_padj_covid <- ancom_padj[-c(1),]
ancom_pval_covid <- as.data.frame(ancom_pval_covid)
ancom_padj_covid <- as.data.frame(ancom_padj_covid)
rownames(ancom_pval_covid) <- NULL
rownames(ancom_padj_covid) <- NULL
```

```{r}
# create taxTableDF object from physeq_gen
taxTableDF <- tax_table(physeq_gen)
taxTableDF <- cbind(newColName = rownames(taxTableDF), taxTableDF)
rownames(taxTableDF) <- 1:nrow(taxTableDF)
taxTableDF <- as.data.frame(taxTableDF)
setnames(taxTableDF, c("newColName"), c("featureID"))
```

```{r}
# write files from ANCOM results, merged on featureID. 
ancom_pval_tax_covid <- merge(ancom_pval_covid, taxTableDF, by = 'featureID')
write.table(ancom_pval_tax_covid, file = "ancom_pval_tax_covid_positive.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

ancom_padj_tax_covid <- merge(ancom_padj_covid, taxTableDF, by = 'featureID')
write.table(ancom_padj_tax_covid, file = "ancom_padj_tax_covid_positive.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
```

# ---- ANCOM-BC severity class script ---- ####
```{r}
# following the selection of covid pos/neg, no pre_antibiotic_use, and no transfers, use this phyloseq object for ANCOM
ancom_motu <- physeq_gen

# optional check to see what is removed during ANCOM-BC data processing
length(taxa_sums(ancom_motu))
ancom_motu <- prune_taxa(taxa_sums(ancom_motu) > 1, ancom_motu) #remove low read taxa
length(taxa_sums(ancom_motu))

#State what to test 
compare_group = "severity" # group
tax_level = "Genus" # level of taxonomy
used_formula = "severity" # for covariate correction add a second group
sample_data(ancom_motu)$severity <- factor(sample_data(ancom_motu)$severity)
sample_data(ancom_motu)$severity <- relevel(sample_data(ancom_motu)$severity, ref = "moderate") # change reference for each analysis of negative, mild, moderate, severe, and death

#group sublevels
comb1 = c("negative" , "moderate")
comb2 = c("mild" , "moderate")
comb3 = c("severe" , "moderate")
comb4 = c("death" , "moderate")

# loop run ANCOMBC: I give featureID instead of OTUID to the column name so that we can combine it with the taxa_table 
# make sure that the taxa table does not have a space between feature and ID

ancom_pval <- matrix(c("empty"),ncol=4,byrow=TRUE)
colnames(ancom_pval) <- c("featureID","grp1", "grp2", "pval") #
ancom_padj <- matrix(c("empty"),ncol=4,byrow=TRUE)
colnames(ancom_padj) <- c("featureID","grp1", "grp2", "padj")

group_list <- list(comb1, comb2, comb3, comb4)
for (comb in group_list){
  print(comb)
  
  #subset table and aggregate to level
  ancom_motu_temp <- subset_samples(ancom_motu, severity %in% c(comb[1], comb[2], comb[3], comb[4]))
  ancom_motu_temp = tax_glom(ancom_motu_temp, tax_level)
  
  #Run ancombc function
  #notice that we have to change the lib_cut to 0 read cutoff - change it according to the dataset. The cutoff was already performed when filtering the phyloseq. 
  out = ancombc(phyloseq = ancom_motu_temp, formula = used_formula,
                p_adj_method = "fdr", zero_cut = 0.90, lib_cut = 0,
                group = compare_group, struc_zero = TRUE, neg_lb = FALSE,
                tol = 1e-5, max_iter = 100, conserve = TRUE,
                alpha = 0.05, global = TRUE)
  
  #select data for table
  res = out$res$p_val %>% as.data.frame()
  res$OTUID = row.names(res)
  res = dplyr::filter(res, res[1] <0.05)
  res_pval = res
  colnames(res_pval)[colnames(res_pval) == paste(compare_group,comb[1], sep="")] <- "pval"
  ifelse(nrow(res_pval) == 0, print("none"), {
    res_pval$grp1 <- comb[1]
    res_pval$grp2 <- comb[2]
    res_pval<-dplyr::select(res_pval, OTUID, grp1, grp2, pval) %>% as.matrix()
    ancom_pval <-  rbind(ancom_pval, res_pval)})
  rm(res)
  
  res = out$res$q_val %>% as.data.frame()
  res$OTUID = row.names(res)
  res_zero = cbind(res,(out$zero_ind%>% as.data.frame()))
  res_padj = subset(res_zero, (out$res$diff_abn)=='TRUE') 
  res_padj=  dplyr::select (res_padj,
                            as.name(paste(compare_group,comb[1], sep="")), OTUID)
  #(dim(df)[1] == 0)
  colnames(res_padj)[colnames(res_padj) == paste(compare_group,comb[1], sep="")] <- "padj"
  ifelse(nrow(res_padj) == 0, print("none"), {
    res_padj$grp1 <- comb[1]
    res_padj$grp2 <- comb[2]
    res_padj<-dplyr::select(res_padj, OTUID, grp1, grp2, padj) %>% as.matrix()
    ancom_padj <- rbind(ancom_padj, res_padj)})
  rm(out)|rm(res)|rm(res_zero)|rm(res_pval)|rm(res_padj)|rm(comb)
}

{rm(compare_group)|rm(used_formula)|rm(group_list)|rm(comb1)} #clean
rm(ancom_motu_temp)
```

```{r}
# create objects which list the p-values and adjusted p-vlues of relevant taxa across severity
ancom_pval_severity <- ancom_pval[-c(1),]
ancom_padj_severity <- ancom_padj[-c(1),]
ancom_pval_severity <- as.data.frame(ancom_pval_severity)
ancom_padj_severity <- as.data.frame(ancom_padj_severity)
rownames(ancom_pval_severity) <- NULL
rownames(ancom_padj_severity) <- NULL
```

```{r}
# create taxTableDF object from physeq_gen
taxTableDF <- tax_table(physeq_gen)
taxTableDF <- cbind(newColName = rownames(taxTableDF), taxTableDF)
rownames(taxTableDF) <- 1:nrow(taxTableDF)
taxTableDF <- as.data.frame(taxTableDF)
setnames(taxTableDF, c("newColName"), c("featureID"))
```

```{r}
# write files from ANCOM-BC results, merged on featureID. 
ancom_pval_tax_severity <- merge(ancom_pval_severity, taxTableDF, by = 'featureID')
write.table(ancom_pval_tax_severity, file = "ancom_pval_tax_severity_moderate.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

ancom_padj_tax_severity <- merge(ancom_padj_severity, taxTableDF, by = 'featureID')
write.table(ancom_padj_tax_severity, file = "ancom_padj_tax_severity_moderate.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
# severe and moderate returned 0 
```


# ---- heatmap on sars-cov-2 positive, negative, and severity classes ---- ####
```{r}
pseq_heat <- physeq_gen_2 # use physeq_gen_2 for top 20 plots 

# pick the core (>0.1% relative abundance in >10% of the samples)
pseq_heat <- aggregate_rare(pseq_heat, level = 'Genus',
                            detection = 0.1/100, prevalence = 10/100)

#find relative abundance across taxa by making compositional
pseq_heat_RA <- microbiome::transform(pseq_heat, "compositional")
rev(sort(taxa_sums(pseq_heat_RA)))
genus_grouped_ra <- c("Streptococcus", "Prevotella", "Corynebacterium", "Veillonella", "Rothia", "Actinomyces", "Staphylococcus", "Haemophilus", "Moraxella", "Neisseria", "Dolosigranulum", "Leptotrichia", "Atopobium", "Gemella", "Granulicatella", "Lactobacillus", "Fusobacterium", "Porphyromonas", "Megasphaera", "Alloprevotella", "Other") # ensure other is at end

# make compositional and ensure taxa with highest RA is at the top 
length(taxa_sums(pseq_heat))
pseq_heat <- pseq.long(pseq_heat, transform.counts = "compositional")
pseq_heat <- pseq_heat[order(pseq_heat$abundance,decreasing = F),]

# order severity
pseq_heat$severity <- factor(pseq_heat$severity, levels = c("negative", "mild", "moderate", "severe", "death"))

# level taxa by relative abundnace and send other to end
pseq_heat$Genus <- factor(pseq_heat$Genus, levels = rev(genus_grouped_ra))

# covid heatmap plot  
heatmap_covid <- ggplot(pseq_heat, aes(x = sample, y = Genus))+
  geom_tile(aes(fill = abundance), color = "white")+
  facet_grid(~covid, scales = "free") +
  scale_x_discrete(label = NULL) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  theme(legend.position = "None", text = element_text(size = 15)) +
  #ggtitle("Genera heatmap log10 of relative abundance") +
  xlab("Sample") +
  #guides(fill = guide_colourbar(title = "Relative Abundance"))+
  rotate_x_text()
#ggsave("pseqheatPosNegGenus1.tiff", height = 8, width = 14)


# severity heatmap plot
heatmap_severity <- ggplot(pseq_heat, aes(x = sample, y = Genus))+
  geom_tile(aes(fill = abundance), color = "white")+
  facet_grid(~severity, scales = "free") +
  scale_x_discrete(label = NULL) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  #ggtitle("Genera heatmap for severity and log10 of relative abundance") +
  xlab("Sample") +
  guides(fill = guide_colourbar(title = "Relative Abundance"))+
  theme(text = element_text(size = 15)) +
  rotate_x_text()
#ggsave("pseqheatSeverityGenus1.tiff", height=8, width=14)

# merge heatmap plots
heatmapFigures <- heatmap_covid + heatmap_severity + 
  plot_annotation(theme = theme(plot.caption = element_text(hjust = 0)), tag_levels = c("a", "b"))

heatmapFigures
#ggsave("heatmapFigure1.tiff", height=8, width=15)
```

# ---- pca and rda on sars-cov-2 positive, negative, and severity classes ---- ####
```{r}
# pca redundancy analysis on CLR transformed taxa using euclidean distance while taking abundance and absence into account on sars-cov-2 positive and negative
pseq_pca <- physeq_gen
pseq_pca <- tax_glom(pseq_pca, "Genus")
ps_clr <- microbiome::transform(pseq_pca, "clr")
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")

# screen plot for component distribution
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

# continue plotting first two components
sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)

# adonis
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 
adon_test <- vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$covid)

pca_covid <- phyloseq::plot_ordination(pseq_pca, ord_clr, type="samples", color="covid") +
  geom_point(size = 3) +theme_bw()+
  stat_ellipse(aes(group = factor(covid)), linetype = 1)+
  scale_color_manual("covid", values = PaletteChoice4[2:1])
#labs("", subtitle = paste0("Pr(>F) = ", round(adon_test$aov.tab$`Pr(>F)`, 3)))
pca_covid

# permanova
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps_clr)$covid)
permutest(dispr, pairwise = TRUE)

ordiYbar(ord_clr, model = "CA")

rm(pseq_pca)|rm(ps_clr)|rm(ord_clr)|rm(clr1)|rm(clr2)|
  rm(clr_dist_matrix)|rm(dispr)|rm(adon_test)

#ggsave("~/covid-microbiome/plots/covid_pca.tiff", width=8, height=4)
```

```{r}
# pca redundancy analysis on CLR transformed taxa using euclidean distance while taking abundance and absence into account on sars-cov-2 severity
pseq_pca <- physeq_gen_1
sample_data(pseq_pca)$severity <- factor(sample_data(pseq_pca)$severity, levels = c("negative", "mild", "moderate", "severe", "death"))
pseq_pca <- tax_glom(pseq_pca, "Genus")
ps_clr <- microbiome::transform(pseq_pca, "clr")
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")

# screen plot for component distribution
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

# continue plotting first two components
sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)

# adonis
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 
adon_test <- vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$severity)

# permanova
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps_clr)$severity)
permutest(dispr, pairwise = TRUE)


pca_severity <- phyloseq::plot_ordination(pseq_pca, ord_clr, type="samples", color="severity") +
  geom_point(size = 3) +theme_bw()+
  stat_ellipse(aes(group = factor(severity)), linetype = 1)+
  scale_color_manual(values = PaletteChoice4[5:1])
#labs("", subtitle = paste0("Pr(>F) = ", round(adon_test$aov.tab$`Pr(>F)`, 3)))
pca_severity

rm(pseq_pca)|rm(ps_clr)|rm(ord_clr)|rm(clr1)|rm(clr2)|
  rm(clr_dist_matrix)|rm(dispr)|rm(adon_test)

#ggsave("~/covid-microbiome/plots/severity_pca.tiff", width=8, height=4)
```

```{r}
# merge sars-cov-2 negative and postive pca with severity pca into one figure 
pca_figures <- pca_covid + pca_severity  + 
  plot_layout(ncol = 1) +
  plot_annotation(theme = theme(plot.caption = element_text(hjust = 0)), tag_levels = c("a", "b"))
pca_figures

#ggsave("pca_figures.tiff", height=6, width=8)
```



#=================================================================
#             Session info
#=================================================================
sessionInfo()

#=================================================================
#             Clean-up
#=================================================================
rm(list = ls(all = T)) 


