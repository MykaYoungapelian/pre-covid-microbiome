#---- load packages and set seed ---- ####
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



# ---- functions ---- ####
#Function change N/A to zero
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

# create phyloseq objects 
p_OTU = otu_table(mOTU, taxa_are_rows = TRUE)
p_TAX = tax_table(mTAX)
p_MAP = sample_data(MAP)
p_MAP2 = sample_data(MAP2)
phyloseq <- merge_phyloseq(p_OTU, p_TAX, p_MAP, p_MAP2)




# ---- filter phyloseq object to only contain controls ---- #### ---- 
rank_names(phyloseq)
get_taxa_unique(phyloseq, "Kingdom")
get_taxa_unique(phyloseq)
phyloseq_filtered <- phyloseq
phyloseq_filtered <- subset_samples(phyloseq_filtered, !sample_type == "Throat_swab")
phyloseq_filtered <- subset_samples(phyloseq_filtered, !sample_source == "Human")
phyloseq_filtered <- subset_samples(phyloseq_filtered, !sample_ID_provider == "Swab_Fernanda")
phyloseq_filtered <- subset_samples(phyloseq_filtered, !sample_ID_provider == "Swab_Shu")
phyloseq_filtered <- subset_samples(phyloseq_filtered, !sample_ID_provider == "Swab_Paul")
phyloseq_filtered <- subset_samples(phyloseq_filtered, !sample_ID_provider == "Swab_Neg")
phyloseq_filtered <- subset_samples(phyloseq_filtered, !sample_ID_provider == "Negative_control")
phyloseq_filtered <- subset_samples(phyloseq_filtered, !sample_ID_provider == "NTC")
#save files to produce otuTableControl and taxTableControl
#write.table(otu_table(phyloseq_filtered), file = "phyloseq_filtered_control_otu.txt", sep = "\t") 
#write.table(tax_table(phyloseq_filtered), file = "phyloseq_filtered_control_tax.txt", sep = "\t")

# change file name to be able to run scripts directly. decide which physeq cutoff file to use.
phyl_utr <- phyloseq_filtered
sam_mat = as(sample_data(phyl_utr), "data.frame")
sort(sample_sums(phyl_utr))
#define colors
cols_group <- c('yes' = rgb(0, 115, 189, maxColorValue = 255),
                'no' = rgb(191, 0, 58, maxColorValue = 255))
# check the number of control samples
length(unique(rownames(sam_mat)))


# ---- filtered phyloseq object to only contain negative templates ---- ####
phyloseq_negative_template <- phyloseq

phyloseq_negative_template <- subset_samples(phyloseq_negative_template, !sample_type == "Throat_swab")
phyloseq_negative_template <- subset_samples(phyloseq_negative_template, !sample_ID_provider == "Swab_Fernanda")
phyloseq_negative_template <- subset_samples(phyloseq_negative_template, !sample_ID_provider == "Swab_Shu")
phyloseq_negative_template <- subset_samples(phyloseq_negative_template, !sample_ID_provider == "Swab_Paul")
phyloseq_negative_template <- subset_samples(phyloseq_negative_template, !sample_ID_provider == "Swab_Neg")
phyloseq_negative_template <- subset_samples(phyloseq_negative_template, !sample_ID_provider == "Mock_Standard_community")
phyloseq_negative_template <- subset_samples(phyloseq_negative_template, !sample_ID_provider == "Mock_Log_community")
phyloseq_negative_template <- subset_samples(phyloseq_negative_template, !sample_ID_provider == "Zymobiomics_gDNA")
# save files to produce otuTableNegativeControl and taxTableNegativeControl
#write.table(otu_table(phyloseq_negative_template), file = "phyloseq_negative_template_otu.txt", sep = "\t") 
#write.table(tax_table(phyloseq_negative_template), file = "phyloseq_negative_template_tax.txt", sep = "\t")

# change file name to be able to run scripts directly. decide which physeq cutoff file to use.
phyl_utr <- phyloseq_negative_template
sam_mat = as(sample_data(phyl_utr), "data.frame")
sort(sample_sums(phyl_utr))

# define colors
cols_group <- c('yes' = rgb(0, 115, 189, maxColorValue = 255),
                'no' = rgb(191, 0, 58, maxColorValue = 255))

# check the number of control samples
length(unique(rownames(sam_mat)))

# make overview control data and adjust 
sam_mat = as(sample_data(phyl_utr), "data.frame")
colnames(sam_mat)

# ---- control top 20 plots ---- ####
# subset table and aggregate to genus level
physeq_gen <- phyl_utr
#write.table(otu_table(physeq_gen), file = "physeq_gen_otu_control.txt", sep = "\t")
#write.table(tax_table(physeq_gen), file = "physeq_gen_tax_control.txt", sep = "\t")


physeq_gen <- phyloseq::tax_glom(physeq_gen, "Genus", NArm = FALSE)
physeq_gen_1 <- microbiome::aggregate_taxa(physeq_gen, "Genus")
physeq_gen_2 <- aggregate_top_taxa2(physeq_gen, top = 10,  "Genus") #reduce to genus
# changed to aggregate_top_taxa2 which is another version of aggregate_top_taxa in microbiome package

# aggregate on species level 
physeq_spec <- phyl_utr
physeq_spec <- phyloseq::tax_glom(physeq_spec, "Species", NArm = FALSE)
physeq_spec_1 <- microbiome::aggregate_taxa(physeq_spec, "Species")
physeq_spec_2 <- aggregate_top_taxa2(physeq_spec, top = 10,  "Species")

topSpecies <- names(sort(taxa_sums(physeq_spec), TRUE)[1:10])
Species <- prune_taxa(topSpecies, physeq_spec)
Species


TopGenus <- names(sort(taxa_sums(physeq_gen), TRUE)[1:10])
Genus <- prune_taxa(TopGenus, physeq_gen)
Genus

# RA
physeq_gen_RA <- microbiome::transform(physeq_gen_2, "compositional") # compositional
taxa_sums(physeq_gen_RA)
sample_data(physeq_gen_RA)$sample_ID_provider <- factor(sample_data(physeq_gen_RA)$sample_ID_provider)
#write.table(taxa_sums(physeq_gen_RA), file = "taxa_sums.txt", sep = "\t")

# ---- palette choices ---- ####
PaletteChoice = rev(rainbow(length(taxa_sums(physeq_gen_RA))))   #pick colors
PaletteChoice

PaletteChoice4 = pal_npg('nrc')(10)
PaletteChoice4 = c(PaletteChoice4, "#8c6bb1", "#bebada", "#fb8072", "#1f78b4", "#fb9a99", "#01665e", "#8e0152", "#d1e5f0", "#a50f15", "#4d4d4d", "#fed976")
length(PaletteChoice4)



# ---- plots of controls ---- ####
#plot ungrouped taxa control - npj
{plot.composition.relAbun <- plot_composition(physeq_gen_RA, 
                                              sample.sort = NULL,
                                              otu.sort = names(sort(
                                                taxa_sums(physeq_gen_RA), F)[1:21]),
                                              x.label = "sample",
                                              group_by = "sample_ID_provider")
#add settings plot
plot.composition.relAbun + 
  scale_fill_manual(values=PaletteChoice4) + 
  theme_bw() +
  ggtitle("Control - relative abundance on genus level") +
  ylab("Relative Abundance")+
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
#ggsave(filename = "plotTaxaControlUngrouped_npj.pdf", width = 55, height = 15, units = "cm")

# plot the x-axis in the order of choice
gg_barplot <- plot_bar(physeq_gen_RA, x="sample_ID_provider", fill = "Genus", 
                       facet_grid = ~Genus) + geom_bar(aes(color=Genus, fill=Genus), 
                                                       stat="identity", position="stack") + theme(legend.position = "none")

newtab1 <- data.table(gg_barplot$data)
newtab1$sample_ID_provider <- factor(newtab1$sample_ID_provider)

ggplot(newtab1, aes(y = Abundance, x = sample_ID_provider, fill = Genus)) + facet_grid(.~ Genus) + geom_bar(aes(color=Genus, fill=Genus), stat = 'identity', position="stack") + theme(legend.position = "right") + theme(axis.text.x = element_text(angle = 90))

# create relative abundance without "composition" function
gg_barplot_RA <- plot_bar(physeq_gen_2, x="sample_ID_provider", fill = "Genus", 
                          facet_grid = ~Genus) + geom_bar(aes(color=Genus, fill=Genus), 
                                                          stat="identity", position="stack") + theme(legend.position = "none")
tabRA <- data.table(gg_barplot_RA$data) #create data table from ggplot data of class data.tabel or data.frame

# add relative abundance column 
tabRA[, abundance_sum := sum(Abundance), by = Sample] # sum of abundance per sample
tabRA[, relative_abundance := Abundance/abundance_sum]
tabRA <- tabRA %>% relocate(relative_abundance, .before = Abundance)

# non stacked bargraph with RA between 0 and 1 
ggplot(newtab1, aes(y = Abundance, x = sample_ID_provider, fill = Genus)) + 
  facet_grid(.~ Genus) + geom_col(aes(color=Genus, fill=Genus, position = "fill"), stat = 'identity', position="stack") + 
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle = 90))

# creating stacked bar graph on abundance and sample_type for negative controls 
covid_stacked <- ggplot(tabRA, aes(x = sample_ID_provider, y = relative_abundance, fill = Genus)) +
  geom_col(position = "fill") + 
  scale_fill_manual(values=PaletteChoice4)+
  xlab("sample_ID_provider") + ylab("Relative Abundance")
covid_stacked
#ggsave(filename = "plotControlRAStackedFilled_npj.pdf", width = 30, height = 13, units = "cm")




#=================================================================
#             Session info
#=================================================================
sessionInfo()
rm(list = ls(all = T))

#####################
