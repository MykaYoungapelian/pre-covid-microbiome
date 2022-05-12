# title: "ancomAddMeanRA"
# author: "Youngapelian Myka Jaap"
# date: "2/21/2022"

  
# script to add mean relative abundance to ancom file
  
#libraries & Palettes
library(data.table)
library(tidyr)
library(dplyr)
library(DT)
library(readr)
library(ggplot2)
library(ggsci)
library(kableExtra)
library(DT)
PaletteChoice = pal_npg('nrc')(10)

# this wd
setwd("/Users/mykajaapyoungapelian/covid-microbiome/fullPhyloseq/ancomResults")

# function to add relative abundances to ANCOM result file and create table from output. Check for proper file names. 
add_RA_ANCOM <- function(ancomResultFile.txt, phyloseq_tax_table.txt, phyloseq_otu_table.txt) {
  # load ANCOM result file
  ancomResultData <- read_delim(ancomResultFile.txt, "\t", escape_double = FALSE, trim_ws = TRUE)
  setDT(ancomResultData, key = "featureID")
  
  # load phyloseq taxonomy table 
  phyloseqTax <- read_delim(phyloseq_tax_table.txt, "\t", escape_double = FALSE, trim_ws = TRUE)
  setDT(phyloseqTax, key = "featureID")
  
  # load phyloseq otu table 
  phyloseqOTU <- read_delim(phyloseq_otu_table.txt, "\t", escape_double = FALSE, trim_ws = TRUE)
  setDT(phyloseqOTU, key = "featureID")
  
  #create combined data.table
  otuTaxDT <- merge(phyloseqOTU, phyloseqTax, by = "featureID")
  
  # melt datatable to only include sample, reads, Genus, Species. Drop all other values.
  idVars <- names(otuTaxDT[, 2:(ncol(otuTaxDT)-7)]) # obtain all sample numbers by excluding taxonomy
  
  otuTaxDTMelted <- melt(otuTaxDT, id.vars = c("featureID"), measure.vars = c(idVars))
  
  setnames(otuTaxDTMelted, c("variable", "value"), c("sample", "reads"))
  
  setDT(otuTaxDTMelted)
  
  # subset datatable to remove rows with an abundance of 0. Calculate relative abundance into new column.
  otuTaxDTFull <- otuTaxDTMelted[!reads == 0]
  otuTaxDTRA <- otuTaxDTFull[, relative_abundance := round(reads / sum(reads) * 100, 1), by = sample]
  
  # write files from ANCOM results, merged on featureID with appropriate file name. 
  ancomResultRA <- merge(ancomResultData, otuTaxDTRA, by = 'featureID')
  outputFileName <- unlist(strsplit(ancomResultFile.txt, ".txt"))
  write.table(ancomResultRA, file = paste(outputFileName,"_RA.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE)
}

# testing function streamline: 
# load ANCOM result file
ancomResultFile.txt <- "ancom_padj_tax_severity.txt"

ancomResultData <- read_delim(ancomResultFile.txt, "\t", escape_double = FALSE, trim_ws = TRUE)
setDT(ancomResultData, key = "featureID")

# load phyloseq taxonomy table 
phyloseq_tax_table.txt <- "physeq_gen_tax.txt"
phyloseqTax <- read_delim(phyloseq_tax_table.txt, "\t", escape_double = FALSE, trim_ws = TRUE)
setDT(phyloseqTax, key = "featureID")

# load phyloseq otu table 
phyloseq_otu_table.txt <- "physeq_gen_otu.txt"
phyloseqOTU <- read_delim(phyloseq_otu_table.txt, "\t", escape_double = FALSE, trim_ws = TRUE)
setDT(phyloseqOTU, key = "featureID")

#create combined data.table
otuTaxDT <- merge(phyloseqOTU, phyloseqTax, by = "featureID")



# melt datatable to include negative template values within datatable. Drop all values aside from Genus. 
idVars <- names(otuTaxDT[, 2:(ncol(otuTaxDT)-7)]) # obtain all sample numbers by excluding taxonomy

otuTaxDTMelted <- melt(otuTaxDT, id.vars = c("featureID"), measure.vars = c(idVars))

setnames(otuTaxDTMelted, c("variable", "value"), c("sample", "reads"))

setDT(otuTaxDTMelted)

# subset datatable to remove rows with an abundance of 0. Calculate relative abundance into new column.
otuTaxDTFull <- otuTaxDTMelted[!reads == 0]
otuTaxDTRA <- otuTaxDTFull[, relative_abundance := round(reads / sum(reads) * 100, 1), by = sample]

# merge otuTaxDTRA with the ancom result file to include everything into one table. 
# write files from ANCOM results, merged on featureID with appropriate file name. 
ancomResultRA <- merge(ancomResultData, otuTaxDTRA, by = 'featureID')
outputFileName <- unlist(strsplit(ancomResultFile.txt, ".txt"))
write.table(ancomResultRA, file = paste(outputFileName,"_RA.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE)


# creating tables from ancom analysis with combined relative abundance 
ancom_RA_Table <- function(ancomResultRAFile.txt, supplementaryTableCaptio) {
  ancomResultRAData <- read_delim(ancomResultRAFile.txt, "\t", escape_double = FALSE, trim_ws = TRUE)
  # create data.table
  setDT(ancomResultRAData, key = "featureID")
  # change names of file
  setnames(ancomResultRAData, c("featureID", "grp1", "grp2", "pval", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "sample", "reads", "relative_abundance"), c("OTU", "Group_1", "Group_2", "p_value", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Sample", "16S_Reads", "Relative_Abundance"))
  # change order of columns
  setcolorder(ancomResultRAData, c("Sample", "Group_1", "Group_2", "p_value", "16S_Reads", "Relative_Abundance", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU"))
  
  # create a data.frame
  ancomResultRADF <- as.data.frame(ancomResultRAData)
  # create table with specified column names
  supplementaryTable <- ancomResultRADF %>%
    kbl(col.names = c("Sample", "Group 1", "Group 2", "p_value", "16S Reads", "Relative Abundance", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU"), align = "lccccccccccccc", table.attr = "style =\"color: black; background-color: white;\"", caption = supplementaryTableCaption) %>%
    kable_styling(bootstrap_options = "striped", full_width = FALSE, position = "left")
  # create relevant filename
  outputFileName <- unlist(strsplit(ancomResultRAFile.txt, ".txt"))
  save_kable(supplementaryTable, file = paste(outputFileName,"_Table.pdf", sep = ""))
}


# supplementary table captions
supplementaryTableCaption <- "Supplementary Table 3: P-Adjusted ANCOM analysis output across severity levels with relative abundance per sample."

supplementaryTableCaption <- "Supplementary Table 4: P-Adjusted ANCOM analysis output across total severity levels with relative abundance per sample."

supplementaryTableCaption <- "Supplementary Table 5: p-value unadjusted ANCOM analysis output across COVID positive and negative levels with relative abundance per sample."


rm(list = ls(all = T))


