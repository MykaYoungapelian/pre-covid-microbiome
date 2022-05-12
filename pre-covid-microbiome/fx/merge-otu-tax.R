---
  title: "merge taxonomy and otu file"
author: "Youngapelian"
date: "2/17/2022"
output: html_document
---
  
  # libraries & Palettes
library(data.table)
library(tidyr)
library(dplyr)
library(DT)
library(readr)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(kableExtra)
library(patchwork)
PaletteChoice = pal_npg('nrc')(10)

mockStdTheoretical <- read_delim("/Users/mykajaapyoungapelian/covid-microbiome/controlPhyloseq/mockStdZymobiomicsTheoretical.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
mockStdLogTheoretical <- read_delim("/Users/mykajaapyoungapelian/covid-microbiome/controlPhyloseq/mockStdLogTheoretical.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# read in taxonomy and otu table. Be sure the column name of all OTUs is listed as 'featureID'. 
otu_table <- read_delim("/Users/mykajaapyoungapelian/covid-microbiome/controlPhyloseq/phyloseq_filtered_control_otu.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
tax_table <- read_delim("/Users/mykajaapyoungapelian/covid-microbiome/controlPhyloseq/phyloseq_filtered_control_tax.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# create data.table objects and merge on OTU as key
otuDT <- data.table(otu_table, key = "featureID")
taxDT <- data.table(tax_table, key = "featureID")
otuTaxDT <- merge(otuDT, taxDT, by = "featureID")

# filter on relevant rows where taxa are present in samples
otuTaxDTFiltered <- otuTaxDT[!`2020.05.COVID.S1000` == 0 | !`2020.05.COVID.S1001` == 0 | !`2020.05.COVID.S1002` == 0 | !`2020.05.COVID.S1003` == 0 | !`2020.05.COVID.S1004` == 0 | !`2020.05.COVID.S1005` == 0 | !`2020.05.COVID.S1045` == 0 | !`2020.05.COVID.S1046` == 0]

# melt datatable to include control values within datatable. Drop all values aside from Genus. 
idVars <- c("2020.05.COVID.S1000", "2020.05.COVID.S1001", "2020.05.COVID.S1002", "2020.05.COVID.S1003", "2020.05.COVID.S1004", "2020.05.COVID.S1005", "2020.05.COVID.S1045", "2020.05.COVID.S1046")

otuTaxDTMelted <- melt(otuTaxDTFiltered, id.vars = c("featureID", "Genus" ), measure.vars = c(idVars))

setnames(otuTaxDTMelted, c("variable", "value"), c("sample", "reads"))

setDT(otuTaxDTMelted)

#Include relevant names for zymobiomics standards into datatable by creating new column "control_type"
otuTaxDTMelted[sample == "2020.05.COVID.S1000" | sample == "2020.05.COVID.S1002" | sample == "2020.05.COVID.S1004", control_type := "Mock_Standard_community"][sample == "2020.05.COVID.S1001" | sample == "2020.05.COVID.S1003" | sample == "2020.05.COVID.S1005", control_type := "Mock_Log_community"][sample == "2020.05.COVID.S1045" | sample == "2020.05.COVID.S1046", control_type := "Zymobiomics_gDNA"]

# subset datatable to remove rows with an abundance of 0. Calculate relative abundance into new column.
otuTaxDTFull <- otuTaxDTMelted[!reads == 0]
otuTaxDTRA <- otuTaxDTFull[, relative_abundance := round(reads / sum(reads) * 100, 1), by = sample]

# recode `Escherichia-Shigella` to `Escherichia` for all samples
otuTaxDTRA[ , Genus := c("Bacillus", "Enterococcus", "Escherichia/Salmonella", "Lactobacillus", "Listeria", "Pseudomonas", "Staphylococcus")[match(Genus, c("Bacillus", "Enterococcus", "Escherichia-Shigella", "Lactobacillus", "Listeria", "Pseudomonas", "Staphylococcus"))]]

# separate into separate datatables. Separate mock log from mock standard and zymobiomics data.tables
# combinded mock standard and zymobiomics 
otuTaxMockZymobiomics <- otuTaxDTRA[!control_type == "Mock_Log_community"]

# mock log alone
otuTaxMockLog <- otuTaxDTRA[!control_type == "Mock_Standard_community"][!control_type == "Zymobiomics_gDNA"]

# mock standard alone 
otuTaxMockStandard <- otuTaxDTRA[!control_type == "Mock_Log_community"][!control_type == "Zymobiomics_gDNA"]

# zymobiomics alone 
otuTaxZymobiomics <- otuTaxDTRA[!control_type == "Mock_Log_community"][!control_type == "Mock_Standard_community"]

# melt mockStdTheoretical and mockStdLogTheoretical and create individual theoretical datatables.
# combined datatable of mock standard and zymobiomics 
setDT(mockStdTheoretical)

mockStdTheoreticalDT <- melt(mockStdTheoretical, id.vars = c("genus"), measure.vars = c("mock_std_16s", "zymobiomics_16s"))
setnames(mockStdTheoreticalDT, c("genus", "variable", "value"), c("Genus", "control_type", "relative_abundance"))

# mock log theoretical datatable
setDT(mockStdLogTheoretical)
setnames(mockStdLogTheoretical, c("16s"), c("mock_log_16s"))
mockStdLogTheoreticalDT <- melt(mockStdLogTheoretical, id.vars = c("genus"), measure.vars = c("mock_log_16s"))
mockStdLogTheoreticalDT <- setnames(mockStdLogTheoreticalDT, c("genus", "variable", "value"), c("Genus", "control_type", "relative_abundance"))

# obtain mock standard theoretical alone 
setDT(mockStdTheoreticalDT)
mockTheoreticalDT <- mockStdTheoreticalDT[!control_type == "zymobiomics_16s"]

# obtain zymobiomics theoretical alone 
zymobiomicsTheoreticalDT <- mockStdTheoreticalDT[!control_type == "mock_std_16s"]


# merge individual datatables for each observed and theoretical set. Merge on... (so far this doesnt work)
# mock standard theoretical and observed 
mockStdObsTheoryDT <- merge(mockTheoreticalDT, otuTaxMockStandard, all = TRUE)
# level and determine order of observed and theoretical
mockStdObsTheoryDT$control_type = factor(mockStdObsTheoryDT$control_type, levels = c("Mock_Standard_community", "mock_std_16s"))

# zymobiomics theoretical and observed 
zymobiomicsObsTheoryDT <- merge(zymobiomicsTheoreticalDT, otuTaxZymobiomics, all = TRUE)
# level and determine order of observed and theoretical 
zymobiomicsObsTheoryDT$control_type = factor(zymobiomicsObsTheoryDT$control_type, levels = c("Zymobiomics_gDNA", "zymobiomics_16s"))

# mock log theoretical and observed 
mockLogObsTheoryDT <- merge(mockStdLogTheoreticalDT, otuTaxMockLog, all = TRUE)
# level and determine order of observed and theoretical
mockLogObsTheoryDT$control_type = factor(mockLogObsTheoryDT$control_type, levels = c("Mock_Log_community", "mock_log_16s"))

# plot individual mock standard observed and theoretical 
# mock standard observed and theoretical  
stackedMockObs <- ggplot(mockStdObsTheoryDT, aes(x = control_type, y = relative_abundance, fill = Genus)) +
  geom_col(position = "fill") + 
  scale_fill_manual(values=PaletteChoice)+
  scale_x_discrete(labels = c("Observed", "Theoretical"))+
  xlab("Mock Standard") + ylab("Relative Abundance") +
  theme(plot.caption = element_text(hjust =0), legend.position = "none", text = element_text(size = 15))
stackedMockObs


# plot individual zymobiomics PCR observed and theoretical 
# PCR mock observed and theoretical 

stackedZymobiomicsObs <- ggplot(zymobiomicsObsTheoryDT, aes(x = control_type, y = relative_abundance, fill = Genus)) +
  geom_col(position = "fill") + 
  scale_fill_manual(values=PaletteChoice)+
  scale_x_discrete(labels = c("Observed", "Theoretical"))+
  xlab("Mock 16S PCR") + ylab(NULL) +
  theme(plot.caption = element_text(hjust =0), legend.position = "none", text = element_text(size = 15))
stackedZymobiomicsObs
# plot individual mock log observed and theoretical 
# mock log observed and theoretical
stackedMockLogObs <- ggplot(mockLogObsTheoryDT, aes(x = control_type, y = relative_abundance, fill = Genus)) +
  geom_col(position = "fill") + 
  scale_fill_manual(values=PaletteChoice)+
  scale_x_discrete(labels = c("Observed", "Theoretical"))+
  xlab("Mock Log") + ylab(NULL) +
  theme(plot.caption = element_text(hjust =0), text = element_text(size = 15))
stackedMockLogObs

# combine control figures into one figure. 
supplementaryFigure1.1 <- stackedMockObs + stackedZymobiomicsObs + stackedMockLogObs + plot_annotation(
  theme = theme(plot.caption = element_text(hjust = 0)), tag_levels = c('a', 'b', 'c'))

supplementaryFigure1.1
#ggsave("supplementaryFigure1.tiff", supplementaryFigure1.1, height = 8, width = 12)

# create table with all data
otuTaxDTFullTable <- otuTaxDTFull
setcolorder(otuTaxDTFullTable, neworder = c("sample", "control_type", "Genus", "reads", "relative_abundance", "featureID"))

#options(knitr.table.format = "html")
otuTaxDF <- as.data.frame(otuTaxDTFullTable)
supplementaryTable2 <- otuTaxDF %>%
  kbl(col.names = c("Sample",  "Sample Type", "Genus Present", "16S Reads", "Relative Abundance", "OTU"), align = "lccccc", table.attr = "style =\"color: black; background-color: white;\"", caption = "Supplementary Table 2: Relative abundance of taxa per sample of controls from the DNA extraction of Nasopharyngeal PCR Samples. Mock_standard_community refers to the standard zymobiomics community, Mock_log_community refers to the standard log zymobiomics community, and Zymobiomics_gdna refers to the control for the PCR.") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE, position = "left")
supplementaryTable2
#save_kable(supplementaryTable2, "controlRelativeAbundanceTable.pdf")

# clean up
rm(list = ls(all = T))


