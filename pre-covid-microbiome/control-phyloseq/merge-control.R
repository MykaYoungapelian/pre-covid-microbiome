


# libraries & Palettes
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


# read in taxonomy and otu table. Be sure the column name of all OTUs is listed as 'featureID'. 
otu_table <- read_delim("/Users/mykajaapyoungapelian/covid-microbiome/controlPhyloseq/phyloseq_negative_template_otu.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
tax_table <- read_delim("/Users/mykajaapyoungapelian/covid-microbiome/controlPhyloseq/phyloseq_negative_template_tax.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# create data.table objects and merge on OTU as key
otuDT <- data.table(otu_table, key = "featureID")
taxDT <- data.table(tax_table, key = "featureID")
otuTaxDT <- merge(otuDT, taxDT, by = "featureID")


# melt datatable to include negative template values within datatable. Drop all values aside from Genus. 
idVars <- names(otuTaxDT[, 2:21])

otuTaxDTMelted <- melt(otuTaxDT, id.vars = c("featureID", "Genus" ), measure.vars = c(idVars))

setnames(otuTaxDTMelted, c("variable", "value"), c("sample", "reads"))

setDT(otuTaxDTMelted)

#Include relevant names for negative control and negative template into datatable by creating new column 'control_type'.
otuTaxDTMelted[sample == "2020.05.COVID.S1010" | sample == "2020.05.COVID.S1011" | sample == "2020.05.COVID.S1013" | sample == "2020.05.COVID.S1014" | sample == "2020.05.COVID.S1015" | sample == "2020.05.COVID.S1016" | sample == "2020.05.COVID.S1017",  control_type := "Negative_control"][sample == "2020.05.COVID.S1021" | sample == "2020.05.COVID.S1022" | sample == "2020.05.COVID.S1023" | sample == "2020.05.COVID.S1026" | sample == "2020.05.COVID.S1028" | sample == "2020.05.COVID.S1029"| sample == "2020.05.COVID.S1037" | sample == "2020.05.COVID.S1038" | sample == "2020.05.COVID.S1040" | sample == "2020.05.COVID.S1041" | sample == "2020.05.COVID.S1042" | sample == "2020.05.COVID.S1043" | sample == "2020.05.COVID.S1044", control_type := "Negative_template"]

# subset datatable to remove rows with an abundance of 0. Calculate relative abundance into new column.
otuTaxDTFull <- otuTaxDTMelted[!reads == 0]
otuTaxDTRA <- otuTaxDTFull[, relative_abundance := round(reads / sum(reads) * 100, 1), by = sample]

# rename and reorder columns 
```{r}
setnames(otuTaxDTRA, c("featureID", "Genus", "sample", "reads", "control_type", "relative_abundance"), c("OTU", "Genus", "Sample", "16S_Reads", "Sample_Type", "Relative_Abundance"))

setcolorder(otuTaxDTRA, c("Sample", "Sample_Type", "Genus", "16S_Reads", "Relative_Abundance", "OTU"))
```

# create table of negative control and negative template 16S read data to export as PDF
#library(knitr)
options(knitr.table.format = "html")
otuTaxDFRA <- as.data.frame(otuTaxDTRA)
supplementaryTable1 <- otuTaxDFRA %>%
  kbl(col.names = c("Sample",  "Sample Type", "Genus Present", "16S Reads", "Relative Abundance", "OTU"), align = "lccccc", table.attr = "style =\"color: black; background-color: white;\"", caption = "Supplementary Table 1: Bacterial presence in negative controls and negative templates from DNA extraction of Nasopharyngeal PCR Samples.") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE, position = "left")

#save_kable(supplementaryTable1, "negativeControlTemplateTable.pdf")

# create a DT table of negative control and negative template 16S read data for browsing
supplementaryTable1.2 <- datatable(otuTaxDTRA, colnames = c("Sample",  "Sample Type", "Genus Present", "OTU", "16S Reads", "Relative_Abundance"), caption = "Supplementary Table 1: Bacterial presence in negative controls and negative templates from DNA extraction of Nasopharyngeal PCR Samples.", filter = "none")
supplementaryTable1.2

rm(list = ls(all = T))

