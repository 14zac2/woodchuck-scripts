# Script created by Zoe Clarke, May 2020

# This script takes orthogroups from orthofinder and gets the woodchuck and human
# gene names into a table

# Set working directory
setwd("/home/zclarke/Documents/AliotoAlignment/ncbi_data")

# Load in required packages
library(dplyr)
library(biomaRt)
#library(biomartr)

# Read in woodchuck orthologues from OrthoFinder
orthogroups <- read.delim(file = "../OrthoFinder/WoodchuckData_proteinSeqs/Results_May03/Orthologues/Orthologues_GCA_901343595.1_MONAX5_protein/GCA_901343595.1_MONAX5_protein__v__GRCh38_latest_protein.tsv",
                          header = TRUE,
                          sep = "\t")
# Table has 3 columns; each orthologue is required to exist in woodchuck
# Multiple orthologues in an orthogroup are separated by commas

# Looking at .gpff file I will need "VERSION" and "locus_tag"
# Read in grepped file
woodchuckGA <- read.delim(file = "proteinAccession_locusTag_woodchuck.txt",
                          header = FALSE)
# Filter protein accession number and locus tag into separate columns
accessionNumber <- grep("VERSION", woodchuckGA$V1, value = TRUE)
accessionNumber <- gsub("VERSION| ", "", accessionNumber) # removes VERSION and whitespace
locusTag <- grep("locus_tag", woodchuckGA$V1, value = TRUE)
locusTag <- gsub(" |/locus_tag=", "", locusTag) # Removes everything except locus tag
woodchuckGA <- cbind.data.frame(accessionNumber, locusTag)

# Save conversion charts
saveRDS(woodchuckGA, "woodchuckConversionChart.RDS")
#saveRDS(humanData, "humanConversionChart.RDS")

# Read in if neccessary
woodchuckGA <- readRDS("woodchuckConversionChart.RDS")
#humanData <- readRDS("humanConversionChart.RDS")

# Now to convert orthogroups to gene names...
# Woodchuck first
woodchuckGeneOrthos <- as.character(orthogroups$GCA_901343595.1_MONAX5_protein)
for (i in 1:nrow(woodchuckGA)) {
  protein <- as.character(woodchuckGA$accessionNumber[i])
  gene <- as.character(woodchuckGA$locusTag[i])
  woodchuckGeneOrthos <- gsub(protein, gene, woodchuckGeneOrthos)
}

# Add this to dataframe
orthogroups$woodchuckGenes <- woodchuckGeneOrthos

# Get "humanLocus" just from the orthogroups file
humanProteins <- as.character(orthogroups$GRCh38_latest_protein)
# Divide orthogroups into single characters
humanProteins <- strsplit(humanProteins, split = ",", fixed = TRUE)
# Remove whitespace from any protein names
humanProteins <- lapply(humanProteins, function(x) gsub(" ", "", as.character(x)))
# Remove decimal
humanProteins <- lapply(humanProteins, function(x) gsub("\\..*", "", as.character(x)))
# Now I have a list of human proteins!!
saveRDS(humanProteins, "humanProteins.RDS")

# Get human data from biomart
useast_ensembl <- useEnsembl(biomart = "ensembl", mirror = "useast")
useast_mart <- useDataset("hsapiens_gene_ensembl", useast_ensembl)
hgncSymbol <- getBM(attributes = c("hgnc_symbol", "refseq_peptide"),
                    filters = "refseq_peptide",
                    values = unlist(humanProteins),
                    mart = useast_mart)

# Gotta get version again
humanVersion <- as.character(orthogroups$GRCh38_latest_protein)
# Divide orthogroups into single characters
humanVersion <- strsplit(humanVersion, split = ",", fixed = TRUE)
# Remove whitespace from any protein names
humanVersion <- lapply(humanVersion, function(x) gsub(" ", "", as.character(x)))
humanVersion <- unlist(humanVersion)

# Make data
humanProteins <- unlist(humanProteins)
humanData <- cbind.data.frame(humanProteins, humanVersion)
colnames(humanData) <- c("refseq_peptide", "humanVersion")
humanData <- left_join(hgncSymbol, humanData, by = "refseq_peptide") # Use left join if want no NAs
saveRDS(humanData, "humanConversionChart.RDS")
# Read in if neccessary
humanData <- readRDS("humanConversionChart.RDS")

# Re-add this to orthogroups
humanGeneOrthos <- as.character(orthogroups$GRCh38_latest_protein)
for  (i in 1:nrow(humanData)) {
  protein <- as.character(humanData$humanVersion[i])
  gene <- as.character(humanData$hgnc_symbol[i])
  humanGeneOrthos <- gsub(protein, gene, humanGeneOrthos)
}
# Much better!
orthogroups$humanGenes <- humanGeneOrthos

# Try to make only unique human info
for (j in 1:nrow(orthogroups)) {
  human <- lapply(strsplit(orthogroups$humanGenes[j], split = ","), function(x) gsub(" ", "", as.character(x)))
  human <- unique(unlist(human))
  orthogroups$uniqueHuman[j] <- paste(human, collapse = ", ")
}

# Reorder columns
orthogroups <- orthogroups %>% dplyr::select(woodchuckGenes, uniqueHuman, humanGenes,
                                             GCA_901343595.1_MONAX5_protein,
                                             GRCh38_latest_protein)

saveRDS(orthogroups, "orthogroupConversion.RDS")
write.csv(orthogroups, "orthogroupConversion.csv")

# Duplicate rows with multiple woodchuck genes
woodchuckSingles <- list()
uniqueHumanGenes <- list()
for (h in 1:nrow(orthogroups)) {
  woodchuck <- lapply(strsplit(orthogroups$woodchuckGenes[h], split = ","), function(x) gsub(" ", "", as.character(x)))
  if (length(woodchuck[[1]]) == 1) {
    woodchuckSingles <- append(woodchuckSingles, woodchuck, after = length(woodchuckSingles))
    uniqueHumanGenes <- append(uniqueHumanGenes, orthogroups$uniqueHuman[h], after = length(uniqueHumanGenes))
  }
  else if (length(woodchuck[[1]]) > 1) {
    for (m in 1:length(woodchuck[[1]])) {
      woodchuckSingles <- append(woodchuckSingles, woodchuck[[1]][m], after = length(woodchuckSingles))
      uniqueHumanGenes <- append(uniqueHumanGenes, orthogroups$uniqueHuman[h], after = length(uniqueHumanGenes))
    }
  }
  else {
    print("Something went wrong!")
  }
}

# Make new dataframe with unique human genes and single woodchuck genes
woodchuckSingles <- data.frame(unlist(woodchuckSingles))
uniqueHumanGenes <- data.frame(unlist(uniqueHumanGenes))
woodchuckHumanGeneList <- cbind.data.frame(woodchuckSingles, uniqueHumanGenes)
colnames(woodchuckHumanGeneList) <- c("woodchuckGenes", "humanGenes")
saveRDS(woodchuckHumanGeneList, "woodchuckHumanGeneList.RDS")
write.csv(woodchuckHumanGeneList, "woodchuckHumanGeneList.csv")


# Do the rest manually