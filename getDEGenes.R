#!/usr/bin/env Rscript
# Script created by Zoe Clarke, June 2020
# Last updated June 7, 2020
# run chmod +x getDEGenes.R to make executable

options(warn = -1)

require(docopt)
'Usage:
  getDEGenes.R [-i <input> -r <resolution>]
  getDEGenes.R (-h | --help)

Options:
  -h --help Show this screen.
  -i Input file base name.
  -r Resolution (as a number, ex. 0.2).
' -> doc

opts <- docopt(doc)

# Load libraries
library(scClustViz)
library(dplyr)

dir <- getwd() # Make everything relative to working directory

# Load in the data
print("Loading in data...")
load(file = file.path(dir, paste(opts["-i"], "_scClustViz.RData", sep = "")))

# Convert resolutino to number
resNum <- as.numeric(opts["r"])*5

print(paste("Resolution gathered is ", opts["r"], " creating a variable of ", resNum, sep = ""))

# Choose resolution
object <- sCVdata_list[[resNum]]

# Get DEvsRest genes
print("Gathering genes...")
genes <- DEvsRest(object)

allClustGenes <- list() # Create empty list to collect marker genes

maxGenes <- 0 # Set counter for max number of genes

print("Ordering gene lists...")

# Loop through clusters to get marker genes for each cluster ordered by p-value
for (i in names(genes)) {
  clustGenes <- as.data.frame(genes[i]) # Get genes for one cluster
  colnames(clustGenes) <- c("logGER", "Wstat", "pVal", "FDR") # Rename columns for that cluster
  clustGenes[paste("cluster", i, "Genes", sep = "")] <- rownames(clustGenes) # Assign genes to a variable specific for cluster in DF
  # If there is a new max number of marker genes, save this value
  if (length(rownames(clustGenes)) > maxGenes) {
    maxGenes <- length(rownames(clustGenes))
  }
  ranked <- arrange(clustGenes, pVal) # Arrange genes by ascending pVal
  assign(paste("cluster", i, "Genes", sep = ""), value = ranked[paste("cluster", i, "Genes", sep = "")]) # Assign ranked genes to new variable
  allClustGenes <- append(allClustGenes, get(paste("cluster", i, "Genes", sep = "")), after = length(allClustGenes)) # Add to list of marker genes for each cluster
}

print("Converting to gene output...")

# Make all genes same length as maxGenes (the other values will be NAs)
sameLengthGenes <- lapply(allClustGenes, `length<-`, maxGenes)

# Make this into a dataframe
clustGenesDF <- rlist::list.cbind(sameLengthGenes)

# Save DF as CSV
write.csv(clustGenesDF, file = file.path(dir, paste(opts["-i"], "_allClustGenes.csv", sep = "")),
          row.names = FALSE)
