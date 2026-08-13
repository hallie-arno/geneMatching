## Created by Hallie Arno on August 13 2026
## Input: data.frame or data.table of SNPs (requires Chromosome, SNP position, and SNP ID)
## For SNP Id, any unique identifier is fine-- I usually just make a column with CHROM_POSITION as ID
## Chromosomes MUST be numeric. To convert chromosomes to numeric, download the table of chromosome names from NCBI, merge with that table, and use the numeric column.
## Output is a data.table with Chromosome, SNP Position, SNP ID, Gene Name/ Symbol, Gene Start, and Gene end in BP from both NCBI and ensembl.
## The included .csv files can be updated as needed from NCBI or ensembl for new alignments or other species.
## Only data for Salmo salar and Salmo trutta are included here.

library(data.table)
library(this.path)
library(dplyr)

findGenes <- function(snpsDf, pos = "POS", chrom = "CHROM", ID = "ID", species) {
  
  print(paste("Files must be stored in this directory: ", (this.dir())))
  ### Load species genes (downloaded from NCBI and ensembl)
  ### If there are files not fiound errors, manually adjust filepaths here
  if (species == "trout") {
    ensGenes <- fread(paste0(this.dir(), "/genes/SalTrutta.genes.ensembl.csv"))
    ncbiGenes <- fread(paste0(this.dir(), "/genes/SalTrutta.genes.symbol.csv"))
  } else if (species == "salmon") {
    ensGenes <- fread(paste0(this.dir(), ("/genes/SalSalar.genes.ensembl.csv")))
    ncbiGenes <- fread(paste0(this.dir(), ("/genes/SalSalar.genes.symbol.csv")))
  } else {
    print("species must be 'salmon' for Salmo salar or 'trout' for Salmo trutta.")
  } 
  
  ### Rename columns for consistency
  snpDfnames <- snpsDf %>% 
    dplyr::rename("POS" = {{pos}}, 
           "CHROM" = {{chrom}},
           "ID" = {{ID}})
  
  ### Convert to data.table instead of data.frame because it is much faster
  snpsdt <- as.data.table(snpDfnames)
  
  ### Make chromosome and POS columns numeric
  snpsdt[, CHROM := as.numeric(CHROM)]
  snpsdt[, POS := as.numeric(POS)]
  
  ### Check chromosome column to ensure format matches .csv files (numeric)
  if (sum(is.na(snpsdt$CHROM)) > 0) {
    print(paste("WARNING:", sum(is.na(snpsdt$CHROM)), "non-numeric chromosome values detected. The SNPs on non-numeric chromosomes will be dropped"))
  }
  
  ### Merge with both sets of genes as two separate data.tables
  
  Ens <- ensGenes[snpsdt, on = .(Chrom = CHROM, Start <= POS, End >= POS), 
                  .(ID, geneStartEns = x.Start, geneEndEns = x.End, Chrom, EnsID)]
  
  ncbi <- ncbiGenes[snpsdt, on = .(Chrom = CHROM, Start <= POS, End >= POS),
                    .(ID, geneStartNcbi = x.Start, geneEndNcbi = x.End, Chrom, GeneSymbol, POS)]
  
  ### Merge NCBI and ensembl data.tables, select columns in an order that makes sense for easier understanding
  
  m <- merge(Ens, ncbi, all.x = TRUE, all.y = TRUE)[, c("ID", "Chrom", "POS", "EnsID", "GeneSymbol", "geneStartEns", "geneEndEns", "geneStartNcbi", "geneEndNcbi")]
  
  return(m)
  
}




