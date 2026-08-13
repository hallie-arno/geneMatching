
findGenes <- function(snpsDf, pos = "POS", chrom = "CHROM", species) {
  
  ### Load species genes
  if (species == "trout") {
    ensGenes <- fread("SalTrutta.genes.ensembl.csv")
    ncbiGenes <- fread("SalTrutta.genes.symbol.csv")
  } else if (species == "salmon") {
    ensGenes <- fread("SalSalar.genes.ensembl.csv")
    ncbiGenes <- fread("SalSalar.genes.symbol.csv")
  } else {
    print("species must be 'salmon' for Salmo salar or 'trout' for Salmo trutta.")
  } 
  
  ### Rename columns for consistency
  snpDfnames <- snpsDf %>% 
    rename("POS" = {{pos}}, 
           "CHROM" = {{chrom}})
  
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