

This function replaces bedtools intersect for finding which genes SNPs
are in for Atlantic salmon and Brown trout and is used directly in R. 
It uses .csv files of the position (start and end) of each gene. It then merges 
the input SNP dataframe with the list of all genes, to return which genes 
contain the SNPs of interest.

Below is an example ofhow to use the function.

First need to source in file, make sure the entire folder is kept. One
issue may be that the entire folder is not in the working directory.
This can be fixed easily by making a copy of the entire folder
(geneMatching) and putting it as a subdirectory in your working
directory. The function requires three packages to be installed: dplyr,
data.table, and this.path.

``` r
library(dplyr) 

source("~/Desktop/Projects/geneMatching/findGenes.R")
```

To test the function, I made a dataframe with some made-up SNP
positions. This is a similar format to how PLINK exports results (such
as FST or GWAS), but column names may be different.

``` r
chroms <- c(1, 5, 16, 24, 26)
position <- c(10000, 750000, 1300, 40000, 88600)
mySNPs <- data.frame(cbind(chroms, position)) %>% 
  mutate(ID = paste(chroms, position, sep = "_"))

print(mySNPs)
```

      chroms position       ID
    1      1    10000  1_10000
    2      5   750000 5_750000
    3     16     1300  16_1300
    4     24    40000 24_40000
    5     26    88600 26_88600

Now to use the function to see if any of the made up SNPs are in genes.
Note that the column names and species (“salmon” or “trout”) must be
specified:

``` r
myGenes <- findGenes(snpsDf = mySNPs, pos = "position", chrom = "chroms", ID = "ID", species = "salmon")
```

    [1] "Files must be stored in this directory:  /Users/hallie/Desktop/Projects/geneMatching"

``` r
print(myGenes)
```

    Key: <ID, Chrom>
             ID Chrom    POS              EnsID   GeneSymbol geneStartEns
         <char> <int>  <int>             <char>       <char>        <int>
    1:  16_1300    16   1300               <NA>         <NA>           NA
    2:  1_10000     1  10000               <NA>         <NA>           NA
    3: 24_40000    24  40000 ENSSSAG00000108255 LOC106584882        28711
    4: 26_88600    26  88600 ENSSSAG00000009116        utp15        87770
    5: 5_750000     5 750000 ENSSSAG00000001494 LOC106604082       738331
       geneEndEns geneStartNcbi geneEndNcbi
            <int>         <int>       <int>
    1:         NA            NA          NA
    2:         NA            NA          NA
    3:     151351          1365      151185
    4:     115464         88254      115464
    5:     804262        738336      795054

Now these can be merged back with results (such as FST or GWAS), plotted
as a manhattan plot, or used for a Gene Ontology Enrichment analysis.
