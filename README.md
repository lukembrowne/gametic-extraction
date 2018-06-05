Modified Twogener gametic extraction
================
Luke Browne
6/2/2018

Overview
--------

This R script provides the code for a modified Twogener gametic extraction used in Browne & Karubian (*in review*). The goal is to determine whether an allele in a biparentally inherited genotype of a seedling came from the maternal or parternal tree. This is useful for determing the separate and combined impacts of pollen and seed dispersal on genetic diversity and structure. We can accomplish this if we have genotype data from biparentally inherited leaf tissue and maternally inherited seedling tissue from the same individual plant. We also need to identify and handle the ambiguous cases where it is unclear whether an allele comes from the paternal or maternal tree.

Reading in data
---------------

-   The script assumes that the genetic data is formatted in Genalex format (<http://biology-assets.anu.edu.au/GenAlEx/Welcome.html>), but the code could be modified to accept a variety of file formats.

-   There must be two separate files for the genotypes of biparentally inherited leaf tissue, and maternally inherited seed tissue. Importantly, the samples in each file must be in the same order!

-   Example data files are available in github repo.

-   Code will also have to be modified accordingly if column names are not are different from the example files.

``` r
# If readGenalex is not installed
    # install.packages("devtools")
    # devtools::install_github("douglasgscofield/readGenalex")
  
  library(readGenalex) 

# Load in a single Genalex format file for the seed genotypes
    seed <- readGenalex("./Seed genotypes BILSA GENALEX 2017_02_17.txt")
    head(seed)
```

    ##   FIELD_NUMBER      PLOT Ob03_1 Ob03_1.2 Ob10_1 Ob10_1.2 Ob04_1 Ob04_1.2
    ## 1      07-0198 Random_01     94      108    150      158    137      141
    ## 2      07-0199 Random_01    108      108    148      154    137      137
    ## 3      07-0200 Random_01    106      110    148      152    137      137
    ## 4      07-0201 Random_01     NA       NA     NA       NA     NA       NA
    ## 5      07-0202 Random_01     NA       NA     NA       NA     NA       NA
    ## 6      07-0203 Random_01     NA       NA     NA       NA     NA       NA
    ##   Ob16_1 Ob16_1.2 Ob12_1 Ob12_1.2 Ob22_1 Ob22_1.2 Ob06_1 Ob06_1.2 Ob07_1
    ## 1    130      138    231      231    199      201    191      191    173
    ## 2    128      130    227      231    199      201    191      191    171
    ## 3    130      138    231      231    199      201    191      191    171
    ## 4     NA       NA     NA       NA     NA       NA     NA       NA     NA
    ## 5     NA       NA     NA       NA     NA       NA     NA       NA     NA
    ## 6     NA       NA     NA       NA     NA       NA     NA       NA     NA
    ##   Ob07_1.2 Ob23_1 Ob23_1.2
    ## 1      177    225      225
    ## 2      181    225      225
    ## 3      173    225      225
    ## 4       NA     NA       NA
    ## 5       NA     NA       NA
    ## 6       NA     NA       NA

``` r
  # Another Genalex format file for the leaf genotypes
    leaf <- readGenalex("./Leaf genotypes BILSA GENALEX 2017_03_11.txt")
    head(leaf)
```

    ##   FIELD_NUMBER      PLOT Ob03_1 Ob03_1.2 Ob10_1 Ob10_1.2 Ob04_1 Ob04_1.2
    ## 1      07-0198 Random_01     94      108    150      158    137      137
    ## 2      07-0199 Random_01    108      108    154      154    137      141
    ## 3      07-0200 Random_01    108      110    152      154    137      137
    ## 4      07-0201 Random_01     94      104    152      154    137      137
    ## 5      07-0202 Random_01    100      108    154      156    137      137
    ## 6      07-0203 Random_01     98      108    156      160    137      139
    ##   Ob16_1 Ob16_1.2 Ob12_1 Ob12_1.2 Ob22_1 Ob22_1.2 Ob06_1 Ob06_1.2 Ob07_1
    ## 1    138      138     NA       NA    167      199    191      191    173
    ## 2    128      128    227      231    175      201    191      195    177
    ## 3    130      138    231      233    175      199    191      191    173
    ## 4    130      138    233      243    163      167    191      193    177
    ## 5    128      138    231      231    193      193    191      191    171
    ## 6    130      130    231      231    175      197    191      193    171
    ##   Ob07_1.2 Ob23_1 Ob23_1.2
    ## 1      173    225      225
    ## 2      181    225      225
    ## 3      177    225      225
    ## 4      177    225      238
    ## 5      173    225      225
    ## 6      181    225      225

``` r
  # Make sure leaf and seed files are in the same order
    # Should return TRUE if they are in the same order
    identical(seed$FIELD_NUMBER, leaf$FIELD_NUMBER)
```

    ## [1] TRUE

``` r
  # Get UTM coordinates / spatial locations of seedlings
    # These data were provided in the Genalex files, as columns appended to the allele calls
    seedling_utms <- attributes(seed)$extra.columns  
```

Extract gametes function
========================

-   The goal is to take the seed and leaf genotypes and extract which allele in the biparentally inherited leaf tissue came from the maternal tree and which came from the paternal tree.

-   The function returns a series of data objects (and optionally prints them to tab-delimited files):
    -   Haploid genetic data of maternal gametes
    -   Haploid genetic data of paternal gametes
    -   Diploid leaf and seed genotypes
    -   List of loci and alleles where gametic assignment is ambiguous
-   In cases where the alleles are ambiguous (heterozygous at the same alleles for both leaf and seed tissue), we can treat them by:
    -   Convert them to missing data (mode = "NA")
    -   Assign partial allele frequencies (the default mode), with a 50% probablity of coming from either maternal or paternal tree (mode = "partial"). This outputs a separate file that lists the loci and alleles of the ambiguous cases.
    -   The information in this file can then be used modify allele frequency or allele count tables to account for these ambiguous cases
-   Importantly, loci with ambiguous gametic assignment are assigned as missing data (NA) in the haploid genetic data output for maternal and paternal gametes! Further manipulation of these datasets are necessary to fully account for these ambiguous cases.

``` r
# Define function 
extract_gametes <- function(leaf, 
                            seed, 
                            seedling_utms, 
                            mode = "partial", 
                            label, 
                            write_files = FALSE){  

  ## Initialize matrix to hold parental genotypes
    paternal <- matrix(NA, nrow = nrow(seed), 
                     ncol = attributes(seed)$n.loci)
    
    maternal <- matrix(NA, nrow = nrow(seed), 
                     ncol = attributes(seed)$n.loci)  
    
  # Initialize counter for ambiguous heterozygotes
    heterozygotes <- 0
    locus_total <- 0
    
  ## Initialize dataframe to hold partial frequencies
    # Column names here will be unique / specific to your study!
    if(mode == "partial"){
      partial <- data.frame(FIELD_NUMBER = NA, 
                                     PLOT = NA, locus = NA, 
                                     allele = NA, freq = NA)
    }
    

    
  ## Loop through each seedling and assign maternal and paternal gamete contributions 
    
  for(i in 1:nrow(seed)){ # Loop through rows
   
  # Print status update   
    cat("\n Working on indv:", seed$FIELD_NUMBER[i], "...\n")
      cat("---------------------------------------------\n")
      
    # Initialize index for locus column
      locus_index <- 1 
    
    ## Loop through loci
      for(locus in attributes(seed)$locus.columns){ 
    
      # Initialize paternal and maternal genotypes  
        paternal_gen <- NA
        maternal_gen <- NA
        
      # Pull out genotypes at seed and leaf tissue
        seed_gen <- as.numeric(seed[i, c(locus, locus + 1)])
        leaf_gen <- as.numeric(leaf[i, c(locus, locus + 1)])
        
        # If interested in testing out test cases:
         
            # Missing data at seed tissue
              # seed_gen <- NA
              # leaf_gen <- c(1,2)
              
            # Missing data at leaf tissue
              # seed_gen <- c(1,2)
              # leaf_gen <- NA
              
            # Mismatch between seed and leaf tissue
              # seed_gen <- c(1,2)
              # leaf_gen <- c(3,4)
              
            # Homozygous at leaf and seed
              # seed_gen <- c(1,1)
              # leaf_gen <- c(1,1)
              
            # Single match between seed tissue and leaf tissue
              # seed_gen <- c(1,2)
              # leaf_gen <- c(2,3)
            
            # Heterozygous at both leaf and seed tissue at same alleles
              # seed_gen <- c(1,2)
              # leaf_gen <- c(1,2)
        
        
        # Find matches between alleles in seed tissue and leaf tissue, 
          # which could be gametes coming from the maternal tree
          # Can be length of 1 or 2
       
           maternal_match_TF <- leaf_gen %in% seed_gen
  
        # If all matches false..
          # Means that there is either missing data at maternal locus
          # Or mismatch between seed and leaf genotypes
          # Either way.. assign as missing data
        if(sum(maternal_match_TF) == 0){
          maternal_gen <- NA
          paternal_gen <- NA
          
          ## Check and print out information about a mismatch..
          if(sum(!is.na(seed_gen)) == 2 & sum(!is.na(leaf_gen)) == 2){
            cat("Leaf and seed mismatch at:", seed$FIELD_NUMBER[i], ": Locus", colnames(seed)[locus], "\n")
          }
        }
          
        ## If homozygous at leaf, one must come from the father.. 
        ## and the same allele from the mother
        ## Must not be just missing data though
          if(leaf_gen[1] == leaf_gen[2] & any(!is.na(leaf_gen))){
            paternal_gen <- leaf_gen[1]
            maternal_gen <- leaf_gen[1]
          }

        ## If there is one match between seed and leaf
        ## That match must come from the maternal tree, and other allele is paternal
          if(sum(maternal_match_TF) == 1){
            maternal_gen <- leaf_gen[maternal_match_TF]
            paternal_gen <- leaf_gen[!maternal_match_TF]
          }
        
        ## If heterozygote at leaf and seed
          if(identical(seed_gen, leaf_gen) & 
             length(unique(c(seed_gen, leaf_gen))) == 2){
  
          # Increment count of heterozygotes at both leaf and seed
            heterozygotes <- heterozygotes + 1 
            
          # If mode is NA, set these ambiguous cases as missing data
            if(mode == "NA"){
              paternal_gen <- NA
              maternal_gen <- NA
            }
          
        
         if(mode == "partial"){
           
           ## Make a DF that has the partial allele frequency for each potential allele
           partial_temp <- data.frame(FIELD_NUMBER = seed$FIELD_NUMBER[i],
                                               PLOT = seed$PLOT[i],
                                               locus =  names(seed)[locus],
                                               allele = leaf_gen,
                                               freq = 0.5, 
                                               stringsAsFactors = FALSE)

  
           ## Add to output
            partial <- rbind(partial, partial_temp)

           ## Save genotype as NA to maintain haploid format and because we don't want to double count alleles..
             paternal_gen <- NA
             maternal_gen <- NA
         }
        
        } # End heterozygote / heterozygote if
      
    # Print status       
       cat("Seed: ", seed_gen, "|| Leaf: ", leaf_gen, "|| Maternal: ", maternal_gen, "|| Paternal:", paternal_gen, "\n")
    
    # Assign alleles into dataframe 
      maternal[i, locus_index] <- maternal_gen
      paternal[i, locus_index] <- paternal_gen
      
    # Increment 
      locus_index <- locus_index + 1
      locus_total <- locus_total + 1
    
     } # End locus loop
    
  } # End row loop
    
    
    
  ### Print out ercentage of locus-ind combinations that were heterozygotes
    
    cat("Percentage of locus-ind combinations that were heterozygotes at same alleles:",
      heterozygotes / locus_total, "\n")
    
   
    
  ### Maternal and paternal allele data frames that are returned from function or written to file
    
    ## Maternal
      maternal_out <- data.frame(FIELD_NUMBER = seed$FIELD_NUMBER, 
                                 PLOT = seed$PLOT,
                                 maternal,
                                 UTM1 = seedling_utms$UTM1,
                                 UTM2 = seedling_utms$UTM2)
      
      colnames(maternal_out)[3:(3+length(attributes(seed)$locus.names)-1)] <- attributes(seed)$locus.names
      
      maternal_out
    
    
    ## Paternal
      paternal_out <- data.frame(FIELD_NUMBER = seed$FIELD_NUMBER, 
                                 PLOT = seed$PLOT,
                                 paternal,
                                 UTM1 = seedling_utms$UTM1,
                                 UTM2 = seedling_utms$UTM2)
      
      colnames(paternal_out)[3:(3+length(attributes(seed)$locus.names)-1)] <- attributes(seed)$locus.names
      
    

    partial <- partial[-1, ] # Remove first NA row in partial DF

    
    # If write_files is true - then write genotype files and partial allele frequencies to file in Genalex format
    if(write_files == TRUE){
      
      
      ########################
      ## WRITING OUTPUT TO FILE
      
      cat("Writing files...\n")
      
      # Suffix for file names
      suffix = paste(Sys.Date(), ".txt", sep = "")
   
      ## Maternal gametes - Write to genalex format
        sink(paste("./maternal_gametes", label, suffix))
          
          ## First line - Number of loci, number of individuals, and total samples per population
          cat(length(attributes(seed)$locus.names), 
              nrow(seed), length(table(seed$PLOT)), as.numeric(table(seed$PLOT)),"\n", sep = "\t")
          
          ## Second line is 3 tabs followed by population names
          cat("","", "", names(table(seed$PLOT)), "\n", sep = "\t")
          
          write.table(maternal_out, sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
          
        sink()
      
      
      ## Paternal gametes - Write to genalex format
        sink(paste("./paternal_gametes", label, suffix))
        
          ## First line - Number of loci, number of individuals, and total samples per population
          cat("8", nrow(leaf), length(table(leaf$PLOT)), as.numeric(table(leaf$PLOT)),"\n", sep = "\t")
          ## Second line is 3 tabs followed by population names
          cat("","", "", names(table(leaf$PLOT)), "\n", sep = "\t")
          
          write.table(paternal_out, sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
          
        sink()
      
      
     ## Write partial allele frequencies to file
      write.table(partial, 
                  paste("./partial_frequencies", label, suffix),
                  sep = "\t", row.names = FALSE)

    } ## End write_files IF
    
    
    ## Return output
    return(list(maternal = maternal_out,
                paternal = paternal_out, 
                leaf = leaf,
                seed = seed,
                partial = partial))

} ## End extract gametes function
```

Run extract gametes function
----------------------------

``` r
gametes_out <- extract_gametes(leaf = leaf, 
                               seed = seed, 
                               seedling_utms =  seedling_utms, 
                               label = "BBS", # Label to include in file names if writing to file
                               mode = "partial", 
                               write_files = TRUE)
```

    ## 
    ##  Working on indv: 07-0198 ...
    ## ---------------------------------------------
    ## Seed:  94 108 || Leaf:  94 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  150 158 || Leaf:  150 158 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 141 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  138 138 || Maternal:  138 || Paternal: 138 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 201 || Leaf:  167 199 || Maternal:  199 || Paternal: 167 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 177 || Leaf:  173 173 || Maternal:  173 || Paternal: 173 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0199 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  148 154 || Leaf:  154 154 || Maternal:  154 || Paternal: 154 
    ## Seed:  137 137 || Leaf:  137 141 || Maternal:  137 || Paternal: 141 
    ## Seed:  128 130 || Leaf:  128 128 || Maternal:  128 || Paternal: 128 
    ## Seed:  227 231 || Leaf:  227 231 || Maternal:  NA || Paternal: NA 
    ## Seed:  199 201 || Leaf:  175 201 || Maternal:  201 || Paternal: 175 
    ## Seed:  191 191 || Leaf:  191 195 || Maternal:  191 || Paternal: 195 
    ## Seed:  171 181 || Leaf:  177 181 || Maternal:  181 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0200 ...
    ## ---------------------------------------------
    ## Seed:  106 110 || Leaf:  108 110 || Maternal:  110 || Paternal: 108 
    ## Seed:  148 152 || Leaf:  152 154 || Maternal:  152 || Paternal: 154 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  231 233 || Maternal:  231 || Paternal: 233 
    ## Seed:  199 201 || Leaf:  175 199 || Maternal:  199 || Paternal: 175 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 173 || Leaf:  173 177 || Maternal:  173 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0201 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  94 104 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  152 154 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  NA NA || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  233 243 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  163 167 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  NA NA || Leaf:  225 238 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0202 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  100 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  154 156 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  NA NA || Leaf:  128 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  193 193 || Maternal:  193 || Paternal: 193 
    ## Seed:  NA NA || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  NA NA || Leaf:  171 173 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0203 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  98 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  156 160 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  137 139 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  NA NA || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  175 197 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  171 181 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0204 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  NA NA || Leaf:  154 158 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  NA NA || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  227 231 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  199 199 || Maternal:  199 || Paternal: 199 
    ## Seed:  NA NA || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  NA NA || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  NA NA || Leaf:  225 238 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0205 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  100 108 || Maternal:  108 || Paternal: 100 
    ## Seed:  148 156 || Leaf:  148 152 || Maternal:  148 || Paternal: 152 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 138 || Leaf:  128 130 || Maternal:  128 || Paternal: 130 
    ## Seed:  NA NA || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  165 201 || Maternal:  NA || Paternal: NA 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 181 || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0206 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  106 108 || Maternal:  108 || Paternal: 106 
    ## Seed:  148 154 || Leaf:  154 156 || Maternal:  154 || Paternal: 156 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  227 231 || Leaf:  227 231 || Maternal:  NA || Paternal: NA 
    ## Seed:  199 201 || Leaf:  197 201 || Maternal:  201 || Paternal: 197 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 181 || Leaf:  171 177 || Maternal:  171 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0207 ...
    ## ---------------------------------------------
    ## Seed:  104 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  148 154 || Leaf:  148 148 || Maternal:  148 || Paternal: 148 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  130 138 || Maternal:  130 || Paternal: 138 
    ## Seed:  231 231 || Leaf:  227 231 || Maternal:  231 || Paternal: 227 
    ## Seed:  199 201 || Leaf:  167 201 || Maternal:  201 || Paternal: 167 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 181 || Leaf:  177 181 || Maternal:  181 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 238 || Maternal:  225 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0208 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  148 154 || Leaf:  148 154 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  227 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 201 || Leaf:  175 201 || Maternal:  201 || Paternal: 175 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 181 || Leaf:  171 177 || Maternal:  171 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0209 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  148 154 || Leaf:  148 154 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  128 128 || Maternal:  128 || Paternal: 128 
    ## Seed:  227 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 201 || Leaf:  175 199 || Maternal:  199 || Paternal: 175 
    ## Seed:  191 191 || Leaf:  191 195 || Maternal:  191 || Paternal: 195 
    ## Seed:  171 181 || Leaf:  177 181 || Maternal:  181 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0001 ...
    ## ---------------------------------------------
    ## Seed:  98 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  150 154 || Leaf:  150 154 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 139 || Leaf:  137 139 || Maternal:  NA || Paternal: NA 
    ## Seed:  128 128 || Leaf:  128 130 || Maternal:  128 || Paternal: 130 
    ## Seed:  231 243 || Leaf:  231 243 || Maternal:  NA || Paternal: NA 
    ## Seed:  193 199 || Leaf:  199 201 || Maternal:  199 || Paternal: 201 
    ## Seed:  185 191 || Leaf:  191 195 || Maternal:  191 || Paternal: 195 
    ## Seed:  173 181 || Leaf:  173 181 || Maternal:  NA || Paternal: NA 
    ## Seed:  225 238 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0004 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  148 156 || Leaf:  156 160 || Maternal:  156 || Paternal: 160 
    ## Seed:  137 137 || Leaf:  137 141 || Maternal:  137 || Paternal: 141 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 243 || Leaf:  231 243 || Maternal:  NA || Paternal: NA 
    ## Seed:  175 199 || Leaf:  195 199 || Maternal:  199 || Paternal: 195 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 181 || Leaf:  171 173 || Maternal:  171 || Paternal: 173 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0005 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  94 108 || Maternal:  108 || Paternal: 94 
    ## Seed:  148 156 || Leaf:  148 154 || Maternal:  148 || Paternal: 154 
    ## Seed:  137 137 || Leaf:  137 139 || Maternal:  137 || Paternal: 139 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 243 || Leaf:  231 243 || Maternal:  NA || Paternal: NA 
    ## Seed:  175 199 || Leaf:  197 199 || Maternal:  199 || Paternal: 197 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 181 || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0006 ...
    ## ---------------------------------------------
    ## Seed:  108 110 || Leaf:  94 110 || Maternal:  110 || Paternal: 94 
    ## Seed:  156 160 || Leaf:  156 156 || Maternal:  156 || Paternal: 156 
    ## Seed:  137 139 || Leaf:  139 139 || Maternal:  139 || Paternal: 139 
    ## Seed:  128 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  223 231 || Maternal:  231 || Paternal: 223 
    ## Seed:  167 199 || Leaf:  197 199 || Maternal:  199 || Paternal: 197 
    ## Seed:  185 191 || Leaf:  185 191 || Maternal:  NA || Paternal: NA 
    ## Seed:  171 177 || Leaf:  171 177 || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0007 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  100 108 || Maternal:  108 || Paternal: 100 
    ## Seed:  148 156 || Leaf:  156 156 || Maternal:  156 || Paternal: 156 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 243 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  175 199 || Leaf:  167 199 || Maternal:  199 || Paternal: 167 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 181 || Leaf:  173 181 || Maternal:  181 || Paternal: 173 
    ## Seed:  225 225 || Leaf:  225 238 || Maternal:  225 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0008 ...
    ## ---------------------------------------------
    ## Seed:  100 108 || Leaf:  100 104 || Maternal:  100 || Paternal: 104 
    ## Seed:  156 156 || Leaf:  148 156 || Maternal:  156 || Paternal: 148 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  128 138 || Maternal:  138 || Paternal: 128 
    ## Seed:  231 245 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 193 || Leaf:  167 199 || Maternal:  167 || Paternal: 199 
    ## Seed:  191 191 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  171 173 || Leaf:  171 179 || Maternal:  171 || Paternal: 179 
    ## Seed:  225 238 || Leaf:  225 236 || Maternal:  225 || Paternal: 236 
    ## 
    ##  Working on indv: 07-0010 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  148 156 || Leaf:  148 156 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 243 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  175 199 || Leaf:  175 203 || Maternal:  175 || Paternal: 203 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 181 || Leaf:  171 177 || Maternal:  171 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 238 || Maternal:  225 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0018 ...
    ## ---------------------------------------------
    ## Seed:  98 108 || Leaf:  98 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  150 154 || Leaf:  150 154 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 139 || Leaf:  137 139 || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  231 231 || Leaf:  220 231 || Maternal:  231 || Paternal: 220 
    ## Seed:  201 203 || Leaf:  201 201 || Maternal:  201 || Paternal: 201 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 171 || Leaf:  171 175 || Maternal:  171 || Paternal: 175 
    ## Seed:  225 238 || Leaf:  225 238 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0020 ...
    ## ---------------------------------------------
    ## Seed:  98 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  150 154 || Leaf:  150 154 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 139 || Leaf:  137 139 || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  231 231 || Leaf:  227 231 || Maternal:  231 || Paternal: 227 
    ## Seed:  201 203 || Leaf:  201 201 || Maternal:  201 || Paternal: 201 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 171 || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  225 238 || Leaf:  225 238 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0021 ...
    ## ---------------------------------------------
    ## Seed:  98 104 || Leaf:  98 104 || Maternal:  NA || Paternal: NA 
    ## Seed:  154 160 || Leaf:  154 156 || Maternal:  154 || Paternal: 156 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 243 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 205 || Leaf:  167 201 || Maternal:  167 || Paternal: 201 
    ## Seed:  191 191 || Leaf:  185 191 || Maternal:  191 || Paternal: 185 
    ## Seed:  173 179 || Leaf:  173 173 || Maternal:  173 || Paternal: 173 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0023 ...
    ## ---------------------------------------------
    ## Seed:  98 108 || Leaf:  108 112 || Maternal:  108 || Paternal: 112 
    ## Seed:  150 154 || Leaf:  154 154 || Maternal:  154 || Paternal: 154 
    ## Seed:  137 139 || Leaf:  137 139 || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  201 203 || Leaf:  203 205 || Maternal:  203 || Paternal: 205 
    ## Seed:  191 193 || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  171 171 || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  225 238 || Leaf:  225 238 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0025 ...
    ## ---------------------------------------------
    ## Seed:  98 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  150 154 || Leaf:  148 154 || Maternal:  154 || Paternal: 148 
    ## Seed:  137 139 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  203 205 || Leaf:  167 203 || Maternal:  203 || Paternal: 167 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 171 || Leaf:  171 181 || Maternal:  171 || Paternal: 181 
    ## Seed:  225 238 || Leaf:  225 238 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0026 ...
    ## ---------------------------------------------
    ## Seed:  108 112 || Leaf:  98 108 || Maternal:  108 || Paternal: 98 
    ## Seed:  154 156 || Leaf:  148 156 || Maternal:  156 || Paternal: 148 
    ## Seed:  129 137 || Leaf:  129 137 || Maternal:  NA || Paternal: NA 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 205 || Leaf:  167 205 || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 171 || Leaf:  171 179 || Maternal:  171 || Paternal: 179 
    ## Seed:  225 238 || Leaf:  225 238 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0027 ...
    ## ---------------------------------------------
    ## Seed:  94 108 || Leaf:  108 110 || Maternal:  108 || Paternal: 110 
    ## Seed:  154 156 || Leaf:  154 154 || Maternal:  154 || Paternal: 154 
    ## Seed:  137 137 || Leaf:  137 141 || Maternal:  137 || Paternal: 141 
    ## Seed:  130 138 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 235 || Maternal:  231 || Paternal: 235 
    ## Seed:  165 199 || Leaf:  165 199 || Maternal:  NA || Paternal: NA 
    ## Seed:  191 193 || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  165 177 || Leaf:  165 177 || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 238 || Maternal:  225 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0028 ...
    ## ---------------------------------------------
    ## Seed:  98 110 || Leaf:  98 104 || Maternal:  98 || Paternal: 104 
    ## Seed:  150 160 || Leaf:  154 160 || Maternal:  160 || Paternal: 154 
    ## Seed:  137 139 || Leaf:  137 141 || Maternal:  137 || Paternal: 141 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 243 || Maternal:  231 || Paternal: 243 
    ## Seed:  165 197 || Leaf:  165 167 || Maternal:  165 || Paternal: 167 
    ## Seed:  185 193 || Leaf:  185 191 || Maternal:  185 || Paternal: 191 
    ## Seed:  171 171 || Leaf:  171 173 || Maternal:  171 || Paternal: 173 
    ## Seed:  225 238 || Leaf:  225 238 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0029 ...
    ## ---------------------------------------------
    ## Seed:  106 110 || Leaf:  108 110 || Maternal:  110 || Paternal: 108 
    ## Seed:  150 158 || Leaf:  150 154 || Maternal:  150 || Paternal: 154 
    ## Seed:  137 143 || Leaf:  141 143 || Maternal:  143 || Paternal: 141 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  165 197 || Leaf:  165 199 || Maternal:  165 || Paternal: 199 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 175 || Leaf:  171 177 || Maternal:  171 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0030 ...
    ## ---------------------------------------------
    ## Seed:  98 110 || Leaf:  98 106 || Maternal:  98 || Paternal: 106 
    ## Seed:  150 160 || Leaf:  150 158 || Maternal:  150 || Paternal: 158 
    ## Seed:  137 139 || Leaf:  139 139 || Maternal:  139 || Paternal: 139 
    ## Seed:  130 130 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  231 231 || Leaf:  231 243 || Maternal:  231 || Paternal: 243 
    ## Seed:  165 197 || Leaf:  165 167 || Maternal:  165 || Paternal: 167 
    ## Seed:  185 193 || Leaf:  185 191 || Maternal:  185 || Paternal: 191 
    ## Seed:  171 171 || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0031 ...
    ## ---------------------------------------------
    ## Seed:  108 110 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  156 160 || Leaf:  154 156 || Maternal:  156 || Paternal: 154 
    ## Seed:  137 139 || Leaf:  139 141 || Maternal:  139 || Paternal: 141 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 199 || Leaf:  175 199 || Maternal:  199 || Paternal: 175 
    ## Seed:  185 191 || Leaf:  185 193 || Maternal:  185 || Paternal: 193 
    ## Seed:  171 177 || Leaf:  173 177 || Maternal:  177 || Paternal: 173 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0037 ...
    ## ---------------------------------------------
    ## Seed:  94 106 || Leaf:  94 108 || Maternal:  94 || Paternal: 108 
    ## Seed:  148 156 || Leaf:  148 156 || Maternal:  NA || Paternal: NA 
    ## Seed:  129 143 || Leaf:  137 143 || Maternal:  143 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  231 231 || Leaf:  227 231 || Maternal:  231 || Paternal: 227 
    ## Seed:  165 167 || Leaf:  165 199 || Maternal:  165 || Paternal: 199 
    ## Seed:  185 193 || Leaf:  185 191 || Maternal:  185 || Paternal: 191 
    ## Seed:  175 179 || Leaf:  171 179 || Maternal:  179 || Paternal: 171 
    ## Seed:  225 236 || Leaf:  236 236 || Maternal:  236 || Paternal: 236 
    ## 
    ##  Working on indv: 07-0043 ...
    ## ---------------------------------------------
    ## Seed:  100 108 || Leaf:  100 100 || Maternal:  100 || Paternal: 100 
    ## Seed:  156 156 || Leaf:  156 156 || Maternal:  156 || Paternal: 156 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 245 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 193 || Leaf:  167 199 || Maternal:  167 || Paternal: 199 
    ## Seed:  191 191 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  171 173 || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  225 238 || Leaf:  238 238 || Maternal:  238 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0044 ...
    ## ---------------------------------------------
    ## Seed:  98 110 || Leaf:  108 110 || Maternal:  110 || Paternal: 108 
    ## Seed:  150 150 || Leaf:  150 154 || Maternal:  150 || Paternal: 154 
    ## Seed:  137 139 || Leaf:  137 139 || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  231 238 || Leaf:  231 243 || Maternal:  231 || Paternal: 243 
    ## Seed:  165 197 || Leaf:  197 199 || Maternal:  197 || Paternal: 199 
    ## Seed:  185 193 || Leaf:  191 193 || Maternal:  193 || Paternal: 191 
    ## Seed:  171 171 || Leaf:  171 173 || Maternal:  171 || Paternal: 173 
    ## Seed:  225 238 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0045 ...
    ## ---------------------------------------------
    ## Seed:  108 110 || Leaf:  106 110 || Maternal:  110 || Paternal: 106 
    ## Seed:  148 154 || Leaf:  154 154 || Maternal:  154 || Paternal: 154 
    ## Seed:  137 141 || Leaf:  137 141 || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  177 179 || Leaf:  171 177 || Maternal:  177 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0046 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  98 104 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  148 154 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  139 141 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  231 242 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  199 205 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  NA NA || Leaf:  173 173 || Maternal:  173 || Paternal: 173 
    ## Seed:  NA NA || Leaf:  225 236 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0047 ...
    ## ---------------------------------------------
    ## Seed:  98 108 || Leaf:  98 100 || Maternal:  98 || Paternal: 100 
    ## Seed:  150 154 || Leaf:  154 160 || Maternal:  154 || Paternal: 160 
    ## Seed:  137 139 || Leaf:  137 139 || Maternal:  NA || Paternal: NA 
    ## Seed:  128 128 || Leaf:  128 130 || Maternal:  128 || Paternal: 130 
    ## Seed:  231 243 || Leaf:  231 242 || Maternal:  231 || Paternal: 242 
    ## Seed:  193 199 || Leaf:  193 205 || Maternal:  193 || Paternal: 205 
    ## Seed:  185 191 || Leaf:  185 193 || Maternal:  185 || Paternal: 193 
    ## Seed:  173 181 || Leaf:  165 173 || Maternal:  173 || Paternal: 165 
    ## Seed:  225 236 || Leaf:  225 236 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0048 ...
    ## ---------------------------------------------
    ## Seed:  106 110 || Leaf:  106 108 || Maternal:  106 || Paternal: 108 
    ## Seed:  150 158 || Leaf:  154 158 || Maternal:  158 || Paternal: 154 
    ## Seed:  137 143 || Leaf:  137 141 || Maternal:  137 || Paternal: 141 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  165 197 || Leaf:  165 199 || Maternal:  165 || Paternal: 199 
    ## Seed:  191 193 || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  171 175 || Leaf:  175 177 || Maternal:  175 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0049 ...
    ## ---------------------------------------------
    ## Seed:  100 108 || Leaf:  98 108 || Maternal:  108 || Paternal: 98 
    ## Seed:  152 160 || Leaf:  156 160 || Maternal:  160 || Paternal: 156 
    ## Seed:  139 141 || Leaf:  141 143 || Maternal:  141 || Paternal: 143 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  NA NA || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  165 165 || Leaf:  165 175 || Maternal:  165 || Paternal: 175 
    ## Seed:  185 193 || Leaf:  185 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  173 177 || Leaf:  173 177 || Maternal:  NA || Paternal: NA 
    ## Seed:  225 238 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0050 ...
    ## ---------------------------------------------
    ## Seed:  98 108 || Leaf:  98 104 || Maternal:  98 || Paternal: 104 
    ## Seed:  150 154 || Leaf:  154 160 || Maternal:  154 || Paternal: 160 
    ## Seed:  137 139 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 128 || Leaf:  128 138 || Maternal:  128 || Paternal: 138 
    ## Seed:  231 243 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  193 199 || Leaf:  199 205 || Maternal:  199 || Paternal: 205 
    ## Seed:  185 191 || Leaf:  185 191 || Maternal:  NA || Paternal: NA 
    ## Seed:  173 181 || Leaf:  173 173 || Maternal:  173 || Paternal: 173 
    ## Seed:  225 236 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0051 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  NA NA || Leaf:  148 152 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  NA NA || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  167 199 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  NA NA || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0052 ...
    ## ---------------------------------------------
    ## Seed:  94 108 || Leaf:  94 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  148 158 || Leaf:  148 156 || Maternal:  148 || Paternal: 156 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  227 231 || Leaf:  231 233 || Maternal:  231 || Paternal: 233 
    ## Seed:  167 199 || Leaf:  175 199 || Maternal:  199 || Paternal: 175 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 177 || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  225 238 || Leaf:  225 238 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0053 ...
    ## ---------------------------------------------
    ## Seed:  94 108 || Leaf:  94 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  148 158 || Leaf:  156 158 || Maternal:  158 || Paternal: 156 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  227 231 || Leaf:  231 242 || Maternal:  231 || Paternal: 242 
    ## Seed:  167 199 || Leaf:  167 167 || Maternal:  167 || Paternal: 167 
    ## Seed:  191 191 || Leaf:  185 191 || Maternal:  191 || Paternal: 185 
    ## Seed:  171 177 || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  225 238 || Leaf:  225 238 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0054 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  106 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  150 154 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  137 139 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  199 201 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  177 181 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0055 ...
    ## ---------------------------------------------
    ## Seed:  94 106 || Leaf:  94 108 || Maternal:  94 || Paternal: 108 
    ## Seed:  148 156 || Leaf:  148 148 || Maternal:  148 || Paternal: 148 
    ## Seed:  129 143 || Leaf:  129 141 || Maternal:  129 || Paternal: 141 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  165 167 || Leaf:  165 199 || Maternal:  165 || Paternal: 199 
    ## Seed:  185 193 || Leaf:  185 191 || Maternal:  185 || Paternal: 191 
    ## Seed:  175 179 || Leaf:  175 177 || Maternal:  175 || Paternal: 177 
    ## Seed:  225 236 || Leaf:  225 236 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0056 ...
    ## ---------------------------------------------
    ## Seed:  108 110 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 139 || Leaf:  137 139 || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  185 191 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  171 177 || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0057 ...
    ## ---------------------------------------------
    ## Seed:  94 108 || Leaf:  94 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  148 156 || Leaf:  148 156 || Maternal:  NA || Paternal: NA 
    ## Seed:  129 143 || Leaf:  129 141 || Maternal:  129 || Paternal: 141 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  165 167 || Leaf:  165 199 || Maternal:  165 || Paternal: 199 
    ## Seed:  185 193 || Leaf:  185 185 || Maternal:  185 || Paternal: 185 
    ## Seed:  175 181 || Leaf:  179 181 || Maternal:  181 || Paternal: 179 
    ## Seed:  225 236 || Leaf:  225 236 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0058 ...
    ## ---------------------------------------------
    ## Seed:  94 110 || Leaf:  106 110 || Maternal:  110 || Paternal: 106 
    ## Seed:  148 156 || Leaf:  154 156 || Maternal:  156 || Paternal: 154 
    ## Seed:  137 139 || Leaf:  137 139 || Maternal:  NA || Paternal: NA 
    ## Seed:  128 128 || Leaf:  128 130 || Maternal:  128 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 199 || Leaf:  165 199 || Maternal:  199 || Paternal: 165 
    ## Seed:  185 191 || Leaf:  185 193 || Maternal:  185 || Paternal: 193 
    ## Seed:  171 171 || Leaf:  165 171 || Maternal:  171 || Paternal: 165 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0060 ...
    ## ---------------------------------------------
    ## Seed:  94 110 || Leaf:  110 112 || Maternal:  110 || Paternal: 112 
    ## Seed:  148 156 || Leaf:  156 156 || Maternal:  156 || Paternal: 156 
    ## Seed:  137 139 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  227 231 || Maternal:  231 || Paternal: 227 
    ## Seed:  199 199 || Leaf:  165 199 || Maternal:  199 || Paternal: 165 
    ## Seed:  185 191 || Leaf:  185 193 || Maternal:  185 || Paternal: 193 
    ## Seed:  171 171 || Leaf:  171 177 || Maternal:  171 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 238 || Maternal:  225 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0061 ...
    ## ---------------------------------------------
    ## Seed:  106 106 || Leaf:  106 112 || Maternal:  106 || Paternal: 112 
    ## Seed:  150 160 || Leaf:  148 150 || Maternal:  150 || Paternal: 148 
    ## Seed:  137 141 || Leaf:  139 141 || Maternal:  141 || Paternal: 139 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  197 203 || Leaf:  165 203 || Maternal:  203 || Paternal: 165 
    ## Seed:  191 191 || Leaf:  185 191 || Maternal:  191 || Paternal: 185 
    ## Seed:  173 181 || Leaf:  173 177 || Maternal:  173 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0062 ...
    ## ---------------------------------------------
    ## Seed:  108 112 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  154 156 || Leaf:  154 154 || Maternal:  154 || Paternal: 154 
    ## Seed:  129 137 || Leaf:  129 137 || Maternal:  NA || Paternal: NA 
    ## Seed:  128 130 || Leaf:  128 138 || Maternal:  128 || Paternal: 138 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  167 205 || Leaf:  167 167 || Maternal:  167 || Paternal: 167 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 171 || Leaf:  171 177 || Maternal:  171 || Paternal: 177 
    ## Seed:  225 238 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0064 ...
    ## ---------------------------------------------
    ## Seed:  94 106 || Leaf:  94 108 || Maternal:  94 || Paternal: 108 
    ## Seed:  148 156 || Leaf:  148 156 || Maternal:  NA || Paternal: NA 
    ## Seed:  129 143 || Leaf:  141 143 || Maternal:  143 || Paternal: 141 
    ## Seed:  130 130 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  229 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  165 201 || Leaf:  165 199 || Maternal:  165 || Paternal: 199 
    ## Seed:  185 193 || Leaf:  191 193 || Maternal:  193 || Paternal: 191 
    ## Seed:  175 179 || Leaf:  177 179 || Maternal:  179 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 236 || Maternal:  225 || Paternal: 236 
    ## 
    ##  Working on indv: 07-0065 ...
    ## ---------------------------------------------
    ## Seed:  100 108 || Leaf:  100 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  156 160 || Leaf:  154 156 || Maternal:  156 || Paternal: 154 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 243 || Leaf:  231 243 || Maternal:  NA || Paternal: NA 
    ## Seed:  165 197 || Leaf:  197 199 || Maternal:  197 || Paternal: 199 
    ## Seed:  185 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  177 179 || Leaf:  177 179 || Maternal:  NA || Paternal: NA 
    ## Seed:  238 240 || Leaf:  225 238 || Maternal:  238 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0068 ...
    ## ---------------------------------------------
    ## Seed:  100 108 || Leaf:  94 108 || Maternal:  108 || Paternal: 94 
    ## Seed:  160 160 || Leaf:  154 160 || Maternal:  160 || Paternal: 154 
    ## Seed:  141 141 || Leaf:  137 141 || Maternal:  141 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  231 245 || Maternal:  231 || Paternal: 245 
    ## Seed:  165 165 || Leaf:  165 175 || Maternal:  165 || Paternal: 175 
    ## Seed:  185 191 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  173 175 || Leaf:  171 173 || Maternal:  173 || Paternal: 171 
    ## Seed:  225 238 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0069 ...
    ## ---------------------------------------------
    ## Seed:  94 106 || Leaf:  94 108 || Maternal:  94 || Paternal: 108 
    ## Seed:  148 156 || Leaf:  148 152 || Maternal:  148 || Paternal: 152 
    ## Seed:  129 143 || Leaf:  129 137 || Maternal:  129 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  229 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  165 201 || Leaf:  165 165 || Maternal:  165 || Paternal: 165 
    ## Seed:  185 193 || Leaf:  191 193 || Maternal:  193 || Paternal: 191 
    ## Seed:  175 179 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 236 || Maternal:  225 || Paternal: 236 
    ## 
    ##  Working on indv: 07-0073 ...
    ## ---------------------------------------------
    ## Seed:  94 106 || Leaf:  106 108 || Maternal:  106 || Paternal: 108 
    ## Seed:  148 156 || Leaf:  154 156 || Maternal:  156 || Paternal: 154 
    ## Seed:  129 143 || Leaf:  129 141 || Maternal:  129 || Paternal: 141 
    ## Seed:  130 130 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  229 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  185 193 || Leaf:  185 191 || Maternal:  185 || Paternal: 191 
    ## Seed:  175 179 || Leaf:  179 181 || Maternal:  179 || Paternal: 181 
    ## Seed:  225 225 || Leaf:  225 236 || Maternal:  225 || Paternal: 236 
    ## 
    ##  Working on indv: 07-0085 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  154 154 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 141 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  128 128 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  197 205 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 193 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  173 181 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0086 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  106 108 || Maternal:  108 || Paternal: 106 
    ## Seed:  154 154 || Leaf:  154 156 || Maternal:  154 || Paternal: 156 
    ## Seed:  137 141 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 128 || Leaf:  128 130 || Maternal:  128 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  197 205 || Leaf:  197 197 || Maternal:  197 || Paternal: 197 
    ## Seed:  191 193 || Leaf:  193 193 || Maternal:  193 || Paternal: 193 
    ## Seed:  173 181 || Leaf:  173 177 || Maternal:  173 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0088 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  98 108 || Maternal:  108 || Paternal: 98 
    ## Seed:  154 154 || Leaf:  148 154 || Maternal:  154 || Paternal: 148 
    ## Seed:  137 141 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 128 || Leaf:  128 138 || Maternal:  128 || Paternal: 138 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 201 || Leaf:  197 199 || Maternal:  199 || Paternal: 197 
    ## Seed:  191 193 || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  173 181 || Leaf:  177 181 || Maternal:  181 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0089 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  154 154 || Leaf:  154 158 || Maternal:  154 || Paternal: 158 
    ## Seed:  137 141 || Leaf:  137 141 || Maternal:  NA || Paternal: NA 
    ## Seed:  128 128 || Leaf:  128 130 || Maternal:  128 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  197 205 || Leaf:  199 205 || Maternal:  205 || Paternal: 199 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 181 || Leaf:  177 181 || Maternal:  181 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0091 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  94 94 || Maternal:  94 || Paternal: 94 
    ## Seed:  NA NA || Leaf:  148 158 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  NA NA || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  NA NA || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  199 199 || Maternal:  199 || Paternal: 199 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  NA NA || Leaf:  225 238 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0092 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  104 108 || Maternal:  108 || Paternal: 104 
    ## Seed:  154 154 || Leaf:  154 156 || Maternal:  154 || Paternal: 156 
    ## Seed:  137 141 || Leaf:  137 139 || Maternal:  137 || Paternal: 139 
    ## Seed:  128 128 || Leaf:  128 130 || Maternal:  128 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  197 205 || Leaf:  167 197 || Maternal:  197 || Paternal: 167 
    ## Seed:  191 193 || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  173 181 || Leaf:  177 181 || Maternal:  181 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0093 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  154 154 || Leaf:  154 154 || Maternal:  154 || Paternal: 154 
    ## Seed:  137 141 || Leaf:  137 141 || Maternal:  NA || Paternal: NA 
    ## Seed:  128 128 || Leaf:  128 130 || Maternal:  128 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  197 205 || Leaf:  167 197 || Maternal:  197 || Paternal: 167 
    ## Seed:  191 193 || Leaf:  193 193 || Maternal:  193 || Paternal: 193 
    ## Seed:  173 181 || Leaf:  173 177 || Maternal:  173 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 238 || Maternal:  225 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0094 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  106 108 || Maternal:  108 || Paternal: 106 
    ## Seed:  154 154 || Leaf:  154 156 || Maternal:  154 || Paternal: 156 
    ## Seed:  137 141 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 128 || Leaf:  128 128 || Maternal:  128 || Paternal: 128 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  197 205 || Leaf:  199 205 || Maternal:  205 || Paternal: 199 
    ## Seed:  191 193 || Leaf:  185 191 || Maternal:  191 || Paternal: 185 
    ## Seed:  173 181 || Leaf:  177 181 || Maternal:  181 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0095 ...
    ## ---------------------------------------------
    ## Seed:  98 104 || Leaf:  98 98 || Maternal:  98 || Paternal: 98 
    ## Seed:  156 156 || Leaf:  150 156 || Maternal:  156 || Paternal: 150 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 201 || Leaf:  167 201 || Maternal:  201 || Paternal: 167 
    ## Seed:  191 191 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  173 181 || Leaf:  173 179 || Maternal:  173 || Paternal: 179 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0099 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  148 154 || Leaf:  148 154 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 141 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  128 128 || Leaf:  128 138 || Maternal:  128 || Paternal: 138 
    ## Seed:  231 231 || Leaf:  231 238 || Maternal:  231 || Paternal: 238 
    ## Seed:  167 197 || Leaf:  167 197 || Maternal:  NA || Paternal: NA 
    ## Seed:  191 193 || Leaf:  193 193 || Maternal:  193 || Paternal: 193 
    ## Seed:  173 181 || Leaf:  177 181 || Maternal:  181 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 238 || Maternal:  225 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0100 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  106 108 || Maternal:  108 || Paternal: 106 
    ## Seed:  154 154 || Leaf:  154 156 || Maternal:  154 || Paternal: 156 
    ## Seed:  137 141 || Leaf:  137 139 || Maternal:  137 || Paternal: 139 
    ## Seed:  128 128 || Leaf:  128 130 || Maternal:  128 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 243 || Maternal:  231 || Paternal: 243 
    ## Seed:  197 205 || Leaf:  197 203 || Maternal:  197 || Paternal: 203 
    ## Seed:  191 193 || Leaf:  185 193 || Maternal:  193 || Paternal: 185 
    ## Seed:  173 181 || Leaf:  173 173 || Maternal:  173 || Paternal: 173 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0102 ...
    ## ---------------------------------------------
    ## Seed:  98 104 || Leaf:  98 106 || Maternal:  98 || Paternal: 106 
    ## Seed:  154 160 || Leaf:  154 156 || Maternal:  154 || Paternal: 156 
    ## Seed:  137 137 || Leaf:  137 141 || Maternal:  137 || Paternal: 141 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 243 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 205 || Leaf:  203 205 || Maternal:  205 || Paternal: 203 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 179 || Leaf:  171 179 || Maternal:  179 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 240 || Maternal:  225 || Paternal: 240 
    ## 
    ##  Working on indv: 07-0103 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  154 154 || Leaf:  148 154 || Maternal:  154 || Paternal: 148 
    ## Seed:  137 141 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 128 || Leaf:  128 138 || Maternal:  128 || Paternal: 138 
    ## Seed:  231 231 || Leaf:  231 243 || Maternal:  231 || Paternal: 243 
    ## Seed:  197 205 || Leaf:  197 197 || Maternal:  197 || Paternal: 197 
    ## Seed:  191 193 || Leaf:  193 193 || Maternal:  193 || Paternal: 193 
    ## Seed:  173 181 || Leaf:  177 181 || Maternal:  181 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0109 ...
    ## ---------------------------------------------
    ## Seed:  106 108 || Leaf:  100 106 || Maternal:  106 || Paternal: 100 
    ## Seed:  148 150 || Leaf:  150 154 || Maternal:  150 || Paternal: 154 
    ## Seed:  137 139 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  242 242 || Leaf:  233 242 || Maternal:  242 || Paternal: 233 
    ## Seed:  165 205 || Leaf:  165 197 || Maternal:  165 || Paternal: 197 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 177 || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 238 || Maternal:  225 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0110 ...
    ## ---------------------------------------------
    ## Seed:  100 112 || Leaf:  100 112 || Maternal:  NA || Paternal: NA 
    ## Seed:  148 154 || Leaf:  148 154 || Maternal:  NA || Paternal: NA 
    ## Seed:  139 139 || Leaf:  137 139 || Maternal:  139 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  185 193 || Leaf:  193 193 || Maternal:  193 || Paternal: 193 
    ## Seed:  171 171 || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0111 ...
    ## ---------------------------------------------
    ## Seed:  104 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  154 154 || Leaf:  148 154 || Maternal:  154 || Paternal: 148 
    ## Seed:  137 141 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 138 || Maternal:  130 || Paternal: 138 
    ## Seed:  231 243 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 199 || Leaf:  199 201 || Maternal:  199 || Paternal: 201 
    ## Seed:  191 191 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  173 177 || Leaf:  173 177 || Maternal:  NA || Paternal: NA 
    ## Seed:  225 238 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0112 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  94 98 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  148 154 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  NA NA || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  NA NA || Leaf:  243 243 || Maternal:  243 || Paternal: 243 
    ## Seed:  NA NA || Leaf:  167 167 || Maternal:  167 || Paternal: 167 
    ## Seed:  NA NA || Leaf:  185 191 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  173 177 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  225 236 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0113 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  154 154 || Leaf:  148 154 || Maternal:  154 || Paternal: 148 
    ## Seed:  137 141 || Leaf:  137 141 || Maternal:  NA || Paternal: NA 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  243 243 || Maternal:  243 || Paternal: 243 
    ## Seed:  167 199 || Leaf:  165 199 || Maternal:  199 || Paternal: 165 
    ## Seed:  191 191 || Leaf:  185 191 || Maternal:  191 || Paternal: 185 
    ## Seed:  177 181 || Leaf:  177 179 || Maternal:  177 || Paternal: 179 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0114 ...
    ## ---------------------------------------------
    ## Seed:  94 94 || Leaf:  94 100 || Maternal:  94 || Paternal: 100 
    ## Seed:  154 156 || Leaf:  154 156 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  137 139 || Maternal:  137 || Paternal: 139 
    ## Seed:  128 128 || Leaf:  128 130 || Maternal:  128 || Paternal: 130 
    ## Seed:  231 238 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  165 199 || Leaf:  165 165 || Maternal:  165 || Paternal: 165 
    ## Seed:  191 191 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  177 181 || Leaf:  171 177 || Maternal:  177 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0115 ...
    ## ---------------------------------------------
    ## Seed:  94 108 || Leaf:  94 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  154 156 || Leaf:  154 156 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  227 231 || Maternal:  231 || Paternal: 227 
    ## Seed:  199 201 || Leaf:  167 199 || Maternal:  199 || Paternal: 167 
    ## Seed:  185 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 179 || Leaf:  177 179 || Maternal:  179 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0116 ...
    ## ---------------------------------------------
    ## Seed:  94 104 || Leaf:  94 94 || Maternal:  94 || Paternal: 94 
    ## Seed:  154 160 || Leaf:  154 156 || Maternal:  154 || Paternal: 156 
    ## Seed:  137 141 || Leaf:  137 141 || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  231 243 || Leaf:  231 243 || Maternal:  NA || Paternal: NA 
    ## Seed:  167 199 || Leaf:  167 167 || Maternal:  167 || Paternal: 167 
    ## Seed:  185 191 || Leaf:  185 191 || Maternal:  NA || Paternal: NA 
    ## Seed:  177 177 || Leaf:  171 177 || Maternal:  177 || Paternal: 171 
    ## Seed:  225 236 || Leaf:  225 236 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0117 ...
    ## ---------------------------------------------
    ## Seed:  100 106 || Leaf:  106 106 || Maternal:  106 || Paternal: 106 
    ## Seed:  148 156 || Leaf:  148 158 || Maternal:  148 || Paternal: 158 
    ## Seed:  137 139 || Leaf:  137 139 || Maternal:  NA || Paternal: NA 
    ## Seed:  128 130 || Leaf:  130 138 || Maternal:  130 || Paternal: 138 
    ## Seed:  231 233 || Leaf:  231 233 || Maternal:  NA || Paternal: NA 
    ## Seed:  199 205 || Leaf:  199 199 || Maternal:  199 || Paternal: 199 
    ## Seed:  191 193 || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  173 177 || Leaf:  171 173 || Maternal:  173 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0118 ...
    ## ---------------------------------------------
    ## Seed:  94 106 || Leaf:  106 108 || Maternal:  106 || Paternal: 108 
    ## Seed:  156 158 || Leaf:  154 158 || Maternal:  158 || Paternal: 154 
    ## Seed:  137 139 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 175 || Leaf:  167 167 || Maternal:  167 || Paternal: 167 
    ## Seed:  185 191 || Leaf:  185 191 || Maternal:  NA || Paternal: NA 
    ## Seed:  173 179 || Leaf:  173 181 || Maternal:  173 || Paternal: 181 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0119 ...
    ## ---------------------------------------------
    ## Seed:  100 106 || Leaf:  106 108 || Maternal:  106 || Paternal: 108 
    ## Seed:  148 156 || Leaf:  148 158 || Maternal:  148 || Paternal: 158 
    ## Seed:  137 139 || Leaf:  137 139 || Maternal:  NA || Paternal: NA 
    ## Seed:  128 130 || Leaf:  128 138 || Maternal:  128 || Paternal: 138 
    ## Seed:  231 233 || Leaf:  227 233 || Maternal:  233 || Paternal: 227 
    ## Seed:  199 205 || Leaf:  199 199 || Maternal:  199 || Paternal: 199 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 177 || Leaf:  171 177 || Maternal:  177 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0121 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  148 154 || Leaf:  154 160 || Maternal:  154 || Paternal: 160 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 205 || Leaf:  199 205 || Maternal:  205 || Paternal: 199 
    ## Seed:  191 193 || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  177 177 || Leaf:  171 177 || Maternal:  177 || Paternal: 171 
    ## Seed:  225 238 || Leaf:  238 238 || Maternal:  238 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0122 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  104 108 || Maternal:  108 || Paternal: 104 
    ## Seed:  154 158 || Leaf:  152 154 || Maternal:  154 || Paternal: 152 
    ## Seed:  139 139 || Leaf:  137 139 || Maternal:  139 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 243 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 199 || Leaf:  167 199 || Maternal:  199 || Paternal: 167 
    ## Seed:  191 191 || Leaf:  185 191 || Maternal:  191 || Paternal: 185 
    ## Seed:  171 181 || Leaf:  171 177 || Maternal:  171 || Paternal: 177 
    ## Seed:  236 240 || Leaf:  225 240 || Maternal:  240 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0124 ...
    ## ---------------------------------------------
    ## Seed:  94 106 || Leaf:  106 108 || Maternal:  106 || Paternal: 108 
    ## Seed:  156 158 || Leaf:  148 156 || Maternal:  156 || Paternal: 148 
    ## Seed:  137 139 || Leaf:  139 139 || Maternal:  139 || Paternal: 139 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 175 || Leaf:  167 197 || Maternal:  167 || Paternal: 197 
    ## Seed:  185 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 179 || Leaf:  177 179 || Maternal:  179 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0125 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  94 108 || Maternal:  108 || Paternal: 94 
    ## Seed:  154 158 || Leaf:  154 160 || Maternal:  154 || Paternal: 160 
    ## Seed:  139 139 || Leaf:  137 139 || Maternal:  139 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 243 || Leaf:  231 243 || Maternal:  NA || Paternal: NA 
    ## Seed:  199 199 || Leaf:  167 199 || Maternal:  199 || Paternal: 167 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 181 || Leaf:  171 177 || Maternal:  171 || Paternal: 177 
    ## Seed:  236 240 || Leaf:  225 236 || Maternal:  236 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0126 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  98 108 || Maternal:  108 || Paternal: 98 
    ## Seed:  148 154 || Leaf:  148 148 || Maternal:  148 || Paternal: 148 
    ## Seed:  137 137 || Leaf:  137 139 || Maternal:  137 || Paternal: 139 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 205 || Leaf:  167 199 || Maternal:  167 || Paternal: 199 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  177 177 || Leaf:  177 179 || Maternal:  177 || Paternal: 179 
    ## Seed:  225 238 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0127 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  106 108 || Maternal:  108 || Paternal: 106 
    ## Seed:  148 154 || Leaf:  148 156 || Maternal:  148 || Paternal: 156 
    ## Seed:  137 137 || Leaf:  137 139 || Maternal:  137 || Paternal: 139 
    ## Seed:  130 138 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 205 || Leaf:  167 175 || Maternal:  167 || Paternal: 175 
    ## Seed:  191 193 || Leaf:  185 193 || Maternal:  193 || Paternal: 185 
    ## Seed:  177 177 || Leaf:  177 179 || Maternal:  177 || Paternal: 179 
    ## Seed:  225 238 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0128 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  154 154 || Leaf:  154 158 || Maternal:  154 || Paternal: 158 
    ## Seed:  135 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  231 233 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  165 199 || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 177 || Leaf:  173 177 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  225 236 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0129 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  104 108 || Maternal:  108 || Paternal: 104 
    ## Seed:  154 158 || Leaf:  152 154 || Maternal:  154 || Paternal: 152 
    ## Seed:  139 139 || Leaf:  137 139 || Maternal:  139 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 243 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 199 || Leaf:  167 199 || Maternal:  199 || Paternal: 167 
    ## Seed:  191 191 || Leaf:  185 191 || Maternal:  191 || Paternal: 185 
    ## Seed:  171 181 || Leaf:  171 177 || Maternal:  171 || Paternal: 177 
    ## Seed:  236 240 || Leaf:  225 236 || Maternal:  236 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0130 ...
    ## ---------------------------------------------
    ## Seed:  94 106 || Leaf:  106 108 || Maternal:  106 || Paternal: 108 
    ## Seed:  154 156 || Leaf:  148 156 || Maternal:  156 || Paternal: 148 
    ## Seed:  137 137 || Leaf:  137 141 || Maternal:  137 || Paternal: 141 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  175 197 || Leaf:  175 195 || Maternal:  175 || Paternal: 195 
    ## Seed:  191 193 || Leaf:  193 193 || Maternal:  193 || Paternal: 193 
    ## Seed:  177 177 || Leaf:  177 181 || Maternal:  177 || Paternal: 181 
    ## Seed:  225 225 || Leaf:  225 242 || Maternal:  225 || Paternal: 242 
    ## 
    ##  Working on indv: 07-0131 ...
    ## ---------------------------------------------
    ## Seed:  104 108 || Leaf:  104 104 || Maternal:  104 || Paternal: 104 
    ## Seed:  154 154 || Leaf:  152 154 || Maternal:  154 || Paternal: 152 
    ## Seed:  137 141 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 138 || Maternal:  130 || Paternal: 138 
    ## Seed:  231 243 || Leaf:  231 243 || Maternal:  NA || Paternal: NA 
    ## Seed:  167 199 || Leaf:  167 195 || Maternal:  167 || Paternal: 195 
    ## Seed:  191 191 || Leaf:  185 191 || Maternal:  191 || Paternal: 185 
    ## Seed:  173 177 || Leaf:  173 177 || Maternal:  NA || Paternal: NA 
    ## Seed:  225 238 || Leaf:  225 240 || Maternal:  225 || Paternal: 240 
    ## 
    ##  Working on indv: 07-0132 ...
    ## ---------------------------------------------
    ## Seed:  104 108 || Leaf:  104 106 || Maternal:  104 || Paternal: 106 
    ## Seed:  154 154 || Leaf:  154 154 || Maternal:  154 || Paternal: 154 
    ## Seed:  137 141 || Leaf:  129 141 || Maternal:  141 || Paternal: 129 
    ## Seed:  130 130 || Leaf:  130 138 || Maternal:  130 || Paternal: 138 
    ## Seed:  231 243 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  167 199 || Leaf:  167 201 || Maternal:  167 || Paternal: 201 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 177 || Leaf:  165 177 || Maternal:  177 || Paternal: 165 
    ## Seed:  225 238 || Leaf:  225 238 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0133 ...
    ## ---------------------------------------------
    ## Seed:  104 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  154 154 || Leaf:  154 154 || Maternal:  154 || Paternal: 154 
    ## Seed:  137 141 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  231 243 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  167 199 || Leaf:  197 199 || Maternal:  199 || Paternal: 197 
    ## Seed:  191 191 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  173 177 || Leaf:  173 177 || Maternal:  NA || Paternal: NA 
    ## Seed:  225 238 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0134 ...
    ## ---------------------------------------------
    ## Seed:  100 112 || Leaf:  108 112 || Maternal:  112 || Paternal: 108 
    ## Seed:  148 154 || Leaf:  148 152 || Maternal:  148 || Paternal: 152 
    ## Seed:  137 139 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  165 167 || Leaf:  167 199 || Maternal:  167 || Paternal: 199 
    ## Seed:  185 193 || Leaf:  191 193 || Maternal:  193 || Paternal: 191 
    ## Seed:  171 177 || Leaf:  171 177 || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 236 || Maternal:  225 || Paternal: 236 
    ## 
    ##  Working on indv: 07-0138 ...
    ## ---------------------------------------------
    ## Seed:  104 108 || Leaf:  94 108 || Maternal:  108 || Paternal: 94 
    ## Seed:  154 156 || Leaf:  156 158 || Maternal:  156 || Paternal: 158 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  229 231 || Leaf:  229 229 || Maternal:  229 || Paternal: 229 
    ## Seed:  165 201 || Leaf:  201 205 || Maternal:  201 || Paternal: 205 
    ## Seed:  185 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  165 173 || Leaf:  165 177 || Maternal:  165 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0139 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  NA NA || Leaf:  150 150 || Maternal:  150 || Paternal: 150 
    ## Seed:  NA NA || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  NA NA || Leaf:  128 128 || Maternal:  128 || Paternal: 128 
    ## Seed:  NA NA || Leaf:  231 242 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  193 199 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  NA NA || Leaf:  181 181 || Maternal:  181 || Paternal: 181 
    ## Seed:  NA NA || Leaf:  225 236 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0140 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  106 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  154 156 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  NA NA || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  205 205 || Maternal:  205 || Paternal: 205 
    ## Seed:  NA NA || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  NA NA || Leaf:  171 181 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0141 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 110 || Maternal:  108 || Paternal: 110 
    ## Seed:  154 156 || Leaf:  154 158 || Maternal:  154 || Paternal: 158 
    ## Seed:  137 139 || Leaf:  137 141 || Maternal:  137 || Paternal: 141 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 243 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 201 || Leaf:  199 199 || Maternal:  199 || Paternal: 199 
    ## Seed:  191 195 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  173 173 || Leaf:  173 173 || Maternal:  173 || Paternal: 173 
    ## Seed:  225 225 || Leaf:  225 238 || Maternal:  225 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0143 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  94 106 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  148 158 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  137 141 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  NA NA || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  197 205 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  NA NA || Leaf:  173 177 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0145 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  98 108 || Maternal:  108 || Paternal: 98 
    ## Seed:  154 156 || Leaf:  154 154 || Maternal:  154 || Paternal: 154 
    ## Seed:  137 139 || Leaf:  137 139 || Maternal:  NA || Paternal: NA 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 243 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 201 || Leaf:  199 201 || Maternal:  NA || Paternal: NA 
    ## Seed:  191 195 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 173 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0147 ...
    ## ---------------------------------------------
    ## Seed:  94 106 || Leaf:  94 108 || Maternal:  94 || Paternal: 108 
    ## Seed:  156 158 || Leaf:  150 156 || Maternal:  156 || Paternal: 150 
    ## Seed:  137 139 || Leaf:  139 139 || Maternal:  139 || Paternal: 139 
    ## Seed:  128 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 175 || Leaf:  167 195 || Maternal:  167 || Paternal: 195 
    ## Seed:  185 191 || Leaf:  185 193 || Maternal:  185 || Paternal: 193 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 236 || Maternal:  225 || Paternal: 236 
    ## 
    ##  Working on indv: 07-0151 ...
    ## ---------------------------------------------
    ## Seed:  104 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 139 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  128 138 || Leaf:  128 130 || Maternal:  128 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  165 199 || Leaf:  165 165 || Maternal:  165 || Paternal: 165 
    ## Seed:  185 191 || Leaf:  185 188 || Maternal:  185 || Paternal: 188 
    ## Seed:  165 173 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0155 ...
    ## ---------------------------------------------
    ## Seed:  98 108 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  148 158 || Leaf:  148 158 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 139 || Leaf:  137 139 || Maternal:  NA || Paternal: NA 
    ## Seed:  130 138 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 243 || Maternal:  231 || Paternal: 243 
    ## Seed:  199 205 || Leaf:  199 199 || Maternal:  199 || Paternal: 199 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 240 || Maternal:  225 || Paternal: 240 
    ## 
    ##  Working on indv: 07-0159 ...
    ## ---------------------------------------------
    ## Seed:  100 108 || Leaf:  100 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  152 160 || Leaf:  148 152 || Maternal:  152 || Paternal: 148 
    ## Seed:  139 141 || Leaf:  139 141 || Maternal:  NA || Paternal: NA 
    ## Seed:  130 138 || Leaf:  138 138 || Maternal:  138 || Paternal: 138 
    ## Seed:  231 245 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  165 165 || Leaf:  165 165 || Maternal:  165 || Paternal: 165 
    ## Seed:  185 193 || Leaf:  185 191 || Maternal:  185 || Paternal: 191 
    ## Seed:  173 175 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 238 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0160 ...
    ## ---------------------------------------------
    ## Seed:  98 108 || Leaf:  98 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  148 154 || Leaf:  150 154 || Maternal:  154 || Paternal: 150 
    ## Seed:  137 139 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 233 || Maternal:  231 || Paternal: 233 
    ## Seed:  197 205 || Leaf:  195 205 || Maternal:  205 || Paternal: 195 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  175 177 || Leaf:  173 177 || Maternal:  177 || Paternal: 173 
    ## Seed:  225 238 || Leaf:  238 238 || Maternal:  238 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0166 ...
    ## ---------------------------------------------
    ## Seed:  106 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  148 156 || Leaf:  148 156 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 139 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 243 || Leaf:  227 243 || Maternal:  243 || Paternal: 227 
    ## Seed:  199 203 || Leaf:  167 203 || Maternal:  203 || Paternal: 167 
    ## Seed:  185 195 || Leaf:  185 191 || Maternal:  185 || Paternal: 191 
    ## Seed:  173 177 || Leaf:  171 173 || Maternal:  173 || Paternal: 171 
    ## Seed:  225 238 || Leaf:  225 238 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0167 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  100 108 || Maternal:  108 || Paternal: 100 
    ## Seed:  154 158 || Leaf:  148 154 || Maternal:  154 || Paternal: 148 
    ## Seed:  139 139 || Leaf:  139 139 || Maternal:  139 || Paternal: 139 
    ## Seed:  128 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 243 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 199 || Leaf:  199 201 || Maternal:  199 || Paternal: 201 
    ## Seed:  191 191 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  171 181 || Leaf:  171 179 || Maternal:  171 || Paternal: 179 
    ## Seed:  236 240 || Leaf:  225 236 || Maternal:  236 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0168 ...
    ## ---------------------------------------------
    ## Seed:  106 108 || Leaf:  94 106 || Maternal:  106 || Paternal: 94 
    ## Seed:  148 156 || Leaf:  148 156 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 139 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  128 138 || Maternal:  128 || Paternal: 138 
    ## Seed:  231 243 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 203 || Leaf:  175 203 || Maternal:  203 || Paternal: 175 
    ## Seed:  185 195 || Leaf:  185 195 || Maternal:  NA || Paternal: NA 
    ## Seed:  173 177 || Leaf:  171 173 || Maternal:  173 || Paternal: 171 
    ## Seed:  225 238 || Leaf:  225 238 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0169 ...
    ## ---------------------------------------------
    ## Seed:  94 94 || Leaf:  94 108 || Maternal:  94 || Paternal: 108 
    ## Seed:  152 158 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  137 141 || Maternal:  137 || Paternal: 141 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  227 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  175 205 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  177 177 || Leaf:  171 177 || Maternal:  177 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0172 ...
    ## ---------------------------------------------
    ## Seed:  100 112 || Leaf:  100 108 || Maternal:  100 || Paternal: 108 
    ## Seed:  148 154 || Leaf:  148 154 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 139 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  165 167 || Leaf:  167 199 || Maternal:  167 || Paternal: 199 
    ## Seed:  185 193 || Leaf:  185 191 || Maternal:  185 || Paternal: 191 
    ## Seed:  171 177 || Leaf:  171 181 || Maternal:  171 || Paternal: 181 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0173 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  154 160 || Leaf:  148 154 || Maternal:  154 || Paternal: 148 
    ## Seed:  137 141 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  130 138 || Maternal:  130 || Paternal: 138 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  165 171 || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  225 242 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0174 ...
    ## ---------------------------------------------
    ## Seed:  106 108 || Leaf:  100 106 || Maternal:  106 || Paternal: 100 
    ## Seed:  154 160 || Leaf:  152 154 || Maternal:  154 || Paternal: 152 
    ## Seed:  137 141 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 138 || Maternal:  130 || Paternal: 138 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  197 205 || Leaf:  195 197 || Maternal:  197 || Paternal: 195 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  165 171 || Leaf:  165 171 || Maternal:  NA || Paternal: NA 
    ## Seed:  225 242 || Leaf:  225 240 || Maternal:  225 || Paternal: 240 
    ## 
    ##  Working on indv: 07-0177 ...
    ## ---------------------------------------------
    ## Seed:  94 94 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  152 158 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  128 130 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  227 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  175 205 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  177 177 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0186 ...
    ## ---------------------------------------------
    ## Seed:  108 110 || Leaf:  110 110 || Maternal:  110 || Paternal: 110 
    ## Seed:  148 154 || Leaf:  148 148 || Maternal:  148 || Paternal: 148 
    ## Seed:  137 141 || Leaf:  137 141 || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 205 || Leaf:  167 167 || Maternal:  167 || Paternal: 167 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  177 179 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0188 ...
    ## ---------------------------------------------
    ## Seed:  108 110 || Leaf:  108 110 || Maternal:  NA || Paternal: NA 
    ## Seed:  148 154 || Leaf:  150 154 || Maternal:  154 || Paternal: 150 
    ## Seed:  137 141 || Leaf:  137 141 || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  167 205 || Leaf:  167 197 || Maternal:  167 || Paternal: 197 
    ## Seed:  191 193 || Leaf:  193 193 || Maternal:  193 || Paternal: 193 
    ## Seed:  177 179 || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0189 ...
    ## ---------------------------------------------
    ## Seed:  100 106 || Leaf:  100 104 || Maternal:  100 || Paternal: 104 
    ## Seed:  148 156 || Leaf:  154 156 || Maternal:  156 || Paternal: 154 
    ## Seed:  137 139 || Leaf:  139 141 || Maternal:  139 || Paternal: 141 
    ## Seed:  128 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 233 || Leaf:  233 233 || Maternal:  233 || Paternal: 233 
    ## Seed:  199 205 || Leaf:  199 199 || Maternal:  199 || Paternal: 199 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 177 || Leaf:  173 173 || Maternal:  173 || Paternal: 173 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0190 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  148 148 || Leaf:  148 160 || Maternal:  148 || Paternal: 160 
    ## Seed:  137 141 || Leaf:  137 139 || Maternal:  137 || Paternal: 139 
    ## Seed:  130 138 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  231 233 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  177 181 || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0191 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  148 148 || Leaf:  148 154 || Maternal:  148 || Paternal: 154 
    ## Seed:  137 141 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  130 138 || Leaf:  124 130 || Maternal:  130 || Paternal: 124 
    ## Seed:  231 243 || Leaf:  233 243 || Maternal:  243 || Paternal: 233 
    ## Seed:  165 199 || Leaf:  193 199 || Maternal:  199 || Paternal: 193 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 171 || Leaf:  171 181 || Maternal:  171 || Paternal: 181 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0192 ...
    ## ---------------------------------------------
    ## Seed:  100 106 || Leaf:  106 108 || Maternal:  106 || Paternal: 108 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 139 || Leaf:  139 139 || Maternal:  139 || Paternal: 139 
    ## Seed:  128 130 || Leaf:  128 128 || Maternal:  128 || Paternal: 128 
    ## Seed:  231 233 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 205 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 193 || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  173 177 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0193 ...
    ## ---------------------------------------------
    ## Seed:  108 110 || Leaf:  94 108 || Maternal:  108 || Paternal: 94 
    ## Seed:  148 156 || Leaf:  154 156 || Maternal:  156 || Paternal: 154 
    ## Seed:  137 141 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 238 || Maternal:  231 || Paternal: 238 
    ## Seed:  165 167 || Leaf:  167 167 || Maternal:  167 || Paternal: 167 
    ## Seed:  193 193 || Leaf:  191 193 || Maternal:  193 || Paternal: 191 
    ## Seed:  177 179 || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 238 || Maternal:  225 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0195 ...
    ## ---------------------------------------------
    ## Seed:  108 110 || Leaf:  108 110 || Maternal:  NA || Paternal: NA 
    ## Seed:  148 154 || Leaf:  148 154 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 141 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 238 || Maternal:  231 || Paternal: 238 
    ## Seed:  167 205 || Leaf:  175 205 || Maternal:  205 || Paternal: 175 
    ## Seed:  191 193 || Leaf:  185 193 || Maternal:  193 || Paternal: 185 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  221 225 || Maternal:  225 || Paternal: 221 
    ## 
    ##  Working on indv: 07-0218 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  148 160 || Leaf:  148 156 || Maternal:  148 || Paternal: 156 
    ## Seed:  137 139 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  165 197 || Leaf:  165 199 || Maternal:  165 || Paternal: 199 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  171 171 || Leaf:  171 177 || Maternal:  171 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0219 ...
    ## ---------------------------------------------
    ## Seed:  98 106 || Leaf:  98 106 || Maternal:  NA || Paternal: NA 
    ## Seed:  148 158 || Leaf:  154 158 || Maternal:  158 || Paternal: 154 
    ## Seed:  137 139 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  130 138 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  165 201 || Leaf:  165 165 || Maternal:  165 || Paternal: 165 
    ## Seed:  191 191 || Leaf:  185 191 || Maternal:  191 || Paternal: 185 
    ## Seed:  177 179 || Leaf:  179 179 || Maternal:  179 || Paternal: 179 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0223 ...
    ## ---------------------------------------------
    ## Seed:  104 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 139 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  128 128 || Maternal:  128 || Paternal: 128 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 203 || Leaf:  199 199 || Maternal:  199 || Paternal: 199 
    ## Seed:  193 193 || Leaf:  193 193 || Maternal:  193 || Paternal: 193 
    ## Seed:  177 179 || Leaf:  171 179 || Maternal:  179 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0227 ...
    ## ---------------------------------------------
    ## Seed:  94 108 || Leaf:  94 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  154 156 || Leaf:  154 160 || Maternal:  154 || Paternal: 160 
    ## Seed:  137 137 || Leaf:  137 141 || Maternal:  137 || Paternal: 141 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 235 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  193 199 || Leaf:  165 199 || Maternal:  199 || Paternal: 165 
    ## Seed:  193 193 || Leaf:  193 193 || Maternal:  193 || Paternal: 193 
    ## Seed:  165 177 || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  225 238 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0230 ...
    ## ---------------------------------------------
    ## Seed:  104 108 || Leaf:  94 108 || Maternal:  108 || Paternal: 94 
    ## Seed:  148 148 || Leaf:  148 158 || Maternal:  148 || Paternal: 158 
    ## Seed:  137 139 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 203 || Leaf:  199 205 || Maternal:  199 || Paternal: 205 
    ## Seed:  193 193 || Leaf:  193 193 || Maternal:  193 || Paternal: 193 
    ## Seed:  177 179 || Leaf:  165 179 || Maternal:  179 || Paternal: 165 
    ## Seed:  225 225 || Leaf:  225 238 || Maternal:  225 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0233 ...
    ## ---------------------------------------------
    ## Seed:  104 108 || Leaf:  94 104 || Maternal:  104 || Paternal: 94 
    ## Seed:  148 148 || Leaf:  148 152 || Maternal:  148 || Paternal: 152 
    ## Seed:  137 139 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  203 203 || Leaf:  199 203 || Maternal:  203 || Paternal: 199 
    ## Seed:  193 193 || Leaf:  193 193 || Maternal:  193 || Paternal: 193 
    ## Seed:  177 179 || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0235 ...
    ## ---------------------------------------------
    ## Seed:  104 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  148 148 || Leaf:  148 148 || Maternal:  148 || Paternal: 148 
    ## Seed:  137 139 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  128 130 || Leaf:  128 128 || Maternal:  128 || Paternal: 128 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 203 || Leaf:  199 199 || Maternal:  199 || Paternal: 199 
    ## Seed:  193 193 || Leaf:  193 193 || Maternal:  193 || Paternal: 193 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0240 ...
    ## ---------------------------------------------
    ## Seed:  104 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  148 148 || Leaf:  148 148 || Maternal:  148 || Paternal: 148 
    ## Seed:  137 139 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 203 || Leaf:  167 199 || Maternal:  199 || Paternal: 167 
    ## Seed:  193 193 || Leaf:  185 193 || Maternal:  193 || Paternal: 185 
    ## Seed:  177 179 || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  225 236 || Leaf:  225 236 || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0242 ...
    ## ---------------------------------------------
    ## Seed:  98 104 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  154 160 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  130 138 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  231 243 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  167 205 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  173 179 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0244 ...
    ## ---------------------------------------------
    ## Seed:  104 108 || Leaf:  104 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  154 154 || Leaf:  154 154 || Maternal:  154 || Paternal: 154 
    ## Seed:  137 141 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  231 243 || Leaf:  243 243 || Maternal:  243 || Paternal: 243 
    ## Seed:  167 199 || Leaf:  199 199 || Maternal:  199 || Paternal: 199 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 177 || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  225 238 || Leaf:  238 238 || Maternal:  238 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0245 ...
    ## ---------------------------------------------
    ## Seed:  98 104 || Leaf:  98 110 || Maternal:  98 || Paternal: 110 
    ## Seed:  154 160 || Leaf:  148 154 || Maternal:  154 || Paternal: 148 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 243 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  167 205 || Leaf:  167 167 || Maternal:  167 || Paternal: 167 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 179 || Leaf:  173 177 || Maternal:  173 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0247 ...
    ## ---------------------------------------------
    ## Seed:  94 94 || Leaf:  94 108 || Maternal:  94 || Paternal: 108 
    ## Seed:  152 158 || Leaf:  152 158 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  227 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  175 205 || Leaf:  175 199 || Maternal:  175 || Paternal: 199 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  177 177 || Leaf:  171 177 || Maternal:  177 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0250 ...
    ## ---------------------------------------------
    ## Seed:  94 110 || Leaf:  108 110 || Maternal:  110 || Paternal: 108 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 141 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  193 193 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  173 181 || Leaf:  181 181 || Maternal:  181 || Paternal: 181 
    ## Seed:  225 238 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0251 ...
    ## ---------------------------------------------
    ## Seed:  98 104 || Leaf:  98 106 || Maternal:  98 || Paternal: 106 
    ## Seed:  154 160 || Leaf:  154 154 || Maternal:  154 || Paternal: 154 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 243 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 205 || Leaf:  167 201 || Maternal:  167 || Paternal: 201 
    ## Seed:  191 191 || Leaf:  185 191 || Maternal:  191 || Paternal: 185 
    ## Seed:  173 179 || Leaf:  173 181 || Maternal:  173 || Paternal: 181 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0253 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  154 160 || Leaf:  148 154 || Maternal:  154 || Paternal: 148 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  167 205 || Leaf:  167 199 || Maternal:  167 || Paternal: 199 
    ## Seed:  191 191 || Leaf:  185 191 || Maternal:  191 || Paternal: 185 
    ## Seed:  173 179 || Leaf:  173 173 || Maternal:  173 || Paternal: 173 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0255 ...
    ## ---------------------------------------------
    ## Seed:  98 104 || Leaf:  104 108 || Maternal:  104 || Paternal: 108 
    ## Seed:  154 160 || Leaf:  156 160 || Maternal:  160 || Paternal: 156 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 243 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  167 205 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 179 || Leaf:  171 173 || Maternal:  173 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0257 ...
    ## ---------------------------------------------
    ## Seed:  94 110 || Leaf:  108 110 || Maternal:  110 || Paternal: 108 
    ## Seed:  158 158 || Leaf:  152 158 || Maternal:  158 || Paternal: 152 
    ## Seed:  137 141 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 181 || Leaf:  173 173 || Maternal:  173 || Paternal: 173 
    ## Seed:  225 238 || Leaf:  225 236 || Maternal:  225 || Paternal: 236 
    ## 
    ##  Working on indv: 07-0260 ...
    ## ---------------------------------------------
    ## Seed:  94 108 || Leaf:  94 110 || Maternal:  94 || Paternal: 110 
    ## Seed:  148 158 || Leaf:  150 158 || Maternal:  158 || Paternal: 150 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 240 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  167 199 || Leaf:  175 199 || Maternal:  199 || Paternal: 175 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  171 177 || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  225 238 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0261 ...
    ## ---------------------------------------------
    ## Seed:  98 104 || Leaf:  98 98 || Maternal:  98 || Paternal: 98 
    ## Seed:  154 160 || Leaf:  160 160 || Maternal:  160 || Paternal: 160 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  128 138 || Maternal:  138 || Paternal: 128 
    ## Seed:  231 243 || Leaf:  231 243 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  173 179 || Leaf:  179 181 || Maternal:  179 || Paternal: 181 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0263 ...
    ## ---------------------------------------------
    ## Seed:  98 104 || Leaf:  104 108 || Maternal:  104 || Paternal: 108 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 243 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  167 205 || Leaf:  167 167 || Maternal:  167 || Paternal: 167 
    ## Seed:  191 191 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  173 179 || Leaf:  173 177 || Maternal:  173 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0267 ...
    ## ---------------------------------------------
    ## Seed:  106 110 || Leaf:  104 106 || Maternal:  106 || Paternal: 104 
    ## Seed:  148 152 || Leaf:  148 156 || Maternal:  148 || Paternal: 156 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 201 || Leaf:  201 201 || Maternal:  201 || Paternal: 201 
    ## Seed:  191 191 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0268 ...
    ## ---------------------------------------------
    ## Seed:  108 112 || Leaf:  106 112 || Maternal:  112 || Paternal: 106 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  129 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 171 || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  225 238 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0269 ...
    ## ---------------------------------------------
    ## Seed:  100 108 || Leaf:  100 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  137 141 || Maternal:  137 || Paternal: 141 
    ## Seed:  130 130 || Leaf:  130 138 || Maternal:  130 || Paternal: 138 
    ## Seed:  231 238 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  165 197 || Leaf:  165 199 || Maternal:  165 || Paternal: 199 
    ## Seed:  185 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  177 179 || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  238 240 || Leaf:  225 238 || Maternal:  238 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0270 ...
    ## ---------------------------------------------
    ## Seed:  106 110 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  148 152 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  130 138 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 201 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  171 173 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0271 ...
    ## ---------------------------------------------
    ## Seed:  98 104 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  154 160 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  130 138 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  231 243 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  167 205 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  173 179 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0272 ...
    ## ---------------------------------------------
    ## Seed:  98 108 || Leaf:  94 98 || Maternal:  98 || Paternal: 94 
    ## Seed:  150 158 || Leaf:  158 158 || Maternal:  158 || Paternal: 158 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  231 235 || Maternal:  231 || Paternal: 235 
    ## Seed:  167 197 || Leaf:  167 175 || Maternal:  167 || Paternal: 175 
    ## Seed:  191 193 || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0273 ...
    ## ---------------------------------------------
    ## Seed:  94 98 || Leaf:  98 108 || Maternal:  98 || Paternal: 108 
    ## Seed:  148 154 || Leaf:  148 154 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 141 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 138 || Maternal:  130 || Paternal: 138 
    ## Seed:  231 235 || Leaf:  231 235 || Maternal:  NA || Paternal: NA 
    ## Seed:  175 197 || Leaf:  175 205 || Maternal:  175 || Paternal: 205 
    ## Seed:  191 193 || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  165 165 || Leaf:  165 177 || Maternal:  165 || Paternal: 177 
    ## Seed:  225 238 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0274 ...
    ## ---------------------------------------------
    ## Seed:  98 100 || Leaf:  98 108 || Maternal:  98 || Paternal: 108 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 141 || Leaf:  141 141 || Maternal:  141 || Paternal: 141 
    ## Seed:  138 138 || Leaf:  130 138 || Maternal:  138 || Paternal: 130 
    ## Seed:  227 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  197 199 || Leaf:  167 197 || Maternal:  197 || Paternal: 167 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 177 || Leaf:  171 173 || Maternal:  173 || Paternal: 171 
    ## Seed:  NA NA || Leaf:  238 238 || Maternal:  238 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0275 ...
    ## ---------------------------------------------
    ## Seed:  106 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  154 160 || Leaf:  160 160 || Maternal:  160 || Paternal: 160 
    ## Seed:  137 141 || Leaf:  141 143 || Maternal:  141 || Paternal: 143 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  197 205 || Leaf:  197 197 || Maternal:  197 || Paternal: 197 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  165 171 || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  225 242 || Leaf:  242 242 || Maternal:  242 || Paternal: 242 
    ## 
    ##  Working on indv: 07-0276 ...
    ## ---------------------------------------------
    ## Seed:  108 112 || Leaf:  100 112 || Maternal:  112 || Paternal: 100 
    ## Seed:  154 156 || Leaf:  152 154 || Maternal:  154 || Paternal: 152 
    ## Seed:  129 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 205 || Leaf:  167 197 || Maternal:  167 || Paternal: 197 
    ## Seed:  191 191 || Leaf:  185 191 || Maternal:  191 || Paternal: 185 
    ## Seed:  171 171 || Leaf:  171 177 || Maternal:  171 || Paternal: 177 
    ## Seed:  225 238 || Leaf:  238 240 || Maternal:  238 || Paternal: 240 
    ## 
    ##  Working on indv: 07-0277 ...
    ## ---------------------------------------------
    ## Seed:  106 110 || Leaf:  94 110 || Maternal:  110 || Paternal: 94 
    ## Seed:  148 152 || Leaf:  152 154 || Maternal:  152 || Paternal: 154 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 201 || Leaf:  165 199 || Maternal:  199 || Paternal: 165 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  171 173 || Leaf:  173 173 || Maternal:  173 || Paternal: 173 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0278 ...
    ## ---------------------------------------------
    ## Seed:  106 110 || Leaf:  94 110 || Maternal:  110 || Paternal: 94 
    ## Seed:  148 152 || Leaf:  148 158 || Maternal:  148 || Paternal: 158 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 201 || Leaf:  167 201 || Maternal:  201 || Paternal: 167 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 171 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0279 ...
    ## ---------------------------------------------
    ## Seed:  106 110 || Leaf:  94 106 || Maternal:  106 || Paternal: 94 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 201 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 173 || Leaf:  173 173 || Maternal:  173 || Paternal: 173 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0280 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  148 152 || Leaf:  148 160 || Maternal:  148 || Paternal: 160 
    ## Seed:  137 137 || Leaf:  137 139 || Maternal:  137 || Paternal: 139 
    ## Seed:  130 138 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 201 || Leaf:  165 199 || Maternal:  199 || Paternal: 165 
    ## Seed:  191 191 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  171 173 || Leaf:  171 173 || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0281 ...
    ## ---------------------------------------------
    ## Seed:  106 110 || Leaf:  106 108 || Maternal:  106 || Paternal: 108 
    ## Seed:  148 152 || Leaf:  152 154 || Maternal:  152 || Paternal: 154 
    ## Seed:  137 137 || Leaf:  137 141 || Maternal:  137 || Paternal: 141 
    ## Seed:  130 138 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 201 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 173 || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0283 ...
    ## ---------------------------------------------
    ## Seed:  94 108 || Leaf:  94 100 || Maternal:  94 || Paternal: 100 
    ## Seed:  150 158 || Leaf:  148 150 || Maternal:  150 || Paternal: 148 
    ## Seed:  137 141 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  128 138 || Maternal:  138 || Paternal: 128 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 199 || Leaf:  199 199 || Maternal:  199 || Paternal: 199 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 177 || Leaf:  173 179 || Maternal:  173 || Paternal: 179 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0284 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  148 152 || Leaf:  148 148 || Maternal:  148 || Paternal: 148 
    ## Seed:  137 137 || Leaf:  137 139 || Maternal:  137 || Paternal: 139 
    ## Seed:  130 138 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 173 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0285 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  104 108 || Maternal:  108 || Paternal: 104 
    ## Seed:  148 154 || Leaf:  148 156 || Maternal:  148 || Paternal: 156 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  130 138 || Maternal:  130 || Paternal: 138 
    ## Seed:  227 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 199 || Leaf:  199 199 || Maternal:  199 || Paternal: 199 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 181 || Leaf:  173 181 || Maternal:  181 || Paternal: 173 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0287 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  106 108 || Maternal:  108 || Paternal: 106 
    ## Seed:  148 154 || Leaf:  154 154 || Maternal:  154 || Paternal: 154 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  128 128 || Maternal:  128 || Paternal: 128 
    ## Seed:  227 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 201 || Leaf:  201 205 || Maternal:  201 || Paternal: 205 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 181 || Leaf:  177 181 || Maternal:  181 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0289 ...
    ## ---------------------------------------------
    ## Seed:  108 114 || Leaf:  106 108 || Maternal:  108 || Paternal: 106 
    ## Seed:  148 150 || Leaf:  148 148 || Maternal:  148 || Paternal: 148 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  185 191 || Maternal:  191 || Paternal: 185 
    ## Seed:  171 181 || Leaf:  177 181 || Maternal:  181 || Paternal: 177 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0292 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  148 154 || Leaf:  154 154 || Maternal:  154 || Paternal: 154 
    ## Seed:  137 137 || Leaf:  137 139 || Maternal:  137 || Paternal: 139 
    ## Seed:  128 130 || Leaf:  128 128 || Maternal:  128 || Paternal: 128 
    ## Seed:  227 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 201 || Leaf:  199 201 || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 181 || Leaf:  171 177 || Maternal:  171 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0293 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  106 108 || Maternal:  108 || Paternal: 106 
    ## Seed:  148 154 || Leaf:  154 154 || Maternal:  154 || Paternal: 154 
    ## Seed:  137 137 || Leaf:  137 139 || Maternal:  137 || Paternal: 139 
    ## Seed:  128 130 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  227 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 201 || Leaf:  199 201 || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0296 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  154 154 || Leaf:  154 158 || Maternal:  154 || Paternal: 158 
    ## Seed:  137 137 || Leaf:  137 139 || Maternal:  137 || Paternal: 139 
    ## Seed:  128 130 || Leaf:  128 128 || Maternal:  128 || Paternal: 128 
    ## Seed:  227 227 || Leaf:  227 231 || Maternal:  227 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  163 163 || Maternal:  163 || Paternal: 163 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 181 || Leaf:  171 177 || Maternal:  171 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0297 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  98 108 || Maternal:  108 || Paternal: 98 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  128 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  227 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 201 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 181 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0299 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 110 || Maternal:  108 || Paternal: 110 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  227 231 || Leaf:  227 231 || Maternal:  NA || Paternal: NA 
    ## Seed:  199 201 || Leaf:  165 201 || Maternal:  201 || Paternal: 165 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 181 || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0301 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  104 108 || Maternal:  108 || Paternal: 104 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  130 138 || Maternal:  130 || Paternal: 138 
    ## Seed:  227 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 201 || Leaf:  199 199 || Maternal:  199 || Paternal: 199 
    ## Seed:  191 191 || Leaf:  185 191 || Maternal:  191 || Paternal: 185 
    ## Seed:  171 181 || Leaf:  171 181 || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0303 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  98 108 || Maternal:  108 || Paternal: 98 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  227 231 || Leaf:  227 231 || Maternal:  NA || Paternal: NA 
    ## Seed:  199 201 || Leaf:  201 201 || Maternal:  201 || Paternal: 201 
    ## Seed:  191 191 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  171 181 || Leaf:  171 171 || Maternal:  171 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0304 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  148 154 || Leaf:  148 154 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  137 141 || Maternal:  137 || Paternal: 141 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  227 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 201 || Leaf:  201 205 || Maternal:  201 || Paternal: 205 
    ## Seed:  191 191 || Leaf:  191 195 || Maternal:  191 || Paternal: 195 
    ## Seed:  171 181 || Leaf:  177 181 || Maternal:  181 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0307 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 110 || Maternal:  108 || Paternal: 110 
    ## Seed:  148 154 || Leaf:  148 148 || Maternal:  148 || Paternal: 148 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  227 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 181 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0310 ...
    ## ---------------------------------------------
    ## Seed:  100 108 || Leaf:  100 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  148 156 || Leaf:  148 154 || Maternal:  148 || Paternal: 154 
    ## Seed:  137 139 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  231 243 || Leaf:  231 235 || Maternal:  231 || Paternal: 235 
    ## Seed:  199 201 || Leaf:  165 199 || Maternal:  199 || Paternal: 165 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  179 179 || Leaf:  173 179 || Maternal:  179 || Paternal: 173 
    ## Seed:  225 225 || Leaf:  225 238 || Maternal:  225 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0312 ...
    ## ---------------------------------------------
    ## Seed:  106 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  148 156 || Leaf:  156 156 || Maternal:  156 || Paternal: 156 
    ## Seed:  137 139 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  128 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  227 227 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 203 || Leaf:  199 199 || Maternal:  199 || Paternal: 199 
    ## Seed:  185 195 || Leaf:  185 191 || Maternal:  185 || Paternal: 191 
    ## Seed:  173 177 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0313 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  106 108 || Maternal:  108 || Paternal: 106 
    ## Seed:  148 154 || Leaf:  154 156 || Maternal:  154 || Paternal: 156 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  128 130 || Maternal:  NA || Paternal: NA 
    ## Seed:  227 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  171 181 || Leaf:  177 181 || Maternal:  181 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0317 ...
    ## ---------------------------------------------
    ## Seed:  110 110 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  148 154 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  139 141 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  165 195 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 191 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  171 177 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0318 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  148 154 || Leaf:  148 148 || Maternal:  148 || Paternal: 148 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  227 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 201 || Leaf:  199 205 || Maternal:  199 || Paternal: 205 
    ## Seed:  191 191 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 238 || Maternal:  225 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0319 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  98 108 || Maternal:  108 || Paternal: 98 
    ## Seed:  148 152 || Leaf:  148 152 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 141 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  138 138 || Leaf:  138 149 || Maternal:  138 || Paternal: 149 
    ## Seed:  231 231 || Leaf:  231 243 || Maternal:  231 || Paternal: 243 
    ## Seed:  199 199 || Leaf:  195 199 || Maternal:  199 || Paternal: 195 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 177 || Leaf:  171 179 || Maternal:  171 || Paternal: 179 
    ## Seed:  225 225 || Leaf:  225 240 || Maternal:  225 || Paternal: 240 
    ## 
    ##  Working on indv: 07-0320 ...
    ## ---------------------------------------------
    ## Seed:  94 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  150 158 || Leaf:  156 158 || Maternal:  158 || Paternal: 156 
    ## Seed:  137 141 || Leaf:  139 141 || Maternal:  141 || Paternal: 139 
    ## Seed:  130 138 || Leaf:  128 130 || Maternal:  130 || Paternal: 128 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  197 199 || Leaf:  199 201 || Maternal:  199 || Paternal: 201 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 177 || Leaf:  173 179 || Maternal:  173 || Paternal: 179 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0321 ...
    ## ---------------------------------------------
    ## Seed:  110 110 || Leaf:  106 110 || Maternal:  110 || Paternal: 106 
    ## Seed:  148 154 || Leaf:  154 156 || Maternal:  154 || Paternal: 156 
    ## Seed:  139 141 || Leaf:  137 141 || Maternal:  141 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  165 195 || Leaf:  175 195 || Maternal:  195 || Paternal: 175 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  171 177 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0322 ...
    ## ---------------------------------------------
    ## Seed:  94 94 || Leaf:  94 100 || Maternal:  94 || Paternal: 100 
    ## Seed:  154 154 || Leaf:  148 154 || Maternal:  154 || Paternal: 148 
    ## Seed:  137 141 || Leaf:  137 141 || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  NA NA || Leaf:  233 245 || Maternal:  NA || Paternal: NA 
    ## Seed:  167 201 || Leaf:  163 201 || Maternal:  201 || Paternal: 163 
    ## Seed:  191 193 || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  181 181 || Leaf:  171 181 || Maternal:  181 || Paternal: 171 
    ## Seed:  225 236 || Leaf:  225 238 || Maternal:  225 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0323 ...
    ## ---------------------------------------------
    ## Seed:  106 108 || Leaf:  106 108 || Maternal:  NA || Paternal: NA 
    ## Seed:  148 150 || Leaf:  148 148 || Maternal:  148 || Paternal: 148 
    ## Seed:  137 139 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  167 167 || Maternal:  167 || Paternal: 167 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  173 177 || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 238 || Maternal:  225 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0325 ...
    ## ---------------------------------------------
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  148 154 || Leaf:  152 154 || Maternal:  154 || Paternal: 152 
    ## Seed:  139 141 || Leaf:  141 141 || Maternal:  141 || Paternal: 141 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  165 195 || Leaf:  165 165 || Maternal:  165 || Paternal: 165 
    ## Seed:  191 191 || Leaf:  185 191 || Maternal:  191 || Paternal: 185 
    ## Seed:  171 177 || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0326 ...
    ## ---------------------------------------------
    ## Seed:  106 108 || Leaf:  108 110 || Maternal:  108 || Paternal: 110 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 139 || Leaf:  137 141 || Maternal:  137 || Paternal: 141 
    ## Seed:  130 130 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  165 205 || Leaf:  165 199 || Maternal:  165 || Paternal: 199 
    ## Seed:  191 193 || Leaf:  193 193 || Maternal:  193 || Paternal: 193 
    ## Seed:  173 177 || Leaf:  171 177 || Maternal:  177 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0331 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  100 108 || Maternal:  108 || Paternal: 100 
    ## Seed:  154 156 || Leaf:  154 156 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 139 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  128 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 243 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0333 ...
    ## ---------------------------------------------
    ## Seed:  110 110 || Leaf:  94 110 || Maternal:  110 || Paternal: 94 
    ## Seed:  148 154 || Leaf:  154 156 || Maternal:  154 || Paternal: 156 
    ## Seed:  139 141 || Leaf:  137 141 || Maternal:  141 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  165 195 || Leaf:  165 175 || Maternal:  165 || Paternal: 175 
    ## Seed:  191 191 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  171 177 || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0334 ...
    ## ---------------------------------------------
    ## Seed:  98 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  148 158 || Leaf:  148 156 || Maternal:  148 || Paternal: 156 
    ## Seed:  137 139 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  130 145 || Maternal:  130 || Paternal: 145 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 205 || Leaf:  165 205 || Maternal:  205 || Paternal: 165 
    ## Seed:  191 193 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 240 || Maternal:  225 || Paternal: 240 
    ## 
    ##  Working on indv: 07-0335 ...
    ## ---------------------------------------------
    ## Seed:  94 94 || Leaf:  94 108 || Maternal:  94 || Paternal: 108 
    ## Seed:  152 158 || Leaf:  152 158 || Maternal:  NA || Paternal: NA 
    ## Seed:  137 137 || Leaf:  137 139 || Maternal:  137 || Paternal: 139 
    ## Seed:  128 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  175 205 || Leaf:  199 205 || Maternal:  205 || Paternal: 199 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  177 177 || Leaf:  171 177 || Maternal:  177 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 238 || Maternal:  225 || Paternal: 238 
    ## 
    ##  Working on indv: 07-0336 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  154 156 || Leaf:  150 154 || Maternal:  154 || Paternal: 150 
    ## Seed:  137 139 || Leaf:  137 139 || Maternal:  NA || Paternal: NA 
    ## Seed:  128 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  199 201 || Leaf:  199 199 || Maternal:  199 || Paternal: 199 
    ## Seed:  191 195 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  173 173 || Leaf:  173 177 || Maternal:  173 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0337 ...
    ## ---------------------------------------------
    ## Seed:  108 108 || Leaf:  106 108 || Maternal:  108 || Paternal: 106 
    ## Seed:  154 156 || Leaf:  152 156 || Maternal:  156 || Paternal: 152 
    ## Seed:  137 141 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  128 138 || Leaf:  138 138 || Maternal:  138 || Paternal: 138 
    ## Seed:  231 231 || Leaf:  231 245 || Maternal:  231 || Paternal: 245 
    ## Seed:  175 201 || Leaf:  165 175 || Maternal:  175 || Paternal: 165 
    ## Seed:  191 191 || Leaf:  185 191 || Maternal:  191 || Paternal: 185 
    ## Seed:  177 181 || Leaf:  175 177 || Maternal:  177 || Paternal: 175 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0338 ...
    ## ---------------------------------------------
    ## Seed:  100 108 || Leaf:  108 108 || Maternal:  108 || Paternal: 108 
    ## Seed:  156 160 || Leaf:  156 156 || Maternal:  156 || Paternal: 156 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  165 197 || Leaf:  197 199 || Maternal:  197 || Paternal: 199 
    ## Seed:  189 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  177 179 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  238 240 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0339 ...
    ## ---------------------------------------------
    ## Seed:  106 108 || Leaf:  94 108 || Maternal:  108 || Paternal: 94 
    ## Seed:  148 154 || Leaf:  148 148 || Maternal:  148 || Paternal: 148 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  185 193 || Leaf:  185 185 || Maternal:  185 || Paternal: 185 
    ## Seed:  173 173 || Leaf:  173 173 || Maternal:  173 || Paternal: 173 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## 
    ##  Working on indv: 07-0340 ...
    ## ---------------------------------------------
    ## Seed:  94 108 || Leaf:  94 94 || Maternal:  94 || Paternal: 94 
    ## Seed:  152 152 || Leaf:  150 152 || Maternal:  152 || Paternal: 150 
    ## Seed:  137 139 || Leaf:  137 139 || Maternal:  NA || Paternal: NA 
    ## Seed:  128 138 || Leaf:  130 138 || Maternal:  138 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 243 || Maternal:  231 || Paternal: 243 
    ## Seed:  181 199 || Leaf:  167 199 || Maternal:  199 || Paternal: 167 
    ## Seed:  191 191 || Leaf:  191 193 || Maternal:  191 || Paternal: 193 
    ## Seed:  171 181 || Leaf:  171 181 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0341 ...
    ## ---------------------------------------------
    ## Seed:  108 110 || Leaf:  94 110 || Maternal:  110 || Paternal: 94 
    ## Seed:  148 154 || Leaf:  152 154 || Maternal:  154 || Paternal: 152 
    ## Seed:  137 141 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  167 205 || Leaf:  167 199 || Maternal:  167 || Paternal: 199 
    ## Seed:  191 193 || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  177 179 || Leaf:  171 177 || Maternal:  177 || Paternal: 171 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0343 ...
    ## ---------------------------------------------
    ## Seed:  108 110 || Leaf:  108 110 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  137 141 || Leaf:  141 141 || Maternal:  141 || Paternal: 141 
    ## Seed:  130 130 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  191 193 || Leaf:  191 193 || Maternal:  NA || Paternal: NA 
    ## Seed:  177 179 || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0345 ...
    ## ---------------------------------------------
    ## Seed:  106 108 || Leaf:  108 110 || Maternal:  108 || Paternal: 110 
    ## Seed:  148 154 || Leaf:  154 154 || Maternal:  154 || Paternal: 154 
    ## Seed:  137 137 || Leaf:  137 141 || Maternal:  137 || Paternal: 141 
    ## Seed:  130 138 || Leaf:  130 130 || Maternal:  130 || Paternal: 130 
    ## Seed:  231 231 || Leaf:  229 231 || Maternal:  231 || Paternal: 229 
    ## Seed:  165 201 || Leaf:  201 205 || Maternal:  201 || Paternal: 205 
    ## Seed:  185 193 || Leaf:  193 193 || Maternal:  193 || Paternal: 193 
    ## Seed:  173 173 || Leaf:  173 179 || Maternal:  173 || Paternal: 179 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0346 ...
    ## ---------------------------------------------
    ## Seed:  106 108 || Leaf:  104 108 || Maternal:  108 || Paternal: 104 
    ## Seed:  148 154 || Leaf:  154 154 || Maternal:  154 || Paternal: 154 
    ## Seed:  137 137 || Leaf:  137 137 || Maternal:  137 || Paternal: 137 
    ## Seed:  130 138 || Leaf:  138 138 || Maternal:  138 || Paternal: 138 
    ## Seed:  231 231 || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  165 201 || Leaf:  165 165 || Maternal:  165 || Paternal: 165 
    ## Seed:  185 193 || Leaf:  185 185 || Maternal:  185 || Paternal: 185 
    ## Seed:  NA NA || Leaf:  NA NA || Maternal:  NA || Paternal: NA 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## 
    ##  Working on indv: 07-0348 ...
    ## ---------------------------------------------
    ## Seed:  94 108 || Leaf:  108 110 || Maternal:  108 || Paternal: 110 
    ## Seed:  152 156 || Leaf:  148 152 || Maternal:  152 || Paternal: 148 
    ## Seed:  137 143 || Leaf:  141 143 || Maternal:  143 || Paternal: 141 
    ## Seed:  130 138 || Leaf:  130 138 || Maternal:  NA || Paternal: NA 
    ## Seed:  NA NA || Leaf:  231 231 || Maternal:  231 || Paternal: 231 
    ## Seed:  199 203 || Leaf:  165 199 || Maternal:  199 || Paternal: 165 
    ## Seed:  191 191 || Leaf:  191 191 || Maternal:  191 || Paternal: 191 
    ## Seed:  171 177 || Leaf:  177 177 || Maternal:  177 || Paternal: 177 
    ## Seed:  225 225 || Leaf:  225 225 || Maternal:  225 || Paternal: 225 
    ## Percentage of locus-ind combinations that were heterozygotes at same alleles: 0.09965636 
    ## Writing files...

``` r
summary(gametes_out)
```

    ##          Length Class      Mode
    ## maternal 13     data.frame list
    ## paternal 13     data.frame list
    ## leaf     20     genalex    list
    ## seed     20     genalex    list
    ## partial   5     data.frame list
