# Accompanying R script for the manuscript:
# "Complex urban environments provide Apis mellifera with a richer plant 
# forage than suburban and more rural landscapes. Diet analysis using a pair 
# of metabarcoding markers provides a broad detection spectrum."

# Graeme Fox, Richard F. Preziosi, Loreto Ros, Joshua Sammy, 
# Jennifer K. Rowntree, and Latha R. Vellaniparambil

# Corresponding author: j.rowntree@mmu.ac.uk

# Script generated, run, and tested on Ubuntu 21.04 with R 4.0.4
# in RStudio 2021.09.0 Build 351

# This script and intermediary data files
# are available to download from GitHub:
# www.github.com/graemefox/Manchester_honey/data

# Download all files to a directory on your disk, setwd() to that location
# and file paths should work 

# raw sequence data will be available from the SRA upon publication
# (BioProject: PRJNA767686) - not necessary to download these to run this script

## install required Ubuntu packages
# sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev

###############################
# INSTALL REQUIRED R PACKAGES #
###############################
requiredpackages <- c("DescTools",
                      "ggplot2",
                      "ggrepel",
                      "gridExtra",
                      "vegan",
                      "psych",
                      "reshape2",
                      "multcomp",
                      "dplyr",
                      "grid",
                      "ggpubr")
for (pkg in requiredpackages) {
  if (pkg %in% rownames(installed.packages()) == FALSE)
  {install.packages(pkg)}
  if (pkg %in% rownames(.packages()) == FALSE)
  {library(pkg, character.only=T)}
}
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools", dependencies=T)}
#devtools::install_github("jbisanz/qiime2R")
library("qiime2R")

########### #
# FUNCTIONS #
#############

## 'not in' function
'%!in%' <- Negate('%in%')

## Count reads in each FASTQ file in list of FASTQ files
## specific to these data files as looks for the '@HISEQ' header
count_reads <- function(list_of_FASTQs){
  for (sequence_file in list_of_FASTQs){
    command <- paste("zgrep -c \"@HISEQ\" ", sequence_file)
    print(sequence_file)
    system(command)
  }
}

## take a community matrix detaframe (or more likely a region of a dataframe)
## containing raw sequence reads per sample (columns) / taxon (rows), 
# compute frequencies PER SAMPLE (Ie. in columns) and remove (over-write with 0.00)
# any values below desired threshold (1% recommended by Taberlet et al.)
remove_low_freq_taxa_per_sample <- function(data_frame, column_to_start, column_to_end, cutoff_freq){
  col_counter = column_to_start
  for (row in data_frame[column_to_start:column_to_end]){
    row_counter = 0
    row <- as.data.frame(row)
    for (value in t(row)){
      row_counter = row_counter + 1
      if ((value/rowSums(t(row))*100)<cutoff_freq){
        ## check values match before making replacement
        if (value != data_frame[row_counter, col_counter]){
          print("Something gone wrong. Values <1% to be replaced do not match.")
        }
        else {
          data_frame[row_counter, col_counter] <- 0.00
        }
      }
    }
    col_counter = col_counter + 1
  } 
  return(data_frame)
}

remove_low_freq_taxa_per_ASV <- function(data_frame, column_to_start, column_to_end, cutoff_freq){
  col_counter = column_to_start
  row_counter = 1
  while (row_counter <= nrow(data_frame)){
    row = 0
    for (value in data_frame[row_counter, column_to_start:column_to_end]){
      row = row + 1
      ## if > 0 - there are some empty ASVs and cannot divide by zero
      if (rowSums(data_frame[row_counter, column_to_start:column_to_end]) > 0){
        if (value/rowSums(data_frame[row_counter, column_to_start:column_to_end])*100<cutoff_freq){
          ## check values match before making replacement
          if (value != data_frame[row_counter, column_to_start:column_to_end][row]){
            print("Something gone wrong. Values <1% to be replaced do not match.")
          }
          else {
            data_frame[row_counter, column_to_start:column_to_end][row] <- 0.00
          }
        }
      }
    }
    col_counter = col_counter + 1
    row_counter = row_counter + 1
  } 
  return(data_frame)
}

## find non-common elements in two lists
## I.e. not the intersect but the outersect
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

########
# MAIN #
########

## raw FastQ files previously processed with both Trimmomatic and cutadapt
# Trimmomatic parameters: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:100"

# cutadapt parameters: 
# cutadapt -g forward_primer_seq -G reverse_primer_seq -o output_R1.fq.gz -p output_R2.fq.gz input_R1.fq.gz input_R2.fq.gz -m 1 -j 0

### set to the location of the data from the GitHub link above
setwd()

### get read counts from the FASTQs
## This section is just counting reads. Values have also been hard-coded to
# save downloading sequence files

## raw FastQ files (available from NCBI SRA) have been 
# processed with both Trimmomatic and Cutadapt
# these are referred to as "Reads" from this point onwards

# Trimmomatic parameters: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:100"

# cutadapt parameters: 

# cutadapt -g forward_primer_seq -G reverse_primer_seq -o output_R1.fq.gz -p output_R2.fq.gz input_R1.fq.gz input_R2.fq.gz -m 1 -j 0

## These "reads" were then processed using the dada2 method built into QIIME2
## Reads passing this quality process are from this point referred to as QIIME2_QC_reads

### get raw, and post QIIME sequence counts
### (values hard coded below to save re-running this)
#count_reads(rbcL_reads)
#count_reads(ITS2p_reads)
#count_reads(rbcL_dada2_qiime2_QC_reads)
#count_reads(ITS2p_dada2_qiime2_QC_reads)

hive_names <- c("PW", "MC", "CM", "WC", "AG", "MM", "CC", "PH", "ML", "HH", "JH", "MB", "SK", "SP")

### read in the data generated by QIIME2
rbcL_dada2_SVs <- read_qza("data/rbcL_table-dada2.qza")
ITS2p_dada2_SVs <- read_qza("data/ITS2p_table-dada2.qza")

## there are two extra hives in this data which ultimately were not used 
# in the analysis - hive07 and hive10
# these are removed by subset later when these files are used

## put the raw sequence counts, and reads assigned to ASV together into a stats table (more to be added later)
## these values from from the count_reads function above but commented out
rbcL_dada2_stats_NC <- data.frame(hive_names, 
                                  "Reads"=c(1918943, 1829390, 1252605, 1085788, 1246936, 618607, 798623, 287206, 824897, 670489, 658049, 1212965, 872350, 1044337),
                                  "QIIME2_QC_reads"=c(232537, 236317, 123027, 113504, 135413, 48091, 83634, 47148, 78546, 67577, 98836, 130220, 126038, 169942))

ITS2p_dada2_stats_NC <- data.frame(hive_names, 
                                   "Reads"=c(913410, 893437, 974364, 431841, 474726, 952382, 682414, 592399, 1218990, 845782, 548228, 1026872, 829165, 989794),
                                   "QIIME2_QC_reads"=c(342250, 389051, 448178, 218144, 130053, 300172, 297665, 309897, 530739, 177575, 267257, 542951, 361542, 525675))

## mean reads per hive and SD
mean(rbcL_dada2_stats_NC$Reads)
mean(ITS2p_dada2_stats_NC$Reads)
sd(rbcL_dada2_stats_NC$Reads)
sd(ITS2p_dada2_stats_NC$Reads)

### percentage of raw reads assigned to ASV (ie. 100-proportion removed by dada2)
#sum(as.data.frame(rbcL_dada2_SVs$data))/sum(rbcL_dada2_stats$Reads)*100
#sum(as.data.frame(ITS2p_dada2_SVs$data))/sum(ITS2p_dada2_stats$Reads)*100

### read in the family, genus, and species level identifications - Generated by QIIME2 workflow
rbcL_dada2_fam <- read.table("data/taxa_assignments/rbcL_dada2_family_assignments.csv", header=T, sep=",")
ITS2p_dada2_fam <- read.table("data/taxa_assignments/ITS2p_dada2_family_assignments.csv", header=T, sep=",")
rbcL_dada2_gen <- read.table("data/taxa_assignments/rbcL_dada2_genus_assignments.csv", header=T, sep=",")
ITS2p_dada2_gen <- read.table("data/taxa_assignments/ITS2p_dada2_genus_assignments.csv", header=T, sep=",")
rbcL_dada2_spec <- read.table("data/taxa_assignments/rbcL_dada2_species_assignments.csv", header=T, sep=",")
ITS2p_dada2_spec <- read.table("data/taxa_assignments/ITS2p_dada2_species_assignments.csv", header=T, sep=",")

## 
rbcL_dada2_fam_NC <- read.table("data/taxa_assignments/rbcL_dada2_family_assignments_gen_not_collapsed.csv", header=T, sep=",")
ITS2p_dada2_fam_NC <- read.table("data/taxa_assignments/ITS2p_dada2_family_assignments_gen_not_collapsed.csv", header=T, sep=",")
rbcL_dada2_gen_NC <- read.table("data/taxa_assignments/rbcL_dada2_genus_assignments_spec_not_collapsed.csv", header=T, sep=",")
ITS2p_dada2_gen_NC <- read.table("data/taxa_assignments/ITS2p_dada2_genus_assignments_spec_not_collapsed.csv", header=T, sep=",")

### illustrate how the above files work
## For example if we look at rbcL ASVs which achieved species level assignment in the genus Rubus.
## There are three species found
#rbcL_dada2_spec[grepl("Rubus", rbcL_dada2_spec$Species, fixed = TRUE),]

## The rbcL totals for the the genus Rubus are larger than the sum of the three species above - this is because
## they ALSO include the reads associated with ASVs which were classified in the Rubus genus, but didn't get species level classification
#rbcL_dada2_gen_NC[which(rbcL_dada2_gen_NC$Genus=="Rubus"),]

## For example, in hive PW, there were 98362 reads associated with Rubus armeniacus (and zero associated with the other two species) 
# and 162633 reads associated with the Rubus genus.
# There were therefore 162633-98362=64271 reads assigned in the Rubus genus but which did not get a species assignment.

## how many unique taxa identified?
## and what resolution did each achieve
nrow(rbcL_dada2_spec[which(rowSums(rbcL_dada2_spec[13:26])>0),])  # 36 species actually found in rbcL
nrow(ITS2p_dada2_spec[which(rowSums(ITS2p_dada2_spec[13:26])>0),]) # 31 species actually found in ITS2p

### how many total unique species were found across the 14 samples?
nrow(rbcL_dada2_spec[which(rowSums(rbcL_dada2_spec[13:26])>0),])+nrow(ITS2p_dada2_spec[which(rowSums(ITS2p_dada2_spec[13:26])>0),])-
  length(intersect(rbcL_dada2_spec[which(rowSums(rbcL_dada2_spec[13:26])>0),]$Species, ITS2p_dada2_spec[which(rowSums(ITS2p_dada2_spec[13:26])>0),]$Species))

## as above but for ASVs which were assigned genus but not species
nrow(rbcL_dada2_gen[which(rowSums(rbcL_dada2_gen[11:24])>0),])  # 43 ASVs assigned genera, but not species
nrow(ITS2p_dada2_gen[which(rowSums(ITS2p_dada2_gen[11:24])>0),]) # 19 ASVs assigned genera, but not species

### how many total unique genera, which did not get species resolution, were found across the 14 samples?
nrow(rbcL_dada2_gen[which(rowSums(rbcL_dada2_gen[11:24])>0),])+nrow(ITS2p_dada2_gen[which(rowSums(ITS2p_dada2_gen[11:24])>0),])-
  length(intersect(rbcL_dada2_gen[which(rowSums(rbcL_dada2_gen[11:24])>0),]$genies, ITS2p_dada2_gen[which(rowSums(ITS2p_dada2_gen[11:24])>0),]$genies))

## as above but for ASVs which were assigned family but not genus
nrow(rbcL_dada2_fam[which(rowSums(rbcL_dada2_fam[9:22])>0),])  # 3 ASVs assigned family, but not genus
nrow(ITS2p_dada2_fam[which(rowSums(ITS2p_dada2_fam[9:22])>0),]) # 11 ASVs assigned family, but not genus

### how many total unique family, which did not get genus resolution, were found across the 14 samples?
nrow(rbcL_dada2_fam[which(rowSums(rbcL_dada2_fam[9:22])>0),])+nrow(ITS2p_dada2_fam[which(rowSums(ITS2p_dada2_fam[9:22])>0),])-
  length(intersect(rbcL_dada2_fam[which(rowSums(rbcL_dada2_fam[9:22])>0),]$famies, ITS2p_dada2_fam[which(rowSums(ITS2p_dada2_fam[9:22])>0),]$famies))

# sum the above
64+62+14  # =140

64/140*100 #  45.71 taxa achieved species assignment
## 24% achieved species
## 64+62=126 achieved genus
126/140*100   # 90.0% taxa achieved genus level assignment

## Look at read counts assigned with ASVs getting low resolution
## family
max(colSums(rbcL_dada2_fam[which(rowSums(rbcL_dada2_fam[9:22])>0),][9:22])/rbcL_dada2_stats_NC$QIIME2_QC_reads*100)
mean(colSums(rbcL_dada2_fam[which(rowSums(rbcL_dada2_fam[9:22])>0),][9:22])/rbcL_dada2_stats_NC$QIIME2_QC_reads*100)
sd(colSums(rbcL_dada2_fam[which(rowSums(rbcL_dada2_fam[9:22])>0),][9:22])/rbcL_dada2_stats_NC$QIIME2_QC_reads*100)
max(colSums(ITS2p_dada2_fam[which(rowSums(ITS2p_dada2_fam[9:22])>0),][9:22])/ITS2p_dada2_stats_NC$QIIME2_QC_reads*100)
mean(colSums(ITS2p_dada2_fam[which(rowSums(ITS2p_dada2_fam[9:22])>0),][9:22])/ITS2p_dada2_stats_NC$QIIME2_QC_reads*100)
sd(colSums(ITS2p_dada2_fam[which(rowSums(ITS2p_dada2_fam[9:22])>0),][9:22])/ITS2p_dada2_stats_NC$QIIME2_QC_reads*100)
# genus
max(colSums(rbcL_dada2_gen[which(rowSums(rbcL_dada2_gen[11:24])>0),][11:24])/rbcL_dada2_stats_NC$QIIME2_QC_reads*100)
mean(colSums(rbcL_dada2_gen[which(rowSums(rbcL_dada2_gen[11:24])>0),][11:24])/rbcL_dada2_stats_NC$QIIME2_QC_reads*100)
sd(colSums(rbcL_dada2_gen[which(rowSums(rbcL_dada2_gen[11:24])>0),][11:24])/rbcL_dada2_stats_NC$QIIME2_QC_reads*100)
max(colSums(ITS2p_dada2_gen[which(rowSums(ITS2p_dada2_gen[11:24])>0),][11:24])/ITS2p_dada2_stats_NC$QIIME2_QC_reads*100)
mean(colSums(ITS2p_dada2_gen[which(rowSums(ITS2p_dada2_gen[11:24])>0),][11:24])/ITS2p_dada2_stats_NC$QIIME2_QC_reads*100)
sd(colSums(ITS2p_dada2_gen[which(rowSums(ITS2p_dada2_gen[11:24])>0),][11:24])/ITS2p_dada2_stats_NC$QIIME2_QC_reads*100)

## table of proportion of reads achieving family/genus/species resolution
sum(rowSums(rbcL_dada2_fam[9:22]))/sum(rbcL_dada2_stats_NC$QIIME2_QC_reads)*100
sum(rowSums(rbcL_dada2_gen[11:24]))/sum(rbcL_dada2_stats_NC$QIIME2_QC_reads)*100
sum(rowSums(rbcL_dada2_spec[13:26]))/sum(rbcL_dada2_stats_NC$QIIME2_QC_reads)*100

# 57.6023% assigned species
# 57.6023 + 39.0177 = 96.62% assigned genus

sum(rowSums(ITS2p_dada2_fam[9:22]))/sum(ITS2p_dada2_stats_NC$QIIME2_QC_reads)*100
sum(rowSums(ITS2p_dada2_gen[11:24]))/sum(ITS2p_dada2_stats_NC$QIIME2_QC_reads)*100
sum(rowSums(ITS2p_dada2_spec[13:26]))/sum(ITS2p_dada2_stats_NC$QIIME2_QC_reads)*100

# 47.12782% assigned species
# 47.12782 + 38.53901 = 85.66683% assigned genus

## Filter by the UK plausibility at either species or genera
## proportion reads removed at species filter
sum(rowSums(rbcL_dada2_spec[which(rbcL_dada2_spec$Species.In.UK=="N"),][13:26]))/sum(rowSums(rbcL_dada2_spec[13:26]))*100
sum(rowSums(ITS2p_dada2_spec[which(ITS2p_dada2_spec$Species.In.UK=="N"),][13:26]))/sum(rowSums(ITS2p_dada2_spec[13:26]))*100

## proportion reads removed at genus filter
sum(rowSums(rbcL_dada2_gen[which(rbcL_dada2_gen$Genus.In.UK =="N"),][11:24]))/sum(rowSums(rbcL_dada2_gen[11:24]))*100
sum(rowSums(ITS2p_dada2_gen[which(ITS2p_dada2_gen$Genus.In.UK =="N"),][11:24]))/sum(rowSums(ITS2p_dada2_gen[11:24]))*100

### remove anything not plausible in the UK for further analysis
rbcL_dada2_spec <- rbcL_dada2_spec[which(rbcL_dada2_spec$Species.In.UK=="Y"),]
ITS2p_dada2_spec <- ITS2p_dada2_spec[which(ITS2p_dada2_spec$Species.In.UK=="Y"),]
rbcL_dada2_gen_NC <- rbcL_dada2_gen_NC[which(rbcL_dada2_gen_NC$Genus.In.UK=="Y"),]
ITS2p_dada2_gen_NC <- ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Genus.In.UK=="Y"),]
rbcL_dada2_gen <- rbcL_dada2_gen[which(rbcL_dada2_gen$Genus.In.UK=="Y"),]
ITS2p_dada2_gen <- ITS2p_dada2_gen[which(ITS2p_dada2_gen$Genus.In.UK=="Y"),]

#### filter combined OTUs in two ways
## first remove reads where ratio between OTU abundance in a sample and total OTU abundance across all samples is <0.03% (Taberlet 2018)
rbcL_dada2_spec <- remove_low_freq_taxa_per_ASV(rbcL_dada2_spec, 13, 26, 0.03)
rbcL_dada2_gen_NC <- remove_low_freq_taxa_per_ASV(rbcL_dada2_gen_NC, 12, 25, 0.03)
rbcL_dada2_fam_NC <- remove_low_freq_taxa_per_ASV(rbcL_dada2_fam_NC, 10, 23, 0.03)
rbcL_dada2_gen <- remove_low_freq_taxa_per_ASV(rbcL_dada2_gen, 11, 24, 0.03)
rbcL_dada2_fam <- remove_low_freq_taxa_per_ASV(rbcL_dada2_fam, 9, 22, 0.03)

ITS2p_dada2_spec <- remove_low_freq_taxa_per_ASV(ITS2p_dada2_spec, 13, 26, 0.03)
ITS2p_dada2_gen_NC <- remove_low_freq_taxa_per_ASV(ITS2p_dada2_gen_NC, 12, 25, 0.03)
ITS2p_dada2_fam_NC <- remove_low_freq_taxa_per_ASV(ITS2p_dada2_fam_NC, 10, 23, 0.03)
ITS2p_dada2_gen <- remove_low_freq_taxa_per_ASV(ITS2p_dada2_gen, 11, 24, 0.03)
ITS2p_dada2_fam <- remove_low_freq_taxa_per_ASV(ITS2p_dada2_fam, 9, 22, 0.03)

## secondly, remove species present at low frequency (<1%) OTUs within a sample
rbcL_dada2_spec <- remove_low_freq_taxa_per_sample(rbcL_dada2_spec, 13, 26, 1)
ITS2p_dada2_spec <- remove_low_freq_taxa_per_sample(ITS2p_dada2_spec, 13, 26, 1)

#### remove any ASVs which do not have any data in either marker 
### the removal of some data as part of low frequency ASV can produce empty rows

### remove empty species rows from species level data
spec_to_remove <- cbind(rbcL_dada2_spec$Species, rbcL_dada2_spec[13:26],
                        ITS2p_dada2_spec[13:26])[which(rowSums(cbind(rbcL_dada2_spec$Species, 
                       rbcL_dada2_spec[13:26], ITS2p_dada2_spec[13:26])[2:29])==0),][1]

rbcL_dada2_spec <- rbcL_dada2_spec[which(rbcL_dada2_spec$Species %!in% spec_to_remove$`rbcL_dada2_spec$Species`),]
ITS2p_dada2_spec <- ITS2p_dada2_spec[which(ITS2p_dada2_spec$Species %!in% spec_to_remove$`rbcL_dada2_spec$Species`),]

### remove empty genus rows from genus level data
gen_to_remove <- cbind(rbcL_dada2_gen_NC$Genus, rbcL_dada2_gen_NC[12:25],
                       ITS2p_dada2_gen_NC[12:25])[which(rowSums(cbind(rbcL_dada2_gen_NC$Genus, 
                       rbcL_dada2_gen_NC[12:25], ITS2p_dada2_gen_NC[12:25])[2:29])==0),][1]

rbcL_dada2_gen_NC <- rbcL_dada2_gen_NC[which(rbcL_dada2_gen_NC$Genus %!in% gen_to_remove$`rbcL_dada2_gen_NC$Genus`), ]
ITS2p_dada2_gen_NC <- ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Genus %!in% gen_to_remove$`rbcL_dada2_gen_NC$Genus`), ]

### remove empty family rows from family level data
fam_to_remove <- cbind(rbcL_dada2_fam_NC$Family, rbcL_dada2_fam_NC[10:23],
                       ITS2p_dada2_fam_NC[10:23])[which(rowSums(cbind(rbcL_dada2_fam_NC$Family, 
                       rbcL_dada2_fam_NC[10:23], ITS2p_dada2_fam_NC[10:23])[2:29])==0),][1]

rbcL_dada2_fam_NC <- rbcL_dada2_fam_NC[which(rbcL_dada2_fam_NC$Family %!in% fam_to_remove$`rbcL_dada2_fam_NC$Family`), ]
ITS2p_dada2_fam_NC <- ITS2p_dada2_fam_NC[which(ITS2p_dada2_fam_NC$Family %!in% fam_to_remove$`rbcL_dada2_fam_NC$Family`),]

## how many species/genera survive the filtering?
nrow(rbcL_dada2_spec[which(rowSums(rbcL_dada2_spec[13:26])>0),])  # 16 species survived the filter (rbcL)
nrow(ITS2p_dada2_spec[which(rowSums(ITS2p_dada2_spec[13:26])>0),]) # 10 species actually found in ITS2p

### how many total unique species were found across the 14 samples?
nrow(rbcL_dada2_spec[which(rowSums(rbcL_dada2_spec[13:26])>0),])+nrow(ITS2p_dada2_spec[which(rowSums(ITS2p_dada2_spec[13:26])>0),])-
  length(intersect(rbcL_dada2_spec[which(rowSums(rbcL_dada2_spec[13:26])>0),]$Species, ITS2p_dada2_spec[which(rowSums(ITS2p_dada2_spec[13:26])>0),]$Species))

## how species survived plausbility filter in each hive?
## note these are just species - there are also unresolved genera/families
colSums(rbcL_dada2_spec[13:26]+ITS2p_dada2_spec[13:26] != 0)

### calculate number of reads surviving QC steps
### calculate number of reads associated with an ASV artefect
## some samples removed here due to us not being sure of the manner of their processing. Samples already removed from other data
rbcL_dada2_ASVs <- c(colSums(rbcL_dada2_SVs$data[,1:6]), colSums(rbcL_dada2_SVs$data[,8:9]), colSums(rbcL_dada2_SVs$data[,11:16]))
ITS2p_dada2_ASVs <- c(colSums(ITS2p_dada2_SVs$data[,1:6]), colSums(ITS2p_dada2_SVs$data[,8:9]), colSums(ITS2p_dada2_SVs$data[,11:16]))

rbcL_dada2_stats_NC <- cbind(rbcL_dada2_stats_NC, "ASVs"=rbcL_dada2_ASVs, "Fam"=colSums(rbcL_dada2_fam_NC[10:23]), "Gen"=colSums(rbcL_dada2_gen_NC[12:25]), "Spec"=colSums(rbcL_dada2_spec[13:26]))
ITS2p_dada2_stats_NC <- cbind(ITS2p_dada2_stats_NC, "ASVs"=ITS2p_dada2_ASVs, "Fam"=colSums(ITS2p_dada2_fam_NC[10:23]), "Gen"=colSums(ITS2p_dada2_gen_NC[12:25]), "Spec"=colSums(ITS2p_dada2_spec[13:26]))

## percentage reads achieving genera or species level ID
sum(rbcL_dada2_stats_NC$Spec)/sum(rbcL_dada2_stats_NC$ASVs)*100
sum(rbcL_dada2_stats_NC$Gen)/sum(rbcL_dada2_stats_NC$ASVs)*100

sum(ITS2p_dada2_stats_NC$Spec)/sum(ITS2p_dada2_stats_NC$ASVs)*100
sum(ITS2p_dada2_stats_NC$Gen)/sum(ITS2p_dada2_stats_NC$ASVs)*100

# create data for plots
time <- c(rep("ASV", each=14), rep("Fam", each=14), rep("Gen", each=14), rep("Spec", each=14))
group <- rep(hive_names, times=4)

### calculate % of raw reads surviving each data processing step
value_rbcL_dada2_NC <- c(((rbcL_dada2_stats_NC$ASVs/rbcL_dada2_stats_NC$QIIME2_QC_reads)*100), ((rbcL_dada2_stats_NC$Fam/rbcL_dada2_stats_NC$QIIME2_QC_reads)*100), ((rbcL_dada2_stats_NC$Gen/rbcL_dada2_stats_NC$QIIME2_QC_reads)*100), ((rbcL_dada2_stats_NC$Spec/rbcL_dada2_stats_NC$QIIME2_QC_reads)*100))
value_ITS2p_dada2_NC <- c((ITS2p_dada2_stats_NC$ASVs/ITS2p_dada2_stats_NC$QIIME2_QC_reads)*100, ((ITS2p_dada2_stats_NC$Fam/ITS2p_dada2_stats_NC$QIIME2_QC_reads)*100), ((ITS2p_dada2_stats_NC$Gen/ITS2p_dada2_stats_NC$QIIME2_QC_reads)*100), ((ITS2p_dada2_stats_NC$Spec/ITS2p_dada2_stats_NC$QIIME2_QC_reads)*100))

## generate plot data
plot_rbcL_dada2_NC <- data.frame(time, value_rbcL_dada2_NC, group)
plot_ITS2p_dada2_NC <- data.frame(time, value_ITS2p_dada2_NC, group)

colnames(plot_rbcL_dada2_NC) <- c("Level", "Value", "Hive")
colnames(plot_ITS2p_dada2_NC) <- c("Level", "Value", "Hive")

rbcL_plot_NC <- ggplot(data=plot_rbcL_dada2_NC, aes(x=Level, y=Value, group=Hive, colour=Hive)) + 
  geom_line(size=1) + 
  geom_label_repel(aes(label = Hive), nudge_x = 1, direction="y", 
  segment.size = 0.2, segment.color = "grey50", 
  na.rm = TRUE, data = plot_rbcL_dada2_NC[which(plot_rbcL_dada2_NC$Level=="Spec"),]) + 
  xlab("rbcL") + ylab("% Raw reads")

ITS2p_plot_NC <- ggplot(plot_ITS2p_dada2_NC, aes(x=Level, y=Value, group=Hive, colour=Hive)) +
  geom_line(size=1) + 
  geom_label_repel(aes(label = Hive), nudge_x = 1, direction="y", 
  segment.size = 0.2, segment.color = "grey50", 
  na.rm = TRUE, data = plot_ITS2p_dada2_NC[which(plot_ITS2p_dada2_NC$Level=="Spec"),]) + 
  xlab("ITS2p") + ylab("% Raw reads")

### plot percentages of reads assigned taxonomy at fam, gen, spec
grid.arrange(rbcL_plot_NC + theme(legend.position="none"), ITS2p_plot_NC + theme(legend.position="none"), ncol=2)

#  Rarefaction plots
rbcL_gen <- t(rbcL_dada2_gen_NC[12:25])
ITS2p_gen <- t(ITS2p_dada2_gen_NC[12:25])
rbcL_raremax <- min(rowSums(rbcL_gen))
ITS2p_raremax <- min(rowSums(ITS2p_gen))

col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
lty <- c("solid", "dashed", "longdash", "dotdash")
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)

par(mfrow=c(1,2)) 
rbcL_gen_rf <- with(pars[1:26,],rarecurve(rbcL_gen, step = 20, sample = F, col = col, lty = lty, ylab="Genera", main="rbcL"))
ITS2p_gen_rf <- with(pars[1:26,],rarecurve(ITS2p_gen, step = 20, sample = F, col = col, lty = lty, ylab="Genera", main="ITS2p"))

### two parallel datasets from here
## 1) just those reads which achieved species level assignment
## 2) those which achieved genus level assignment, inc. spec collapsed to gen

### Shannon diversity per hive
rbcL_spp_sh <- as.numeric(lapply(diversity(t(rbcL_dada2_spec[13:26]), index="shannon"), round, 2))
min(rbcL_spp_sh)
max(rbcL_spp_sh)
ITS2p_spp_sh <- as.numeric(lapply(diversity(t(ITS2p_dada2_spec[13:26]), index="shannon"), round, 2))
min(ITS2p_spp_sh)
max(ITS2p_spp_sh)
rbcL_gen_sh <- as.numeric(lapply(diversity(t(rbcL_dada2_gen_NC[12:25]), index="shannon"), round, 2))
min(rbcL_gen_sh)
max(rbcL_gen_sh)
ITS2p_gen_sh <- as.numeric(lapply(diversity(t(ITS2p_dada2_gen_NC[12:25]), index="shannon"), round, 2))
min(ITS2p_gen_sh)
max(ITS2p_gen_sh)

## how do estimates of Shannon diversity compare between data derived from each marker?
div_data <- cbind(diversity(t(rbcL_dada2_spec[13:26]), index="shannon"), diversity(t(ITS2p_dada2_spec[13:26]), index="shannon"))
div_gen_data <- cbind(diversity(t(rbcL_dada2_gen_NC[12:25]), index="shannon"), diversity(t(ITS2p_dada2_gen_NC[12:25]), index="shannon"))
colnames(div_data) <- c("rbcL", "ITS2p")
colnames(div_gen_data) <- c("rbcL", "ITS2p")
t.test(div_data)
t.test(div_gen_data)

### count genera and species in each hive
spec_gen_counts <- rbind(t(apply(rbcL_dada2_gen_NC[12:25], 2, function(c)sum(c!=0))),
                         t(apply(rbcL_dada2_spec[13:26], 2, function(c)sum(c!=0))),
                         t(apply(ITS2p_dada2_gen_NC[12:25], 2, function(c)sum(c!=0))),
                         t(apply(ITS2p_dada2_spec[13:26], 2, function(c)sum(c!=0))))
row.names(spec_gen_counts) <- c("rbcL Genera", "rbcL Species", "ITS2p Genera", "ITS2p Species")

### create total taxa tables for manuscript supp info
#write.table(file="all_rbcL_fam_data.txt", rbcL_dada2_fam[which(rowSums(rbcL_dada2_fam[9:22])>0),][6:22][-c(3)], sep=",", quote=FALSE)
#write.table(file="all_rbcL_gen_data.txt", rbcL_dada2_gen[which(rowSums(rbcL_dada2_gen[11:24])>0),][6:24][-c(4,5)], sep=",", quote=FALSE)
#write.table(file="all_rbcL_spec_data.txt", rbcL_dada2_spec[which(rowSums(rbcL_dada2_spec[13:26])>0),][7:26][-c(5,6)], sep=",", quote=FALSE)
#write.table(file="all_ITS2p_fam_data.txt", ITS2p_dada2_fam[which(rowSums(ITS2p_dada2_fam[9:22])>0),][6:22][-c(3)], sep=",", quote=FALSE)
#write.table(file="all_ITS2p_gen_data.txt", ITS2p_dada2_gen[which(rowSums(ITS2p_dada2_gen[11:24])>0),][6:24][-c(4,5)], sep=",", quote=FALSE)
#write.table(file="all_ITS2p_spec_data.txt", ITS2p_dada2_spec[which(rowSums(ITS2p_dada2_spec[13:26])>0),][7:26][-c(5,6)], sep=",", quote=FALSE)

### Read in the file of land use
## This data is derived from the UK land cover map and this analysis generated in QGIS using the LECOS plugin
Land<-read.csv("data/Land.csv", header=T)

### how similar are the rbcL vs ITS2p descriptions of each sample?
## caluclate the Jaccard and Bray-Curtis genetic distances between each pair of descriptions of each sample
## do at both species and genus level (many more reads assigned to genus than species)
spec_level_jaccard <- list()
genus_level_jaccard <- list()
family_level_jaccard <- list()
spec_level_BC <- list()
genus_level_BC <- list()
family_level_BC <- list()
for (sample in seq(1:14)){
  spec_jaccard <- (vegdist(t(cbind(rbcL_dada2_spec[sample+12], ITS2p_dada2_spec[sample+12])), method="jaccard", binary=TRUE, diag=FALSE, upper=FALSE, na.rm=FALSE))
  gen_jaccard <- (vegdist(t(cbind(rbcL_dada2_gen_NC[sample+11], ITS2p_dada2_gen_NC[sample+11])), method="jaccard", binary=TRUE, diag=FALSE, upper=FALSE, na.rm=FALSE))
  fam_jaccard <- (vegdist(t(cbind(rbcL_dada2_fam_NC[sample+9], ITS2p_dada2_fam_NC[sample+9])), method="jaccard", binary=TRUE, diag=FALSE, upper=FALSE, na.rm=FALSE))
  spec_BC <- (vegdist(t(cbind(rbcL_dada2_spec[sample+12], ITS2p_dada2_spec[sample+12])), method="bray", binary=TRUE, diag=FALSE, upper=FALSE, na.rm=FALSE))
  gen_BC <- (vegdist(t(cbind(rbcL_dada2_gen_NC[sample+11], ITS2p_dada2_gen_NC[sample+11])), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE))
  fam_BC <- (vegdist(t(cbind(rbcL_dada2_fam_NC[sample+9], ITS2p_dada2_fam_NC[sample+9])), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE))
  
  spec_level_jaccard <- cbind(spec_level_jaccard, spec_jaccard)
  genus_level_jaccard <- cbind(genus_level_jaccard, gen_jaccard)
  family_level_jaccard <- cbind(family_level_jaccard, fam_jaccard)
  spec_level_BC <- cbind(spec_level_BC, spec_BC)
  genus_level_BC <- cbind(genus_level_BC, gen_BC)
  family_level_BC <- cbind(family_level_BC, fam_BC)
}
spec_dissim_data <- as.data.frame(cbind("Spec_Jacc"=unlist(spec_level_jaccard), "Spec_BC"=unlist(spec_level_BC)))
dissim_data <- as.data.frame(cbind("Gen_Jacc"=unlist(genus_level_jaccard), "Gen_BC"=unlist(genus_level_BC)))
rownames(dissim_data) <- hive_names   

#### Mantel tests for similarity between two markers
rbcl_species_dist <- vegdist(t(rbcL_dada2_spec[13:26]), method = "bray")
ITS2p_species_dist <- vegdist(t(ITS2p_dada2_spec[13:26]), method = "bray")
rbcL_genus_dist <- vegdist(t(rbcL_dada2_gen_NC[12:25]), method="bray")
ITS2p_genus_dist <- vegdist(t(ITS2p_dada2_gen_NC[12:25]), method="bray")

rbcl_species_dist_jacc <- vegdist(t(rbcL_dada2_spec[13:26]), method = "jaccard")
ITS2p_species_dist_jacc <- vegdist(t(ITS2p_dada2_spec[13:26]), method = "jaccard")
rbcL_genus_dist_jacc <- vegdist(t(rbcL_dada2_gen_NC[12:25]), method="jaccard")
ITS2p_genus_dist_jacc <- vegdist(t(ITS2p_dada2_gen_NC[12:25]), method="jaccard")

rbcl_ITS2p_species_mantel <- mantel(rbcl_species_dist, ITS2p_species_dist, method="pearson", permutations=999)
rbcl_ITS2p_genus_mantel <- mantel(rbcL_genus_dist, ITS2p_genus_dist, method="pearson", permutations=999)
rbcl_ITS2p_species_mantel_jacc <- mantel(rbcl_species_dist_jacc, ITS2p_species_dist_jacc, method="pearson", permutations=999)
rbcl_ITS2p_genus_mantel_jacc <- mantel(rbcL_genus_dist_jacc, ITS2p_genus_dist_jacc, method="pearson", permutations=999)

## Bray-curtis, species level
rbcl_ITS2p_species_mantel$statistic^2
p.adjust(rbcl_ITS2p_species_mantel$signif)

## Bray-curtis, genus level
rbcl_ITS2p_genus_mantel$statistic^2
p.adjust(rbcl_ITS2p_genus_mantel$signif)

## Jaccard, species level
rbcl_ITS2p_species_mantel_jacc$statistic^2
p.adjust(rbcl_ITS2p_species_mantel_jacc$signif)

## Jaccard, genus level
rbcl_ITS2p_genus_mantel_jacc$statistic^2
p.adjust(rbcl_ITS2p_genus_mantel_jacc$signif)

## print mean and SD
mean(dissim_data$Gen_Jacc)
sd(dissim_data$Gen_Jacc)

mean(dissim_data$Gen_BC)
sd(dissim_data$Gen_BC)

mean(spec_dissim_data$Spec_Jacc)
sd(spec_dissim_data$Spec_Jacc)

mean(spec_dissim_data$Spec_BC)
sd(spec_dissim_data$Spec_BC)

##### look for patterns in families/genera that were particularly good/bad at getting low level resolution
## pull out anything unable to be resolved beyond family in either marker
nrow(rbcL_dada2_fam[which(rowSums(rbcL_dada2_fam[9:22])>0),])
nrow(ITS2p_dada2_fam[which(rowSums(ITS2p_dada2_fam[9:22])>0),])

## proportion of rosaceae, oleaceae and arecaceae unable to achieve genus
rowSums(rbcL_dada2_fam[which(rbcL_dada2_fam$Family=="Rosaceae"),][9:22])/sum(rowSums(rbcL_dada2_gen_NC[which(rbcL_dada2_gen_NC$Family=="Rosaceae"),][12:25]))
rowSums(rbcL_dada2_fam[which(rbcL_dada2_fam$Family=="Oleaceae"),][9:22])/sum(rowSums(rbcL_dada2_gen_NC[which(rbcL_dada2_gen_NC$Family=="Oleaceae"),][12:25]))
rowSums(rbcL_dada2_fam[which(rbcL_dada2_fam$Family=="Arecaceae"),][9:22])/sum(rowSums(rbcL_dada2_gen_NC[which(rbcL_dada2_gen_NC$Family=="Arecaceae"),][12:25]))

## proportion of total reads in unresolvable families
rowSums(rbcL_dada2_fam[which(rbcL_dada2_fam$Family=="Rosaceae"),][9:22])/sum(rbcL_dada2_stats_NC$Reads)*100
rowSums(rbcL_dada2_fam[which(rbcL_dada2_fam$Family=="Oleaceae"),][9:22])/sum(rbcL_dada2_stats_NC$Reads)*100
rowSums(rbcL_dada2_fam[which(rbcL_dada2_fam$Family=="Arecaceae"),][9:22])/sum(rbcL_dada2_stats_NC$Reads)*100

## same as above but for those families ITS2p was unable to resolve
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Fabaceae"),][9:22])/sum(rowSums(ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Family=="Fabaceae"),][12:25]))
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Poaceae"),][9:22])/sum(rowSums(ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Family=="Poaceae"),][12:25]))
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Asteraceae"),][9:22])/sum(rowSums(ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Family=="Asteraceae"),][12:25]))
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Brassicaceae"),][9:22])/sum(rowSums(ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Family=="Brassicaceae"),][12:25]))
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Melastomataceae"),][9:22])/sum(rowSums(ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Family=="Melastomataceae"),][12:25]))
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Acanthaceae"),][9:22])/sum(rowSums(ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Family=="Acanthaceae"),][12:25]))
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Nothofagaceae"),][9:22])/sum(rowSums(ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Family=="Nothofagaceae"),][12:25]))
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Boraginaceae"),][9:22])/sum(rowSums(ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Family=="Boraginaceae"),][12:25]))
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Crassulaceae"),][9:22])/sum(rowSums(ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Family=="Crassulaceae"),][12:25]))
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Pentaphylacaceae"),][9:22])/sum(rowSums(ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Family=="Pentaphylacaceae"),][12:25]))
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Anacardiaceae"),][9:22])/sum(rowSums(ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Family=="Anacardiaceae"),][12:25]))

rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Fabaceae"),][9:22])/sum(ITS2p_dada2_stats_NC$Reads)*100
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Poaceae"),][9:22])/sum(ITS2p_dada2_stats_NC$Reads)*100
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Asteraceae"),][9:22])/sum(ITS2p_dada2_stats_NC$Reads)*100
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Brassicaceae"),][9:22])/sum(ITS2p_dada2_stats_NC$Reads)*100
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Melastomataceae"),][9:22])/sum(ITS2p_dada2_stats_NC$Reads)*100
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Acanthaceae"),][9:22])/sum(ITS2p_dada2_stats_NC$Reads)*100
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Nothofagaceae"),][9:22])/sum(ITS2p_dada2_stats_NC$Reads)*100
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Boraginaceae"),][9:22])/sum(ITS2p_dada2_stats_NC$Reads)*100
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Crassulaceae"),][9:22])/sum(ITS2p_dada2_stats_NC$Reads)*100
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Pentaphylacaceae"),][9:22])/sum(ITS2p_dada2_stats_NC$Reads)*100
rowSums(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family=="Anacardiaceae"),][9:22])/sum(ITS2p_dada2_stats_NC$Reads)*100

### Find genera which were COMPLETELY unresolavable to species
## I.e. there is nothing in that genera in the species level dataset

### genera that are completely absent from the list in the species table
rbcL_gen_absent_from_spec <- as.data.frame(cbind(rbcL_dada2_gen[which(rowSums(rbcL_dada2_gen_NC[12:25])>0),]$Genus, rbcL_dada2_gen[which(rowSums(rbcL_dada2_gen_NC[12:25])>0),]$Genus%!in%rbcL_dada2_spec$Genus))
rbcL_gen_absent_from_spec <- rbcL_gen_absent_from_spec[which(rbcL_gen_absent_from_spec$V2=="TRUE"),]
rbcL_gen_absent_from_spec$V2 <- NULL
## add on genera which are in the table, but have no reads assigned to them
spec_empty_rows <- as.data.frame(rbcL_dada2_spec[which(rowSums(rbcL_dada2_spec[13:26])==0),]$Genus)
colnames(spec_empty_rows) <- c("V1")
rbcL_gen_absent_from_spec <- rbind(as.data.frame(rbcL_gen_absent_from_spec), spec_empty_rows)
rbcL_gen_absent_from_spec

## repeat above but for ITS2p
ITS2p_gen_absent_from_spec <- as.data.frame(cbind(ITS2p_dada2_gen[which(rowSums(ITS2p_dada2_gen_NC[12:25])>0),]$Genus, ITS2p_dada2_gen[which(rowSums(ITS2p_dada2_gen_NC[12:25])>0),]$Genus%!in%ITS2p_dada2_spec$Genus))
ITS2p_gen_absent_from_spec <- ITS2p_gen_absent_from_spec[which(ITS2p_gen_absent_from_spec$V2=="TRUE"),]
ITS2p_gen_absent_from_spec$V2 <- NULL

## add on genera which are in the table, but have no reads assigned to them
spec_empty_rows <- as.data.frame(ITS2p_dada2_spec[which(rowSums(ITS2p_dada2_spec[13:26])==0),]$Genus)
colnames(spec_empty_rows) <- c("V1")
ITS2p_gen_absent_from_spec <- rbind(as.data.frame(ITS2p_gen_absent_from_spec), spec_empty_rows)
ITS2p_gen_absent_from_spec

### get proportion of total reads associated with a completely unresolvable genus
for (genus in rbcL_gen_absent_from_spec$V1){
  print(genus)
  print(rowSums(rbcL_dada2_gen[which(rbcL_dada2_gen$Genus==genus),][11:24])/sum(rbcL_dada2_stats_NC$Reads)*100)
}

for (genus in ITS2p_gen_absent_from_spec$V1){
  print(genus)
  print(rowSums(ITS2p_dada2_gen[which(ITS2p_dada2_gen$Genus==genus),][11:24])/sum(ITS2p_dada2_stats_NC$Reads)*100)
}

## Citrus in Pendlebury hive a particularly interesting example of the above
## high proportion of Citrus spp. in Pendlebury (hive HH) in both markers
rbcL_dada2_gen_NC[which(rbcL_dada2_gen_NC$Genus=="Citrus"),][21]/sum(rbcL_dada2_gen_NC[21])*100
ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Genus=="Citrus"),][21]/sum(ITS2p_dada2_gen_NC[21])*100

## Very high Citrus at species level in rbcL but very low in ITS2p, due to differing resolution in each marker
rbcL_dada2_spec[which(rbcL_dada2_spec$Genus=="Citrus"),][22]/sum(rbcL_dada2_spec[22])*100
ITS2p_dada2_spec[which(ITS2p_dada2_spec$Genus=="Citrus"),][22]/sum(ITS2p_dada2_spec[22])*100

## transform to presence absence data
## combine species data only, for most of the following analysis
rbcL_df <- as.data.frame(rbcL_dada2_spec[13:26])
rownames(rbcL_df) <- rbcL_dada2_spec$Species
ITS2p_df <- as.data.frame(ITS2p_dada2_spec[13:26])
rownames(ITS2p_df) <- ITS2p_dada2_spec$Species
rbcL_df[rbcL_df>0] <-1
ITS2p_df[ITS2p_df>0] <-1

### as above but for genus to calculate shannon div
rbcL_gen_NC_df <- as.data.frame(rbcL_dada2_gen_NC[12:25])
rownames(rbcL_gen_NC_df) <- rbcL_dada2_gen_NC$Genus
ITS2p_gen_NC_df <- as.data.frame(ITS2p_dada2_gen_NC[12:24])
rownames(ITS2p_gen_NC_df) <- ITS2p_dada2_gen_NC$Genus
rbcL_gen_NC_df[rbcL_gen_NC_df>0]<-1
ITS2p_gen_NC_df[ITS2p_gen_NC_df>0]<-1

### combine data from the two markers
combined <- rbcL_df+ITS2p_df
combined[combined>0] <- 1
all_pres_abs <- cbind(rbcL_dada2_spec[3:10], combined)

## number of species and genera in the total data
rowSums(rbcL_df)
rowSums(ITS2p_df)

total_gen_counts <- rowSums(rbcL_gen_NC_df)+rowSums(ITS2p_gen_NC_df)
total_gen_counts[total_gen_counts>0] <- 1
length(total_gen_counts)


### what is the difference between the two markers?
rbcL_uniq_spec <- rbcL_dada2_spec[which(rowSums(rbcL_dada2_spec[13:26])>0),]$Species
ITS2p_uniq_spec <- ITS2p_dada2_spec[which(rowSums(ITS2p_dada2_spec[13:26])>0),]$Species
intersect(rbcL_uniq_spec, ITS2p_uniq_spec)

## two species in the overlap. what proportion of reads are these?
((rowSums(rbcL_dada2_spec[which(rbcL_dada2_spec$Species=="Impatiens glandulifera"),][13:26])+
    rowSums(rbcL_dada2_spec[which(rbcL_dada2_spec$Species=="Trifolium repens"),][13:26]))/sum(rbcL_dada2_spec[13:26]))*100

((rowSums(ITS2p_dada2_spec[which(ITS2p_dada2_spec$Species=="Impatiens glandulifera"),][13:26])+
    rowSums(ITS2p_dada2_spec[which(ITS2p_dada2_spec$Species=="Trifolium repens"),][13:26]))/sum(ITS2p_dada2_spec[13:26]))*100


### same as above but at genera level (also taking into account species reads collapsed into their genus)
rbcL_uniq_gen <- rbcL_dada2_gen_NC[which(rowSums(rbcL_dada2_gen_NC[11:24])>0),]$Genus
ITS2p_uniq_gen <- ITS2p_dada2_gen_NC[which(rowSums(ITS2p_dada2_gen_NC[11:24])>0),]$Genus

## which are common?
intersect_gen <- intersect(rbcL_uniq_gen, ITS2p_uniq_gen)
length(intersect_gen)

# which are not-common?
outersect_gen <- outersect(rbcL_uniq_gen, ITS2p_uniq_gen)

### what total prop of genus assigned reads fall into the overlap?
sum(rowSums(rbcL_dada2_gen_NC[which(rbcL_dada2_gen_NC$Genus%in%intersect_gen),][12:25]))/sum(rbcL_dada2_gen_NC[12:25])*100
sum(rowSums(ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Genus%in%intersect_gen),][12:25]))/sum(ITS2p_dada2_gen_NC[12:25])*100

### what total prop of genus assigned reads fall into the outerlap?
sum(rowSums(rbcL_dada2_gen_NC[which(rbcL_dada2_gen_NC$Genus%in%outersect_gen),][12:25]))/sum(rbcL_dada2_gen_NC[12:25])*100
sum(rowSums(ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Genus%in%outersect_gen),][12:25]))/sum(ITS2p_dada2_gen_NC[12:25])*100

## count spec per gen, and gen per fam
rbcL_spec_per_gen <- data.frame()
for (gen in rbcL_dada2_gen_NC$Genus){
  rbcL_spec_per_gen <- rbind(rbcL_spec_per_gen, nrow(rbcL_dada2_spec[which(rbcL_dada2_spec$Genus==gen & rowSums(rbcL_dada2_spec[13:26])>0),]))
}
rbcL_gen_per_fam <- data.frame()
for (fam in rbcL_dada2_fam_NC$Family){
  rbcL_gen_per_fam <- rbind(rbcL_gen_per_fam, nrow(rbcL_dada2_gen_NC[which(rbcL_dada2_gen_NC$Family==fam & rowSums(rbcL_dada2_gen_NC[12:25])>0),]))
}
ITS2p_gen_per_fam <- data.frame()
for (fam in ITS2p_dada2_fam_NC$Family){
  ITS2p_gen_per_fam <- rbind(ITS2p_gen_per_fam, nrow(ITS2p_dada2_gen_NC[which(ITS2p_dada2_gen_NC$Family==fam & rowSums(ITS2p_dada2_gen_NC[12:25])>0),]))
}
ITS2p_spec_per_gen <- data.frame()
for (gen in ITS2p_dada2_gen_NC$Genus){
  ITS2p_spec_per_gen <- rbind(ITS2p_spec_per_gen, nrow(ITS2p_dada2_spec[which(ITS2p_dada2_spec$Genus==gen & rowSums(ITS2p_dada2_spec[13:26])>0),]))
}
rbcL_spec_per_gen <- cbind(rbcL_spec_per_gen, "Genus"=rbcL_dada2_gen_NC$Genus)
rbcL_gen_per_fam <- cbind(rbcL_gen_per_fam, "Family"=rbcL_dada2_fam_NC$Family)

### generate gen per family, and species per genera heatmaps for comparision between markers
gen_per_fam_data <- cbind(ITS2p_gen_per_fam, rbcL_gen_per_fam)
spec_per_gen_data <- cbind(ITS2p_spec_per_gen, rbcL_spec_per_gen)
colnames(gen_per_fam_data) <- c("ITS2p", "rbcL", "Family")
colnames(spec_per_gen_data) <- c("ITS2p", "rbcL", "Genus")

## remove any empty rows (families with no genera resolved)
gen_per_fam_data <- gen_per_fam_data[which(rowSums(gen_per_fam_data[1:2])>0),]
spec_per_gen_data <- spec_per_gen_data[which(rowSums(spec_per_gen_data[1:2])>0),]
row.names(gen_per_fam_data) <- gen_per_fam_data$Family
row.names(spec_per_gen_data) <- spec_per_gen_data$Genus
gen_per_fam_data$Family <- NULL
spec_per_gen_data$Genus <- NULL

### are the counts significantly different between markers?
fisher.test(gen_per_fam_data)
fisher.test(spec_per_gen_data)

## read in plant metadata info
Taxa <- read.csv("data/plant_metadata.csv", header=T)

## check for species not present in the plant metadata, and add them in if needs be 
# (this should produce 0 rows)
all_pres_abs[which(all_pres_abs$Species %!in% Taxa$Species && is.na(all_pres_abs$Species)==F),]
all_pres_abs[which(all_pres_abs$Genus %!in% Taxa$Genus && is.na(all_pres_abs$Genus)==F),]

### add status, habit, habitat, and flowering period to our raw data
species_level_raw_data <- list()
for (spec in all_pres_abs$Species){
  if (spec %in% Taxa$Species){
    species_level_raw_data <- rbind(species_level_raw_data, cbind(all_pres_abs[which(all_pres_abs$Species==spec),], Taxa[which(Taxa$Species==spec),][5:13], Taxa[which(Taxa$Species==spec),][17]))
  }
  else{
    if (spec!="NA"){
      print("MISSING")
    }
  }
}
species_level_raw_data

#### calculate proportions of native/non-native and neophyte in each hive  - SPECIES LEVEL
Status_summary <- rbind(colSums(species_level_raw_data[which(species_level_raw_data$Status=="Native"),][9:22]), 
                        colSums(species_level_raw_data[which(species_level_raw_data$Status=="Non_native"),][9:22]), 
                        colSums(species_level_raw_data[which(species_level_raw_data$Status=="Neophyte"),][9:22]))

## total species found in native, non-native and neophyte
nrow(species_level_raw_data[which(species_level_raw_data$Status=="Native"),][9:22])
nrow(species_level_raw_data[which(species_level_raw_data$Status=="Non_native"),][9:22])
nrow(species_level_raw_data[which(species_level_raw_data$Status=="Neophyte"),][9:22])

row.names(Status_summary) <- c("Native", "Non_native", "Neophyte")

#Are proportions of native/non native/neophytes independent of hive?
## Fisher's test more appropriate than chi^2 here as the measurements are relatively low values
fisher.test(Status_summary)
## p = 0.9984, so proportions independent of hive

Stat_long <- melt(Status_summary)
colnames(Stat_long) <- c("Status", "Hive", "Count")
#Stat_long <- Stat_long[order(Stat_long$Hive),]

## get total for each hive, to work out proportions later
temp_tots <- aggregate(Stat_long[, 3], list(Stat_long$Hive), sum)
## add to df
Stat_long$Hive_tot <- rep(temp_tots$x, each=3)
## calc props
Stat_long$Prop <- Stat_long$Count/Stat_long$Hive_tot

## do a one-way ANOVA (aov) to see if differences between categories
## check the categories
levels(Stat_long$Status)

## check means and SDs
group_by(Stat_long, Status) %>%
  summarise(
    count = n(),
    mean = mean(Prop, na.rm = TRUE),
    sd = sd(Prop, na.rm = TRUE)
  )

## ANOVA
res.aov <- aov(Prop~Status, data=Stat_long)
summary(res.aov)
## significant p-value indicates that there are differences
## between the group means
## but we don't know which pairs of groups are different

## to do all the pairwise tests, because our ANOVA gave a sig
## result, we can to TUKEY HONEST SIGNIFICANT DIFFERENT (Tukey HSD)
TukeyHSD(res.aov)

## two pairs significant different
## only pair not is neophyte-native

pairwise.t.test(Stat_long$Prop, Stat_long$Status, p.adjust.method="BY")

## are the proportions of native, non native and neophyte correlated with
# the proportion of Urban land use in the immediate 500m buffer?

Land$UB_Prop_500[1:6]
Land$UB_Prop_500[8:15]

Prop_500_UB <- c("0.9921875", "0.9852025", "0.1632812", "0.4813120", 
                 "1.0000000", "0.8465732", "0.549653580", "0.006245121", 
                 "0.188405797", "0.145638629", "0.008500773", "0.059190031", 
                 "0.002321981", "0.014797508")

Stat_long$Prop_500_UB <- rep(Prop_500_UB, each=3)
Stat_long$Prop_500_UB <- factor(Stat_long$Prop_500_UB)

print(corr.test(cbind("Prop_UB"=unique(as.numeric(as.character(Stat_long$Prop_500_UB))), 
                      "Native"=Stat_long[which(Stat_long$Status=="Native"),]$Prop, 
                      "Non_native"=Stat_long[which(Stat_long$Status=="Non_native"),]$Prop, 
                      "Neophyte"=Stat_long[which(Stat_long$Status=="Neophyte"),]$Prop), method="spearman", 
                use="pairwise", alpha=.05, adjust="BY"),short=F)

## remove heaton park from land use data - sample not used
Land <- Land[-c(7), ]

### pull out just the "Prop" columns, at various buffer sizes
buffer_500 <- Land[, grep("Prop_500$", colnames(Land))]
buffer_1000 <- Land[, grep("Prop_1000", colnames(Land))]
buffer_2500 <- Land[, grep("Prop_2500", colnames(Land))]
buffer_5000 <- Land[, grep("Prop_5000", colnames(Land))]
LandDiv <- Land[, grep("DIV", colnames(Land))]

### combine data from the two markers
all_rbcL_taxa <- as.data.frame(rbind(rbcL_dada2_spec[13:26], rbcL_dada2_gen[11:24], rbcL_dada2_fam[9:22]))
all_ITS2p_taxa <- as.data.frame(rbind(ITS2p_dada2_spec[13:26], ITS2p_dada2_gen[11:24], ITS2p_dada2_fam[9:22]))
all_rbcL_taxa[all_rbcL_taxa>0] <- 1
all_ITS2p_taxa[all_ITS2p_taxa>0] <- 1

all_taxa_combined <- all_rbcL_taxa+all_ITS2p_taxa
all_taxa_combined[all_taxa_combined>0] <- 1

### which hive had the most/least unqiue ASVs (taking into account those which only achieved family or genus level assignment)
max(colSums(all_taxa_combined))
min(colSums(all_taxa_combined))

#Diversity indices
LandDiv <- cbind(t(all_taxa_combined), LandDiv)
#LandDiv <- cbind(all_taxa_combined, LandDiv)

## add Shannon per hive
LandDiv$Hive_Shan <- diversity(LandDiv[1:170], index="shannon")
### calc total number of different taxa in each hive (these include species, and genus not given a species)
LandDiv$Hive_Tots <- rowSums(LandDiv[1:170])

### calculate spearman p for hive div against all land div (for multiple p tests)
## replace a NA in the 500M buffer with zero - breaks calculations otherwise
LandDiv$DIV_SH_500[is.na(LandDiv$DIV_SH_500)] <- 0

### do all the corr.test for Hive shannnon diversity against land use diversity in each buffer. Automagically adjusts p-values
sig_vals <- corr.test(LandDiv$Hive_Shan, LandDiv[, grep("DIV_SH_", colnames(LandDiv))], method="pearson", use="pairwise", adjust="BH", alpha=.05)

## the corrected p-vals are in here
sig_vals$p.adj

#### look at diversity in land use with diversity of plant taxa
hive_shan_vs_5000m_shan <- as.data.frame(cbind(LandDiv$Hive_Shan, LandDiv$DIV_SH_5000))

hive_shan_vs_5000m_shan_plot <- ggscatter(hive_shan_vs_5000m_shan, y="V1", x="V2", add="reg.line", xlab="Land Cover Shannon Diversity (5000m)" ,col.method="pearson", 
                                          ylab="Plant Taxa Shannon Diversity", add.params=list(color="blue", fill="darkgray"), conf.int=T, label=hive_names, repel=T)
sig_vals <- corr.test(hive_shan_vs_5000m_shan)

## add in R^2 values of the Pearson's calculated above, and adjust p-values by the BH method
annotation_text <- grobTree(textGrob(paste("R^2 = ", format(round(sig_vals$r[1,2]^2,3),nsmall=2), "p=",format(round(sig_vals$p.adj,5),nsmall=2)), x=0.4,  y=0.85, hjust=0,
                                     gp=gpar(col="black", fontsize=15, fontface="italic")))

hive_shan_vs_5000m_shan_plot <- hive_shan_vs_5000m_shan_plot + annotation_custom(annotation_text)
hive_shan_vs_5000m_shan_plot

#### dominant taxa
## was anything present in every sample?
# what are the most typical species?
all_taxa_with_info <- cbind(rbind(rbcL_dada2_spec[4:10], 
                                  cbind(rbcL_dada2_gen[3:8], "Species"=rep("NA", each=nrow(rbcL_dada2_gen))), 
                                  cbind(rbcL_dada2_fam[3:7], "Genus"=rep("NA", each=nrow(rbcL_dada2_fam)), "Species"=rep("NA", each=nrow(rbcL_dada2_fam)))),
                            all_taxa_combined)

## Anything species present in every sample?
all_taxa_with_info[which((rowSums(all_taxa_with_info[8:21])/14)*100==100 & all_taxa_with_info$Species!="NA"),]
## species present in >50% sample?
all_taxa_with_info[which((rowSums(all_taxa_with_info[8:21])/14)*100>50 & all_taxa_with_info$Species!="NA"),]

## as above but for taxa achieving genus assignments but not species assignments
all_taxa_with_info[which((rowSums(all_taxa_with_info[8:21])/14)*100>66.6 & all_taxa_with_info$Species=="NA" & all_taxa_with_info$Genus!="NA"),]

## highlight Canabis sativa
all_taxa_with_info[which(all_taxa_with_info$Species=="Cannabis sativa"),]

### UK land cover definitions
## (https://www.ceh.ac.uk/sites/default/files/LCM2015_Dataset_Documentation_22May2017.pdf)
BW <- "Broadleaved woodland"
CW <- "Coniferous woodland"
AR <- "Arable and horticulture"
IG <- "Improved grassland"
NG <- "Neutral grassland"
CG <- "Calcareous grassland"
AG <- "Acid grassland"
FS <- "Fen, marsh and swamp"
HR <- "Heather"
HG <- "Heather grassland"
BO <- "Bog"
IR <- "Inland rock"
SW <- "Saltwater"
FW <- "Freshwater"
SS <- "Surpa-littoral sediment"
LS <- "Littoral sediment"
SM <- "Saltmarsh"
UB <- "Urban"
SB <- "Suburban"

land_cover_abbrv <- c("BW","CW","AR","IG","NG","CG","AG","FS","HR","HG","BO","IR","SW","FW","SS","LS","SM","UB","SB")
restricted_land_cover_abbrvs <- c("UB", "IG")
buffer_sizes = c("500m", "1000m", "1500m", "2000m", "2500m", "3000m", "3500m", "4000m", "4500m", "5000m")

restr_buffers <- cbind(Land[, grep("Prop_500$", colnames(Land))], 
                       Land[, grep("Prop_1000", colnames(Land))],
                       Land[, grep("Prop_2000", colnames(Land))],
                       Land[, grep("Prop_5000", colnames(Land))])
restr_buffers <- cbind(Land[, grep("Prop_500$", colnames(Land))], 
                       Land[, grep("Prop_1000", colnames(Land))],
                       Land[, grep("Prop_2000", colnames(Land))],
                       Land[, grep("Prop_5000", colnames(Land))])

#very_restr_buffers <- cbind(Land[, grep("Prop_1000", colnames(Land))])
buffer_sizes = c("500m", "1000m", "2500m", "5000m")

#### get correlations between taxa diversity and land use types
## the p-values have been corrected for multiple testing

### generate the plots for IG (improved grassland) and UB (urban) in 5000m against platn diverity)
for (landtype in restricted_land_cover_abbrvs){
  plots <- list()
  count=0
  for (buffer in restr_buffers[, grep(paste(landtype,"_Prop", sep=""), colnames(restr_buffers))]){
    count = count+1
    temp <- as.data.frame(cbind(LandDiv$Hive_Shan, buffer))
    sig_vals <- capture.output(print(corr.test(temp$V1, temp$buffer, method="pearson", adjust="BY", alpha=.05), short=F))
    p_val <- as.numeric(capture.output(cat(corr.test(temp$V1, temp$buffer, method="pearson", adjust="BY", alpha=.05)$p.adj)))
    R <- as.numeric(strsplit(sig_vals[4], " ")[[1]][2])  # get the pearson r value
    add_text <- grobTree(textGrob(paste("R^2=", R*R, 
                                        " p=",format(round(p_val,3),nsmall=2)), 
                                  x=0.3,  y=0.8, hjust=0,
                                  gp=gpar(col="black", fontsize=15, fontface="italic")))
    plots[[count]] <- ggscatter(temp, y="V1", x="buffer", add="reg.line", 
                                xlab=paste("Proportion of ", landtype," (", buffer_sizes[count], " buffer)", sep=""),
                                ylab="Plant Taxa Shannon Diversity", add.params=list(color="blue", 
                                fill="darkgray"), conf.int=T, label=hive_names, repel=T) 
    plots[[count]] <- plots[[count]] + annotation_custom(add_text)
    stat_cor(aes(label = "blah"), method="spearman", label.x=3)
    
  }
  print(plots)
}

## as above but for all land use types- plot four to a page
## values for BW and NG highlighted in the supplementary material
for (landtype in land_cover_abbrv){
  plots <- list()
  count=0
  for (buffer in restr_buffers[, grep(paste(landtype,"_Prop", sep=""), colnames(restr_buffers))]){
    count = count+1
    temp <- as.data.frame(cbind(LandDiv$Hive_Shan, buffer))
    sig_vals <- capture.output(print(corr.test(temp$V1, temp$buffer, method="pearson", adjust="BY", alpha=.05), short=F))
    p_val <- as.numeric(capture.output(cat(corr.test(temp$V1, temp$buffer, method="pearson", adjust="BY", alpha=.05)$p.adj)))
    R <- as.numeric(strsplit(sig_vals[4], " ")[[1]][2])  # get the pearson r value
    add_text <- grobTree(textGrob(paste("R^2=", R*R, 
                                  " p=",format(round(p_val,4),nsmall=2)), 
                                  x=0.4,  y=0.95, hjust=0,
                                  gp=gpar(col="black", fontsize=15, fontface="italic")))
    plots[[count]] <- ggscatter(temp, y="V1", x="buffer", add="reg.line", 
                                xlab=paste("Proportion of ", landtype," (", buffer_sizes[count], " buffer)", sep=""),
                                ylab="Plant Taxa Shannon Diversity", add.params=list(color="blue", 
                                fill="darkgray"), conf.int=T, label=hive_names, repel=T) 
    plots[[count]] <- plots[[count]] + annotation_custom(add_text)
    stat_cor(aes(label = "blah"), method="spearman", label.x=3)
    
  }
  do.call("grid.arrange", c(plots, ncol=2))
}

### are the land types independent of one another?
colnames(buffer_500) <- c("BW","CW","AR","IG","NG","CG","AG","FS","HR","HG","BO","IR","SW","FW","SS","LS","SM","UB","SB")
colnames(buffer_1000) <- c("BW","CW","AR","IG","NG","CG","AG","FS","HR","HG","BO","IR","SW","FW","SS","LS","SM","UB","SB")
colnames(buffer_2500) <- c("BW","CW","AR","IG","NG","CG","AG","FS","HR","HG","BO","IR","SW","FW","SS","LS","SM","UB","SB")
colnames(buffer_5000) <- c("BW","CW","AR","IG","NG","CG","AG","FS","HR","HG","BO","IR","SW","FW","SS","LS","SM","UB","SB")

pairs.panels(buffer_500[, which(colSums(buffer_500)>0),], method="spearman", stars = T)
pairs.panels(buffer_1000[, which(colSums(buffer_1000)>0),], method="spearman", stars = T)
pairs.panels(buffer_2500[, which(colSums(buffer_2500)>0),], method="spearman", stars = T)
pairs.panels(buffer_5000[, which(colSums(buffer_5000)>0),], method="spearman", stars = T)

print(corr.test(buffer_500[, which(colSums(buffer_500)>0),], method="spearman", use="pairwise", adjust="BY", alpha=.05), short=F)

##### is a particular land type, a good predictor of a particular plant family?
## what are the families with the largest proportions?
rbcL_top_fams <- cbind(as.data.frame(rbcL_dada2_fam_NC$Family),rowSums(rbcL_dada2_fam_NC[10:23]))
rbcL_top_fams <- rbcL_top_fams[order(-rbcL_top_fams$`rowSums(rbcL_dada2_fam_NC[10:23])`),]
ITS2p_top_fams <- cbind(as.data.frame(ITS2p_dada2_fam_NC$Family),rowSums(ITS2p_dada2_fam_NC[10:23]))
ITS2p_top_fams <- ITS2p_top_fams[order(-ITS2p_top_fams$`rowSums(ITS2p_dada2_fam_NC[10:23])`),]

## rbcL top three families Rosaceae, Balsaminaceae, Oleaceae = account for 71.6% of total rbcL reads
sum(head(rbcL_top_fams$`rowSums(rbcL_dada2_fam_NC[10:23])`, n=3))/sum(rbcL_top_fams$`rowSums(rbcL_dada2_fam_NC[10:23])`)*100

## each of the top three individually
paste(rbcL_top_fams$`rbcL_dada2_fam_NC$Family`[1], rbcL_top_fams$`rowSums(rbcL_dada2_fam_NC[10:23])`[1]/sum(rbcL_top_fams$`rowSums(rbcL_dada2_fam_NC[10:23])`)*100)
paste(rbcL_top_fams$`rbcL_dada2_fam_NC$Family`[2],rbcL_top_fams$`rowSums(rbcL_dada2_fam_NC[10:23])`[2]/sum(rbcL_top_fams$`rowSums(rbcL_dada2_fam_NC[10:23])`)*100)
paste(rbcL_top_fams$`rbcL_dada2_fam_NC$Family`[3],rbcL_top_fams$`rowSums(rbcL_dada2_fam_NC[10:23])`[3]/sum(rbcL_top_fams$`rowSums(rbcL_dada2_fam_NC[10:23])`)*100)

## ITS2p top three families: Balsaminaceae, Rosaceae, Rutaceae =account for 91% of total ITS2p reads
sum(head(ITS2p_top_fams$`rowSums(ITS2p_dada2_fam_NC[10:23])`, n=3))/sum(ITS2p_top_fams$`rowSums(ITS2p_dada2_fam_NC[10:23])`)*100

## each of the top three individually
paste(ITS2p_top_fams$`ITS2p_dada2_fam_NC$Family`[1], ITS2p_top_fams$`rowSums(ITS2p_dada2_fam_NC[10:23])`[1]/sum(ITS2p_top_fams$`rowSums(ITS2p_dada2_fam_NC[10:23])`)*100)
paste(ITS2p_top_fams$`ITS2p_dada2_fam_NC$Family`[2],ITS2p_top_fams$`rowSums(ITS2p_dada2_fam_NC[10:23])`[2]/sum(ITS2p_top_fams$`rowSums(ITS2p_dada2_fam_NC[10:23])`)*100)
paste(ITS2p_top_fams$`ITS2p_dada2_fam_NC$Family`[3],ITS2p_top_fams$`rowSums(ITS2p_dada2_fam_NC[10:23])`[3]/sum(ITS2p_top_fams$`rowSums(ITS2p_dada2_fam_NC[10:23])`)*100)

## print total taxa associated with each family
for (family in unique(Taxa$Family)){
  print(paste(family,nrow(rbcL_dada2_spec[which(rbcL_dada2_spec$Family==family),])+
    nrow(rbcL_dada2_gen[which(rbcL_dada2_gen$Family==family),])+
    nrow(rbcL_dada2_fam[which(rbcL_dada2_fam$Family==family),])))
}

for (family in unique(Taxa$Family)){
  print(paste(family,nrow(ITS2p_dada2_spec[which(ITS2p_dada2_spec$Family==family),])+
                nrow(ITS2p_dada2_gen[which(ITS2p_dada2_gen$Family==family),])+
                nrow(ITS2p_dada2_fam[which(ITS2p_dada2_fam$Family==family),])))
}

### More info on the Rosaceae discovered (interesting because this is the most populous family)
for (species in rbcL_dada2_spec[which(rbcL_dada2_spec$Family=="Rosaceae"),]$Species){
  print(paste(species,Taxa[which(Taxa$Species==species),]$Status))
}
for (species in ITS2p_dada2_spec[which(ITS2p_dada2_spec$Family=="Rosaceae"),]$Species){
  print(paste(species,Taxa[which(Taxa$Species==species),]$Status))
}