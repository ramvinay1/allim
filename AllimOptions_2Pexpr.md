

# Introduction #
Allim pipeline can be run one of multiple replicates RNA-seq data. In order to use Allim user has to modify the **AllimOptions\_2Pexpr** configuration file. This file allows user to define all input files and other parameters once and whole pipeline will be running without more interactions. This option is especially useful when user have RNA-seq data for parent1, parent2 and F1 (hybrid).

This configuration file can be downloaded from [AllimOptions\_2Pexpr](http://allim.googlecode.com/files/AllimOptions_2Pexpr)


```
##############################################################################################################
## Allim Options file (AllimOptions_2Pexpr.txt)
##
## Note: This configuration file is used when user has expression data for parent1, parent2 and F1.
##
## This file sets the customizable options for running Allim. The options include the paths to the
## sequencing data, paths to executables and other parameters of the pipeline e.g. mapping parameters.
## This file contains example parameter settings, which can be used to run the test data set
## provided with Allim. 
## For running Allim with your own data, the example parameter values have to be modified.
## Please note that only the input value but not the parameter name itself should be changed.
##
## Note: 1) Lines starting with # are comment lines. They do not need to be changed.
##       2) Parameters in upper case letters have to be adjusted. The parameters can be changed
##          by replacing the values after the “=” sign.
##
##############################################################################################################



# Output directory: the name of the directory can be changed.
# The directory will be generated in <Allim_1.0> 
OUTPUT_DIRECTORY = /Allim_1.0/test_output

## Give the number of replicates your data consists of. 
REPLICATE_COUNT=1

# Give the path to the fasta file of the genomic reference. 
REFERENCEE_FASTA = /Allim_1.0/test_data/reference.fa

# Give the gene annotation file of your genomic reference in GTF format.
REFERENCE_GTF_FILE = /Allim_1.0/test_data/reference.gtf

# Give the number of processors that should be used for read mapping.
THREAD = 15

# Give the minimum base quality required to call a fixed SNP and to use a nucleotide position for
# the assessment of allele specific expression profiles. The base quality ranges between 0-40.
# The sequencer for each nucleotide assigns it in read sequence.
MINIMUM_BASE_QUALITY = 20

#
# GSNAP or any other short read mappers (BWA, Bowtie, TopHat)
# assign a mapping quality to each mapped read. It ranges between 0-40.
# Give the minimum mapping quality of a read required to be included in the analysis. 
MINIMUM_MAPPING_QUALITY = 20

#
# Specify the encoding scheme for the base quality of the sequence reads.
# This can either be “sanger” or “illumina”. This parameter value is case sensitive. 
# Note:
# “sanger” (Illumina 1.3+ Phred+33, raw reads typically (0, 40))
# “illumina” (Illumina 1.5+ Phred+64, raw reads typically (3, 40))
QUALITY_ENCODING= illumina



##############################################################################################################
#
# Fixed SNP OPTIONS
# Step1: Identification of fixed SNPs (fixed differences between parental strains)
#        The fixed SNPs will later be used to estimate allele specific gene expression
#        in simulated and experimental data. 
#
##############################################################################################################

#### Give the sequencing type. Allim can calculate fixed SNP either with RNA-Seq FASTQ files OR
#### with DNA sequencing FASTQ files from parent1 and parent2. For Fixed SNP identification
#### user can specify the Sequence type: SEQUENCE_TYPE=mRNA or SEQUENCE_TYPE=DNA
#### It is case sensitive to give exact word mRNA or DNA
SEQUENCE_TYPE=mRNA


## Provide the full paths (relative to <Allim_1.0>) of the paired-end sequencing files
## of the parents and all present replicates.
#
## Paired-end sequencing files and insert size (fragment size - 2(read length)) should be
## provided in following format:
# read1.fq,read2.fq,78
# The read1.fq file contains all “read 1” reads of the paired end reads. read2.fq contains
# all corresponding “read 2” reads.
# 78 is the insert. Make sure that the order of the paired end reads is identical in both fastq files.

# The insert size is an optional parameter. Fastq files for both paired-end reads are mandatory.

# For each replicate provide 3 files 1) parent1; 2) parent2; and 3) hybrid (F1)
#### Replicate 1: parent1,parent2 and hybrid fastq file data
PARENT1_FASTQ_FILE = /Allim_1.0/test_data/parent1_read1.fastq,/Allim_1.0/test_data/parent1_read2.fastq,68
PARENT2_FASTQ_FILE = /Allim_1.0/test_data/parent2_read1.fastq,/Allim_1.0/test_data/parent2_read2.fastq,78


#### Replicate 2: parent1,parent2 and hybrid fastq file data
#PARENT1_FASTQ_FILE = 
#PARENT2_FASTQ_FILE = 

#### Replicate 2: parent1,parent2 and hybrid fastq file data
#PARENT1_FASTQ_FILE = 
#PARENT2_FASTQ_FILE = 


# Provide parameters for the mapping with GSNAP and the identification of fixed SNPs.
#
# Fixed SNPsbetween both parents are called based on the mapping of the reads from both parents.
# However, the mapping with GSNAP can be improved by integrating the fixed SNP information
# to remap the sequence reads of the parents.

# The following parameter allows to specify the number of iterations that is done to call the final fixed SNPs
# for the subsequent analysis. (The minimum number is 2.) 
CALCULATE_FIXED_SNP_ITERATION = 2


# Give the minimum coverage required to call a fixed SNP.
# This means that at least X mapped reads from each parent have to map to the respective nucleotide
# position in the reference (X = minimum coverage).
MINIMUM_COVERAGE = 2

# Give the minimum frequency of a nucleotide at a genomic position in one parent that is-
# required to call a fixed SNP. This frequency to identify the major allele in each parent-
# has to be between 0.51-1.0.
# Note that the choice of this frequency has an effect on the number of reads that can be-
# used to determine allele specific expression (ASE) (power) AND the accuracy to determine ASE.
# higher frequency: lower throughput, higher accuracy
# lower frequency: higher throughput, lower accuracy
FIXED_ALLELE_FREQURNCY = 1.0








##############################################################################################################
#
# Simulation OPTIONS
# Step2: Simulating RNA-Seq paired-end reads for the whole transcriptome. 
#        This simulation data will be used to estimate the residual mapping bias after -
#        the generation of a polymorphism aware reference genome. 
#
##############################################################################################################

# Give the ASCII letter, which will be assigned to each nucleotide of simulated read in FASTQ file.
# Example: ASCII_QUALITY_LETTER=e for illumina encoding; or ASCII_QUALITY_LETTER=H for sanger encoding.
# For illumina encoding the ASCII QUALITY letters UVWXYZ[\]^_`abcdefgh encode for base qualities from 20-40.
# For sanger encoding the ASCII QUALITY letters 6789:;<=>?@ABCDEFGHI encode for base qualities from 20-40.
# Give a single ASCII letter, which will be assigned for each nucleotide in the simulated read.
ASCII_QUALITY_LETTER = e

# Provide the read length of the simulated reads. For example: READ_LENGTH=100 to simulate 100 bp reads.
READ_LENGTH = 100


# Provide the insert size between read1 and read2 for the paired-end read simulation.
# (fragment size = insert size + 2 * read length)
INSERT_SIZE = 78










##############################################################################################################
#
# ASE OPTIONS
# Step3: Measuring Allele-specific expression (ASE) with simulated as well as experimental data 
#
##############################################################################################################

## Provide the paths (relative to <Allim_1.0>) of the paired-end sequencing files-
#  of the two parents (parent1, parent2) and the hybrid (F1) for all present replicates.
#
## Paired-end sequencing files and insert size (fragment size - 2(read length)) should-
## be provided in following format:
# read1.fq,read2.fq,78
# The read1.fq file contains all “read 1” of the paired end reads. read2.fq contains-
# all corresponding “read 2”. 78 is the insert. Make sure that the order of the paired-
# end reads is identical in both fastq files.
# The insert size is an optional parameter. Fastq files for both paired-end reads are mandatory.
#
# For each replicate provide 3 RNA-seq files 1) parent1; 2) parent2; and 3) hybrid (F1)
# In case no RNA-seq data for the parents is available,
# PARENT1_FASTQ_FILE_EXPR and PARENT2_FASTQ_FILE_EXPR should be empty.

#### Replicate 1: parent1,parent2 and hybrid fastq file data
PARENT1_FASTQ_FILE_EXPR = /Allim_1.0/test_data/parent1_read1.fastq,/Allim_1.0/test_data/parent1_read2.fastq,68
PARENT2_FASTQ_FILE_EXPR = /Allim_1.0/test_data/parent2_read1.fastq,/Allim_1.0/test_data/parent2_read2.fastq,78
HYBRID_FASTQ_FILE_EXPR = /Allim_1.0/test_data/hybrid_read1.fastq,/Allim_1.0/test_data/hybrid_read2.fastq,128
#
#### Replicate 2: parent1,parent2 and hybrid fastq file data
#PARENT1_FASTQ_FILE_EXPR = 
#PARENT2_FASTQ_FILE_EXPR = 
#HYBRID_FASTQ_FILE_EXPR = 

#
#### Replicate 3: parent1,parent2 and hybrid fastq file data
#PARENT1_FASTQ_FILE_EXPR = 
#PARENT2_FASTQ_FILE_EXPR = 
#HYBRID_FASTQ_FILE_EXPR = 

##############################################################################################################
#
# Third party software/tool executables
# Provide the full path to the third party executables to run the Allim pipeline.
#  
#
##############################################################################################################


### Get samtools to manipulate SAM and BAM files from http://samtools.sourceforge.net/
# Provide the full path of the executables relative to the <Allim_1.0>
SAMTOOLS = /Allim_1.0/executables/samtools


### Get intersectBed executables from a collection of useful utilities called bedtools-
### from http://code.google.com/p/bedtools/
# Provide the full paths of the executables relative to the <Allim_1.0>
INTERSECTBED = /Allim_1.0/executables/intersectBed
SORTBED = /Allim_1.0/executables/sortBed
BAMTOBED = /Allim_1.0/executables/bamToBed


#### Three PICARD JAR files from http://picard.sourceforge.net/index.shtml
# Provide the full paths of the executables relative to the <Allim_1.0>
SORTSAM = /Allim_1.0/executables/SortSam.jar
MERGESAMFILES = /Allim_1.0/executables/MergeSamFiles.jar
BUILDBAMINDEX = /Allim_1.0/executables/BuildBamIndex.jar

```