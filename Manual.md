


# 1. Introduction #
Allim, Allelic imbalance meter, offers an integrated and user-friendly solution for measuring allele specific gene expression (ASE) within species. Allim estimates allelic imbalance in F1 hybrids. Since mapping bias is the largest problem for reliable estimates of allele specific gene expression using RNA-seq, Allim combines two different measures to account for mapping biases. First, Allim generates a polymorphism aware reference genome that accounts for the sequence variation between the alleles of both parents (or parental lines). Second, Allim includes a sequence specific simulation tool to estimate the remaining mapping bias. This estimated mapping bias is then incorporated in the statistical tests for allelic imbalance.

The pipeline requires whole transcript high throughput RNA sequencing (RNA-seq) reads of F1 hybrids. Additionally, either RNA-seq reads for both parents, genomic sequencing reads for both parents or two private genomes for both parents have to be provided in separate files. Allim was tested on Illumina paired-end RNA-seq reads but it can also handle FASTQ files from other NGS sequencing platforms. The provided parental RNA-seq libraries can be from homozygous as well as heterozygous parents (however, all heterozygous loci will be excluded during the analysis).
Allim has five modules that can be run by a single command. All input parameters can be specified in the AllimOptions file. These parameters are then used to run the complete pipeline. Allim provides two different input options:

**AllimOptions\_2Pexpr:** This configuration file can be used when RNA-seq data for parent1, parent2 and F1 (hybrid) is provided. In this case Allim uses parent1 and parent2 RNA-seq reads to identify the fixed SNPs. The fixed SNPs are subsequently used to generate two “private -genomes”, one for each parent. Note, rather than RNA-seq data, also DNA reads can be provided for the two parents.<br>
With Allim if user does not have reference genome and corresponding gene annotation in GTF format. It is possible to provide a single reference transcriptome along with RNA-seq data of both parents. The use of a reference transcriptome (or contigs from a RNA-Seq de novo assembly) instead of a reference genome is accompanied by the following differences:<br>
<ul><li>In a transcriptome assembly different isoforms are typically presented by different contigs. Furthermore, assemblies often contain additional redundancies between contigs due to sequencing errors and polymorphisms in the RNA-seq reads used for de novo transcriptome assembly.<br>
</li><li>As Allim creates two parental references (genomes or here transcriptomes), and maps RNA-seq reads of the F1 individual to the “diploid genome”, only reads that map non-ambiguously can be used to determine ASE profiles. It is therefore recommended to remove redundancy from the transcriptome assembly before Allim usage.<br>
</li><li>As transcriptome assemblies typically do not have a gtf file containing gene features (needed as Allim input). We have therefore added an additional short script to obtain a simple gtf file based on contig ids.</li></ul>


<b>AllimOptions_2Pgenomes:</b> This configuration file can be used when the user provides RNA-seq data for the F1 hybrid and a genomic sequence for each parental line. If genomic sequences of both parents are provided, both FASTA files need to have the same size and identical FASTA IDs.<br>
<br>
<br>
<br>
<h2>Five Modules of Allim</h2>

<ul><li>Identification of fixed SNPs and creating two parental genomes<br>
</li><li>Computer simulation of RNA-seq reads with fixed SNPs<br>
</li><li>Estimation of the remaining mapping bias with simulated data<br>
</li><li>Estimation of allele specific expression for experimental data<br>
</li><li>Statistical test of significant allelic imbalance</li></ul>

The source code and the user manual of Allim are available from <a href='http://code.google.com/p/allim/'>http://code.google.com/p/allim/</a>


<h1>2. System Requirements</h1>

To use Allim pipeline, user needs to meet these requirements.<br>
<ul><li>Linux or Macintosh OSX system or any other Unix 64 bit system with at least 4 GB of RAM and 2 CPU (processors)<br>
</li><li>Python 2.7.3 <a href='http://www.python.org/download/'>http://www.python.org/download/</a>
</li><li>Biopython 1.59 or higher <a href='http://biopython.org/wiki/Download'>http://biopython.org/wiki/Download</a>
</li><li>PICARD (BuildBamIndex.jar, MergeSamFiles.jar, SortSam.jar) <a href='http://sourceforge.net/projects/picard/files/picard-tools/'>http://sourceforge.net/projects/picard/files/picard-tools/</a>
</li><li>SAMTOOLS <a href='http://sourceforge.net/projects/samtools/files/samtools/'>http://sourceforge.net/projects/samtools/files/samtools/</a>
</li><li>BedTools (bamToBed, intersectBed, sortBed) <a href='http://code.google.com/p/bedtools/'>http://code.google.com/p/bedtools/</a>
</li><li>GSNAP <a href='http://research-pub.gene.com/gmap/'>http://research-pub.gene.com/gmap/</a>
</li><li>R version 2.15.0 <a href='http://cran.r-project.org/'>http://cran.r-project.org/</a>
</li><li>NumPy <a href='http://numpy.scipy.org/'>http://numpy.scipy.org/</a>
</li><li>RPy 2.2.2 <a href='http://rpy.sourceforge.net/rpy2.html'>http://rpy.sourceforge.net/rpy2.html</a>
</li><li>R package car version 2.0-12 <a href='http://cran.r-project.org/web/packages/car/index.html'>http://cran.r-project.org/web/packages/car/index.html</a>
</li><li>R package multcomp version 1.2-12 <a href='http://cran.r-project.org/web/packages/multcomp/index.html'>http://cran.r-project.org/web/packages/multcomp/index.html</a>
</li><li>Bioconductor package edgeR <a href='http://www.bioconductor.org/packages/2.10/bioc/html/edgeR.html'>http://www.bioconductor.org/packages/2.10/bioc/html/edgeR.html</a>
</li><li>Bioconductor package limma <a href='http://www.bioconductor.org/packages/release/bioc/html/limma.html'>http://www.bioconductor.org/packages/release/bioc/html/limma.html</a></li></ul>

<blockquote>Note: To install R packages listed above 11-14, R >=2.15.0 should be installed.<br>
<h2>2.1 Operating system</h2>
The Allim package is designed to work with a 64-bit Unix operating system with at least 4 GB of RAM and 2 CPUs (processors).</blockquote>

<h2>2.2 Python installation</h2>
Allim is developed on Python version 2.7.3. Python 2.7.3 for your operating system can be obtained from <a href='http://www.python.org/download/'>http://www.python.org/download/</a>. After the download follow the instructions given on the web page to complete the python installation and configuration.<br>
<br>
<h2>2.3 Biopython installation</h2>
Allim requires biopython for efficient fasta sequence reading, writing and other sequence manipulations. Once python is installed, biopython can be obtained and installed with the following steps:<br>
<br>
<b>Step1:</b> Download the biopython source code from <a href='http://biopython.org/DIST/biopython-1.60.tar.gz'>http://biopython.org/DIST/biopython-1.60.tar.gz</a>

<b>Step2:</b> Uncompress the downloaded file with the following command:<br>
<br>
<pre><code>tar –zxvf biopython-1.60.tar.gz<br>
</code></pre>

This command will return the folder/directory biopython-1.6.0<br>
<br>
<b>Step3:</b> Enter the uncompressed folder with following command:<br>
<br>
<pre><code>cd biopython-1.60<br>
</code></pre>

<b>Step4:</b> To install the biopython package run these two commands:<br>
<br>
<pre><code>python setup.py build<br>
python setup.py install<br>
</code></pre>


<h2>2.4 NumPy installation</h2>
NumPy is the fundamental package for scientific computing with Python. In the Allim pipeline NumPy is used to integrate sequencing errors during the simulation of RNA-seq reads. In order to install NumPy on any Unix system the following steps are required:<br>
<br>
<b>Step1:</b> Download the NumPy package source code from <a href='http://sourceforge.net/projects/numpy/files/NumPy/1.6.0/'>http://sourceforge.net/projects/numpy/files/NumPy/1.6.0/</a>


<b>Step2:</b> Uncompress the downloaded file with following command:<br>
<br>
<pre><code>tar –zxvf numpy-1.6.0.tar.gz<br>
</code></pre>

This command will return the folder/directory numpy-1.6.0<br>
<br>
<b>Step3:</b> Enter the uncompressed folder with following command:<br>
<br>
<pre><code>cd numpy-1.6.0<br>
</code></pre>


<b>Step4:</b> To install this python package run these two commands:<br>
<br>
<pre><code>python setup.py build<br>
python setup.py install<br>
</code></pre>



<h2>2.5 PICARD installation</h2>
Picard comprises Java-based command-line utilities that manipulate SAM/BAM files. It is a collection of many JAVA jar files, which can be downloaded and used directly without prior installation.<br>
<br>
<br>
<b>Step1:</b> Download the latest version of PICARD from <a href='http://sourceforge.net/projects/picard/files/picard-tools/1.75/'>http://sourceforge.net/projects/picard/files/picard-tools/1.75/</a>

<b>Step2:</b> Uncompress the downloaded file with the following command:<br>
<br>
<pre><code>unzip picard-tools-1.75.zip<br>
</code></pre>

The Allim pipeline uses three jar files: 1) <i>BuildBamIndex.jar</i>, 2) <i>MergeSamFiles.jar</i>, 3) <i>SortSam.jar</i>. Provide the full paths of these three jar files in AllimOptions run file to run the Allim pipeline.<br>
<br>
<h2>2.6 SAMTOOLS installation</h2>
SAMTOOLS provide various utilities for the manipulation of alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.<br>
<br>
<br>
<b>Step1:</b> Download latest SAMTOOLS source code from <a href='http://sourceforge.net/projects/samtools/files/samtools/'>http://sourceforge.net/projects/samtools/files/samtools/</a>


<b>Step2:</b> Uncompress the downloaded file samtools-0.1.18.tar.bz2 with following command:<br>
<br>
<pre><code>tar -xvjf samtools-0.1.18.tar.bz2<br>
</code></pre>

This command will return the folder/directory samtools-0.1.18<br>
<br>
<b>Step3:</b> Enter the uncompressed folder samtools-0.1.18 with following command:<br>
<br>
<pre><code>cd samtools-0.1.18<br>
</code></pre>


<b>Step4:</b> To make samtools executable run the make command:<br>
<br>
<pre><code>make<br>
</code></pre>

After running the make command an executable called samtools will be created. Provide the full path of the samtools executable in the AllimOptions run file to run the Allim pipeline.<br>
<br>
<br>
<br>
<h2>2.7 BedTools installation</h2>
BedTools provide a set of functions to sort and intersect various genomic formats including bam files and files that contain genomic annotation such as gene position and SNP information.<br>
<br>
<b>Step1:</b> Download BedTools from <a href='http://bedtools.googlecode.com/files/BEDTools.v2.16.2.tar.gz'>http://bedtools.googlecode.com/files/BEDTools.v2.16.2.tar.gz</a>

<b>Step2:</b> Uncompress the downloaded file BEDTools.v2.16.2.tar.gz with the following command:<br>
<br>
<pre><code>tar -xzvf BEDTools.v2.16.2.tar.gz<br>
</code></pre>

This command will return the folder/directory BEDTools-Version-2.16.2<br>
<br>
<br>
<b>Step3:</b> Enter the uncompressed folder BEDTools-Version-2.16.2 with following command:<br>
<br>
<pre><code>cd BEDTools-Version-2.16.2<br>
</code></pre>


<b>Step4:</b> To make BedTools executable run the make all command:<br>
<br>
<pre><code>make all<br>
</code></pre>

After running the “make all” command, BEDTools-Version-2.16.2/bin sub-folder/sub-directory will be created, which contains all BedTool executables.<br>
The Allim pipeline uses three executables: 1) <i>intersectBed</i>, 2) <i>sortBed</i>, 3) <i>bamToBed</i>. Provide the full paths of these three executables in the AllimOptions run file to run the Allim pipeline.<br>
<br>
<br>
<br>
<br>
<h2>2.8 GSNAP installation</h2>
GSNAP (Genomic Short-read Nucleotide Alignment Program) is a mapper for RNA-seq data. It has the advantages that it can detect splicing events and is capable of SNP tolerant alignments (Wu & Nacu 2010).<br>
<br>
<b>Step1:</b> Download the GSNAP source code from <a href='http://research-pub.gene.com/gmap/src/gmap-gsnap-2012-07-20.tar.gz'>http://research-pub.gene.com/gmap/src/gmap-gsnap-2012-07-20.tar.gz</a>

<b>Step2:</b> Uncompress the downloaded file “gmap-gsnap-2012-07-20.tar.gz” with following command:<br>
<pre><code>tar -xzvf gmap-gsnap-2012-07-20.tar.gz<br>
</code></pre>


<b>Step3:</b> Enter the uncompressed folder “gmap-2012-07-20” with following command:<br>
<br>
<pre><code>cd gmap-2012-07-20<br>
</code></pre>


<b>Step4:</b> To install the gmap package run these four commands:<br>
<br>
<pre><code>./configure<br>
make<br>
make check   (optional)<br>
make install<br>
</code></pre>

The above four commands will build GSNAP and other executables in “/usr/bin”<br>
<br>
<b>Step5:</b> To check the GSNAP installations run this command:<br>
<br>
<pre><code>gsnap –help<br>
</code></pre>

This command will return the detailed help manual for various GSNAP parameters if the installation was successful.<br>
<br>
<b>Note:</b> For more detailed information how to install GSNAP please read the gmap-2012-07-20/README file.<br>
<br>
<br>
<br>
<h2>2.9 R installation</h2>
The R source code specific to your operating system can be obtained from <a href='http://cran.r-project.org/'>http://cran.r-project.org/</a>. To install R, please follow the instructions provided with R.<br>
<br>
<br>
<h2>2.10 RPy2 installation</h2>
RPy2 is a python class designed to do statistics programming using R in python. In order to install RPy2 on any Unix system the following steps are required:<br>
<br>
<b>Step1:</b> Download the RPy2 package source code from <a href='http://sourceforge.net/projects/rpy/files/rpy2/2.2.x/'>http://sourceforge.net/projects/rpy/files/rpy2/2.2.x/</a>

<b>Step2:</b> Uncompress the downloaded file with following command:<br>
<br>
<pre><code>tar –zxvf rpy2-2.2.2.tar.gz<br>
</code></pre>

This command will return the folder/directory rpy2-2.2.2<br>
<br>
<b>Step3:</b> Enter the uncompressed folder with following command:<br>
<br>
<pre><code>cd rpy2-2.2.2<br>
</code></pre>


<b>Step4:</b> To install the RPy2 package run the following two commands:<br>
<br>
<pre><code>python setup.py build<br>
python setup.py install<br>
</code></pre>



<h2>2.11 Installation of the R package <i>car</i></h2>

<b>Step1:</b> Download the R package car source code from <a href='http://cran.r-project.org/src/contrib/car_2.0-12.tar.gz'>http://cran.r-project.org/src/contrib/car_2.0-12.tar.gz</a>

<b>Step2:</b> Open the terminal and run the following command:<br>
<br>
<pre><code>R CMD INSTALL car_2.0-12.tar.gz<br>
</code></pre>



<h2>2.12 Installation of the R package <i>multcomp</i></h2>

<b>Step1:</b> Download R package multcomp source code from <a href='http://cran.r-project.org/src/contrib/multcomp_1.2-12.tar.gz'>http://cran.r-project.org/src/contrib/multcomp_1.2-12.tar.gz</a>

<b>Step2:</b> Open the terminal and run the following command:<br>
<br>
<pre><code>R CMD INSTALL multcomp_1.2-12.tar.gz<br>
</code></pre>



<h2>2.13 Installation of the R package <i>edgeR</i></h2>

<b>Step1:</b> Download R package edgeR source code from www.bioconductor.org/packages/2.3/bioc/src/contrib/edgeR_1.0.4.tar.gz<br>
<br>
<b>Step2:</b> Open the terminal and run following command:<br>
<br>
<pre><code>R CMD INSTALL edgeR_1.0.4.tar.gz<br>
</code></pre>



<h2>2.14 Installation of the R package <i>limma</i></h2>

<b>Step1:</b> Download R package limma source code from <a href='http://phase.hpcc.jp/mirrors/stat/R/CRAN/src/contrib/limma_2.0.2.tar.gz'>http://phase.hpcc.jp/mirrors/stat/R/CRAN/src/contrib/limma_2.0.2.tar.gz</a>

<b>Step2:</b> Open the terminal and run the following command:<br>
<br>
<pre><code>R CMD INSTALL limma_2.0.2.tar.gz<br>
</code></pre>

<h1>3. Allim Installation</h1>
The user can download the latest version of Allim from <a href='http://code.google.com/p/allim/'>http://code.google.com/p/allim/</a>. The file to download is called Allim_1.0.tar.gz. Move the file to an appropriate directory and run the following command to uncompress the file:<br>
<pre><code>tar –zxvf Allim_1.0.tar.gz<br>
</code></pre>

Note that after uncompressing the tar.gz file, a new folder will be created named Allim_1.0. This directory contains the following files:<br>
<br>
<br>
<img src='http://allim.googlecode.com/files/Allim_folder.png' />

<h1>4. The Validation of the Installation & Sample Input Files</h1>
To validate the installation of the Allim pipeline it can be run with a small test data set. The test data set and the corresponding Allim configuration files can be obtained from the following URLs:<br>
<br>
<b>Parent1 RNA-seq:</b> <a href='http://allim.googlecode.com/files/parent1_RNAseq_fastq.tar.gz'>http://allim.googlecode.com/files/parent1_RNAseq_fastq.tar.gz</a>

<b>Parent2 RNA-seq:</b> <a href='http://allim.googlecode.com/files/parent2_RNAseq_fastq.tar.gz'>http://allim.googlecode.com/files/parent2_RNAseq_fastq.tar.gz</a>

<b>F1 RNA-seq:</b> <a href='http://allim.googlecode.com/files/F1_RNAseq_fastq.tar.gz'>http://allim.googlecode.com/files/F1_RNAseq_fastq.tar.gz</a>

<b>Parent1 genome:</b> <a href='http://allim.googlecode.com/files/parent1_genome.fa'>http://allim.googlecode.com/files/parent1_genome.fa</a>

<b>Parent2 genome:</b> <a href='http://allim.googlecode.com/files/parent2_genome.fa'>http://allim.googlecode.com/files/parent2_genome.fa</a>

<b>Reference fasta file:</b> <a href='http://allim.googlecode.com/files/reference.fa'>http://allim.googlecode.com/files/reference.fa</a>

<b>Gene annotation file:</b> <a href='http://allim.googlecode.com/files/reference.gtf'>http://allim.googlecode.com/files/reference.gtf</a>

<b>Allim configuration file:</b> <a href='http://allim.googlecode.com/files/AllimOptions_2Pexpr'>http://allim.googlecode.com/files/AllimOptions_2Pexpr</a>

<b>Allim configuration file:</b> <a href='http://allim.googlecode.com/files/AllimOptions_2Pgenomes'>http://allim.googlecode.com/files/AllimOptions_2Pgenomes</a>



<b>Note:</b> You either need to use the parental RNA-seq or the parental genomic data. The option you choose depends on the type of input data you have for your own analysis.<br>
<br>
Create a folder named <b>“test_data”</b> in the Allim_1.0 directory and download all the above test data files into this folder, unzip and extract the compressed files. Move the <b>AllimOptions_2Pexpr</b>  and <b>AllimOptions_2Pgenomes</b> files into the directory <b>Allim_1.0</b> directly. The Allim pipeline can then be run on the test data set with the following two steps:<br>
<br>
<br>
<b>1.</b> Open the terminal and enter the Allim_1.0 directory.<br>
cd Allim_1.0<br>
<br>
<b>2.</b> Run the Allim pipeline with the test data set via the following command (the type of the AllimOptions file that is used depends on the type of your input data):<br>
<pre><code>python allim.py --option-file AllimOptions_2Pexpr<br>
OR<br>
python allim.py --option-file AllimOptions_2Pgenomes<br>
</code></pre>

To get help on how to run Allim and required parameters enter:<br>
<pre><code>	python allim.py –help<br>
</code></pre>



<h1>5. Allim Input Description</h1>
Allim can be run with the following command, which should be run under the Allim_1.0 directory:<br>
<pre><code>	python allim.py --option-file &lt;path to AllimOptions file&gt;<br>
</code></pre>

However, before running Allim with your own dataset all parameters have to be specified in the appropriate Allim Options files (AllimOptions_2Pexpr OR AllimOptions_2Pgenomes).<br>
<br>
<h2>Global input parameters</h2>

Figure 1 shows global input parameters, which are essential for all modules of Allim. Figure 2 shows how the paths to third party tools that are used by Allim can be specified in the AllimOptions file.<br>
<br>
<br>
<img src='http://allim.googlecode.com/files/global-parameters.png' />

<b>Figure 1: Global input parameters of the Allim pipeline.</b> The figure shows an extract from the AllimOptions_2Pexpr file. In AlimOptions_2Pgenomes, the parameter REFERENCEE_FASTA is replaced by the two parameters: PARENT1_REFERENCEE_FASTA AND PARENT2_REFERENCEE_FASTA.<br>
<br>
<br>
<br>
<img src='http://allim.googlecode.com/files/third-party-softwares.png' />

<b>Figure 2: Specifications of full paths of external software/tools used in Allim.</b>



The remaining input parameters are described in the following chapter along with the Allim modules description.<br>
<br>
<br>
<h1>6. Allim Modules Description</h1>
<h2>Module 1: Identification of fixed SNPs</h2>
This module is only used when the input option, “AllimOptions_2Pexpr” is chosen (RNA-seq or genomic DNA reads of both parents are provided). With the input option, “AllimOptions_2Pgenomes” (two genomes, one for each parent) is used this module will be skipped.<br>
<br>
The module determines fixed SNPs between the parental genotypes (lines) based on the user provided RNA-seq libraries. It accepts multiple replicates for the calculation of fixed SNPs to increase the power and accuracy of SNP detection. For each parental genotype and replicate RNA-seq data (alternatively genomic DNA reads) can be provided for paired-end sequence data (pairs of fastq files, read1.fq and read2.fq, for each condition). Besides the RNA-seq data (genomic reads) further parameters can be specified as shown in Figure 3. The type of data that is used for calling fixed SNPs, either RNA or genomic sequence is reads is specified via the parameter SEQUENCE_TYPE (“mRNA” or “DNA”).<br>
<br>
In Allim the identification of fixed SNPs is based on the alignments of the RNA-seq reads (alternatively genomic DNA reads) performed with the GSNAP mapper (Wu & Nacu 2010). GSNAP is a unique mapper that can integrate given SNP information into the mapping algorithm in order to improve read mapping of allelic variants. Allim makes use of this functionality by integrating the fixed SNP information into the mapping process to improve the quality of the alignments. However, in the beginning no fixed SNP information is available. Therefore, multiple cycles of read mapping and subsequent calling of fixed SNPs are required to improve the quality of the alignment and the identified fixed SNPs. The user can specify how often this procedure should be repeated (minimum & default: 2 cycles). The number of cycles given as input provides the opportunity to fine-tune SNP calling accuracy, as this is dependent on the accuracy of the alignment. As output a fixed SNP table for each cycle is saved in the given output directory. The fixed SNP information determined in the last cycle is used in subsequent modules.<br>
<br>
The parameters MINIMUM_COVERAGE, FIXED_ALLELE_FREQUENCY, MINIMUM_MAPPING_QUALITY and MINIMUM_BASE_QUALITY are used in each cycle to determine fixed SNPs between both parental genotypes (Figure 3).<br>
<br>
<br>
<img src='http://allim.googlecode.com/files/fixed-snp.jpg' />


<b>Figure 3: Input options for the identification of fixed SNPs between the two parental genotypes.</b> The figure shows an extract from the AllimOptions_2Pexpr file.<br>
<br>
<br>
<br>
<h2>Module 2: Computer simulation of RNA-seq reads</h2>
This module simulates paired-end Illumina reads in fastq format. It internally uses the fixed SNP information to create two “private genomes”, one for each parent or alternatively the user provided reference genomes. Further it used the information from the provided GTF file (provided as a global input parameter).<br>
<br>
For the read simulation, first two parent specific genomes (“private genomes”), which only differ with respect to the identified fixed SNPs, are generated. For each private genome all possible paired-end reads that cover at least one fixed SNP position are simulated once. Therefore, the expression ratio been reads simulated for each parent should equal 1 for each gene. In contrast, the deviation from an expression ratio of 1 in the simulated data for a gene indicates a remaining mapping bias, which is caused by the mapper itself.<br>
<br>
<b>Algorithm:</b>

<b>1.</b>        Creation of two parent specific (“private”) genomes via the substitution of the base at a fixed SNP position in initial reference genome (this step is skipped if two parental genomes are provided initially).<br>
<br>
<b>2.</b>	Construction of the longest transcript for each gene for both private genomes separately. Due to step 1, the transcripts for both private genomes only differ at the fixed SNP positions (AllimOptions_2Pexpr). If two parental genomes are provided initially the both genomes must have identical genome length and identical gene annotation. Each transcript consists of the 5’UTR, the CDS and the 3’UTR.<br>
<br>
<b>3.</b>	Simulation of the reads for both “transcriptomes” (as defined in step 2) separately. This results in the same number of reads from the identical genomic locations for both parents. The simulation is implemented via a sliding window approach where the “simulation window” slides one base pair at a time from 5’ towards 3’ of transcript. Only the reads that span a fixed SNP are kept as the remaining ones are not informative for the subsequent analysis (they cannot be used to determine ASE). This approach results in equal coverage of SNP position in the middle of a transcript and decreasing coverage of SNPs towards both ends of the transcript (end of a transcript: 2 x read length + insert size).<br>
<br>
Additional parameters that define the properties of the simulated reads are shown in Figure 4.<br>
<br>
<br>
<img src='http://allim.googlecode.com/files/simulation.png' />

<b>Figure 4: Input options for the computer simulation of RNA-seq reads from two parent genotypes.</b> The figure shows an extract from the AllimOptions file.<br>
<br>
<br>
<br>
<h2>Modules 3 & 4: Estimation of the remaining mapping bias with simulated data & Estimation of allele specific expression for experimental data</h2>

Module 3 is designed to estimate the remaining mapping bias for each gene (exon) via the simulated RNA-seq reads. In the simulated data the number of reads for each transcript is identical for both parental genotypes (parental lines). Resulting expression ratios between counts for both genotypes (parental lines) should therefore be 1 for each gene (exon) in the absence of a remaining mapping bias due to the mapper. Ratios deviating from 1 indicate a remaining mapping bias and are used as a correction factor to determine statistical significance of allelic imbalance in module 5.<br>
<br>
Module 4 is designed to estimate the total gene expression in both parental lines and allele specific gene expression in the F1 generation. The module generates an expression table for all genes and all exons in two separate files. This expression table will be subsequently used in module 5 to test for significant allelic imbalance.<br>
<br>
Modules 3 & 4 make use of a two-genome approach to determine allele specific expression. In this approach the two “private” genomes (one for each parental genotype / parental line) that where created for the RNA-seq read simulation in module 2 (AllimOptions_2Pexpr) or provided by the user (AllimOptions_2Pgenomes) are used as a “combined reference” (genome of parent 1 plus genome of parent 2) to map RNA-seq reads. Only the reads that map unambiguously (no multiple equally good mapping positions) are used to determine allele specific expression. This approach excludes reads that do not span fixed SNPs as these have at least two equally good mapping locations in the “combined reference”.<br>
<br>
The input files and other parameters for modules 3& 4 can be specified in the AllimOptions file as shown in figure 5.<br>
<br>
<br>
<img src='http://allim.googlecode.com/files/ASE.png' />

<b>Figure 5: Input options for estimating total allele specific expression.</b> The figure shows an extract from the AllimOptions_2Pexpr file. In AlimOptions_2Pgenomes only fastq files for the F1 hybrid have to be provided.<br>
<br>
<br>
<br>
<h2>Module 5: Statistical test of significant allelic imbalance</h2>
G-tests of allelic imbalance for each gene/exon are provided for samples without biological replication. If biological replicates are available, Allim provides an alternative approach to determine allelic imbalance across replicates. In an ANOVA framework it is tested whether differences are present in the expression strength between both alleles in the F1 generation. Module 5 returns p-values for each gene (exon) including multiple testing correction via the false discovery rate (FDR).<br>
<br>
For both approaches (G-test, ANOVA) three different correction factors are implemented in Allim, namely 1) Library size normalization, 2) Mapping bias correction and 3) Rescaling of the corrected read count. These different normalization steps are explained in detail below together with examples (tables 1-4).<br>
<br>
<br>
<b>1) Library size normalization</b>

Libraries are normalized by the TMM factor (trimmed mean of M-values normalization method; Robinson and Oshlack 2010). The factor is calculated using the bioconductor package edgeR (Robinson et al. 2010) for each parental library and the F1 library for each replicate. Individual gene counts are then corrected by the TMM factor of the respective library. Note that for each hybrid library only one TMM factor has to be calculated. This factor is used for the expression profile of both parental alleles. Table 1 provides read counts before and table 2 after the first normalization step.<br>
<br>
<img src='http://allim.googlecode.com/files/Table1.png' />



<img src='http://allim.googlecode.com/files/Table2.png' />


<b>2) Mapping bias correction</b>

The residual mapping bias after accounting for fixed SNPs between the two parental lines is accounted for via the mapping bias coefficient. This coefficient is defined as the expression ratio between both parental alleles obtained from the simulated data (module 2).<br>
<br>
<pre><code>Mapping bias coefficient = readcount_simuldata(parent1) / readcount_simuldata(parent2)<br>
</code></pre>

The mapping bias coefficient is calculated for each gene. All read counts normalized by library size (table 2) are then multiplied by the mapping bias coefficient of the respective gene (table 3).<br>
<br>
<br>
<img src='http://allim.googlecode.com/files/Table3.png' />



<b>3) Rescaling of the corrected read count</b>

After applying the above two normalization factors, read counts for the different genes across libraries can be inflated or decreased depending on the values of the respective factors used. As inflation results in an unjustified increase of power in the statistical test and a lowered count in an unjustified decrease, the values have to be rescaled on a gene-wise level. This is done via the following rescaling factor.<br>
<br>
<pre><code>Rescaling factor = expression value before normalization (sum over all samples) /after normalization (sum over all samples)<br>
</code></pre>


All corrected read counts by the two previous correction factors (table 3) are multiplied by the rescaling factor of the respective gene (table 4).<br>
<br>
<br>
<img src='http://allim.googlecode.com/files/Table4.png' />



<h1>7. Allim Output / Results Files</h1>
After running Allim, the results of the Allim pipeline can be found in the output directory specified in the Allim Options file. Results are provided in the following format.<br>
<br>
<img src='http://allim.googlecode.com/files/Output.png' />


<b>Note:</b>

<b>‘</b>1’<b>this directory is only provided if the input option “AllimOptions_2Pexpr” was chosen.<br></b>‘<b>2’</b> indicates the final output files in above structure; the remaining files are intermediate output files constructed by the pipeline.<br>
<br>
<br>
<br>
<br>
<01_fixed_SNP><br>
<br>
:<b><br>
This folder contains the identified fixed SNPs between the parental genotypes (parental lines). The number at the beginning of the file indicates in which cycle of remapping the fixed SNPs were called. “n” is a parameter that has to be provided by the user.</b>

<br>
<br>
<02_simulation_data_ase><br>
<br>
:<b><br>
This folder contains the expression tables with the raw read counts of the simulated data for both parental genotypes (lines). Profiles are obtained for genes as well as for exons as units for which gene expression is measured.</b>


<br>
<br>
<03_experimental_data_ase><br>
<br>
:<b><br>
This folder contains the expression tables with the raw read counts of the experimental data for all replicates of both parental genotypes (lines) and the F1 generation. Profiles are obtained for genes as well as for exons as units of expression.</b>

<b><04_statistical test>:</b>

This folder contains the read counts of the experimental data, which were normalized and corrected for the residual mapping bias. Profiles for all replicates of both parental genotypes (lines) and the F1 generation are obtained for genes as well as for exons as units of expression.<br>
<br>
The remaining two files provide the results of the significance testing for allelic imbalance for the F1 generation.<br>
<br>
<br>
<h1>8. Benchmark of the Allim Pipeline</h1>

We have run the Allim pipeline with the following data set of RNA-seq reads from <i>Drosophila pseudoobscura</i> to benchmark the runtime of Allim. The reference genome release 2.23 was taken from Flybase and an improved annotation file of <i>D. pseudoobscura</i> was taken from (<i>Palmieri et al 2012</i>). The GTF file with the gene annotation contained 17,112 gene models.<br>
<br>
<br>
<b>RNA-seq reads (Illumina, GA II):</b>

Parent1 reads: 2 million 100 bp paired-end reads with 68 bp average insert size.<br>
<br>
Parent2 reads: 2 million 100 bp paired-end reads with 78 bp average insert size.<br>
<br>
F1 reads: 2 million 100 bp paired-end reads with 98 bp average insert size<br>
<br>
<br>
<br>
<b>The Allim run simulated the following amount of sequence reads:</b>

Parent1 reads: 10.26 million 100 bp paired-end reads with 78 bp insert size<br>
<br>
Parent2 reads: 10.26 million 100 bp paired-end reads with 78 bp insert size<br>
<br>
<br>
<b>Table 5: Benchmarks for processing time for each module of Allim</b>

<table><thead><th><b>S.No.</b></th><th><b>Module</b></th><th><b>Time (min:sec)</b></th></thead><tbody>
<tr><td>1.</td><td>Finding fixed SNP</td><td>95:00</td></tr>
<tr><td>2.</td><td>Computer simulation of RNA-Seq reads</td><td>20:00</td></tr>
<tr><td>3.</td><td>Estimating Allele Specific Expression with experimental data</td><td>24:00</td></tr>
<tr><td>4.</td><td>Estimating Allele Specific Expression with simulated data</td><td>34:00</td></tr>
<tr><td>5 </td><td>Statistical test of significant allelic imbalance</td><td>3:00</td></tr></tbody></table>


<b>Note:</b> All benchmarks have been done on a Mac OS X 10.6.8, 2x2.8 GHz Quad-Core Intel Xeon, with 4 GB of RAM using 4 processors (CPU)<br>
<br>
<br>
<br>
<br>
<h1>9. Validation</h1>
We validated the performance of our pipeline with experimental as well as simulated RNA-seq reads. The experimental data consisted of paired-end RNA-seq reads from males and females of two different isofemale lines of <i>Drosophila pseudoobscura</i> (Table 6; <i>Palmieri et al 2012</i>).<br>
<br>
<br>
<br>
<br>
<b>Table 6: Number of paired-end RNA-seq reads of D. pseudoobscura, which were used for validation. (ps88 and ps94 stand for two different isofemale lines.)</b>

<table><thead><th>Samples</th><th>Reads (in million)</th><th>Insert Size (bp)</th><th>Read length (bp)</th></thead><tbody>
<tr><td>ps88 males</td><td>79.21</td><td>78</td><td>100</td></tr>
<tr><td>ps88 females</td><td>80.00</td><td>78</td><td>100</td></tr>
<tr><td>ps94 males</td><td>79.21</td><td>128</td><td>100</td></tr>
<tr><td>ps94 females</td><td>80.00</td><td>68</td><td>100</td></tr></tbody></table>

The four libraries were given as input to Allim to call fixed differences between both isofemale lines (different sexes were used as biological replicates). Module 1 of the pipeline identified154,920 fixed SNPs present in a total of 9,626 genes. The following Allim parameters were used to call the fixed SNPs:<br>
CALCULATE_FIXED_SNP_ITERATION = 2<br>
MINIMUM_COVERAGE = 2<br>
FIXED_ALLELE_FREQURNCY = 1.0<br>
MINIMUM_BASE_QUALITY = 20<br>
MINIMUM_MAPPING_QUALITY = 20<br>
The reference genome release 2.23 was taken from Flybase and an improved annotation file of <i>D. pseudoobscura</i> was taken from (<i>Palmieri et al 2012</i>). The GTF file with the gene annotation contained 17,112 gene models.<br>
<br>
<br>
<br>
<h2>Assignment of reads to the correct parental line:</h2>

To test the accuracy of Allim to identify the parental origin of a read, ASE expression profiles were determined for experimental and simulated data:<br>
<br>
<b>1. Experimental data:</b>

a) pooled RNA-seq reads from libraries ps88 males and ps94 males (Table 6)<br>
<br>
b) pooled RNA-seq reads from libraries ps88 females and ps94 females (Table 6)<br>
<br>
<br>
In contrast to RNA-seq data of individuals from the F1 generation, the parental origin of the pooled reads from both lines is known. Therefore, the percentage of reads that were assigned to the correct parental line by Allim could be determined. The results show that between 98.60% and 99.97% of the reads were assigned to the correct parental line by Allim (Table 7). The marginal difference between the success rates of both lines is due to the fact that ps94 is a derivative of the strain for which the reference genome sequence was provided to Allim. For this reason reads originating from ps94 that span polymorphic sites, which did not meet the threshold for calling fixed SNPs, are more easily mapped than reads originating from ps88. Since we used the most extreme case for a reference bias, the small bias observed suggests that under less extreme settings, Allim will provide almost no bias caused by the reference genome.<br>
<br>
<b>2. Simulated data:</b>

RNA-seq reads were simulated with the procedure described in module 2 based on the 154,920 fixed SNPs determined for the experimental data previously. (Approach: For both parental genomes the same number of RNA-seq reads were simulated for the identical genomic positions.) For the simulation the following Allim parameters were used:<br>
ASCII_QUALITY_LETTER = e<br>
READ_LENGTH = 100<br>
INSERT_SIZE = 78<br>
<br>
The results show that 99.99% of the simulated reads were assigned to the correct parental line by Allim (Table 7). As the reads were directly simulated based on the reference genome with only fixed differences edited, this remaining error occurs during the mapping with GSNAP (possibly for reads spanning splice junctions).<br>
<br>
<br>
<br>
<b>Table 7: RNA-seq data sets used to test accuracy of Allim to identify the parental origin of a read.</b>

<table><thead><th>Dataset</th><th># Correctly identified reads, ps88 (%)</th><th># Correctly identified reads, ps94 (%)</th></thead><tbody>
<tr><td>Pooled reads from females of both lines</td><td>99.97</td><td>98.96</td></tr>
<tr><td>Pooled reads from males of both lines</td><td>99.96</td><td>98.60</td></tr>
<tr><td>Simulated reads for both parental genomes</td><td>99.99</td><td>99.99</td></tr></tbody></table>



<h2>Improvement of mapping quality via modification of the reference genome</h2>

We validated the improvement of the mapping quality via reference modification by mapping RNA-seq reads (experimental & simulated data) to the original as well as the modified genomes. Mapping success was then compared between both approaches. The mapping success varies between 85.97% and 92.15% between the different data sets. However, in all data sets an increase in mapping success after editing of the genomes between 0.07% and 0.32% can be detected (Table 8).<br>
<br>
<br>
<img src='http://allim.googlecode.com/files/Table8.png' />


We validated the performance of our pipeline with experimental as well as simulated RNA-seq reads. The experimental data consisted of paired-end RNA-seq reads from males and females of two different isofemale lines of <i>Drosophila pseudoobscura</i> (Table 6, <i>Palmieri et al 2012</i>).<br>
<br>
<h1>10. Contact Information</h1>

<b>Prof. Dr. Christian Schlötterer</b>

christian.schloetterer@vetmeduni.ac.at<br>
<br>
<b>Ram Vinay Pandey</b>

ramvinay.pandey@vetmeduni.ac.at<br>
<br>
<b>Susanne U. Franssen</b>

susanne.franssen@vetmeduni.ac.at<br>
<br>
<br>
<h1>11. References</h1>

1.  Wu TD, Nacu S (2010) Fast and SNP-tolerant detection of complex variants and splicing in short reads. Bioinformatics 26: 873-881.<br>
<br>
<br>
2. Robinson MD, Oshlack A (2010) A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biol. 11(3):<a href='https://code.google.com/p/allim/source/detail?r=25'>R25</a>

3. Robinson MD, McCarthy DJ, Smyth GK (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26(1):139-140.<br>
<br>
4. Palmieri N, Nolte V, Suvorov A, Kosiol C, Schlötterer C (2012) Evaluation of Different Reference Based Annotation Strategies Using RNA-Seq – A Case Study in Drososphila pseudoobscura. PLoS One, 7, e46415.