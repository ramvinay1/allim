

# Introduction #
Allim, Allelic imbalance meter, offers an integrated and user-friendly solution for measuring allele specific gene expression (ASE) within species. Allim estimates allelic imbalance in F1 hybrids. Since mapping bias is the largest problem for reliable estimates of allele specific gene expression using RNA-seq, Allim combines two different measures to account for mapping biases. First, Allim generates a polymorphism aware reference genome that accounts for the sequence variation between the alleles of both parents (or parental lines). Second, Allim includes a sequence specific simulation tool to estimate the remaining mapping bias. This estimated mapping bias is then incorporated in the statistical tests for allelic imbalance.

The pipeline requires whole transcript high throughput RNA sequencing (RNA-seq) reads of F1 hybrids. Additionally, either RNA-seq reads for both parents, genomic sequencing reads for both parents or two parental genomes for both parents have to be provided in separate files. Allim was tested on Illumina paired-end RNA-seq reads but it can also handle FASTQ files from other NGS sequencing platforms. The provided parental RNA-seq libraries can be from homozygous as well as heterozygous parents (however, all heterozygous loci will be excluded during the analysis).
Allim has five modules that can be run by a single command. All input parameters can be specified in the AllimOptions file. These parameters are then used to run the complete pipeline. Allim provides two different input options:

**AllimOptions\_2Pexpr:** This configuration file can be used when RNA-seq data for parent1, parent2 and F1 (hybrid) is provided. In this case Allim uses parent1 and parent2 RNA-seq reads to identify the fixed SNPs. The fixed SNPs are subsequently used to generate two “parental -genomes”, one for each parent. Note, rather than RNA-seq data, also DNA reads can be provided for the two parents.<br>

In addition with Allim if user does not have reference genome and corresponding gene annotation in GTF format. It is possible to provide a single reference transcriptome along with RNA-seq data of both parents. The use of a reference transcriptome (or contigs from a RNA-Seq de novo assembly) instead of a reference genome is accompanied by the following differences:<br>
<br>
<ul><li>In a transcriptome assembly different isoforms are typically presented by different contigs. Furthermore, assemblies often contain additional redundancies between contigs due to sequencing errors and polymorphisms in the RNA-seq reads used for de novo transcriptome assembly.</li></ul>

<ul><li>As Allim creates two parental references (genomes or here transcriptomes), and maps RNA-seq reads of the F1 individual to the “diploid genome”, only reads that map non-ambiguously can be used to determine ASE profiles. It is therefore recommended to remove redundancy from the transcriptome assembly before Allim usage.</li></ul>

<ul><li>As transcriptome assemblies typically do not have a gtf file containing gene features (needed as Allim input). We have therefore added an additional short script to obtain a simple gtf file based on contig ids.</li></ul>

<b>AllimOptions_2Pgenomes:</b> This configuration file can be used when the user provides RNA-seq data for the F1 hybrid and a genomic sequence for each parental line. If genomic sequences of both parents are provided, both FASTA files need to have the same size and identical FASTA IDs.<br>
<br>
<br>
<br>
<h2>Five Modules of Allim</h2>

<ul><li>Identification of fixed SNPs and ceating two parental genomes<br>
</li><li>Computer simulation of RNA-seq reads with fixed SNPs<br>
</li><li>Estimation of the remaining mapping bias with simulated data<br>
</li><li>Estimation of allele specific expression for experimental data<br>
</li><li>Statistical test of significant allelic imbalance</li></ul>

The source code and the user manual of Allim are available from <a href='http://code.google.com/p/allim/'>http://code.google.com/p/allim/</a>

<h1>Usage</h1>

<ul><li>Allim Configuration file when two parents expression data available: <a href='http://code.google.com/p/allim/wiki/AllimOptions_2Pexpr'>http://code.google.com/p/allim/wiki/AllimOptions_2Pexpr</a>
</li><li>Allim Configuration file when two parents expression data not available but two parents genome fasta files available: <a href='http://code.google.com/p/allim/wiki/AllimOptions_2Pgenomes'>http://code.google.com/p/allim/wiki/AllimOptions_2Pgenomes</a>
</li><li>User Manual: <a href='http://code.google.com/p/allim/wiki/Manual'>http://code.google.com/p/allim/wiki/Manual</a><br></li></ul>

<h1>How to cite <code>Allim</code>?</h1>

Please cite the following <code>Allim</code> paper.<br>
<br>
Pandey RV, Franssen SU, Futschik A, Schlötterer C. (2013) Allelic imbalance metre (Allim), a new tool for measuring allele-specific gene expression with RNA-seq data. Mol Ecol Resour. 13(4):740-745.<br>
