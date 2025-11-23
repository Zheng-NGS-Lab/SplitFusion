## Introduction

Gene fusion is a hallmark of cancer. Many gene fusions are effective therapeutic targets such as BCR-ABL in chronic myeloid leukemia and EML4-ALK in lung cancer. Accurate detection of gene fusion plays a pivotal role in precision medicine by matching the right drugs to the right patients.

Challenges in the diagnosis of gene fusions include that there could be many and sometimes unknown fusion partners, low gene expression (e.g. ALK), non fusion-specific protein expression (e.g. ROS1), potential involvments of cryptic splice sites, low sequence diversity at genomic breakpoints and associated mapping difficulty, and poor sample (poor quality, limited amount and low tumor cellularity) in clinical specimens. The anchored multiplex PCR (AMP) is a clinically proven technology that addresses all these issues and has accelerated gene fusion discoveries and supported robust clinical diagnosis ([Zheng Z, et al. Anchored multiplex PCR for targeted next-generation sequencing. Nat Med. 2014](http://www.nature.com/nm/journal/v20/n12/full/nm.3729.html)).

Equally important to a robust wet lab technology is a high-performing computational method for calling gene fusions, for which we develop SplitFusion ([Bian et al. SplitFusion enables ultrasensitive gene fusion detection and reveals fusion variant-associated tumor heterogeneity. Patterns, 2025](https://www.cell.com/patterns/fulltext/S2666-3899(25)00022-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2666389925000224%3Fshowall%3Dtrue)). SplitFusion leverages BWA-MEM split alignments, can detect **cryptic splice-site fusions (e.g., EML4::ALK v3b and ARv7)**, call fusions involving **highly repetitive gene partners (e.g., CIC::DUX4)**, and **infer frame-ness and exon-boundary** alignments for functional prediction and minimizing false positives. Using 1,848 datasets of various sizes, SplitFusion demonstrated superior sensitivity and specificity compared to three other tools. In 1,076 formalin-fixed paraffin-embedded lung cancer samples, SplitFusion identified novel fusions and revealed that EML4::ALK variant 3 was associated with multiple fusion variants coexisting in the same tumor. Additionally, SplitFusion can call targeted splicing variants. Using data from 515 The Cancer Genome Atlas (TCGA) samples, SplitFusion showed the highest sensitivity and **uncovered two cases of SLC34A2::ROS1 that were missed** in previous studies. These capabilities make SplitFusion highly suitable for clinical applications and the study of fusion-defined tumor heterogeneity.

**SplitFusion** can be used for standard RNA-seq data and the Anchored Multiplex PCR (AMP) data.


## How does SplitFusion work?  

The analsyis consists of the following computational steps:

1. Reference alignment and deduplication.

2. CIGAR transformation.

3. Candidate breakpoint calling.

4. Initial breakpoint filtering.

5. Breakpoint gene annotation, frame-ness, exon boundary, further filtering and target reporting.

6. Result reporting and visualization.


## Dependencies

When running SplitFusion, you can specify paths to the tools and genome files you already have. If not, here are the human genome data and tools for installation.

- [human genome](https://data.broadinstitute.org/snowman/hg19/) 

	E.g. I saved my large database files under /home/user1/database/:

		cd /home/user1/database
		wget https://data.broadinstitute.org/snowman/hg19/Homo_sapiens_assembly19.fasta 

- [samtools](http://samtools.sourceforge.net)

	E.g. I installed the above tool in /home/user1/tools/:
	
		cd /home/user1/tools
		wget -O samtools.tar.bz2 https://sourceforge.net/projects/samtools/files/latest/download
			tar -xvf samtools.tar.bz2
			cd samtools-1.10 ## Note the samtools-x.xx version
			mkdir /home/user1/tools/samtools
			./configure --prefix=/home/user1/tools/samtools
			make; make install; cd ..

- [bwa](https://sourceforge.net/projects/bio-bwa/files)

		git clone https://github.com/lh3/bwa.git
		        cd bwa; make; cd ..

- [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)

		wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
			tar -zxvf bedtools-2.29.1.tar.gz
			cd bedtools2; make

- [perl](https://www.perl.org/get.html)


	Currently, SplitFusion by default uses snpEff but it also supports ANNOVA. snpEff is a free annotation software. When you run it for the first time, it will automatically download the genome database (either hg19 or hg38) from the internet. To run with snpEff, you only need to specify your snpEff directory via --annovar /home/user1/tools/snpEff

- [snpEff](https://pcingola.github.io/SnpEff/download/)
- [annovar](http://download.openbioinformatics.org/annovar_download_form.php)	

        If you use ANNOVAR, you need a free registration. Then run SplitFusion with the flag --AnnotationMethod annovar. Note that the annovar sub-directory structure should be maintained, e.g. if you install annovar under /home/user1/tools:
	
		/home/user1/tools/annovar/annotate_variation.pl
		
		/home/user1/tools/annovar/table_annovar.pl
		
		/home/user1/tools/annovar/humandb/hg19_refGeneMrna.fa
		
		/home/user1/tools/annovar/humandb/hg19_refGene.txt

- [R](https://www.r-project.org/) Install requried R packages within R:

	> install.packages(c("Rcpp", "data.table", "plyr"))

- [Python](https://www.python.org/downloads/) Install a Python module by [pip](https://pip.pypa.io/en/stable/installing/):

		pip install future



## Installation

	cd /home/user1/tools/
	git clone https://github.com/Zheng-NGS-Lab/SplitFusion.git


## Run

#### 1. Help

```
	python /home/user1/tools/SplitFusion/exec/SplitFusion.py -h

usage: SplitFusion.py [-h] --refGenome REFGENOME --annovar SNPEFF_PATH --samtools
                      SAMTOOLS --bedtools BEDTOOLS --bwa BWA --R R --perl PERL
                      --output OUTPUT --sample_id SAMPLE_ID
                      [--bam_file BMS_FILE] [--fastq_file1 FASTQ_FILE1] [--fastq_file2 FASTQ_FILE2]
                      [--panel_dir PANEL_DIR] [--panel PANEL] [--steps STEPS]
                      [--AnnotationMethod ANNOTATIONMETHOD] [--thread THREAD]
                      [--minMQ MINMQ] [--minMQ1 MINMQ1]
                      [--minMapLength MINMAPLENGTH]
                      [--minMapLength2 MINMAPLENGTH2]
                      [--maxQueryGap MAXQUERYGAP] [--maxOverlap MAXOVERLAP]
                      [--minExclusive MINEXCLUSIVE]
                      [--FusionMinStartSite FUSIONMINSTARTSITE]
                      [--minPartnerEnds_BothExonJunction MINPARTNERENDS_BOTHEXONJUNCTION]
                      [--minPartnerEnds_OneExonJunction MINPARTNERENDS_ONEEXONJUNCTION]

Split-Fusion is a fast data analysis pipeline detects gene fusion based on
chimeric split-read alignments.

optional arguments:
  -h, --help            show this help message and exit
  --refGenome REFGENOME
                        The reference genome file, with a full path
                        [required].
  --annovar SNPEFF	The snpEff path [required].
  --samtools SAMTOOLS   The samtools executable file with full path [Optional].
  --bedtools BEDTOOLS   The bedtools executable file with full path [Optional].
  --bwa BWA             The bwa executable file with full path [Optional].
  --R R                 The R executable file with full path [Optional].
  --perl PERL           The perl executable file with full path [Optional].
  --output OUTPUT       The directory for output SplitFusion results [required].
  --sample_id SAMPLE_ID
                        The name of sample to be analyzed [required].
  --bam_file BAM_FILE   If the bam_file is specified, the Kickstart mode will
                        be used. [Either fastq_file or bam_file should be
                        specified]
  --fastq_file1 FASTQ_FILE1
                        The fastq file (Read 1 of paired-end) to be analyzed
                        [Either fastq_file or bam_file should be specified].
  --fastq_file2 FASTQ_FILE2
                        Read 2 of paired-end fastq file.
  --panel_dir PANEL_DIR
                        For TARGET mode: the path where known significant
                        fusions or splicing isoforms (whitelist) or unwanted
                        fusions involving homologous genes or recurrent falsed
                        positives (blacklist) are stored. Default='NA'
  --panel PANEL         The prefix name of TARGET gene panel file. E.g.,
                        'LungFusion' for LungFusion.GSP2.bed. Default='NA'
  --steps STEPS         Specify steps to run. Default='1_fastq-bam,2_bam-
                        breakpoint,3_breakpoint-filter,4_breakpoint-anno
                        ,5_breakpoint-anno-post'

```


#### 2. Examples

```java

# Examples of running different modes of SplitFusion

## First, copy example data files from pipeline to local testing directory, e.g.:

	mkdir -p /home/user1/SplitFusion-test/data
	cp /home/user1/tools/SplitFusion/inst/data/example_data/Lib001.* /home/user1/SplitFusion-test/data/

## Before using BWA, the reference genome needs to be indexed, e.g.:

	/home/user1/tools/bwa/bwa index -a bwtsw /home/user1/database/Homo_sapiens_assembly19.fasta

##=========================================================
## Start from FASTQ files, no panel info
## , compatible with RNA-seq whole transcriptome analysis
##=========================================================
python /home/user1/tools/SplitFusion/exec/SplitFusion.py \
        --refGenome /home/user1/database/Homo_sapiens_assembly19.fasta \
        --annovar /home/user1/tools/snpEff \
        --samtools /home/user1/tools/samtools/bin/samtools \
        --bedtools /home/user1/tools/bedtools2/bin/bedtools \
        --bwa /home/user1/tools/bwa/bwa \
        --R /usr/bin/R \
        --perl /usr/bin/perl \
        --output /home/user1/SplitFusion-test/output \
        --sample_id "Lib001" \
        --fastq_file1 /home/user1/SplitFusion-test/data/Lib001.R1.fq \
        --fastq_file2 /home/user1/SplitFusion-test/data/Lib001.R2.fq \
        --thread 4

##=========================================================
## Kickstart mode, no panel info
## , compatible with RNA-seq whole transcriptome analysis
##=========================================================
python /home/user1/tools/SplitFusion/exec/SplitFusion.py \
        --refGenome /home/user1/database/Homo_sapiens_assembly19.fasta \
        --annovar /home/user1/tools/snpEff \
        --samtools /home/user1/tools/samtools/bin/samtools \
        --bedtools /home/user1/tools/bedtools2/bin/bedtools \
        --bwa /home/user1/tools/bwa/bwa \
        --R /usr/bin/R \
        --perl /usr/bin/perl \
        --output /home/user1/SplitFusion-test/kickstart-mode-output \
        --sample_id "Lib001" \
	--bam_file /home/user1/SplitFusion-test/data/Lib001.bam \
        --thread 1

##=========================================================
## TARGET gene panel mode, with panel info
##=========================================================
mkdir -p /home/user1/SplitFusion-test/panel
cp /home/user1/tools/SplitFusion/inst/data/panel/* /home/user1/SplitFusion-test/panel/

python /home/user1/tools/SplitFusion/exec/SplitFusion.py \
        --refGenome /home/user1/database/Homo_sapiens_assembly19.fasta \
        --annovar /home/user1/tools/snpEff \
        --samtools /home/user1/tools/samtools/bin/samtools \
        --bedtools /home/user1/tools/bedtools2/bin/bedtools \
        --bwa /home/user1/tools/bwa/bwa \
        --R /usr/bin/R \
        --perl /usr/bin/perl \
        --output /home/user1/SplitFusion-test/target-mode-output \
        --sample_id "Lib001" \
        --fastq_file1 /home/user1/SplitFusion-test/data/Lib001.R1.fq \
        --fastq_file2 /home/user1/SplitFusion-test/data/Lib001.R2.fq \
	--panel_dir /home/user1/SplitFusion-test/panel \
	--thread 4 \
        --panel LungFusion

##=======================================================================
## Selecting only some steps to run.
##   If you have run Steps 1 and 2, and just want to repeat analyses from
## the Step 3 (e.g. the whitelist/blacklist has been updated), then:
##=======================================================================
python /home/user1/tools/SplitFusion/exec/SplitFusion.py \
        --refGenome /home/user1/database/Homo_sapiens_assembly19.fasta \
        --annovar /home/user1/tools/snpEff \
        --samtools /home/user1/tools/samtools/bin/samtools \
        --bedtools /home/user1/tools/bedtools2/bin/bedtools \
        --bwa /home/user1/tools/bwa/bwa \
        --R /usr/bin/R \
        --perl /usr/bin/perl \
        --output /home/user1/SplitFusion-test/output \
        --sample_id "Lib001" \
        --fastq_file1 /home/user1/SplitFusion-test/data/Lib001.R1.fq \
        --fastq_file2 /home/user1/SplitFusion-test/data/Lib001.R2.fq \
        --panel_dir /home/user1/SplitFusion-test/panel \
        --panel LungFusion \
	--steps "3_breakpoint-filter,4_breakpoint-anno,5_breakpoint-anno-post"

```

## Output 

[An example brief output table:](https://github.com/Zheng-NGS-Lab/SplitFusion/blob/master/inst/data/example_data/output/Lib001/Lib001.brief.summary)

|	SampleID	|	GeneExon5_GeneExon3	|	frame	| num_partner_ends |        num_unique_reads |        exon.junction |   breakpoint |      transcript_5 |    transcript_3 |    function_5 |      function_3 |      gene_5 |  cdna_5 |  gene_3 |  cdna_3 |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| Lib001 | EML4_intronic---ALK_exon20 | in-frame | 9 | 11 | Both | 2_29446396__2_42492091 | NM_019063 | NM_004304 | intronic | exonic | EML4 | 667 | ALK | 3171 |
| Lib001 | EML4_exon6---ALK_exon20 | in-frame | 2 | 1 | Both | 2_29446396__2_42491871 | NM_019063 | NM_004304 | exonic | exonic | EML4 | 667 | ALK | 3171 |



[An example output fastq file for the EML4_intronic---ALK_exon20 fusion of sample Lib001 is:](https://github.com/Zheng-NGS-Lab/SplitFusion/blob/master/inst/data/example_data/output/Lib001/Lib001.EML4_intronic---ALK_exon20.txt)

```java

>NS500673:45:HHK2HAFXX:1:11206:14487:5244:
ATTGCTTGCAGCTCCTGGTGCTTCCGGCGGTACACTTGGCTATTTTTTTCGCGAGTTGACATTTTTGCTTGGTTGATGATGACATCTTTATGCTTGTCTGCAATTTTGGTAACTTTTG
>NS500673:45:HHK2HAFXX:1:11307:26804:17812:
ATGGCTTGCAGCTCCTGGTGCTTCCGGCGGTACACTTGGCTGTTTTTTTCGCGAGTTGACATTTTTGCTTGGTTGATGATGACATCTTTATGCTTG
>NS500673:45:HHK2HAFXX:2:11207:15517:8468:
ATGGCTTGCAGCTCCTGGTGCTTCCGGCGGTACACTTGGCTGTTTTTTTCGCGAGTTGACATTTTTG
>NS500673:45:HHK2HAFXX:2:21212:23872:20046:
ATGGCTTGCAGCTCCTGGTGCTTCCGGCGGTACACTTGGCTGTTTTTTTCGCGAGTTGACATTTTTGCTTGGTTGATGATGACATCTTTATGCTTG
>NS500673:45:HHK2HAFXX:3:21410:6281:16202:
ATGGCTTGCAGCTCCTGGTGCTTCCGGCGGTACACTTGGCTGTTTTTTTCGCGAGTTGACATTTTTGCTTGGTTGATGATGACATCTTTATGCTTG
>NS500673:45:HHK2HAFXX:3:21605:12646:6668:
GTCATCATCAACCAAGCAAAAATGTCAACTCGCGAAAAAAACAGCCAAGTGTTCCGCCGGAAGCACCAGGAGCTGCAAGCCAT
>NS500673:45:HHK2HAFXX:4:21510:10840:13732:
ATGGCTTGCAGCTCCTGGTGCTTCCGGCGGTACACTTGGCTGTTTTTTTCGCGAGTTGACATTTTTGCTTAGTTGATGATGACATC
>NS500673:45:HHK2HAFXX:4:21605:21501:18033:
ATGGCTTGCAGCTCCTGGTGCTTCCGGCGGTACACTTGGCTGTTTTTTTCGCGAGTTGACATTTTTGCTTGGTTGGTGATGACATCTTTATGCTTG

```

## Visualization (on PC or Mac)

```java

Within R, run:

> if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
> BiocManager::install("igvR")
> source("SplitFusion/R/bam2igv.R")
> bam2igv(bamfile = "/home/user1/SplitFusion-test/output/Lib001/Lib001.EML4_intronic---ALK_exon20.bam")

```

[An visualization of the EML4 intronic fusion site:](https://github.com/Zheng-NGS-Lab/SplitFusion/blob/master/inst/data/example_data/output/Lib001/Lib001.EML4_intronic---ALK_exon20.bam.2.42492091.svg)
![Segmentfault](https://github.com/Zheng-NGS-Lab/SplitFusion/blob/master/inst/data/example_data/output/Lib001/Lib001.EML4_intronic---ALK_exon20.bam.2.42492091.svg)

[An visualization of the ALK exon20 fusion site:](https://github.com/Zheng-NGS-Lab/SplitFusion/blob/master/inst/data/example_data/output/Lib001/Lib001.EML4_intronic---ALK_exon20.bam.2.29446361.svg)
![Segmentfault](https://github.com/Zheng-NGS-Lab/SplitFusion/blob/master/inst/data/example_data/output/Lib001/Lib001.EML4_intronic---ALK_exon20.bam.2.29446361.svg)

