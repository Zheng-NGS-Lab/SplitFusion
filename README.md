# SplitFusion - a fast pipeline for detection of gene fusion based on fusion-supporting split alignment.

Gene fusion is a hallmark of cancer. Many gene fusions are effective therapeutic targets such as BCR-ABL in chronic myeloid leukemia, EML4-ALK in lung cancer, and any of a number of partners-ROS1 in lung cancer. Accurate detection of gene fusion plays a pivotal role in precision medicine by matching the right drugs to the right patients.

Challenges in the diagnosis of gene fusions include poor sample quality, limited amount of available clinical specimens, and complicated gene rearrangements. The anchored multiplex PCR (AMP) is a clinically proven technology designed, in one purpose, for robust detection of gene fusions across clinical samples of different types and varied qualities, including RNA extracted from FFPE samples.

**SplitFusion** can be used for RNA-seq data and the Anchored Multiplex PCR (AMP) data, for detecting gene fusions based on BWA-MEM split alignments, i.e. reads crossing fusion breakpoints, with the ability to accurately infer frame-ness and exon-boundary alignments for prediction functional fusions. SplitFusion also outputs example breakpoint-supporting seqeunces in FASTA and BAM format, allowing for further investigations.


## Reference publication
[Zheng Z, et al. Anchored multiplex PCR for targeted next-generation sequencing. Nat Med. 2014](http://www.nature.com/nm/journal/v20/n12/full/nm.3729.html)

## How does SplitFusion work?  


The analysis consists of ## computational steps:

1. Reference alignment and deduplication.

2. CIGAR transformation.

3. Candidate breakpoint calling.

4. Initial breakpoint filtering

5. Breakpoint gene annotation, frame-ness, exon boundary, further filtering and target reporting

6. Result reporting and visualization.


Lastly, outputs a summary table and breakpoint-spanning reads.



## Installation

### 1. Downloading database required by SplitFusion (--database_dir)

```java

1. refGenome:

wget https://data.broadinstitute.org/snowman/hg19/Homo_sapiens_assembly19.fasta
wget https://data.broadinstitute.org/snowman/hg19/Homo_sapiens_assembly19.fasta.amb
wget https://data.broadinstitute.org/snowman/hg19/Homo_sapiens_assembly19.fasta.ann
wget https://data.broadinstitute.org/snowman/hg19/Homo_sapiens_assembly19.fasta.bwt
wget https://data.broadinstitute.org/snowman/hg19/Homo_sapiens_assembly19.fasta.fai
wget https://data.broadinstitute.org/snowman/hg19/Homo_sapiens_assembly19.fasta.pac
wget https://data.broadinstitute.org/snowman/hg19/Homo_sapiens_assembly19.fasta.sa

```

### 2. Downloading third-party tools required by SplitFusion (--tool_dir)
```java

1. R: download from https://www.r-project.org/

2. perl: download from https://www.perl.org/get.html

3. bwa: download from https://sourceforge.net/projects/bio-bwa/files/

4. bedtools: download from https://bedtools.readthedocs.io/en/latest/content/installation.html

5. samtools: download from http://samtools.sourceforge.net/

6. annovar: http://download.openbioinformatics.org/annovar_download_form.php

```

### 3. Dependencies

- R packages ("data.table", "plyr")


```java

> install.packages(c("data.table", "plyr"))

```

### 4. Installing SplitFusion

```java
wget https://github.com/Zheng-NGS-Lab/SplitFusion/archive/master.zip
unzip master.zip
mv SplitFusion-master SplitFusion
R CMD INSTALL SplitFusion
```


## Run

### 1. Help (NOTE: Python 2.7 not version 3)
```java
python ./SplitFusion/exec/SplitFusion.py -h

usage: SplitFusion.py [-h] --refGenome REFGENOME --database_dir DATABASE_DIR
                      --annovar ANNOVAR --samtools SAMTOOLS --bedtools
                      BEDTOOLS --bwa BWA --R R --perl PERL --output OUTPUT
                      --sample_id SAMPLE_ID [--bam_dir BAM_DIR]
                      [--fastq_dir FASTQ_DIR] [--r1filename R1FILENAME]
                      [--r2filename R2FILENAME] [--panel_dir PANEL_DIR]
                      [--panel PANEL] [--steps STEPS]
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
split reads and/or paired-end reads.

optional arguments:
  -h, --help            show this help message and exit
  --refGenome REFGENOME
                        [required]. the path where human genome reference is
                        stored
  --database_dir DATABASE_DIR
                        [required]. the path where large databases e.g.
                        reference genome and annotation databases are stored
  --annovar ANNOVAR     [required]. the path of annovar
  --samtools SAMTOOLS   [required]. the path of samtools
  --bedtools BEDTOOLS   [required]. the path of bedtools
  --bwa BWA             [required]. the path of bwa
  --R R                 [required]. the path of R
  --perl PERL           [required]. the path of perl
  --output OUTPUT       [required]. the path where output is stored
  --sample_id SAMPLE_ID
                        [required]. the sample name of running
  --bam_dir BAM_DIR     the path where bam or fastq file is stored.
                        [Kickstart] The bam file of the sameple_id (xxx.bam or
                        xxx.consolidated.bam) will be used. Either fastq_dir
                        or bam_dir should be specified
  --fastq_dir FASTQ_DIR
                        the path where fastq file is stored. Either fastq_dir
                        or bam_dir should be specified
  --r1filename R1FILENAME
                        Read 1 fastq filename. Can be in gzipped format. If
                        not specified, $fastq_dir/$sample_id.R1.fq will be
                        used.
  --r2filename R2FILENAME
                        Read 2 fastq filename. Can be in gzipped format. If
                        not specified, $fastq_dir/$sample_id.R2.fq will be
                        used.
  --panel_dir PANEL_DIR
                        For Target mode: the path where known significant
                        fusions or splicing isoforms (white.list) or unwanted
                        fusions involving homologous genes or recurrent falsed
                        positives (black.list) are stored. default='NA'
  --panel PANEL         the prefix name of target genes panel file is stored,
                        e.g., LungFusion for LungFusion.GSP2.bed. default='NA'
  --steps STEPS         specify steps to run. default='1_fastq-bam,2_bam-
                        breakpoint,3_breakpoint-filter,4_breakpoint-
                        anno,5_breakpoint-anno-post'
  --AnnotationMethod ANNOTATIONMETHOD
                        the name of annotation tools. default = 'annovar'
  --thread THREAD       number of threads for computing. default=4
  --minMQ MINMQ         minimum mapping quality. default=13
  --minMQ1 MINMQ1       minimum mapping quality of a leftmost of Read1
                        (rightmost of Read2). default=30
  --minMapLength MINMAPLENGTH
                        minimum read mapping length. default=18
  --minMapLength2 MINMAPLENGTH2
                        minimum mapping length of rightmost of Read1 (leftmost
                        of Read2). default=25
  --maxQueryGap MAXQUERYGAP
                        maximum gap length on a query read of split
                        alignments. default=0
  --maxOverlap MAXOVERLAP
                        maximum overlap bases of two split alignments.
                        default=6
  --minExclusive MINEXCLUSIVE
                        minimum exclusive length between two split alignments.
                        default=18
  --FusionMinStartSite FUSIONMINSTARTSITE
                        minimum number of Adaptor Ligation Read Starting Sites
                        to call Structure Variation/Fusion. Should be less or
                        equal minPartnerEnds_BothExonJunction. default=1
  --minPartnerEnds_BothExonJunction MINPARTNERENDS_BOTHEXONJUNCTION
                        minimum number of fusion partner ends (ligation site),
                        when both breakpoints are at exon junctions, to call
                        Structure Variation/Fusion. default=1
  --minPartnerEnds_OneExonJunction MINPARTNERENDS_ONEEXONJUNCTION
                        minimum number of fusion partner ends (ligation site),
                        when one breakpoint is at exon junction, to call
                        Structure Variation/Fusion. default=3

```

### 2. run SplitFusion
[An example command:](https://github.com/Zheng-NGS-Lab/SplitFusion/blob/master/example.command.txt)

[An example data for test (click here):](https://github.com/Zheng-NGS-Lab/SplitFusion/blob/master/inst/data/example_data/)

[An internal fusion library data for downloading as --panel_dir of target mode (click here):](https://github.com/Zheng-NGS-Lab/SplitFusion/blob/master/inst/data/panel/)

```java
#python ./SplitFusion/exec/SplitFusion.py --SplitFusionPath SplitFusionPath --refGenome refGenome --bam_dir bam_dir --sample_id sample_id --output output --R R --perl perl --database_dir database_dir --panel_dir SplitFusionPath/data/panel/
# I installed SplitFusion under
#       /home/zz/repo/
# copy ./SplitFusion/inst/data/example_data/Lib001.R1.fq ./SplitFusion/inst/data/example_data/Lib001.R2.fq /home/zz/repo/test
# Example run:

##=========================================================
## Start from FASTQ files, no panel info
## , compatible with RNA-seq whole transcriptome analysis
##=========================================================
python /home/zz/repo/SplitFusion/exec/SplitFusion.py \
        --SplitFusionPath /home/zz/R/x86_64-pc-linux-gnu-library/3.5/SplitFusion \
        --sample_id "Lib001" \
        --fastq_dir /home/zz/repo/test \
        --database_dir /home/zz/repo/database \
	--tool_dir /home/zz/repo/tool \
        --r1filename "Lib001".R1.fq \
        --r2filename "Lib001".R2.fq \
        --output /home/zz/repo/test \
        --refGenome Homo_sapiens_assembly19.fasta \
        --thread 6 &


##=========================================================
## Kickstart mode, no panel info
## , compatible with RNA-seq whole transcriptome analysis
##=========================================================
python /home/zz/repo/SplitFusion/exec/SplitFusion.py \
        --SplitFusionPath /home/zz/R/x86_64-pc-linux-gnu-library/3.5/SplitFusion \
        --sample_id "Lib001" \
        --bam_dir /home/zz/repo/test \
        --database_dir /home/zz/repo/database \
	--tool_dir /home/zz/repo/tool \
        --output /home/zz/repo/test \
        --refGenome Homo_sapiens_assembly19.fasta \
        --thread 6 &


##===============================
## TARGET mode, with panel info
##===============================
python /home/zz/repo/SplitFusion/exec/SplitFusion.py \
        --SplitFusionPath /home/zz/R/x86_64-pc-linux-gnu-library/3.5/SplitFusion \
        --sample_id Lib001 \
        --fastq_dir /home/zz/repo/test \
        --database_dir /home/zz/repo/database \
	--tool_dir /home/zz/repo/tool \
        --panel_dir /home/zz/repo/panel \
        --panel LungFusion \
        --r1filename "Lib001".R1.fq \
        --r2filename "Lib001".R2.fq \
        --output /home/zz/repo/test \
        --refGenome Homo_sapiens_assembly19.fasta \
        --thread 6 &


##===============================
## Selecting only some steps to run
##===============================
python /home/zz/repo/SplitFusion/exec/SplitFusion.py \
        --SplitFusionPath /home/zz/R/x86_64-pc-linux-gnu-library/3.5/SplitFusion \
        --sample_id "Lib001" \
        --bam_dir /home/zz/repo/test \
        --database_dir /home/zz/repo/database \
	--tool_dir /home/zz/repo/tool \
        --panel_dir /home/zz/repo/panel \
        --panel LungFusion \
        --output /home/zz/repo/test \
        --refGenome Homo_sapiens_assembly19.fasta \
        --steps "3_breakpoint-filter,4_breakpoint-anno,5_breakpoint-anno-post" \
        --thread 6 &
```

## Output 
[An example brief output table:](https://github.com/Zheng-NGS-Lab/SplitFusion/blob/master/inst/data/example_data/target_mode_result/Lib001/Lib001.brief.summary)

|	SampleID	|	GeneExon5_GeneExon3	|	frame	| num_partner_ends |        num_unique_reads |        exon.junction |   breakpoint |      transcript_5 |    transcript_3 |    function_5 |      function_3 |      gene_5 |  cdna_5 |  gene_3 |  cdna_3 |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| Lib001 | EML4_exon6---ALK_exon20 | in-frame | 1 | 1 | Both | 2_29446396__2_42491871 | NM_019063 | NM_004304 | exonic | exonic | EML4 | 667 | ALK | 3171 |
| Lib001 | EML4_intronic---ALK_exon20 | N.A. | 6 | 11 | One | 2_29446396__2_42492091 | NM_019063 | NM_004304 | intronic | exonic | EML4 | - | ALK | 3171 |

[An example output fastq file for the EML4_intronic---ALK_exon20 fusion of sample Lib001 is:](https://github.com/Zheng-NGS-Lab/SplitFusion/blob/master/inst/data/example_data/target_mode_result/Lib001/Lib001.EML4_intronic---ALK_exon20.txt)

```java

>NS500673:45:HHK2HAFXX:1:21106:26233:4628:
ATGGCTTGCAGCTCCTGGTGCTTCCGGCGGTACACTTGGCTGTTTTTTTCGCGAGTTGACATTTTTGCTTGATTAAAGATGTCATCATT
>NS500673:45:HHK2HAFXX:1:21108:4972:6200:
ATGGCTTGCAGCTCCTGGTGCTTCCGGCGGTACACTGGCTGTTTTTTTCGCGAGTTTACATTTTTGCTTGGTTGATT
>NS500673:45:HHK2HAFXX:2:11207:14331:12205:
ATGGCTTGCAGCTCCTGGTGCTTCCGGCGGTACACTGGCTGTTTTTTTCGCGAGTTGACATTTTTGCTTGGTTGATG
>NS500673:45:HHK2HAFXX:2:11301:14903:19850:
ATGGCTTGCAGCTCCTGGTGCTTCCGGCGGTACACTTGGCTGTTTTTTTCGCGAGTTGACATTTTTG
>NS500673:45:HHK2HAFXX:2:21111:24355:8828:
ATGGCTTGCAGCTCCTGGTGCTTCCGGCGGTACACTTGGCTGTTTTTTTCGCGAGTTGACATTTTTG
>NS500673:45:HHK2HAFXX:4:11406:15146:11569:
ATGGCTTGCAGCTCCTGGTGCTTCCGGCGGTACACTTGGCCGTTTTTTTCGCGAGTTGACATTTTTG
>NS500673:45:HHK2HAFXX:4:11606:2779:2081:
ATGGCTTGCAGCTCCTGGTGCTTCCGGCGGTACACTTGGCTGTTTTTTTCGCGAGTTGACATTTTTGCTTGGTTGATGATGACATCTTT
>NS500673:45:HHK2HAFXX:4:21409:22050:11159:
ATGGCTTGCAGCTCCTGGTGCTTCCGGCGGTACACTGGCTGTTATTTTCGCGAGTAGACATTTTTGCTTGGTTGATG
>NS500673:45:HHK2HAFXX:4:21508:24201:16676:
ATGGCTTGCAGCTCCTGGTGCTTCCGGCGGTACACTTGGCTGTTTTTTTCGCGAGTTGACATTTTTG

```

## Visualization (PC or Mac)

```java

> if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

> BiocManager::install("igvR")

> library(SplitFusion)

> bam2igv(bamfile = "Lib001.EML4_intronic---ALK_exon20.bam")


```

[An visualization of the EML4 intronic fusion site:](https://github.com/Zheng-NGS-Lab/SplitFusion/blob/master/inst/data/example_data/target_mode_result/Lib001/Lib001.EML4_intronic---ALK_exon20.bam.2.42492091.svg)
![Segmentfault](https://github.com/Zheng-NGS-Lab/SplitFusion/blob/master/inst/data/example_data/target_mode_result/Lib001/Lib001.EML4_intronic---ALK_exon20.bam.2.42492091.svg)

[An visualization of the ALK exon20 fusion site:](https://github.com/Zheng-NGS-Lab/SplitFusion/blob/master/inst/data/example_data/target_mode_result/Lib001/Lib001.EML4_intronic---ALK_exon20.bam.2.29446361.svg)
![Segmentfault](https://github.com/Zheng-NGS-Lab/SplitFusion/blob/master/inst/data/example_data/target_mode_result/Lib001/Lib001.EML4_intronic---ALK_exon20.bam.2.29446361.svg)

