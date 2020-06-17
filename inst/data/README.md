Big data files that needed to be first downloaded from public webportal and stored under $database_dir.

| Filename                   | Content                                                        |
|:--------------------------:|:--------------------------------------------------------------:|
| Homo_sapiens_assembly19.fasta | Contains a list of human genome reference, please mannually downloaded from ucsc or other official site. |
| hg19_refGeneMrna.fa | File needed by ANNOVAR. |
| hg19_refGene.txt | File needed by ANNOVAR. |
| hg19_refGeneVersion.txt | File needed by ANNOVAR. |
| hg19_refGeneWithVer.txt | File needed by ANNOVAR. |
| hg19_refGeneWithVerMrna.fa | File needed by ANNOVAR. |


For TARGET mode, panel-specific information can be used (under $panel_dir).

| Filename                   | Content                                                        |
|:--------------------------:|:--------------------------------------------------------------:|
| known.gene-gene.txt | A list of known significant fusion between two genes for targeted output even evidence is weak e.g. low number of reads and out-of-frame. |
| known.gene-exon.txt | By default, SplitFusion only outputs fusions that involve two different genes. This file contains clinically relevant isoforms of a gene, e.g. "MET_exon13---MET_exon15" that is an important theraputic target. Many exon skipping/alternative splicing events are normal or of unknown clinical relevance and are thus not output by default. |
| known.partners.txt | Contains a list of known fusion partners of targets. |
| known.3UTR.truncation.txt | Contains a list of known significant genes with 3'UTR trancation, e.g. FGFRs exon18 deletion |
| filter.gene-gene.txt | Contains recurrent breakpoints identified as data accumulates, but are not of interest. |

The above files could be updated periodically as a backend supporting database that facilitates automatc filtering and outputing of fusion candidates.


Two files come with SplitFusion are premade reference gene transcripts to resolve some nomenclatural conflicts

| Filename                   | Content                                                        |
|:--------------------------:|:--------------------------------------------------------------:|
| refGene.mostExon.NM | For genes with multiple isoforms, selecte the one with the most number of exons. |
| refGene0.txt | For genes with multiple isoforms, selecte the one with the smallest NM##, which is often the classic isoform. |

