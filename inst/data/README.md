

For TARGET mode, panel-specific information can be used (under $panel_dir).

* Whitelist

| Filename                   | Content                                                        |
|:--------------------------:|:--------------------------------------------------------------:|
| known.gene-gene.txt | A list of known significant fusion between two genes for targeted output even evidence is weak e.g. low number of reads and out-of-frame. |
| known.gene-exon.txt | By default, SplitFusion only outputs fusions that involve two different genes. This file contains clinically relevant isoforms of a gene, e.g. "MET_exon13---MET_exon15" that is an important theraputic target. Many exon skipping/alternative splicing events are normal or of unknown clinical relevance and are thus not output by default. |
| known.partners.txt | Contains a list of known fusion partners of targets. |
| known.3UTR.truncation.txt | Contains a list of known significant genes with 3'UTR trancation, e.g. FGFRs exon18 deletion |


* Blacklist

| Filename                   | Content                                                        |
|:--------------------------:|:--------------------------------------------------------------:|
| filter.gene-gene.txt | Contains recurrent breakpoints identified as data accumulates, but are not of interest. |

The whitelist and blacklist could be updated periodically as a backend supporting database to improve final reporting.


Two files come with SplitFusion are premade reference gene transcripts to resolve some nomenclatural conflicts

| Filename                   | Content                                                        |
|:--------------------------:|:--------------------------------------------------------------:|
| hg19_refGene.mostExon.NM | For genes with multiple isoforms, selecte the one with the most number of exons. hg19 version |
| hg19_refGene0.txt | For genes with multiple isoforms, selecte the one with the smallest NM##, which is often the classic isoform. hg19 version |
| hg38_refGene.mostExon.NM | For genes with multiple isoforms, selecte the one with the most number of exons. hg38 version |
| hg38_refGene0.txt | For genes with multiple isoforms, selecte the one with the smallest NM##, which is often the classic isoform. hg38 version |

