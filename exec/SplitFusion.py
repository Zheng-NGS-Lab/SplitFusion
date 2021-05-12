#!/usr/bin/python
from __future__ import print_function
import sys
import re
import os
import argparse
from future.utils import viewitems

def parseArgs():

    parser = argparse.ArgumentParser(description='Split-Fusion is a fast data analysis pipeline detects gene fusion based on chimeric split-read alignments.')
    parser.add_argument('--refGenome', required=True
                        , help="The reference genome file, with a full path [required].")
    parser.add_argument('--annovar', required=True
                        , help="The annovar executable file with full path [required].")
    parser.add_argument('--output', required=True
                        , help="The directory for output SplitFusion results [required].")
    parser.add_argument('--sample_id', required=True
                        , help="The name of sample to be analyzed [required].")
    parser.add_argument('--bam_file', required=False
                        , help="If the bam_file is specified, the Kickstart mode will be used. [Either fastq_file or bam_file should be specified]")
    parser.add_argument('--fastq_file1', required=False
                        , help="The fastq file (Read 1 of paired-end) to be analyzed [Either fastq_file or bam_file should be specified].")
    parser.add_argument('--fastq_file2', required=False
                        , help="Read 2 of paired-end fastq file.")
    parser.add_argument('--samtools', type=str
                        , default='samtools'
                        , help="The samtools executable file with full path [Optional].")
    parser.add_argument('--bedtools', type=str
                        , default='bedtools'
                        , help="The bedtools executable file with full path [Optional].")
    parser.add_argument('--bwa', type=str
                        , default='bwa'
                        , help="The bwa executable file with full path [Optional].")
    parser.add_argument('--R', type=str
                        , default='R'
                        , help="The R executable file with full path [Optional].")
    parser.add_argument('--perl', type=str
                        , help="The perl executable file with full path [Optional]."
                        , default='perl')
    parser.add_argument('--panel_dir', required=False
                        , default='NA'
                        , help="For TARGET mode: the path where known significant fusions or splicing isoforms (whitelist) or unwanted fusions involving homologous genes or recurrent falsed positives (blacklist) are stored. Default='NA'")
    parser.add_argument('--panel', required=False
                        , default='NA'
                        , help="The prefix name of TARGET gene panel file. E.g., 'LungFusion' for LungFusion.GSP2.bed. Default='NA'")
    parser.add_argument('--steps', required=False
                        , default='1_fastq-bam,2_bam-breakpoint,3_breakpoint-filter,4_breakpoint-anno,5_breakpoint-anno-post'
			, help="Specify steps to run. Default='1_fastq-bam,2_bam-breakpoint,3_breakpoint-filter,4_breakpoint-anno,5_breakpoint-anno-post'")
    parser.add_argument('--AnnotationMethod', required=False
                        , default = 'annovar'
                        , help="the name of annotation tools. Default = 'annovar'")
    parser.add_argument('--thread', type=int
                        , default=1
                        , help="number of threads for parallel computing. Default=1")
    parser.add_argument('--genomeVer', type=str
                        , default="hg38"
                        , help="genome version for annotation software. Default=hg38")
    parser.add_argument('--minMQ', type=int
                        , default=13
                        , help="minimum mapping quality for all split alignments (both Ligation and Anchored ends). Default=13")
    parser.add_argument('--minMQ1', type=int
                        , default=30
                        , help="minimum mapping quality of the leftmost of Read1 (Ligation end). Default=30")
    parser.add_argument('--minMapLength', type=int
                        , default=18
                        , help="minimum read mapping length for all split alignments (both Ligation and Anchored ends). Default=18")
    parser.add_argument('--minMapLength2', type=int
                        , default=25
                        , help="minimum mapping length of the leftmost of Read1 (Ligation end). Default=25")
    parser.add_argument('--maxQueryGap', type=int
                        , default=0
                        , help="maximum gap length between split alignments on a query read. Default=0")
    parser.add_argument('--maxOverlap', type=int
                        , default=6
                        , help="maximum overlapping bases of two split alignments on a query read. Default=6")
    parser.add_argument('--minExclusive', type=int
                        , default=18
                        , help="minimum exclusive length between two split alignments. Default=18")
    parser.add_argument('--FusionMinStartSite', type=int
                        , default=1
                        , help="minimum number of fusion partner ends (ligation end) to call CANDIDATE structure variation/fusion. Should be less or equal minPartnerEnds_BothExonJunction. Default=1")
    parser.add_argument('--minPartnerEnds_BothExonJunction', type=int
                        , default=1
                        , help="minimum number of fusion partner ends (ligation end), when both breakpoints are at exon boundaries/junctions, in the final call of structure variation/fusion. Default=1")
    parser.add_argument('--minPartnerEnds_OneExonJunction', type=int
                        , default=3
                        , help="minimum number of fusion partner ends (ligation end), when only one breakpoint is at exon boundary/junction, in the final call of structure variation/fusion. Default=3")

    args = vars(parser.parse_args())


    return args

def mkdir(path):
 
	folder = os.path.exists(path)
 
	if not folder:
		os.makedirs(path)

#os.system(cleanup)
if __name__ == '__main__':
    
    args = parseArgs()
    output=[str(k) + "=" + "\"" + str(v)+ "\"" for k,v in viewitems(args) if v != None]
    config_p = args['output'] + "/" + args['sample_id'] + "/" 
    mkdir(config_p)
    config_o = config_p + "config.txt"
    config = open(config_o, "w+")
    print("\n".join(output), file = config, sep="\n")

    config.close()

#design =  args['R'] +  " -e " + "'suppressMessages(library(SplitFusion)); runSplitFusion(configFile = " + "\"" + config_o + "\"" + ")'" #+ "> /dev/null 2>&1"
cmd =  "cd " + config_p + "; " + args['R'] +  " -e " + "'suppressMessages(library(SplitFusion)); runSplitFusion()'"
os.system(cmd)

## END
