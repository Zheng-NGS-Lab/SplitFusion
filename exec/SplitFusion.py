#!/usr/bin/python
from __future__ import print_function
import sys
import re
import os
import argparse

def parseArgs():

    parser = argparse.ArgumentParser(description='Split-Fusion is a fast data analysis pipeline detects gene fusion based on chimeric split-read alignments.')
    parser.add_argument('--refGenome', required=True
                        , help="The reference genome file, with a full path [required].")
    parser.add_argument('--annovar', required=True
                        , help="The annovar executable file [required].")
    parser.add_argument('--samtools', required=True
                        , help="The samtools executable file [required].")
    parser.add_argument('--bedtools', required=True
                        , help="The bedtools executable file [required].")
    parser.add_argument('--bwa', required=True
                        , help="The bwa executable file [required].")
    parser.add_argument('--R', required=True
                        , help="The R executable file [required].")
    parser.add_argument('--perl', required=True
                        , help="The perl executable file [required].")
    parser.add_argument('--output', required=True
                        , help="The directory for output SplitFusion results [required].")
    parser.add_argument('--sample_id', required=True
                        , help="The name of sample to be analyzed [required].")
    parser.add_argument('--bam_dir', required=False
                        , help="The path to the bam file to be analyzed. The Kickstart mode will use the bam file ('$sample_id'.bam or '$sample_id'.consolidated.bam) in this directory [Either fastq_dir or bam_dir should be specified].")
    parser.add_argument('--fastq_dir', required=False
                        , help="The path to the fastq file to be analyzed [Either fastq_dir or bam_dir should be specified].")
    parser.add_argument('--r1filename', required=False
                        , help="Read 1 fastq filename. Can be in gzipped format. If not specified, $fastq_dir/$sample_id.R1.fq will be used [optional]")
    parser.add_argument('--r2filename', required=False
                        , help="Read 2 fastq filename. Can be in gzipped format. If not specified, $fastq_dir/$sample_id.R2.fq will be used [optional].")
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
    output=[str(k) + "=" + "\"" + str(v)+ "\"" for k,v in args.iteritems() if v != None]
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
