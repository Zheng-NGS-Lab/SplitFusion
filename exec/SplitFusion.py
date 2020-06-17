#!/usr/bin/python
from __future__ import print_function
import sys
import re
import os
import argparse

def parseArgs():

    parser = argparse.ArgumentParser(description='Split-Fusion is a \
                    fast data analysis pipeline detects gene fusion based \
                    on split reads and/or paired-end reads.')
    parser.add_argument('--refGenome', required=True
                        , help="[required]. the path where human genome reference is stored")
    parser.add_argument('--database_dir', required=True
                        , help="[required]. the path where large databases e.g. reference genome and annotation databases are stored")
    parser.add_argument('--annovar', required=True
                        , help="[required]. the path of annovar")
    parser.add_argument('--samtools', required=True
                        , help="[required]. the path of samtools")
    parser.add_argument('--bedtools', required=True
                        , help="[required]. the path of bedtools")
    parser.add_argument('--bwa', required=True
                        , help="[required]. the path of bwa")
    parser.add_argument('--R', required=True
                        , help="[required]. the path of R")
    parser.add_argument('--perl', required=True
                        , help="[required]. the path of perl")
    parser.add_argument('--output', required=True
                        , help="[required]. the path where output is stored")
    parser.add_argument('--sample_id', required=True
                        , help="[required]. the sample name of running")
    parser.add_argument('--bam_dir', required=False
                        , help="the path where bam or fastq file is stored. [Kickstart] The bam file of the sameple_id (xxx.bam or xxx.consolidated.bam) will be used. Either fastq_dir or bam_dir should be specified")
    parser.add_argument('--fastq_dir', required=False
                        , help="the path where fastq file is stored. Either fastq_dir or bam_dir should be specified")
    parser.add_argument('--r1filename', required=False
                        , help="Read 1 fastq filename. Can be in gzipped format. If not specified, $fastq_dir/$sample_id.R1.fq will be used.")
    parser.add_argument('--r2filename', required=False
                        , help="Read 2 fastq filename. Can be in gzipped format. If not specified, $fastq_dir/$sample_id.R2.fq will be used.")
    parser.add_argument('--panel_dir', required=False
                        , default='NA'
                        , help="For Target mode: the path where known significant fusions or splicing isoforms (white.list) or unwanted fusions involving homologous genes or recurrent falsed positives (black.list) are stored. default='NA'")
    parser.add_argument('--panel', required=False
                        , default='NA'
                        , help="the prefix name of target genes panel file is stored, e.g., LungFusion for LungFusion.GSP2.bed. default='NA'")
    parser.add_argument('--steps', required=False
                        , default='1_fastq-bam,2_bam-breakpoint,3_breakpoint-filter,4_breakpoint-anno,5_breakpoint-anno-post'
			, help="specify steps to run. default='1_fastq-bam,2_bam-breakpoint,3_breakpoint-filter,4_breakpoint-anno,5_breakpoint-anno-post'")
    parser.add_argument('--AnnotationMethod', required=False
                        , default = 'annovar'
                        , help="the name of annotation tools. default = 'annovar'")
    parser.add_argument('--thread', type=int
                        , default=4
                        , help="number of threads for computing. default=4")
    parser.add_argument('--minMQ', type=int
                        , default=13
                        , help="minimum mapping quality. default=13")
    parser.add_argument('--minMQ1', type=int
                        , default=30
                        , help="minimum mapping quality of a leftmost of Read1 (rightmost of Read2). default=30")
    parser.add_argument('--minMapLength', type=int
                        , default=18
                        , help="minimum read mapping length. default=18")
    parser.add_argument('--minMapLength2', type=int
                        , default=25
                        , help="minimum mapping length of rightmost of Read1 (leftmost of Read2). default=25")
    parser.add_argument('--maxQueryGap', type=int
                        , default=0
                        , help="maximum gap length on a query read of split alignments. default=0")
    parser.add_argument('--maxOverlap', type=int
                        , default=6
                        , help="maximum overlap bases of two split alignments. default=6")
    parser.add_argument('--minExclusive', type=int
                        , default=18
                        , help="minimum exclusive length between two split alignments. default=18")
    parser.add_argument('--FusionMinStartSite', type=int
                        , default=1
                        , help="minimum number of Adaptor Ligation Read Starting Sites to call Structure Variation/Fusion. Should be less or equal minPartnerEnds_BothExonJunction. default=1")
    parser.add_argument('--minPartnerEnds_BothExonJunction', type=int
                        , default=1
                        , help="minimum number of fusion partner ends (ligation site), when both breakpoints are at exon junctions, to call Structure Variation/Fusion. default=1")
    parser.add_argument('--minPartnerEnds_OneExonJunction', type=int
                        , default=3
                        , help="minimum number of fusion partner ends (ligation site), when one breakpoint is at exon junction, to call Structure Variation/Fusion. default=3")

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
