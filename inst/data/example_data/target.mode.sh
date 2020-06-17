python2 ../../../exec/SplitFusion.py \
	--refGenome Homo_sapiens_assembly19.fasta \
	--database_dir /media/storage1/tempdir/zhangbaifeng/database \
        --annovar /media/storage1/tempdir/zhangbaifeng/tool/annovar \
        --samtools /media/storage1/tempdir/zhangbaifeng/tool/samtools \
        --bedtools /media/storage1/tempdir/zhangbaifeng/tool/bedtools \
        --bwa /media/storage1/tempdir/zhangbaifeng/tool/bwa \
        --R /media/storage1/tempdir/zhangbaifeng/tool/R \
        --perl /media/storage1/tempdir/zhangbaifeng/tool/perl \
        --output /media/storage1/tempdir/zhangbaifeng/software/SplitFusion/inst/data/example_data/target_mode_result \
        --sample_id "Lib001" \
        --fastq_dir /media/storage1/tempdir/zhangbaifeng/software/SplitFusion/inst/data/example_data \
        --r1filename "Lib001".R1.fq \
        --r2filename "Lib001".R2.fq \
	--panel_dir /media/storage1/tempdir/zhangbaifeng/software/SplitFusion/inst/data/panel \
        --panel LungFusion \
        --thread 6 &

