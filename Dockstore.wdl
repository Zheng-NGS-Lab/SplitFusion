version 1.0

task splitfusion {
  input {
    String fastq1
    String fastq1
    String sample
    String outdir
    String genome
    String thread
  }

  command {
    python /home/prod/exec/SplitFusion.py --bwa bwa --R R --samtools samtools --bedtools bedtools --perl perl --fastq_file1 ${fastq1} --fastq_file2 ${fastq2} --panel_dir /home/prod/inst/data/panel --annovar ${annovar} --output ${outdir} --sample_id ${sample} --refGenome ${genome} --thread ${thread}
  }
  output {
    File response = stdout()
  }
  runtime {
   docker: 'ymcki/splitfusion:0.4'
  }
}

workflow TCGA {
  call splitfusion
}
