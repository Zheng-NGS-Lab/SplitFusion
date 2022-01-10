#### Pipeline for running SplitFusino for kinase gene fusion detection
### https://github.com/Zheng-NGS-Lab/SplitFusion

task samtofastq {

    File input_bam
    String prefix

    Boolean? err_on_pe_and_unpaired=true

    String? picard_jar="/opt/picard-tools/picard.jar"
    String? base_args="INCLUDE_NON_PF_READS=true INCLUDE_NON_PRIMARY_ALIGNMENTS=false VALIDATION_STRINGENCY=SILENT"
    String? other_args=""

    Int? memory="7"
    Float java_memory_subtract = "1"
    Float java_memory_scaler = "1"
    Int java_memory = floor(memory * java_memory_scaler - java_memory_subtract)
    Int? num_threads="1"
    Int? disk_buffer_gb = "50"
    Int? disk_scaler = "3"
    Int disk_space = ceil(size(input_bam, "GB") * disk_scaler + disk_buffer_gb)
    Int? num_preempt="4"
    String docker="gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"

    File monitor_script

    command {
        set -euo pipefail
        chmod +x ${monitor_script}; ${monitor_script} > monitoring.log &

        java -jar "-Xmx${java_memory}g" ${picard_jar} SamToFastq I=${input_bam} F=${prefix}_1.fastq.gz F2=${prefix}_2.fastq.gz FU=${prefix}_unpaired.fastq.gz ${base_args} ${other_args}

        #Check if BAM was PE/SE by checking for empty files
        hchar_fq1=`gunzip -c ${prefix}_1.fastq.gz | head | wc -c`
        hchar_fq2=`gunzip -c ${prefix}_2.fastq.gz | head | wc -c`
        hchar_fq_s=`gunzip -c ${prefix}_unpaired.fastq.gz | head | wc -c`

        if [[ $hchar_fq1 -gt 0 && $hchar_fq2 -gt 0 ]]; then #PE
          if [[ $hchar_fq_s -gt 0]]; then
            if [[ ${err_on_pe_and_unpaired} = true ]]; then
              echo 'Error: there are unpaired reads in your input BAM!'; exit 1 # terminate and indicate error
            else
              nrow_fastq_s=`gunzip -c ${prefix}_unpaired.fastq.gz | wc -l`
              echo "Warning: you have $((nrow_fastq_s / 4)) unpaired reads in your input BAM being dropped!"
            fi
          fi
        elif [[ $hchar_fq_s -gt 0 ]]; then #SE (put unpaired in fastq1)
          mv ${prefix}_unpaired.fastq.gz ${prefix}_1.fastq.gz
        else
          echo 'Error: no paired or unpaired reads retrieved from input BAM!'; exit 1 # terminate and indicate error
        fi

    }

    output {
      File fastq1="${prefix}_1.fastq.gz"
      File fastq2="${prefix}_2.fastq.gz"
      File monitoring_log="monitoring.log"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Binyamin A. Knisbacher"
    }
}

task bwa {
  File? fastq1
  File? fastq2
  File? fastq_tar_gz
  String prefix
  File reference_fa
  File reference_index_tar_gz

  String? tar_fastq_suffix1="_1.fastq"
  String? tar_fastq_suffix2="_2.fastq"

  String? other_args="-T 18 -M -q -K 10000000"

  String? docker="docker.io/ymcki/splitfusion@sha256:de93d2fd85816a2a402c5384a77cc1abb45abcfd05d87f1261f32f5bfd92012b"
  Int? disk_buffer_gb = "10"
  Int? disk_scaler = "2"
  Int? disk_space = ceil( (size(fastq1, "GB") + size(fastq2, "GB") + size(fastq_tar_gz, "GB")*10 + size(reference_fa, "GB") + size(reference_index_tar_gz, "GB")) * disk_scaler ) + disk_buffer_gb
  Int? memory="12"
  Int? num_threads="16"
  Int? num_preempt="4"

  File monitor_script

  command {
    set -euo pipefail
    chmod +x ${monitor_script}; ${monitor_script} > monitoring.log &

    tar xvf ${reference_index_tar_gz} -C `dirname ${reference_fa}`

    ### Run per tar.gz or specified files
    if [[ -n "${fastq_tar_gz}" ]]; then
      mkdir fastqs; tar xvvf ${fastq_tar_gz} -C fastqs
      bwa mem ${other_args} -t ${num_threads} ${reference_fa} <(zcat -f fastqs/*${tar_fastq_suffix1}*) <(zcat -f fastqs/*${tar_fastq_suffix2}*) | samtools view -b - > ${prefix}.bam 2> ${prefix}.bwa.log
    elif [[ -n "${fastq1}" ]]; then # will work for PE or SE
      fastq_files="${fastq1} ${fastq2}"
      echo "fastq files being used for bwa-mem: $fastq_files"
      bwa mem ${other_args} -t ${num_threads} ${reference_fa} $fastq_files | samtools view -b - > ${prefix}.bam 2> ${prefix}.bwa.log
    else
      echo "Error: must input fastq1 or fastq_tar_gz, exiting!"; exit 1
    fi
  }

  output {
      File bam="${prefix}.bam"
      File bwa_log="${prefix}.bwa.log"
      File monitoring_log="monitoring.log"
  }

  runtime {
      docker: "${docker}"
      memory: "${memory}GB"
      disks: "local-disk ${disk_space} HDD"
      cpu: "${num_threads}"
      preemptible: "${num_preempt}"
  }

  meta {
      author: "Binyamin A. Knisbacher"
  }
}

task splitfusion {

    String sample_id
    File? fastq1
    File? fastq2
    File? fastq_tar_gz
    File? bam
    File reference_fa
    File reference_index_tar_gz
    File annovar_db_tar_gz

    ## Allow to place/override panel files in default panel location in docker
    #see important parameter --panel to be specified in $other_args
    File? panel_tar_gz

    String? other_args=""

    Boolean? save_extended_output=false

    String? docker="docker.io/ymcki/splitfusion@sha256:de93d2fd85816a2a402c5384a77cc1abb45abcfd05d87f1261f32f5bfd92012b"
    Int? disk_buffer_gb="25"
    Float? disk_scaler="15" #bam is converted to sam and multiple copies are made in SplitFusion
    Int? disk_space=ceil( (size(fastq1,"GB") + size(fastq2,"GB") + size(bam,"GB") + size(fastq_tar_gz,"GB")) * disk_scaler + size(reference_fa, "GB") + size(reference_index_tar_gz, "GB") + size(annovar_db_tar_gz, "GB")*5) + disk_buffer_gb
    Int? num_threads="1"
    Float? memory=8*ceil(size(fastq1,"GB") + size(bam,"GB") + size(fastq_tar_gz,"GB"))
    Int? num_preempt="4"

    File monitor_script

    command {
        set -euo pipefail
        chmod +x ${monitor_script}; ${monitor_script} > monitoring.log &

        ## Unpack bwa index
        echo "Extracting reference index ${reference_index_tar_gz}:"
        reference_dir=`dirname ${reference_fa}`
        tar xvf ${reference_index_tar_gz} -C $reference_dir
        rm ${reference_index_tar_gz}

        ## Allow to place/override panel files in default panel location in docker
        #For the'se to work would need to specify --panel in $other_args
        if [[ -n "${panel_tar_gz}" ]]; then
          echo "Extracting panel ${panel_tar_gz}:"
          tar xvf ${panel_tar_gz} -C /SplitFusion/inst/data/panel
          other_args="$other_args --panel_dir /SplitFusion/inst/data/panel"
        fi

        #unpack tar.gz to annovar directory
        echo "Extracting annovar DB ${annovar_db_tar_gz}:"
        mkdir annovar; tar xvf ${annovar_db_tar_gz} -C annovar --strip-components=1; chmod +x annovar/*.pl

        all_other_args=${other_args}
        if [[ -n "${bam}" ]]; then
          all_other_args="$all_other_args --bam_file ${bam}"
        elif [[ -n  "${fastq_tar_gz}" ]]; then
          mkdir fastqs; tar xvf ${fastq_tar_gz} -C fastqs
          all_other_args="$all_other_args --fastq_file1 $PWD/fastqs/*_1.fastq* --fastq_file2 $PWD/fastqs/*_2.fastq*"
        elif [[ -n  "${fastq1}" ]]; then
          all_other_args="$all_other_args --fastq_file1 ${fastq1}"
          if [[ -n  "${fastq2}" ]]; then
            all_other_args="$all_other_args --fastq_file2 ${fastq2}"
          fi
        else
          echo 'Error: must specify bam or fastq1, exiting!'; exit 1
        fi

        python3 /SplitFusion/exec/SplitFusion.py \
                --refGenome ${reference_fa} \
                --annovar $PWD/annovar \
                --samtools /usr/bin/samtools \
                --bedtools /usr/bin/bedtools \
                --bwa /usr/bin/bwa \
                --R /usr/bin/R \
                --perl /usr/bin/perl \
                --output $PWD \
                --sample_id ${sample_id} \
                --thread ${num_threads} \
                --panel_dir /SplitFusion/inst/data/panel \
                $all_other_args

        ### Organize output
        cp "${sample_id}/${sample_id}.brief.summary" "${sample_id}.brief.summary.tsv"

        ##Save output tar.gz
        # 1) save_extended_output=true: for saving all output from splitfusion (for debugging)
        # 2) save_extended_output=false: for saving all per-fusion event files
        if [ ${save_extended_output} = true ]; then
          tar zcvf ${sample_id}.tar.gz ${sample_id}
          echo '### Printing file tree in extended output mode ###'; ls -l **
        else
          tar zcvf ${sample_id}.tar.gz $( ls ${sample_id}/*brief.summary ${sample_id}/*---* ${sample_id}/*fusion* 2>/dev/null )
        fi


    }

    output {
        File summary_tsv = "${sample_id}.brief.summary.tsv"
        File output_tgz = "${sample_id}.tar.gz"
        File monitoring_log = "monitoring.log"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Binyamin A. Knisbacher"
    }
}


workflow splitfusion_workflow {

    String sample_id
    File? input_bam
    File? fastq1
    File? fastq2
    File? fastq_tar_gz
    File reference_fa
    File reference_index_tar_gz
    File annovar_db_tar_gz
    File? panel_tar_gz
    String? docker="docker.io/ymcki/splitfusion@sha256:de93d2fd85816a2a402c5384a77cc1abb45abcfd05d87f1261f32f5bfd92012b"
    String? samtofastq_docker
    Int? num_preempt

    Boolean? kickstart=false
    String? splitfusion_other_args #do not set defaults here
    String? bwa_other_args #do not set defaults here

    Boolean? save_extended_output #debug mode or save additional files

    File? monitor_script='gs://fc-secure-b7402204-58ce-4566-8450-31016f3489e8/monitor_script.sh'


    if( defined(fastq_tar_gz) ){
      call splitfusion as splitfusion1 {
        input:
            sample_id=sample_id,
            fastq_tar_gz=fastq_tar_gz,
            reference_fa=reference_fa,
            reference_index_tar_gz=reference_index_tar_gz,
            annovar_db_tar_gz=annovar_db_tar_gz,
            panel_tar_gz=panel_tar_gz,
            other_args=splitfusion_other_args,
            save_extended_output=save_extended_output,
            docker=docker,
            num_preempt=num_preempt,
            monitor_script=monitor_script
      }
    }

    #Two options:
    ### 1) PE fastq input
    ### 2) BAM input, run splitfusion on input bam as in kickstart mode
    if( defined(fastq1) || ( defined(input_bam) && kickstart ) ){
      call splitfusion as splitfusion2 {
      input:
            sample_id=sample_id,
            fastq1=fastq1,
            fastq2=fastq2,
            bam=input_bam,
            reference_fa=reference_fa,
            reference_index_tar_gz=reference_index_tar_gz,
            annovar_db_tar_gz=annovar_db_tar_gz,
            panel_tar_gz=panel_tar_gz,
            other_args=splitfusion_other_args,
            save_extended_output=save_extended_output,
            docker=docker,
            num_preempt=num_preempt,
            monitor_script=monitor_script
        }
    }

    ### BAM input, no kickstart (samtofastq -> bwa -> splitfusion)
    if( ! defined(fastq1) && defined(input_bam) && ! kickstart){
      call samtofastq {
        input:
          input_bam=input_bam,
          prefix=sample_id,
          docker=samtofastq_docker,
          num_preempt=num_preempt,
          monitor_script=monitor_script
        }

        call bwa as bwa2 {
          input:
            fastq1=samtofastq.fastq1,
            fastq2=samtofastq.fastq2,
            fastq_tar_gz=fastq_tar_gz,
            prefix=sample_id,
            reference_fa=reference_fa,
            reference_index_tar_gz=reference_index_tar_gz,
            other_args=bwa_other_args,
            docker=docker,
            num_preempt=num_preempt,
            monitor_script=monitor_script
        }

        call splitfusion as splitfusion3 {
          input:
            sample_id=sample_id,
            fastq1=fastq1,
            fastq2=fastq2,
            bam=bwa2.bam,
            reference_fa=reference_fa,
            reference_index_tar_gz=reference_index_tar_gz,
            annovar_db_tar_gz=annovar_db_tar_gz,
            panel_tar_gz=panel_tar_gz,
            other_args=splitfusion_other_args,
            save_extended_output=save_extended_output,
            docker=docker,
            num_preempt=num_preempt,
            monitor_script=monitor_script
        }
    }

    output{
      File summary_tsv = select_first([splitfusion1.summary_tsv, splitfusion2.summary_tsv, splitfusion3.summary_tsv])
      File output_tgz = select_first([splitfusion1.output_tgz, splitfusion2.output_tgz, splitfusion3.output_tgz])
      Array[File] monitoring_logs = select_all([splitfusion1.monitoring_log, samtofastq.monitoring_log, bwa2.monitoring_log, splitfusion2.monitoring_log, splitfusion3.monitoring_log])
    }

}

