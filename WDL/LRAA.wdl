version 1.0

workflow lraa_wf {

    input {
        File genome_fa
        File aligned_reads_bam
        String sample_id
        File? gtf_file
        Boolean quant_only

        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        Int cpu = 4
        String memory = "32G"
        Int disk_space = "100"
        Int preemptible = 0
        
    }

    call lraa_task {
        input:
            genome_fa=genome_fa,
            aligned_reads_bam=aligned_reads_bam,
            sample_id=sample_id,
            gtf_file=gtf_file,
            quant_only=quant_only,
            docker=docker,
            cpu=cpu,
            memory=memory,
            disk_space=disk_space,
            preemptible=preemptible
    }


    output {
          File? LRAA_gtf = lraa_task.LRAA_gtf
          File? LRAA_quant_expr = lraa_task.LRAA_quant_expr
          File? LRAA_quant_tracking = lraa_task.LRAA_quant_tracking
    }
        
}



task lraa_task {
    
    input {
        File genome_fa
        File aligned_reads_bam
        String sample_id
        File? gtf_file
        Boolean quant_only

        String docker
        Int cpu
        String memory
        Int disk_space
        Int preemptible
  
    }

    String output_prefix = if defined (gtf_file) then ".refGuided" else ".refFree"

    
    command <<<
        
        set -ex


        out_prefix = ~{output_prefix}
        if [[ quant_only == "true" ]]; then
            out_prefix=""
        fi
        
        /usr/local/src/LRAA/LRAA --genome ~genome_fa \
                                 --bam ~aligned_reads_bam \
                                 --output_prefix ~{sample_id}.LRAA \
                                 ~{true='--quant_only' false='' quant_only} \
                                 ~{"--gtf " + gtf_file} \
                                 --output_prefix ~{sample_id}.LRAA${out_prefix}

      >>>


      output {
          File? LRAA_gtf = "~{sample_id}.LRAA~{output_prefix}.gtf"
          File? LRAA_quant_expr = "~{sample_id}.LRAA.quant.expr"
          File? LRAA_quant_tracking = "~{sample_id}.LRAA.quant.tracking"
      }

      runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + disk_space + " HDD"
        cpu: cpu
        preemptible: preemptible
     }
      
}
