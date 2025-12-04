configfile: "config.yaml"

rule all:
    input:
        f"{config['out_dir']}/consensus_peak_calling/bed_paths.tsv",
        f"{config['out_dir']}/consensus_peak_calling/bw_paths.tsv",
	f"{config['out_dir']}/consensus_peak_calling/macs2_peaks",
	f"{config['out_dir']}/consensus_peak_calling/consensus_peaks",
        config['out_dir'] + "/QC/" + config['modified_tss_bed'],
        config['qc_cmd'],
        config['out_dir'] + "/fragments_dict.pkl",
        expand("{sample}.cbs_for_otsu_thresholds.tsv", sample=config["samples"]),
        expand("{sample}.fragments_stats_per_cb_for_otsu_thresholds.parquet", sample=config["samples"]),
        expand("{sample}.fragments_stats_per_cb_for_otsu_thresholds.tsv", sample=config["samples"]),
        expand("{sample}.fragments_stats_per_cb.parquet", sample=config["samples"]),
        expand("{sample}.fragments_insert_size_dist.parquet", sample=config["samples"]),
        expand("{sample}.otsu_thresholds.tsv", sample=config["samples"]),
        expand("{sample}.tss_norm_matrix_per_cb.parquet", sample=config["samples"]),
        expand("{sample}.tss_norm_matrix_sample.parquet", sample=config["samples"]),
        expand("{sample}.pycistopic_qc.log", sample=config["samples"]),
        config['out_dir'] + "/QC/qc_barcodes_thresholds.pkl",
        config['out_dir'] + "/QC/qc_barcodes_thresholds.pkl", 
        config['out_dir'] + "/cistopic_objects_mm10.pkl", 

rule pseudobulk:
    output:
        bed_paths = f"{config['out_dir']}/consensus_peak_calling/bed_paths.tsv",
        bw_paths  = f"{config['out_dir']}/consensus_peak_calling/bw_paths.tsv"
    shell:
        """
        python pseudobulk.py \
            --ctrl_fragments {config['ctrl_fragments']} \
            --trt_fragments {config['trt_fragments']} \
            --chrom_sizes {config['chrom_sizes']} \
            --ctrl_name {config['ctrl_name']} \
            --trt_name {config['trt_name']} \
            --out_dir {config['out_dir']} \
            --n_cpu {config['n_cpu']} \
            --normalize_bigwig \
            --temp_dir /tmp
        """


rule peak_calling:
    input:
        bed_dir=f"{config['out_dir']}/consensus_peak_calling/pseudobulk_bed_files"
    output:
        directory(f"{config['out_dir']}/consensus_peak_calling/macs2_peaks")
    shell:
        """
        python peak_calling.py \
            -i {input.bed_dir} \
            -o {output} \
            -g {config[g]}
        """


rule consensus_peaks:
    input:
         f"{config['out_dir']}/consensus_peak_calling/macs2_peaks" 
    output:
        directory(f"{config['out_dir']}/consensus_peak_calling/consensus_peaks")
    params:
        combined_bed = config['combined_bed'],
        consensus_bed = config['consensus_bed']
    shell:
        """
        python consensus_peaks.py \
            -i {input} \
            -o {output} \
            -c {params.combined_bed} \
            -p {params.consensus_bed}
        """



rule gtf_to_tss:
    input:
        config['modified_gtf']
    output:
        config['out_dir'] + "/QC/" + config['modified_tss_bed']
    shell:
        """
        mkdir -p $(dirname {output})
        bash gtf_to_tss.sh {input} {output}
        """




rule prepQC: 
    input:
          config['out_dir'] + "/QC/" + config['modified_tss_bed'] 
    params: 
        control = config['ctrl_fragments'],
        treatment = config['trt_fragments'],
        outdir = config['out_dir'],
    output:
         qc_cmd = config['qc_cmd'],
         pkl = config['out_dir'] + "/fragments_dict.pkl",  
    shell: 
      """
      python prepQC.py \
       --out_dir {params.outdir} \
       --consensus_dir  {params.outdir}/consensus_peak_calling/consensus_peaks \
       --tss_bed {input} \
       --th1_fragments {params.control} \
       --th2_fragments {params.treatment} \
       --qc_commands_filename {output.qc_cmd} 
      """

rule runQC:
    input:
        config['qc_cmd']
    params:
        tss_bed=config['out_dir'] + "/QC/" + config['modified_tss_bed']
    output:
       expand("{sample}.cbs_for_otsu_thresholds.tsv", sample=config["samples"]),
       expand("{sample}.fragments_stats_per_cb_for_otsu_thresholds.parquet", sample=config["samples"]),
       expand("{sample}.fragments_stats_per_cb_for_otsu_thresholds.tsv", sample=config["samples"]),
       expand("{sample}.fragments_stats_per_cb.parquet", sample=config["samples"]),
       expand("{sample}.fragments_insert_size_dist.parquet", sample=config["samples"]),
       expand("{sample}.otsu_thresholds.tsv", sample=config["samples"]),
       expand("{sample}.tss_norm_matrix_per_cb.parquet", sample=config["samples"]),
       expand("{sample}.tss_norm_matrix_sample.parquet", sample=config["samples"]),
       expand("{sample}.pycistopic_qc.log", sample=config["samples"])
    shell:
        """
        bash {input}
        """


rule collect_barcodes: 
         input: 
           config['out_dir'] + "/fragments_dict.pkl",
         params: 
           qc_dir= config['out_dir'] + "/QC/", 
           tss_cutoff = config['minTSS'],
           min_frag = config['minFrag'],
           min_frip = config['minFrip'] 
         output: 
            config['out_dir'] + "/QC/qc_barcodes_thresholds.pkl"
         shell:
            """
               python collect_qc_barcodes.py \
               --fragments_dict {input} \
               --qc_output_dir {params.qc_dir} \
               --output_pickle {output} \
               --unique_fragments_threshold {params.tss_cutoff} \
               --tss_enrichment_threshold {params.min_frag} \
               --frip_threshold {params.min_frip}
            """ 




rule plotQC: 
     input:
           config['out_dir'] + "/fragments_dict.pkl",
     params:
           qc_dir= config['out_dir'] + "/QC/",
     output: 
         config['out_dir'] + "/QC/qc_barcodes_thresholds.pkl"
     shell: 
        """  
        python plot_pycistopicQC.py \
           --fragments_dict {input} \
           --qc_output_dir {params.qc_dir} \
           --plots_output_dir {params.qc_dir} \
           --qc_results_pickle {output} \
           --barcode_plots_output_dir {params.qc_dir}
     """


rule create_cisObject: 
     input:
        config['out_dir'] + "/fragments_dict.pkl",
        config['out_dir'] + "/QC/qc_barcodes_thresholds.pkl",
     params: 
         blacklist = config['blacklist'],
         qc_dir= config['out_dir'] + "/QC/",
         regions_bed = config['out_dir'] + "/consensus_peak_calling/consensus_peaks",
         nCPUs = config['n_cpu']
     output: 
         config['out_dir'] + "/cistopic_objects_mm10.pkl"
     shell: 
          """ 
          python create_cistopic_objects.py \
          --fragments_dict {input[0]} \
          --qc_results_pickle {input[1]} \
          --regions_bed {params.regions_bed} \
          --blacklist_bed {params.blacklist} \
          --qc_output_dir {params.qc_dir}\
          --output_pickle {output} \
          --n_cpu {params.nCPUs}
       """ 
