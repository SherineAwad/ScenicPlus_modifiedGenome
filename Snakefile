configfile: "config.yaml"

rule all:
    input:
        f"{config['out_dir']}/consensus_peak_calling/bed_paths.tsv",
        f"{config['out_dir']}/consensus_peak_calling/bw_paths.tsv",
	f"{config['out_dir']}/consensus_peak_calling/macs2_peaks",
        f"{config['out_dir']}/consensus_peak_calling/consensus_peaks/{config['combined_bed']}",
        f"{config['out_dir']}/consensus_peak_calling/consensus_peaks/{config['consensus_bed']}",
        config['out_dir'] + "/QC/" + config['modified_tss_bed'],
        config['qc_cmd'],
        config['out_dir'] + "/fragments_dict.pkl",
        expand("scenicOuts/QC/{sample}.cbs_for_otsu_thresholds.tsv", sample=config["samples"]),
        expand("scenicOuts/QC/{sample}.fragments_stats_per_cb_for_otsu_thresholds.parquet", sample=config["samples"]),
        expand("scenicOuts/QC/{sample}.fragments_stats_per_cb_for_otsu_thresholds.tsv", sample=config["samples"]),
        expand("scenicOuts/QC/{sample}.fragments_stats_per_cb.parquet", sample=config["samples"]),
        expand("scenicOuts/QC/{sample}.fragments_insert_size_dist.parquet", sample=config["samples"]),
        expand("scenicOuts/QC/{sample}.otsu_thresholds.tsv", sample=config["samples"]),
        expand("scenicOuts/QC/{sample}.tss_norm_matrix_per_cb.parquet", sample=config["samples"]),
        expand("scenicOuts/QC/{sample}.tss_norm_matrix_sample.parquet", sample=config["samples"]),
        expand("scenicOuts/QC/{sample}.pycistopic_qc.log", sample=config["samples"]),
        config['out_dir'] + "/QC/qc_barcodes_thresholds.pkl",
        config['out_dir'] + "/cistopic_objects_mm10.pkl", 
        config['out_dir'] + "/merged_cistopic.pkl",
        f"{config['PROJ_NAME']}.h5ad",
        f"clustered_{config['PROJ_NAME']}.h5ad", 
        f"annotated_clustered_{config['PROJ_NAME']}.h5ad", 
        config['RNA_Barcodes'], 
        config['out_dir'] + "/merged_with_meta.pkl", 
        MALLET = config['out_dir'] +"/MALLET", 
rule pseudobulk:
    output:
        bed_paths = f"{config['out_dir']}/consensus_peak_calling/bed_paths.tsv",
        bw_paths  = f"{config['out_dir']}/consensus_peak_calling/bw_paths.tsv"
    params: 
        ctrl_frag = config['ctrl_fragments'], 
        trt_frag = config['trt_fragments'],
        ctrl_name = config['ctrl_name'],
        trt_name = config['trt_name'],
        chrom_sizes = config['chrom_sizes'],
        out_dir = config['out_dir'], 
        n_cpu = config['n_cpu'],
        
    shell:
        """
        python src/pseudobulk.py \
            --ctrl_fragments {params.ctrl_frag} \
            --trt_fragments {params.trt_frag} \
            --chrom_sizes {params.chrom_sizes} \
            --ctrl_name {params.ctrl_name} \
            --trt_name {params.trt_name} \
            --out_dir {params.out_dir} \
            --n_cpu {params.n_cpu} \
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
        python src/peak_calling.py \
            -i {input.bed_dir} \
            -o {output} \
            -g {config[g]}
        """


rule consensus_peaks:
    input:
         f"{config['out_dir']}/consensus_peak_calling/macs2_peaks" 
    output:
        combined_bed = f"{config['out_dir']}/consensus_peak_calling/consensus_peaks/{config['combined_bed']}",
        consensus_bed = f"{config['out_dir']}/consensus_peak_calling/consensus_peaks/{config['consensus_bed']}"
    shell:
        """
        python src/consensus_peaks.py \
            -i {input} \
            -o {output} \
            -c {output.combined_bed} \
            -p {output.consensus_bed}
        """


rule gtf_to_tss:
    input:
        config['modified_gtf']
    output:
        config['out_dir'] + "/QC/" + config['modified_tss_bed']
    shell:
        """
        mkdir -p $(dirname {output})
        bash src/gtf_to_tss.sh {input} {output}
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
      python src/prepQC.py \
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
        tss_bed = config['out_dir'] + "/QC/" + config['modified_tss_bed']
    output:
        expand("scenicOuts/QC/{sample}.cbs_for_otsu_thresholds.tsv", sample=config["samples"]),
        expand("scenicOuts/QC/{sample}.fragments_stats_per_cb_for_otsu_thresholds.parquet", sample=config["samples"]),
        expand("scenicOuts/QC/{sample}.fragments_stats_per_cb_for_otsu_thresholds.tsv", sample=config["samples"]),
        expand("scenicOuts/QC/{sample}.fragments_stats_per_cb.parquet", sample=config["samples"]),
        expand("scenicOuts/QC/{sample}.fragments_insert_size_dist.parquet", sample=config["samples"]),
        expand("scenicOuts/QC/{sample}.otsu_thresholds.tsv", sample=config["samples"]),
        expand("scenicOuts/QC/{sample}.tss_norm_matrix_per_cb.parquet", sample=config["samples"]),
        expand("scenicOuts/QC/{sample}.tss_norm_matrix_sample.parquet", sample=config["samples"]),
        expand("scenicOuts/QC/{sample}.pycistopic_qc.log", sample=config["samples"])
    shell:
        """
        bash {input}
        """


rule collect_barcodes: 
         input: 
           config['out_dir'] + "/fragments_dict.pkl",
         params: 
           qc_dir= config['out_dir'] + "/QC/", 
           uniq_frag = config['uniqFrag'],
           min_tss = config['minTSS'],
           min_frip = config['minFrip'] 
         output: 
            config['out_dir'] + "/QC/qc_barcodes_thresholds.pkl"
         shell:
            """
               python src/collect_qc_barcodes.py \
                   --fragments_dict {input} \
                   --qc_output_dir {params.qc_dir} \
                   --output_pickle {output} \
                   --unique_fragments_threshold {params.uniq_frag} \
                   --tss_enrichment_threshold {params.min_tss} \
                   --frip_threshold {params.min_frip}
              python src/plotQC.py \
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
        consensus_bed = f"{config['out_dir']}/consensus_peak_calling/consensus_peaks/{config['consensus_bed']}"
     params: 
         blacklist = config['blacklist'],
         qc_dir= config['out_dir'] + "/QC/",
         nCPUs = config['n_cpu']
     output: 
         config['out_dir'] + "/cistopic_objects_mm10.pkl"
     shell: 
          """ 
          python src/create_cistopic_objects.py \
          --fragments_dict {input[0]} \
          --qc_results_pickle {input[1]} \
          --regions_bed {input.consensus_bed} \
          --blacklist_bed {params.blacklist} \
          --qc_output_dir {params.qc_dir}\
          --output_pickle {output} \
          --n_cpu {params.nCPUs}
       """


rule mergeCisObj: 
     input: 
         config['out_dir'] + "/cistopic_objects_mm10.pkl"
     output: 
        config['out_dir'] + "/merged_cistopic.pkl"
     shell: 
         """python src/merge_cistopic.py {input} {output}"""

rule prepRNA: 
    input:  
       expand("{sample}_filtered_feature_bc_matrix.h5", sample=config["samples"])
    params: 
       samples_file = config['samples_file'],
       project = config['PROJ_NAME'] 
    output: 
        f"{config['PROJ_NAME']}.h5ad"
    shell: 
         """python src/preprocess_scRNA.py {params.project} {params.samples_file}"""


rule clusterRNA: 
     input: 
        f"{config['PROJ_NAME']}.h5ad"
     params: 
        markers = config['markers']
     output: 
        f"clustered_{config['PROJ_NAME']}.h5ad"
     shell: 
        """
        python src/cluster_scRNA.py {input} {params.markers}
        """  

rule annotateRNA:
     input: 
       f"clustered_{config['PROJ_NAME']}.h5ad"
     output: 
       f"annotated_clustered_{config['PROJ_NAME']}.h5ad"
     params: 
         config['annotations_file'] 
     shell: 
         """ 
         python src/annotate_scRNA.py {input} {params} 
         """ 
rule barcodesRNA: 
     input: 
         f"annotated_clustered_{config['PROJ_NAME']}.h5ad"
     output:
        config['RNA_Barcodes']
     shell: 
       """
       python src/barcodes_scRNA.py {input} {output} 
       """ 


rule metaRNA: 
    input: 
       config['RNA_Barcodes'],
       config['out_dir'] + "/merged_cistopic.pkl"
    output: 
       config['out_dir'] + "/merged_with_meta.pkl"
    shell: 
      """
      python src/add_scrna_metadata.py \
        --cistopic_pickle {input[1]} \
        --scrna_csv {input[0]} \
        --output_pickle {output}
      """


rule runMallet: 
     input: 
         config['out_dir'] + "/merged_with_meta.pkl"
     params: 
        tmp = config['out_dir'] +"/TMP",
        MALLET_PATH = config['MALLET_PATH'],
        nCPU = config['n_cpu'], 
     output:
        MALLET = config['out_dir'] +"/MALLET",
     shell: 
       """
       python src/run_mallet.py \
         --cistopic_obj_pickle {input} \
         --mallet_path {params.MALLET_PATH} \
         --n_topics 5  \
         --n_cpu {params.nCPU} \
         --n_iter 500 \
         --tmp_path {params.tmp} \
         --save_path {output} \
         --mallet_memory 80G \
         --random_state 555 \
         --alpha 0.1 \
         --alpha_by_topic \
         --eta 0.01 \
         --eta_by_topic
     """ 
