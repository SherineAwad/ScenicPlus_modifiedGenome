configfile: "config.yaml"

rule all:
    input:
        f"{config['out_dir']}/consensus_peak_calling/bed_paths.tsv",
        f"{config['out_dir']}/consensus_peak_calling/bw_paths.tsv",
	f"{config['out_dir']}/consensus_peak_calling/macs2_peaks",
	f"{config['out_dir']}/consensus_peak_calling/consensus_peaks",
        config['out_dir'] + "/QC/" + config['modified_tss_bed'],
	config['qc_cmd'], 
        "done.txt"
rule pseudobulk:
    output:
        bed_paths = f"{config['out_dir']}/consensus_peak_calling/bed_paths.tsv",
        bw_paths  = f"{config['out_dir']}/consensus_peak_calling/bw_paths.tsv"
    shell:
        """
        python pseudobulk.py \
            --ctrl_fragments {config[ctrl_fragments]} \
            --trt_fragments {config[trt_fragments]} \
            --chrom_sizes {config[chrom_sizes]} \
            --ctrl_name {config[ctrl_name]} \
            --trt_name {config[trt_name]} \
            --out_dir {config[out_dir]} \
            --n_cpu {config[n_cpu]} \
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
        directory(f"{config['out_dir']}/consensus_peak_calling/macs2_peaks")
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
         config['qc_cmd']
    shell: 
      """
      python prepQC.py \
       --out_dir {params.outdir} \
       --consensus_dir  {params.outdir}/consensus_peak_calling/consensus_peaks \
       --tss_bed {input} \
       --th1_fragments {params.control} \
       --th2_fragments {params.treatment} \
       --qc_commands_filename {output} 
      """

rule runQC:
    input:
        config['qc_cmd']
    params:
        tss_bed=config['out_dir'] + "/QC/" + config['modified_tss_bed']
    output:
        "done.txt"
    shell:
        """
        tmp_file=$(mktemp)
        echo -e "# Chromosome\tstart\tend\tgene_name\tscore\tstrand" | cat - {params.tss_bed} > $tmp_file
        mv $tmp_file {params.tss_bed}
        bash {input}
        touch {output}
        """
 
