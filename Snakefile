configfile: "config.yaml"

rule all:
    input:
        f"{config['out_dir']}/consensus_peak_calling/bed_paths.tsv",
        f"{config['out_dir']}/consensus_peak_calling/bw_paths.tsv"

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

