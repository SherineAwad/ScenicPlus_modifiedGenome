#https://scenicplus.readthedocs.io/en/latest/human_cerebellum_ctx_db.html
####Creating cistarget DB
create_cisTarget_databases:
	git clone https://github.com/aertslab/create_cisTarget_databases

v10nr_clust_public:
	wget https://resources.aertslab.org/cistarget/programs/cbust
	mkdir -p aertslab_motif_colleciton
	wget -O aertslab_motif_colleciton/v10nr_clust_public.zip https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip
	cd aertslab_motif_colleciton; unzip -q v10nr_clust_public.zip


####Preparing the modified reference genome 
neurog2.fa:
	cp /nfs/turbo/umms-thahoang/Thanh/Genome_10x/Mouse_10xM_Neurog29SAmScarlet3/genes.gtf neurog2.gtf
	cp /nfs/turbo/umms-thahoang/Thanh/Genome_10x/Mouse_10xM_Neurog29SAmScarlet3/genome.fa neurog2.fa

neurog2.chrom.sizes: neurog2.gtf
	samtools faidx neurog2.fa
	cut -f1,2 neurog2.fa.fai > neurog2.chrom.sizes

neurog2_1kb_windows.bed: neurog2.chrom.sizes
	bedtools makewindows -g neurog2.chrom.sizes -w 1000 > neurog2_1kb_windows.bed 

neurog2.with_1kb_bg_padding.fa: create_cisTarget_databases neurog2_1kb_windows.bed neurog2.fa neurog2.chrom.sizes
	create_cisTarget_databases/create_fasta_with_padded_bg_from_bed.sh \
	neurog2.fa \
	neurog2.chrom.sizes \
	neurog2_1kb_windows.bed \
	neurog2.with_1kb_bg_padding.fa \
	1000 \
	yes


motifs.txt: aertslab_motif_colleciton/v10nr_clust_public/singletons 
	ls aertslab_motif_colleciton/v10nr_clust_public/singletons > motifs.txt



cbust:
	wget https://github.com/weng-lab/cluster-buster/archive/refs/heads/master.zip -O cluster-buster.zip
	unzip cluster-buster.zip
	make 
	./cbust -h

neurog2_cistarget: neurog2.with_1kb_bg_padding.fa aertslab_motif_colleciton/v10nr_clust_public/singletons
	ls aertslab_motif_colleciton/v10nr_clust_public/singletons > motifs.txt
	python create_cisTarget_databases/create_cistarget_motif_databases.py \
	-f neurog2.with_1kb_bg_padding.fa \
	-M aertslab_motif_colleciton/v10nr_clust_public/singletons \
	-m motifs.txt \
	-o  neurog2_cistarget \
	--bgpadding 1000 \
 	-t 20

#----------------------------------------------------
##Downloadig some needed resource

mm10-blacklist.v2.bed:
	wget https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz
	gunzip ENCFF547MET.bed.gz
	mv ENCFF547MET.bed mm10-blacklist.v2.bed



## Set up snakemake scenicplus init_snakemake --out_dir scplus_pipeline
motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl:
	wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl

mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather:
	wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather


mm10_screen_v10_clust.regions_vs_motifs.rankings.feather:
	wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather

mm10.fa:
	wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
	gunzip mm10.fa.gz


mm10-limited-upstream10000-tss-downstream10000-full-transcript.bed: 
	wget https://resources.aertslab.org/cistarget/regions/mm10-limited-upstream10000-tss-downstream10000-full-transcript.bed


v10nr_clust_public.zip:
# Download ranking DB
	get https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather
# Download scores DB
	wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.scores.feather
# Download motif collection (annotations)
	wget https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip



### Some random checking plots 
plot_bw_tracks: 
	python plot_bw_tracks.py \
  --bw1 scenicOuts/consensus_peak_calling/pseudobulk_bw_files/Control.bw \
  --bw2 scenicOuts/consensus_peak_calling/pseudobulk_bw_files/KO.bw \
  --chrom chr3 \
  --start 10000000 \
  --end 11000000 \
  --output chr3_tracks.png


KO_peaks.png:	
	python plot_narrowpeaks.py scenicOuts/consensus_peak_calling/macs2_peaks/Control_peaks.narrowPeak -o Control_peaks.png 
	python plot_narrowpeaks.py scenicOuts/consensus_peak_calling/macs2_peaks/KO_peaks.narrowPeak -o KO_peaks.png
