#https://scenicplus.readthedocs.io/en/latest/human_cerebellum_ctx_db.html
g:
	awk '/^>/ {next} {n+=length($0)} END {print n}' neurog2.fa
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


cbust:
	wget https://github.com/weng-lab/cluster-buster/archive/refs/heads/master.zip -O cluster-buster.zip
	unzip cluster-buster.zip
	make 
	./cbust -h #export PATH 

neurog2_cistarget: neurog2.with_1kb_bg_padding.fa aertslab_motif_colleciton/v10nr_clust_public/singletons
	ls aertslab_motif_colleciton/v10nr_clust_public/singletons > motifs.txt
	python create_cisTarget_databases/create_cistarget_motif_databases.py \
	-f neurog2.with_1kb_bg_padding.fa \
	-M aertslab_motif_colleciton/v10nr_clust_public/singletons \
	-m motifs.txt \
	-o  neurog2_cistarget \
	--bgpadding 1000 \
	-t 12 \
	-c /nfs/turbo/umms-thahoang/sherine/mScenicP/cluster-buster-master/cbust
	
