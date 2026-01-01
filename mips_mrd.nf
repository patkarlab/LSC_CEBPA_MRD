#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""

process trimming {
	tag "${Sample}" 
	input:
		val (Sample)
	output:
		//tuple val (Sample), file("${Sample}.R1.trimmed.fastq"), file("${Sample}.R2.trimmed.fastq")
		tuple val (Sample), file("*1P.fq.gz"), file("*2P.fq.gz")
	script:
	"""	
	#/home/programs/ea-utils/clipper/fastq-mcf -o ${Sample}.R1.trimmed.fastq -o ${Sample}.R2.trimmed.fastq -l 53 -k 0 -q 0 /home/diagnostics/pipelines/smMIPS_pipeline/code/functions/preprocess_reads_miseq/smmip_adaptors.fa ${params.sequences}/${Sample}_S*_R1_*.fastq.gz  ${params.sequences}/${Sample}_S*_R2_*.fastq.gz
	trimmomatic PE \
	${params.sequences}/${Sample}_*R1.fastq.gz ${params.sequences}/${Sample}_*R2.fastq.gz \
	-baseout ${Sample}.fq.gz \
	ILLUMINACLIP:${params.adaptors}:2:30:10:2:keepBothReads \
	ILLUMINACLIP:${params.nextera_adapters}:2:30:10:2:keepBothReads \
	LEADING:3 SLIDINGWINDOW:4:15 MINLEN:40	
	sleep 5s
	"""
}

process fastqc{
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*'
	input:
		val Sample
	output:
		path "*"
	"""
	${params.fastqc} -o ./ -f fastq ${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R2_*.fastq.gz
	"""
}

process FastqToBam {
	input:
		tuple val (Sample), file(R1_trim_fastq), file(R2_trim_fastq)
	output:
		tuple val (Sample), file ("*.unmapped.bam")
	script:	
	"""
	java -Xmx12g -jar "/home/programs/fgbio/fgbio-2.0.1.jar" --compression 1 --async-io FastqToBam --input ${R1_trim_fastq} ${R2_trim_fastq} --read-structures 4M+T 4M+T --sample ${Sample} --library ${Sample} --output ${Sample}.unmapped.bam
	"""
}

process FastqToBam_NARASIMHA {
	input:
		tuple val (Sample), file(R1_trim_fastq), file(R2_trim_fastq)
	output:
		tuple val (Sample), file ("*.unmapped.bam")
	script:	
	"""
	java -Xmx12g -jar "/home/programs/fgbio/fgbio-2.0.1.jar" --compression 1 --async-io FastqToBam --input ${R1_trim_fastq} ${R2_trim_fastq} --read-structures 8M+T +T --sample ${Sample} --library ${Sample} --output ${Sample}.unmapped.bam
	"""
}

process FastqToBam_U2AF1 {
	tag "${Sample}"
	input:
		tuple val (Sample), file(R1_trim_fastq), file(R2_trim_fastq)
		// val (Sample)
	output:
		tuple val (Sample), file ("*.unmapped.bam")
	script:	
	"""
	java -Xmx12g -jar "/home/programs/fgbio/fgbio-2.0.1.jar" --compression 1 --async-io FastqToBam --input ${R1_trim_fastq} ${R2_trim_fastq} --read-structures 8M+T 8M+T --sample ${Sample} --library ${Sample} --output ${Sample}.unmapped.bam
	#java -Xmx12g -jar "/home/programs/fgbio/fgbio-2.0.1.jar" --compression 1 --async-io FastqToBam --input ${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R2_*.fastq.gz --read-structures 17S8M+T 20S8M+T --sample ${Sample} --library ${Sample} --output ${Sample}.unmapped.bam
	"""
}

process trimming_trimmomatic {
	input:
		val Sample
	output:
		tuple val (Sample), file("*1P.fq.gz"), file("*2P.fq.gz")
	script:
	"""
	trimmomatic PE \
	${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R3_*.fastq.gz \
	-baseout ${Sample}.fq.gz \
	ILLUMINACLIP:${params.adaptors}:2:30:10:2:keepBothReads \
	LEADING:3 SLIDINGWINDOW:4:15 MINLEN:40
	sleep 5s
	"""
}

process pair_assembly_pear {
	tag "${Sample}"
	input:
		tuple val (Sample), file(paired_forward), file(paired_reverse)
	output:
		tuple val (Sample), file("*.assembled.fastq")
	script:
	"""
	${params.pear_path} -f ${paired_forward} -r ${paired_reverse} -o ${Sample} -n 53 -j 25
	"""
}

process mapping_reads {
	input:
		tuple val (Sample), file (pairAssembled)
	output:
		tuple val (Sample), file ("*.bam")
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200" \
	-M -t 20 ${params.genome} ${pairAssembled} | ${params.samtools} sort -@ 8 -o ${Sample}.bam -
	"""
}

process sam_conversion {
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.uncollaps.bam*'
	input:
		tuple val (Sample), file (bamfile)	
	output:
		tuple val(Sample), file ("*.uncollaps.bam"), file ("*.uncollaps.bam.bai")
	script:
	"""
	${params.samtools} sort ${bamfile} > ${Sample}.uncollaps.bam
	${params.samtools} index ${Sample}.uncollaps.bam > ${Sample}.uncollaps.bam.bai
	"""
}

process MapBam {
	tag "${Sample}"	
	//publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.mapped.bam*'
	input:
		tuple val (Sample), file(unmapped_bam)
	output:
		tuple val (Sample), file("*.mapped.bam")
	script:
	"""
	# Align the data with bwa and recover headers and tags from the unmapped BAM
	${params.samtools} fastq ${unmapped_bam} | bwa mem -t 16 -p -K 150000000 -Y ${params.genome} - | java -Xmx4g -jar ${params.fgbio_path} --compression 1 --async-io ZipperBams --unmapped ${unmapped_bam} --ref ${params.genome} --output ${Sample}.mapped.bam
	#${params.samtools} sort ${Sample}.mapped.bam -o ${Sample}.Uncollaps.bam
	#${params.samtools} index ${Sample}.Uncollaps.bam > ${Sample}.Uncollaps.bam.bai
	"""
}

process GroupReadsByUmi {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: "${Sample}_family.png"
	input:
		tuple val (Sample), file(mapped_bam)
	output: 
		tuple val (Sample), file ("*.grouped.bam"), emit : grouped_bam_ch
		path ("${Sample}_family.png")
	script:
	"""
	java -Xmx20g -jar ${params.fgbio_path} --compression 1 --async-io GroupReadsByUmi --input ${mapped_bam} --strategy Adjacency --edits 1 --output ${Sample}.grouped.bam --family-size-histogram ${Sample}.tag-family-sizes.txt
	plot_family_sizes.py ${Sample}.tag-family-sizes.txt ${Sample}_family
	"""
}

process GroupReadsByUmi_amplicon {
	tag "${Sample}"
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: "${Sample}_family.png"
	input:
		tuple val (Sample), file(mapped_bam)
	output: 
		tuple val (Sample), file ("*.grouped.bam"), emit : filtered_files_ch
		path ("${Sample}_family.png")
	script:
	"""
	java -Xmx25g -jar ${params.fgbio_new} --compression 1 --async-io GroupReadsByUmi --input ${mapped_bam} --strategy Identity --edits 0 --output ${Sample}.grouped.bam --family-size-histogram ${Sample}.tag-family-sizes.txt
	${params.plot_family_size} ${Sample}.tag-family-sizes.txt ${Sample}_family
	"""
}

process CallMolecularConsensusReads {
	tag "${Sample}"
	input:
		tuple val (Sample), file (grouped_bam)
	output:
		tuple val (Sample), file ("*.cons.unmapped.bam")
	script:
	"""
	# This step generates unmapped consensus reads from the grouped reads and immediately filters them
	java -Xmx25g -jar ${params.fgbio_path} --compression 0 CallMolecularConsensusReads --input ${grouped_bam} --output /dev/stdout --min-reads 4 --threads 4 | java -Xmx4g -jar ${params.fgbio_path} --compression 1 FilterConsensusReads --input /dev/stdin --output ${Sample}.cons.unmapped.bam --ref ${params.genome} --min-reads 4 --min-base-quality 20 --max-base-error-rate 0.25
	#java -Xmx4g -jar ${params.fgbio_path} --compression 0 CallMolecularConsensusReads --input ${grouped_bam} --output /dev/stdout --min-reads 2 --min-input-base-quality 20 --threads 4 | java -Xmx4g -jar ${params.fgbio_path} --compression 1 FilterConsensusReads --input /dev/stdin --output ${Sample}.cons.unmapped.bam --ref ${params.genome} --min-reads 2 --min-base-quality 20 --max-base-error-rate 0.35
	#java -Xmx4g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 1 CallMolecularConsensusReads --input NPM1-078.grouped.bam --output /dev/stdout --min-reads 2 --min-input-base-quality 20 --error-rate-post-umi=30 --threads 30 | java -Xmx4g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 1 FilterConsensusReads --input /dev/stdin --output NPM1-078.cons.unmapped.bam --ref /home/reference_genomes/hg19_broad/hg19_all.fasta --min-reads 2 --min-base-quality 20 --max-base-error-rate 0.2
	#java -Xmx4g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 0 CallMolecularConsensusReads --input NPM1-078.grouped.bam --output /dev/stdout --min-reads 2 --min-input-base-quality 20 --threads 10 | java -Xmx4g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 0 FilterConsensusReads --input /dev/stdin --output NPM1-078.cons.unmapped.bam --ref /home/reference_genomes/hg19_broad/hg19_all.fasta --min-reads 2 --min-base-quality 20 
	"""
}

process CallMolecularConsensusReadsAmplicon {
	tag "${Sample}"
	input:
		tuple val (Sample), file (grouped_bam)
	output:
		tuple val (Sample), file ("*.cons.unmapped.bam")
	script:
	"""
	# This step generates unmapped consensus reads from the grouped reads and immediately filters them
	java -Xmx30g -jar ${params.fgbio_new} --compression 0 CallMolecularConsensusReads --input ${grouped_bam} --output /dev/stdout --min-reads 20 --threads $task.cpus | java -Xmx25g -jar ${params.fgbio_new} --compression 1 FilterConsensusReads --input /dev/stdin --output ${Sample}.cons.unmapped.bam --ref ${params.genome} --min-reads 20 --min-base-quality 20 --max-base-error-rate 0.2
	#java -Xmx30g -jar ${params.fgbio_new} --compression 0 CallDuplexConsensusReads --input ${grouped_bam} --output /dev/stdout --min-reads 4 --threads $task.cpus | java -Xmx25g -jar ${params.fgbio_new} --compression 1 FilterConsensusReads --input /dev/stdin --output ${Sample}.cons.unmapped.bam --ref ${params.genome} --min-reads 4 --min-base-quality 20 --max-base-error-rate 0.2
	"""
}

process FilterConsBam {
	tag "${Sample}"
	input:
		tuple val (Sample), file (cons_unmapped_bam)
	output:
		tuple val (Sample), file ("*.cons.filtered.bam"), file ("*.cons.filtered.bam.bai")
	script:
	"""
	${params.samtools} fastq ${cons_unmapped_bam} | bwa mem -t 16 -p -K 150000000 -Y ${params.genome} - | java -Xmx12g -jar ${params.fgbio_path} --compression 0 --async-io ZipperBams --unmapped ${cons_unmapped_bam} --ref ${params.genome} --tags-to-reverse Consensus --tags-to-revcomp Consensus | ${params.samtools} sort --threads 8 -o ${Sample}.cons.filtered.bam

	${params.samtools} index ${Sample}.cons.filtered.bam > ${Sample}.cons.filtered.bam.bai
	"""
}

process RealignerTargetCreator {
	input:
		tuple val (Sample), file (bamFile), file(bamBai)
	output:
		tuple val(Sample), file ("*.intervals")
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T RealignerTargetCreator -R ${params.genome} -nt 10 -I ${bamFile} --known ${params.site1} -o ${Sample}.intervals
	"""
}

process IndelRealigner{
	input:
		tuple val(Sample), file(bamFile), file(bamBai)
	output:
		tuple val(Sample), file ("*.realigned.bam")
	script:
	"""
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T RealignerTargetCreator -R ${params.genome} -nt 10 -I ${bamFile} --known ${params.site1} -o ${Sample}.intervals
	${params.java_path}/java -Xmx8G -jar ${params.GATK38_path} -T IndelRealigner -R ${params.genome} -I ${bamFile} -known ${params.site1} --targetIntervals ${Sample}.intervals -o ${Sample}.realigned.bam
	"""
}

process SyntheticFastq {
	input:
		tuple val (Sample), file (cons_filt_bam), file (cons_filt_bam_bai)
	output:
		tuple val (Sample), file ("*.synreads.bam"), file ("*.synreads.bam.bai")
	script:
	"""
	# Making a synthetic fastq and bam based on the number of reads in the sample bam file
	${params.PosControlScript} ${cons_filt_bam} ${Sample}_inreads
	${params.fastq_bam} ${Sample}_inreads		# This script will output a file named *_inreads.fxd_sorted.bam

	# Merging the Sample bam file with positive control bam file
	mv ${cons_filt_bam} ${Sample}.cons.filtered_merge.bam
	${params.samtools} merge -f ${Sample}.synreads.bam ${Sample}.cons.filtered_merge.bam ${Sample}_inreads.fxd_sorted.bam
	${params.samtools} index ${Sample}.synreads.bam > ${Sample}.synreads.bam.bai
	"""
}

process ABRA2_realign {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.ErrorCorrectd.bam*'
	input:
		tuple val (Sample), file (filt_bam), file (filt_bai)
	output:
		tuple val (Sample), file ("*.ErrorCorrectd.bam"), file ("*.ErrorCorrectd.bam.bai")
	script:
	"""
	${params.java_path}/java -Xmx16G -jar ${params.abra2_path}/abra2-2.23.jar --in ${filt_bam} --out ${Sample}.ErrorCorrectd.bam --ref ${params.genome} --threads 8 --targets ${params.bedfile}.bed --tmpdir ./ > abra.log
	${params.samtools} index ${Sample}.ErrorCorrectd.bam > ${Sample}.ErrorCorrectd.bam.bai
	"""
}

process CNS_filegen {
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.final.cns'
	input:
		tuple val (Sample), file (abra_bam), file (abra_bai)
	output:
		tuple val(Sample), file ("*.final.cns")
	script:
	"""
	#${params.samtools} mpileup -s -x -BQ 0 -q 1 -d 1000000 --skip-indels -f ${params.genome} ${abra_bam} > ${Sample}.mpileup

	#${params.samtools} mpileup -A -a --skip-indels -l ${params.bedfile}.bed -f ${params.genome} ${abra_bam} > ${Sample}.mpileup
	#${params.samtools} mpileup -d 1000000 -A -a --skip-indels -l ${params.bedfile}.bed -f ${params.genome} ${abra_bam} > ${Sample}.mpileup

	${params.samtools} mpileup -s -x -BQ 0 -q 1 -d 1000000 -A -a --skip-indels -l ${params.bedfile}.bed -f ${params.genome} ${abra_bam} > ${Sample}.mpileup
	# These parmeters were taken from the Espresso paper(Abelson et al., Sci. Adv. 2020; 6 : eabe3722)
	# A stringent (lower) p value removes all the snps hence, p-value is 1
	# ${params.java_path}/java -jar ${params.varscan_path} pileup2cns ${Sample}".mpileup" --variants SNP --min-coverage 10 --min-reads2 1 --min-avg-qual 30 --min-var-freq 0.0001 --p-value 1 --strand-filter 0 > ${Sample}".cns"

	#${params.java_path}/java -jar ${params.varscan_path} pileup2cns ${Sample}".mpileup" --variants SNP --min-coverage 2 --min-reads2 1 --min-avg-qual 15 --min-var-freq 0.0001 --p-value 1 --strand-filter 0 > ${Sample}".cns"

	${params.java_path}/java -jar ${params.varscan_path} pileup2cns ${Sample}".mpileup" --variants SNP --min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --strand-filter 0 > ${Sample}".cns"

	grep -v '^chrM' ${Sample}".cns" > ${Sample}".nochrM.cns"
	${params.filter_cns} ${Sample}".nochrM.cns" ${Sample}".final.cns"
	"""
}

process espresso {
	input:
		val (tuple)
	script:	
	"""
	for i in `cat ${params.input}`
	do
		ln -s ${params.outdir}/\${i}/\${i}.final.cns ./
	done
	${params.umi_error_model} ${params.input} ${params.bedfile}.bed ${params.outdir}
	"""
}

process coverage_bedtools {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_umi_cov.regions_bedtools.bed'
	input:
		tuple val (Sample), file (abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("${Sample}_umi_cov.regions_bedtools.bed")
	script:
	"""
	${params.bedtools} bamtobed -i ${abra_bam} > ${Sample}_ErrorCorrected.bed
	${params.bedtools} coverage -counts -a ${params.bedfile}.bed -b ${Sample}_ErrorCorrected.bed > "${Sample}_umi_cov.regions_bedtools.bed"
	
	"""	
}

process coverage_mosdepth {
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_umi_cov.regions.bed'
	input:
		tuple val (Sample), file (abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("${Sample}_umi_cov.regions.bed")
	script:
	"""
	${params.mosdepth_script} ${abra_bam} ${Sample}_umi_cov ${params.bedfile}.bed
	#${params.mosdepth_script} ${abra_bam} ${Sample}_exon_umi_cov ${params.bedfile2}.bed
	"""
}

process coverage_bedtools_uncoll {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_uncollaps_bedtools.bed'
	//publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_exon_uncollaps_bedtools.bed'
	input:
		tuple val (Sample), file (uncollaps_bam), file (uncollaps_bam_bai)
	output:
		tuple val (Sample), file ("${Sample}_uncollaps_bedtools.bed")
	script:
	"""
	${params.bedtools} bamtobed -i ${uncollaps_bam} > ${Sample}_Uncollaps.bed
	${params.bedtools} coverage -counts -a ${params.bedfile}.bed -b ${Sample}_Uncollaps.bed > "${Sample}_uncollaps_bedtools.bed"
	#${params.bedtools} coverage -counts -a ${params.bedfile2}.bed -b ${Sample}_Uncollaps.bed > "${Sample}_exon_uncollaps_bedtools.bed"
	"""	
}

process coverage_mosdepth_uncoll {
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.regions.bed'
	input:
		tuple val (Sample), file (abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("${Sample}_uncollaps.regions.bed"), file ("${Sample}_exon_uncollaps.regions.bed")
	script:
	"""
	${params.mosdepth_script} ${abra_bam} ${Sample}_uncollaps ${params.bedfile}.bed
	${params.mosdepth_script} ${abra_bam} ${Sample}_exon_uncollaps ${params.bedfile2}.bed
	"""
}

process hsmetrics_run {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_hsmetrics.txt'
	input:
		tuple val(Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*_hsmetrics.txt")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${abra_bam} O= ${Sample}_hsmetrics.txt BAIT_INTERVALS= ${params.bedfile}.interval_list TARGET_INTERVALS= ${params.bedfile}.interval_list R= ${params.genome} VALIDATION_STRINGENCY=LENIENT
	"""	
}

process hsmetrics_run_uncoll {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_uncoll_hsmetrics.txt'
	input:
		tuple val(Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*_uncoll_hsmetrics.txt")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${abra_bam} O= ${Sample}_uncoll_hsmetrics.txt BAIT_INTERVALS= ${params.bedfile}.interval_list TARGET_INTERVALS= ${params.bedfile}.interval_list R= ${params.genome} VALIDATION_STRINGENCY=LENIENT
	"""
}

process minimap_getitd {
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		path "*_getitd"
	script:
	"""
	${params.samtools} sort ${abra_bam} -o ${Sample}.sorted.bam
	${params.samtools} index ${Sample}.sorted.bam
	${params.samtools} view ${Sample}.sorted.bam -b -h chr13 > ${Sample}.chr13.bam
	${params.bedtools} bamtofastq -i ${Sample}.chr13.bam -fq ${Sample}_chr13.fastq
	python ${params.get_itd_path}/getitd.py -reference ${params.get_itd_path}/anno/amplicon.txt -anno ${params.get_itd_path}/anno/amplicon_kayser.tsv -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern 8 ${Sample} ${Sample}_chr13.fastq
	"""
}

process mutect2_run {
	maxForks 10
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("${Sample}.mutect2.vcf")
	script:
	"""
	${params.java_path}/java -Xmx10G -jar ${params.GATK42_path} Mutect2 -R ${params.genome} -I:tumor ${abra_bam} -O ${Sample}.mutect2.vcf -L ${params.bedfile}.bed -mbq 20
	${params.java_path}/java -Xmx10G -jar ${params.GATK42_path} FilterMutectCalls -V ${Sample}.mutect2.vcf --stats ${Sample}.mutect2.vcf.stats -O ${Sample}_filtered.vcf -R ${params.genome} --unique-alt-read-count 3 --min-median-base-quality 20 --min-median-mapping-quality 30
	#perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.mutect2.vcf --outfile ${Sample}.mutect2.avinput -allsample -withfreq --includeinfo
	#perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.mutect2.avinput --out ${Sample}.mutect2.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
	#
	#if [ -s ${Sample}.mutect2.somaticseq.hg19_multianno.csv ]; then
	#	python3 ${PWD}/scripts/somaticseqoutput-format_v2_mutect2.py ${Sample}.mutect2.somaticseq.hg19_multianno.csv mutect2_.csv
	#else
	#	touch mutect2_.csv
	#fi	
	"""
}

process mutect2_run_uncoll {
	tag "${Sample}"
	maxForks 10
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("${Sample}.mutect2.vcf")
	script:
	"""
	${params.java_path}/java -Xmx80G -jar ${params.GATK42_path} Mutect2 -R ${params.genome} -I:tumor ${abra_bam} -O ${Sample}.mutect2.vcf -L ${params.bedfile}.bed -mbq 20
	${params.java_path}/java -Xmx80G -jar ${params.GATK42_path} FilterMutectCalls -V ${Sample}.mutect2.vcf --stats ${Sample}.mutect2.vcf.stats -O ${Sample}_filtered.vcf -R ${params.genome} --unique-alt-read-count 3 --min-median-base-quality 20 --min-median-mapping-quality 30
	#perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.mutect2.vcf --outfile ${Sample}.mutect2.avinput -allsample -withfreq --includeinfo
	#perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.mutect2.avinput --out ${Sample}.mutect2.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
	#
	#if [ -s ${Sample}.mutect2.somaticseq.hg19_multianno.csv ]; then
	#	python3 ${PWD}/scripts/somaticseqoutput-format_v2_mutect2.py ${Sample}.mutect2.somaticseq.hg19_multianno.csv mutect2_.csv
	#else
	#	touch mutect2_.csv		
	#fi	
	"""
}

process vardict {
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.vardict.vcf")
	script:
	"""
	#VarDict -G ${params.genome} -f 0.0001 -N ${Sample} -b ${abra_bam} -c 1 -S 2 -E 3 -g 4 ${params.bedfile}.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.0001 > ${Sample}.vardict.vcf
	VarDict -G ${params.genome} -f 0.0001 -r 8 -N ${Sample} -b ${abra_bam} -c 1 -S 2 -E 3 -g 4 ${params.bedfile}.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.0001 > ${Sample}.vardict.vcf
	"""
}

process vardict_uncoll {
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.vardict.vcf")
	script:
	"""
	VarDict -G ${params.genome} -f 0.0001 -r 8 -N ${Sample} -b ${abra_bam} -c 1 -S 2 -E 3 -g 4 ${params.bedfile}.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.0001 > ${Sample}.vardict.vcf
	"""
}

process lofreq {
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val(Sample), file ("*.lofreq.filtered.vcf")
	script:
	"""
	${params.lofreq_path} viterbi -f ${params.genome} -o ${Sample}.lofreq.pre.bam ${abra_bam}
	${params.samtools} sort ${Sample}.lofreq.pre.bam > ${Sample}.lofreq.bam
	#${params.lofreq_path} call -b dynamic -C 2 -a 0.00005 -q 20 -Q 20 -m 50 -f ${params.genome} -l ${params.bedfile}.bed -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam --no-default-filter
	${params.lofreq_path} call -f ${params.genome} -l ${params.bedfile}.bed -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam --no-default-filter
	#${params.lofreq_path} filter --no-defaults -a 0.0001 -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf
	${params.lofreq_path} filter -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf
	"""
}

process lofreq_uncoll {
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val(Sample), file ("*.lofreq.filtered.vcf")
	script:
	"""
	${params.lofreq_path} viterbi -f ${params.genome} -o ${Sample}.lofreq.pre.bam ${abra_bam}
	${params.samtools} sort ${Sample}.lofreq.pre.bam > ${Sample}.lofreq.bam
	${params.lofreq_path} call -b dynamic -C 2 -a 0.00005 -q 20 -Q 20 -m 50 -f ${params.genome} -l ${params.bedfile}.bed -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam --no-default-filter
	${params.lofreq_path} filter --no-defaults -a 0.0001 -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf
	"""
}

process varscan {
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.varscan.vcf")
	script:
	"""
	#${params.samtools} mpileup -d 1000000 -f ${params.genome} ${abra_bam} > ${Sample}.mpileup
	${params.samtools} mpileup -d 1000000 -A -a -l ${params.bedfile}.bed -f ${params.genome} ${abra_bam} > ${Sample}.mpileup
	#${params.java_path}/java -jar ${params.varscan_path} mpileup2snp ${Sample}.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 15 --min-var-freq 0.0001 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_snp.vcf

	if [ -s ${Sample}.mpileup ]; then
		${params.java_path}/java -jar ${params.varscan_path} mpileup2snp ${Sample}.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_snp.vcf
	else 
		echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > ${Sample}.varscan_snp.vcf
	fi

	#${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 15 --min-var-freq 0.0001 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_indel.vcf
	if [ -s ${Sample}.mpileup ]; then
		${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_indel.vcf
	else
		echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > ${Sample}.varscan_indel.vcf
	fi

	bgzip -c ${Sample}.varscan_snp.vcf > ${Sample}.varscan_snp.vcf.gz
	bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_snp.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} concat -a ${Sample}.varscan_snp.vcf.gz ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf
	"""
}

process strelka {
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.strelka*.vcf'
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.strelka*.vcf")
	script:
	"""
	${params.strelka_path}/configureStrelkaGermlineWorkflow.py --bam ${abra_bam} --referenceFasta ${params.genome} --callRegions ${params.bedfile}.bed.gz --targeted --runDir ./
	./runWorkflow.py -m local -j 20
	gunzip -f ./results/variants/variants.vcf.gz
	mv ./results/variants/variants.vcf ${Sample}.strelka.vcf

	${params.strelka_path}/configureStrelkaSomaticWorkflow.py --normalBam ${params.NA12878_bam} --tumorBam ${abra_bam} --referenceFasta ${params.genome} --callRegions ${params.bedfile}.bed.gz --targeted --runDir ./strelka-somatic
	./strelka-somatic/runWorkflow.py -m local -j 20
	${params.bcftools_path} concat -a strelka-somatic/results/variants/somatic.indels.vcf.gz strelka-somatic/results/variants/somatic.snvs.vcf.gz -o ./${Sample}.strelka-somatic.vcf
	"""
}

process varscan_uncoll {
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.varscan.vcf")
	script:
	"""
	${params.samtools} mpileup -d 1000000 -A -a -l ${params.bedfile}.bed -f ${params.genome} ${abra_bam} > ${Sample}.mpileup
	${params.java_path}/java -jar ${params.varscan_path} mpileup2snp ${Sample}.mpileup -min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_snp.vcf
	${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup -min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_indel.vcf
	bgzip -c ${Sample}.varscan_snp.vcf > ${Sample}.varscan_snp.vcf.gz
	bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_snp.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} concat -a ${Sample}.varscan_snp.vcf.gz ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf
	"""
}

process somaticSeq_run {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.ErrorCorrectd.vcf'
	//publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.ErrorCorrectd.csv'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_ErrorCorrectd.xlsx'
	input:
		tuple val (Sample), file (Coverage), file (mutectVcf), file (vardictVcf), file (varscanVcf), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("${Sample}.ErrorCorrectd.vcf"), emit: correctd_vcf_files
		tuple val (Sample), file ("${Sample}_ErrorCorrectd.xlsx"), emit: correctd_excels
	script:
	"""
	${params.vcf_filter} ${mutectVcf} ${Sample}.mutect2_sort.vcf
	${params.vcf_filter} ${vardictVcf} ${Sample}.vardict_sort.vcf
	${params.vcf_filter} ${varscanVcf} ${Sample}.varscan_sort.vcf
	echo "inside the somaticSeq"
	
	somaticseq_parallel.py --output-directory ${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads 25 --pass-threshold 0 --lowqual-threshold 0 --algorithm xgboost -minMQ 0 -minBQ 0 -mincaller 0 --dbsnp-vcf /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${finalBam} --mutect2-vcf ${Sample}.mutect2_sort.vcf --vardict-vcf ${Sample}.vardict_sort.vcf --varscan-vcf ${Sample}.varscan_sort.vcf --sample-name ${Sample}

	${params.vcf_sorter_path} ${Sample}.somaticseq/Consensus.sSNV.vcf ${Sample}.somaticseq/somaticseq_snv.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_snv.vcf > ${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_snv.vcf.gz

	${params.vcf_sorter_path} ${Sample}.somaticseq/Consensus.sINDEL.vcf ${Sample}.somaticseq/somaticseq_indel.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_indel.vcf > ${Sample}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_indel.vcf.gz

	${params.bcftools_path} concat -a ${Sample}.somaticseq/somaticseq_snv.vcf.gz ${Sample}.somaticseq/somaticseq_indel.vcf.gz -o ${Sample}.ErrorCorrectd.vcf

	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.ErrorCorrectd.vcf --outfile ${Sample}.ErrorCorrectd.avinput -allsample -withfreq --includeinfo
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.ErrorCorrectd.avinput --out ${Sample}.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt

	if [ -s ${Sample}.somaticseq.hg19_multianno.csv ]; then
		python3 ${params.format_somaticseq_script} ${Sample}.somaticseq.hg19_multianno.csv ${Sample}.ErrorCorrectd.csv
	else
		touch ${Sample}.ErrorCorrectd.csv
	fi
	python3 ${PWD}/scripts/final_output.py ${Sample}_ErrorCorrectd.xlsx ${Sample}.ErrorCorrectd.csv ${Coverage}
	"""
}

process somaticSeq_run_uncoll {
	tag "${Sample}"
	//publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.UncollapsVariant.vcf'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_uncollapsed.xlsx'
	input:
		tuple val (Sample), file (Coverage_uncollaps), file (mutectVcf), file (vardictVcf), file (varscanVcf), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	${params.vcf_filter} ${mutectVcf} ${Sample}.mutect2_sort.vcf
	${params.vcf_filter} ${vardictVcf} ${Sample}.vardict_sort.vcf
	${params.vcf_filter} ${varscanVcf} ${Sample}.varscan_sort.vcf

	somaticseq_parallel.py --output-directory ${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads 25 --pass-threshold 0 --lowqual-threshold 0 --algorithm xgboost -minMQ 0 -minBQ 0 -mincaller 0 --dbsnp-vcf /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${finalBam} --mutect2-vcf ${Sample}.mutect2_sort.vcf --vardict-vcf ${Sample}.vardict_sort.vcf --varscan-vcf ${Sample}.varscan_sort.vcf --sample-name ${Sample}

	${params.vcf_sorter_path} ${Sample}.somaticseq/Consensus.sSNV.vcf ${Sample}.somaticseq/somaticseq_snv.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_snv.vcf > ${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_snv.vcf.gz

	${params.vcf_sorter_path} ${Sample}.somaticseq/Consensus.sINDEL.vcf ${Sample}.somaticseq/somaticseq_indel.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_indel.vcf > ${Sample}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_indel.vcf.gz

	${params.bcftools_path} concat -a ${Sample}.somaticseq/somaticseq_snv.vcf.gz ${Sample}.somaticseq/somaticseq_indel.vcf.gz -o ${Sample}.ErrorCorrectd.vcf

	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.ErrorCorrectd.vcf --outfile ${Sample}.ErrorCorrectd.avinput -allsample -withfreq --includeinfo
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.ErrorCorrectd.avinput --out ${Sample}.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt

	if [ -s ${Sample}.somaticseq.hg19_multianno.csv ]; then
		python3 ${params.format_somaticseq_script} ${Sample}.somaticseq.hg19_multianno.csv ${Sample}.ErrorCorrectd.csv
	else
		touch ${Sample}.ErrorCorrectd.csv
	fi
	python3 ${PWD}/scripts/final_output.py ${Sample}_uncollapsed.xlsx ${Sample}.ErrorCorrectd.csv ${Coverage_uncollaps}
	"""
}

include { FASTP } from './scripts/fastp.nf'
include { CUTADAPT } from './scripts/cutadapt.nf'
include { BAMTOFASTQ } from './scripts/cutadapt.nf'

workflow MRD {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.set { samples_ch }

	main:
		trimming(samples_ch)
		FastqToBam(trimming.out) 
		MapBam(FastqToBam.out) 
		GroupReadsByUmi(MapBam.out) 
		CallMolecularConsensusReads(GroupReadsByUmi.out.grouped_bam_ch) 
		FilterConsBam(CallMolecularConsensusReads.out) 
		SyntheticFastq(FilterConsBam.out)
		ABRA2_realign(SyntheticFastq.out)
		CNS_filegen(ABRA2_realign.out)
		espresso(CNS_filegen.out.collect())

		pair_assembly_pear(trimming.out) | mapping_reads | sam_conversion

		coverage_bedtools(ABRA2_realign.out)
		coverage_bedtools_uncoll(sam_conversion.out)

		hsmetrics_run(ABRA2_realign.out)
		hsmetrics_run_uncoll(sam_conversion.out)
		minimap_getitd(sam_conversion.out)

		mutect2_run(ABRA2_realign.out)
		vardict(ABRA2_realign.out)
		//lofreq(ABRA2_realign.out)
		varscan(ABRA2_realign.out)
		//strelka(ABRA2_realign.out)

		mutect2_run_uncoll(sam_conversion.out)
		vardict_uncoll(sam_conversion.out)
		//lofreq_uncoll(MapBam.out)
		varscan_uncoll(sam_conversion.out)

		somaticSeq_run(coverage_bedtools.out.join(mutect2_run.out.join(vardict.out.join(varscan.out.join(ABRA2_realign.out)))))
		somaticSeq_run_uncoll(coverage_bedtools_uncoll.out.join(mutect2_run_uncoll.out.join(vardict_uncoll.out.join(varscan_uncoll.out.join(sam_conversion.out)))))
}

workflow NARASIMHA_MRD {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.set { samples_ch }

	main:
		trimming(samples_ch)
		fastqc(samples_ch)
		FastqToBam_NARASIMHA(trimming.out) 
		MapBam(FastqToBam_NARASIMHA.out) 
		GroupReadsByUmi(MapBam.out) 
		CallMolecularConsensusReads(GroupReadsByUmi.out.grouped_bam_ch) 
		FilterConsBam(CallMolecularConsensusReads.out) 
		//IndelRealigner(FilterConsBam.out)
		SyntheticFastq(FilterConsBam.out)
		ABRA2_realign(SyntheticFastq.out)
		CNS_filegen(ABRA2_realign.out)
		espresso(CNS_filegen.out.collect())

		pair_assembly_pear(trimming.out) | mapping_reads | sam_conversion

		coverage_bedtools(ABRA2_realign.out)
		coverage_bedtools_uncoll(sam_conversion.out)

		hsmetrics_run(ABRA2_realign.out)
		hsmetrics_run_uncoll(sam_conversion.out)
		//minimap_getitd(sam_conversion.out)

		mutect2_run(ABRA2_realign.out)
		vardict(ABRA2_realign.out)
		//lofreq(ABRA2_realign.out)
		varscan(ABRA2_realign.out)
		//strelka(ABRA2_realign.out)

		mutect2_run_uncoll(sam_conversion.out)
		vardict_uncoll(sam_conversion.out)
		//lofreq_uncoll(MapBam.out)
		varscan_uncoll(sam_conversion.out)

		somaticSeq_run(coverage_bedtools.out.join(mutect2_run.out.join(vardict.out.join(varscan.out.join(ABRA2_realign.out)))))
		somaticSeq_run_uncoll(coverage_bedtools_uncoll.out.join(mutect2_run_uncoll.out.join(vardict_uncoll.out.join(varscan_uncoll.out.join(sam_conversion.out)))))
}

workflow NARASIMHA_MRD_VALIDATION {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.set { samples_ch }

	main:
		trimming(samples_ch)
		fastqc(samples_ch)
		FastqToBam_NARASIMHA(trimming.out) 
		MapBam(FastqToBam_NARASIMHA.out) 
		GroupReadsByUmi(MapBam.out) 
		CallMolecularConsensusReads(GroupReadsByUmi.out.grouped_bam_ch) 
		FilterConsBam(CallMolecularConsensusReads.out) 
		//SyntheticFastq(FilterConsBam.out)
		ABRA2_realign(FilterConsBam.out)
		CNS_filegen(ABRA2_realign.out)
		//espresso(CNS_filegen.out.collect())

		pair_assembly_pear(trimming.out) | mapping_reads | sam_conversion

		coverage_mosdepth(ABRA2_realign.out)
		coverage_mosdepth_uncoll(sam_conversion.out)

		hsmetrics_run(ABRA2_realign.out)
		hsmetrics_run_uncoll(sam_conversion.out)
		minimap_getitd(sam_conversion.out)

		mutect2_run(ABRA2_realign.out)
		vardict(ABRA2_realign.out)
		//lofreq(ABRA2_realign.out)
		varscan(ABRA2_realign.out)
		//strelka(ABRA2_realign.out)

		mutect2_run_uncoll(sam_conversion.out)
		vardict_uncoll(sam_conversion.out)
		//lofreq_uncoll(MapBam.out)
		varscan_uncoll(sam_conversion.out)

		somaticSeq_run(coverage_mosdepth.out.join(mutect2_run.out.join(vardict.out.join(varscan.out.join(ABRA2_realign.out)))))
		somaticSeq_run_uncoll(coverage_mosdepth_uncoll.out.join(mutect2_run_uncoll.out.join(vardict_uncoll.out.join(varscan_uncoll.out.join(sam_conversion.out)))))
}

workflow U2AF1MRD {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.set { samples_ch }

	main:
		trimming(samples_ch)
		// FASTP(trimming.out)
		// fastqc(samples_ch)
		CUTADAPT (trimming.out)
		FastqToBam_U2AF1(CUTADAPT.out) 
		MapBam(FastqToBam_U2AF1.out) 
		GroupReadsByUmi_amplicon(MapBam.out) 
		CallMolecularConsensusReadsAmplicon(GroupReadsByUmi_amplicon.out.filtered_files_ch)
		// FilterConsBam(CallMolecularConsensusReadsAmplicon.out) 
		// ABRA2_realign(FilterConsBam.out)

		// Temporary arrangement
		// BAMTOFASTQ (FilterConsBam.out)

		// pair_assembly_pear(trimming.out) | mapping_reads | sam_conversion
		// pair_assembly_pear(BAMTOFASTQ.out) | mapping_reads | sam_conversion

		// coverage_bedtools(ABRA2_realign.out)
		// coverage_bedtools_uncoll(sam_conversion.out)

		// hsmetrics_run(ABRA2_realign.out)
		// hsmetrics_run_uncoll(sam_conversion.out)

		// mutect2_run(ABRA2_realign.out)
		// vardict(ABRA2_realign.out)
		// varscan(ABRA2_realign.out)

		// mutect2_run_uncoll(sam_conversion.out)
		// vardict_uncoll(sam_conversion.out)
		// varscan_uncoll(sam_conversion.out)

		// somaticSeq_run(coverage_bedtools.out.join(mutect2_run.out.join(vardict.out.join(varscan.out.join(ABRA2_realign.out)))))
		// somaticSeq_run_uncoll(coverage_bedtools_uncoll.out.join(mutect2_run_uncoll.out.join(vardict_uncoll.out.join(varscan_uncoll.out.join(sam_conversion.out)))))
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
	println "Completed at: ${workflow.complete}"
	println "Total time taken: ${workflow.duration}"
}

