#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process TRIM {
	tag "${Sample}" 
    label 'process_medium'
	input:
		tuple val (Sample), file (read1), file (read2)
		file (truseq_adapters)
		file (nextera_adapters)
	output:
		tuple val (Sample), file("*1P.fq.gz"), file("*2P.fq.gz")
	script:
	"""	
	trimmomatic PE \
	${read1} ${read2} \
	-baseout ${Sample}.fq.gz -threads ${task.cpus}\
	ILLUMINACLIP:${truseq_adapters}:2:30:10:2:keepBothReads \
	ILLUMINACLIP:${nextera_adapters}:2:30:10:2:keepBothReads \
	LEADING:3 SLIDINGWINDOW:4:15 MINLEN:40	
	sleep 5s
	"""
}

process ADD_UMI {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), file(trim1), file(trim2)
	output:
		tuple val (Sample), file ("${Sample}.unmapped.bam")
	script:	
	"""
	fgbio -Xmx${task.memory.toGiga()}g --tmp-dir=. --async-io=true --compression 1 FastqToBam --input ${trim1} ${trim2} --read-structures 8M+T +T --sample ${Sample} --library ${Sample} --output ${Sample}.unmapped.bam
	"""
}

process MAPBAM {
	tag "${Sample}"
	label 'process_medium'	
	input:
		tuple val (Sample), file(unmapped_bam)
		path (GenFile)
		path (GenDir)
	output:
		tuple val (Sample), path("${unmapped_bam}"), path("${Sample}_mapped.bam")
	script:
	"""
	#  Align the data with bwa and recover headers and tags from the unmapped BAM
	samtools fastq ${unmapped_bam} | bwa mem -t ${task.cpus} -p -K 150000000 -Y ${GenFile} - | samtools view -b -@ ${task.cpus} -o ${Sample}_mapped.bam -
	"""
} 

process ZIPPERBAM {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val(Sample), path(unmappedbam), path(mappedbam)
		path (GenFile)
		path (GenDir)
	output:
		tuple val(Sample), path("${Sample}.mapped.bam")
	script:
	"""
	fgbio -Xmx${task.memory.toGiga()}g --compression 1 --async-io ZipperBams \
	--input ${mappedbam} --unmapped ${unmappedbam} --ref ${GenFile} --output ${Sample}.mapped.bam
	"""
}

process SORT_INDEX {
	tag "${Sample}"
	label 'process_low'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
		tuple val(Sample), path(unsortedbam)
	output:
		tuple val(Sample), path("${Sample}.uncollaps.bam"), path("${Sample}.uncollaps.bam.bai")
	script:
	"""
	samtools sort -@ ${task.cpus} ${unsortedbam} -o ${Sample}.uncollaps.bam
	samtools index ${Sample}.uncollaps.bam > ${Sample}.uncollaps.bam.bai
	"""
}

process GROUPREADSBYUMI {
	tag "${Sample}"
	label 'process_medium'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: "${Sample}_family.png"
	input:
		tuple val (Sample), file(mapped_bam), file(mapped_bamBai)
	output: 
		tuple val (Sample), file ("*.grouped.bam"), emit : grouped_bam_ch
	script:
	"""
	fgbio -Xmx${task.memory.toGiga()}g --async-io --compression 1 GroupReadsByUmi --input ${mapped_bam} --strategy Adjacency --edits 1 --output ${Sample}.grouped.bam --family-size-histogram ${Sample}.tag-family-sizes.txt
	"""
}

process CALLMOLCONSREADS {
	tag "${Sample}"
	label 'process_high'
	input:
		tuple val (Sample), file (grouped_bam)
		path (GenFile)
		path (GenDir)
	output:
		tuple val (Sample), file ("*.cons.unmapped.bam")
	script:
	"""
	# This step generates unmapped consensus reads from the grouped reads and immediately filters them
	fgbio -Xmx${task.memory.toGiga()}g --compression 0 CallMolecularConsensusReads --input ${grouped_bam} --output /dev/stdout --min-reads 2 --threads ${task.cpus} | \
	fgbio -Xmx${task.memory.toGiga()}g --compression 1 FilterConsensusReads --input /dev/stdin --output ${Sample}.cons.unmapped.bam --ref ${GenFile} --min-reads 2 --min-base-quality 20 --max-base-error-rate 0.25
	"""
}

process MAPBAM_CONS {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), path(unmapped_bam)
		path (GenFile)
		path (GenDir)
	output:
		tuple val(Sample), path("${unmapped_bam}"), path("${Sample}_mapped.bam")
	script:
	"""
	# Align the data with bwa and recover headers and tags from the unmapped BAM
	samtools fastq ${unmapped_bam} | \
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200\\tPU:UNIT_1" -t ${task.cpus} -p -K 150000000 -Y ${GenFile} - | \
	samtools view -b -@ ${task.cpus} -o ${Sample}_mapped.bam -
	""" 
}

process FILTERCONSBAM {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), path (cons_unmapped_bam), path(cons_mapped_bam)
		path (GenFile)
		path (GenDir)		
	output:
		tuple val (Sample), file ("${Sample}_cons.bam")
	script:
	"""
	fgbio -Xmx${task.memory.toGiga()}g --compression 0 --async-io ZipperBams --input ${cons_mapped_bam} --unmapped ${cons_unmapped_bam} --ref ${GenFile} \
	--tags-to-reverse Consensus --tags-to-revcomp Consensus --output ${Sample}_cons.bam
	"""
}

process ADDGROUPS {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), path(cons_mapped_bam)
	output:
		tuple val (Sample), path("${Sample}_filt.bam")
	script:
	"""
	gatk AddOrReplaceReadGroups I=${cons_mapped_bam} O=${Sample}_filt.bam RGID=AML RGLB=LIB-MIPS RGPU=UNIT_1 RGPL=ILLUMINA RGSM=${Sample}
	"""
}

process SORT_INDEX_CONS {
	tag "${Sample}"
	label 'process_low'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
		tuple val(Sample), path(unsortedbam)
	output:
		tuple val(Sample), path("${Sample}_cons_sortd.bam"), path("${Sample}_cons_sortd.bam.bai")
	script:
	"""
	samtools sort -@ ${task.cpus} ${unsortedbam} -o ${Sample}_cons_sortd.bam
	samtools index ${Sample}_cons_sortd.bam > ${Sample}_cons_sortd.bam.bai
	"""
}

process COVERAGE {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), file (bam), file (bai)
		file (bedfile)
		val (bamtype)
	output:
		tuple val (Sample), file ("${Sample}_${bamtype}_coverage.bed")
	script:
	"""
	bedtools coverage -counts -a ${bedfile} -b ${bam} > "${Sample}_${bamtype}_coverage.bed"
	"""	
}

process DICT_GEN {
	label 'process_low'
	input:
		path(genfile)
	output:
		path("*.dict")
	script:
	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}g" CreateSequenceDictionary -R ${genfile}
	"""
}

process MUTECT2 {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), path(bam), path(bai)
		path (bedfile)
		path (GenFile)
		path (dict)
		path (GenDir)
		val (bamtype)
	output:
		tuple val(Sample), path("${Sample}_${bamtype}_mutect2.vcf")
	script:
	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}g" Mutect2 -R ${GenFile} -I ${bam} \
	-O ${Sample}_${bamtype}_mutect2.vcf -L ${bedfile} --native-pair-hmm-threads ${task.cpus} -mbq 20 \
	--max-reads-per-alignment-start 0 --af-of-alleles-not-in-resource 1e-6
	"""
}

process VARDICT {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), path(bam), path(bai)
		path (bedfile)
		path (GenFile)
		path (dict)
		path (GenDir)
		val (bamtype)
	output:
		tuple val(Sample), path("${Sample}.vardict.vcf")
	script:
	"""
	vardict-java -G ${GenFile} -th ${task.cpus} -f 0.0001 -r 8 -N ${Sample} -b ${bam} -c 1 -S 2 -E 3 -g 4 ${bedfile} | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.0001 > ${Sample}.vardict.vcf
	"""
}

process MPILEUP {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val(Sample), path(bam), path(bai)
		path (bedfile)
		path (GenFile)
		path (dict)
		path (GenDir)
		val (bamtype)
	output:
		tuple val(Sample), path("${Sample}.mpileup")
	script:
	"""
	samtools mpileup -d 1000000 -A -l ${bedfile} -f ${GenFile} ${bam} > ${Sample}.mpileup
	"""
}

process VARSCAN {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), path(bam), path(bai), path(mpileup)
		path (bedfile)
		path (GenFile)
		path (dict)
		path (GenDir)
		val (bamtype)
	output:
		tuple val(Sample), path("${Sample}_varscan.vcf")
	script:
	"""
	java -jar /usr/local/share/varscan-2.4.4-0/VarScan.jar mpileup2cns ${mpileup} --min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1  --variants 1 > ${Sample}_varscan.vcf
	"""
}

process ANNOVAR {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), path(Vcf)
		val(variant_caller)
	output:
		tuple val (Sample), path("${Sample}_${variant_caller}.out.hg19_multianno.csv")
	script:
	"""
	convert2annovar.pl -format vcf4 ${Vcf} --outfile ${Sample}_${variant_caller}.avinput --withzyg --includeinfo
	table_annovar.pl ${Sample}_${variant_caller}.avinput --out ${Sample}_${variant_caller}.out --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout \
	--thread ${task.cpus} /databases/humandb -xreffile /databases/gene_fullxref.txt

	touch ${Sample}_${variant_caller}.out.hg19_multianno.csv
	"""
}


process COMBINE_CALLERS {
	tag "${Sample}"
	label 'process_single'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
			tuple val (Sample), path(mutect2), path(vardict), path(varscan), path(coverage)
			val (bamtype)
	output:
			tuple val (Sample), path ("${Sample}_${bamtype}.xlsx")
	script:
	"""
	format_v2_mutect2.py ${mutect2} mutect2_.csv
	format_v2_varscan.py ${varscan} varscan_.csv
	format_v2_vardict.py ${vardict} vardict_.csv
	combine_callers.py ${Sample}_${bamtype}.xlsx mutect2_.csv varscan_.csv vardict_.csv ${coverage}
	"""
}

process ERROR_CORRECTN {
	tag "${Sample}"
	label 'process_single'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:	
		tuple val (Sample), path (collapsed_excel)
	output:
		tuple val (Sample), path ("${Sample}_ErrorCorrected.xlsx")
	script:
	"""
	background_error.py --input_collapsed_excel ${collapsed_excel} --output_excel ${Sample}_ErrorCorrected.xlsx
	"""		
}
