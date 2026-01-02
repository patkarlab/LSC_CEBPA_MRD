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