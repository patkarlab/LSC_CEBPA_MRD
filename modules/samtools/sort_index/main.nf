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