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