process GROUPREADSBYUMI {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), file(mapped_bam), file(mapped_bamBai)
	output: 
		tuple val (Sample), file ("*.grouped.bam"), emit : grouped_bam_ch
		tuple val (Sample), file("${Sample}.tag-family-sizes.txt"), emit : family_sizes
	script:
	"""
	fgbio -Xmx${task.memory.toGiga()}g --async-io --compression 1 GroupReadsByUmi --input ${mapped_bam} --strategy Adjacency --edits 1 -@ ${task.cpus} --output ${Sample}.grouped.bam --family-size-histogram ${Sample}.tag-family-sizes.txt
	"""
}
