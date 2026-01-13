process COVERAGE {
	tag "${Sample}"
	label 'process_inter'
	input:
		tuple val (Sample), file (bam), file (bai)
		file (genome_index)
		file (bedfile)
		val (bamtype)
	output:
		tuple val (Sample), file ("${Sample}_${bamtype}_coverage.bed")
	script:
	"""
	bedtools coverage -counts -sorted -g ${genome_index} -a ${bedfile} -b ${bam} > "${Sample}_${bamtype}_coverage.bed"
	"""	
}

