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

