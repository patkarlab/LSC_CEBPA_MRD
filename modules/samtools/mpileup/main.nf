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