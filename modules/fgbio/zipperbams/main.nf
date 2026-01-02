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