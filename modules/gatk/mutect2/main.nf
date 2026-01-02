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