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
