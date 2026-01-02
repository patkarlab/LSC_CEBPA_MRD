process CALLMOLCONSREADS {
	tag "${Sample}"
	label 'process_high'
	input:
		tuple val (Sample), file (grouped_bam)
		path (GenFile)
		path (GenDir)
	output:
		tuple val (Sample), file ("*.cons.unmapped.bam")
	script:
	"""
	# This step generates unmapped consensus reads from the grouped reads and immediately filters them
	fgbio -Xmx${task.memory.toGiga()}g --compression 0 CallMolecularConsensusReads --input ${grouped_bam} --output /dev/stdout --min-reads 2 --threads ${task.cpus} | \
	fgbio -Xmx${task.memory.toGiga()}g --compression 1 FilterConsensusReads --input /dev/stdin --output ${Sample}.cons.unmapped.bam --ref ${GenFile} --min-reads 2 --min-base-quality 20 --max-base-error-rate 0.25
	"""
}