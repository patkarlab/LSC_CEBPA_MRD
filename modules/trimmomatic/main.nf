process TRIM {
	tag "${Sample}" 
	label 'process_low'
	input:
		tuple val (Sample), file (read1), file (read2)
		file (truseq_adapters)
		file (nextera_adapters)
	output:
		tuple val (Sample), file("*1P.fq.gz"), file("*2P.fq.gz")
	script:
	"""	
	java -Xmx${task.memory.toGiga()}g -jar /usr/local/share/trimmomatic-0.39-2/trimmomatic.jar PE \
	${read1} ${read2} \
	-baseout ${Sample}.fq.gz -threads ${task.cpus}\
	ILLUMINACLIP:${truseq_adapters}:2:30:10:2:keepBothReads \
	ILLUMINACLIP:${nextera_adapters}:2:30:10:2:keepBothReads \
	LEADING:3 SLIDINGWINDOW:4:15 MINLEN:40
	"""
}