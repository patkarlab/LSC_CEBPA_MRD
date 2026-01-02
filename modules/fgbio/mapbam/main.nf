process MAPBAM {
	tag "${Sample}"
	label 'process_medium'	
	input:
		tuple val (Sample), file(unmapped_bam)
		path (GenFile)
		path (GenDir)
	output:
		tuple val (Sample), path("${unmapped_bam}"), path("${Sample}_mapped.bam")
	script:
	"""
	#  Align the data with bwa and recover headers and tags from the unmapped BAM
	samtools fastq ${unmapped_bam} | bwa mem -t ${task.cpus} -p -K 150000000 -Y ${GenFile} - | samtools view -b -@ ${task.cpus} -o ${Sample}_mapped.bam -
	"""
} 

process MAPBAM_CONS {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), path(unmapped_bam)
		path (GenFile)
		path (GenDir)
	output:
		tuple val(Sample), path("${unmapped_bam}"), path("${Sample}_mapped.bam")
	script:
	"""
	# Align the data with bwa and recover headers and tags from the unmapped BAM
	samtools fastq ${unmapped_bam} | \
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200\\tPU:UNIT_1" -t ${task.cpus} -p -K 150000000 -Y ${GenFile} - | \
	samtools view -b -@ ${task.cpus} -o ${Sample}_mapped.bam -
	""" 
}