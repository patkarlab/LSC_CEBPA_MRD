process MUTECT2_UNCOLL {
	tag "${Sample}"
	label 'process_mutect'
	input:
		tuple val(Sample), path(bam), path(bai)
		path (bedfile)
		path (GenFile)
		path (dict)
		path (GenDir)
		val (bamtype)
	output:
		tuple val(Sample), path("${Sample}_${bamtype}_mutect2.vcf")
	script:
	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}g" Mutect2 -R ${GenFile} -I ${bam} \
	-O ${Sample}_${bamtype}_mutect2.vcf -L ${bedfile} --native-pair-hmm-threads 1 -mbq 20 \
	--af-of-alleles-not-in-resource 1e-6
	"""
}

process MUTECT2_COLL {
	tag "${Sample}"
	label 'process_mutect'
	input:
		tuple val(Sample), path(bam), path(bai)
		path (bedfile)
		path (GenFile)
		path (dict)
		path (GenDir)
		val (bamtype)
	output:
		tuple val(Sample), path("${Sample}_${bamtype}_mutect2.vcf")
	script:
	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}g" Mutect2 -R ${GenFile} -I ${bam} \
	-O ${Sample}_${bamtype}_mutect2.vcf -L ${bedfile} --native-pair-hmm-threads 1 -mbq 20 \
	--max-reads-per-alignment-start 0 --af-of-alleles-not-in-resource 1e-6
	"""
}

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