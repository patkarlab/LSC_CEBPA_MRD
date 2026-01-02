process VARSCAN {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), path(bam), path(bai), path(mpileup)
		path (bedfile)
		path (GenFile)
		path (dict)
		path (GenDir)
		val (bamtype)
	output:
		tuple val(Sample), path("${Sample}_varscan.vcf")
	script:
	"""
	java -jar /usr/local/share/varscan-2.4.4-0/VarScan.jar mpileup2cns ${mpileup} --min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1  --variants 1 > ${Sample}_varscan.vcf
	"""
}