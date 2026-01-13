process VARDICT {
	tag "${Sample}"
	label 'process_inter'
	input:
		tuple val(Sample), path(bam), path(bai)
		path (bedfile)
		path (GenFile)
		path (dict)
		path (GenDir)
		val (bamtype)
	output:
		tuple val(Sample), path("${Sample}.vardict.vcf")
	script:
	"""
	java -Xmx${task.memory.toGiga()}g -jar /usr/local/share/vardict-java-1.8.3-0/lib/VarDict-1.8.3.jar -G ${GenFile} -th ${task.cpus} -f 0.0001 -r 8 -N ${Sample} -b ${bam} -c 1 -S 2 -E 3 -g 4 ${bedfile} | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.0001 > ${Sample}.vardict.vcf
	"""
}
