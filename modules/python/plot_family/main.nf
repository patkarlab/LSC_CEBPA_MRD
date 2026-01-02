process PLOT_FAMILY_SIZE {
	tag "${Sample}"
	label 'process_single'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: "${Sample}_family.png"
	input:
		tuple val (Sample), file(family_sizes)
	output: 
		tuple val (Sample), file ("${Sample}_family.png")
	script:
	"""
	plot_family_sizes.py ${Sample}.tag-family-sizes.txt ${Sample}_family
	"""
}