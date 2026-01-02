process ERROR_CORRECTN {
	tag "${Sample}"
	label 'process_single'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:	
		tuple val (Sample), path (collapsed_excel)
	output:
		tuple val (Sample), path ("${Sample}_ErrorCorrected.xlsx")
	script:
	"""
	background_error.py --input_collapsed_excel ${collapsed_excel} --output_excel ${Sample}_ErrorCorrected.xlsx
	"""		
}