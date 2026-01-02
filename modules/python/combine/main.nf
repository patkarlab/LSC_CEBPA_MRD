process COMBINE_CALLERS {
	tag "${Sample}"
	label 'process_single'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
			tuple val (Sample), val(mutect2_vartype), path(mutect2_anno), val(vardict_vartype), path(vardict_anno), val(varscan_vartype), path(varscan_anno), path(coverage)
			val (bamtype)
	output:
			tuple val (Sample), path ("${Sample}_${bamtype}.xlsx")
	script:
	"""
	format_v2_mutect2.py ${mutect2_anno} mutect2_.csv
	format_v2_varscan.py ${varscan_anno} varscan_.csv
	format_v2_vardict.py ${vardict_anno} vardict_.csv
	combine_callers.py ${Sample}_${bamtype}.xlsx mutect2_.csv varscan_.csv vardict_.csv ${coverage}
	"""
}