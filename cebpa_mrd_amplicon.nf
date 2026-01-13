#!/usr/bin/env nextflow
nextflow.enable.dsl=2


genome_file = file("${params.genome}", checkIfExists: true)
index_files = file("${params.genome_dir}/${params.ind_files}.*")
fai_index_file = file("${params.genome_dir}/${params.ind_files}.fai")
bedfile = file("${params.bedfile}", checkIfExists: true)
truseq_adapters = file("${params.truseq_adapters}", checkIfExists: true)
nextera_adapters = file("${params.nextera_adapters}", checkIfExists: true)
collapsed = params.collapsed
uncollapsed = params.uncollapsed
mutect2 = params.mutect2
vardict = params.vardict
varscan = params.varscan


include { TRIM } from './modules/trimmomatic/main.nf'
include { ADD_UMI } from './modules/fgbio/fastqtobam/main.nf'
include { MAPBAM; MAPBAM_CONS } from './modules/fgbio/mapbam/main.nf'
include { ZIPPERBAM } from './modules/fgbio/zipperbams/main.nf'
include { SORT_INDEX; SORT_INDEX_CONS } from './modules/samtools/sort_index/main.nf'
include { GROUPREADSBYUMI } from './modules/fgbio/groupreads/main.nf'
include { PLOT_FAMILY_SIZE } from './modules/python/plot_family/main.nf'
include { CALLMOLCONSREADS } from './modules/fgbio/callmolecularconsensus/main.nf'
include { FILTERCONSBAM } from './modules/fgbio/fiterconsbam/main.nf'
include { ADDGROUPS } from './modules/gatk/addreplacegroups/main.nf'
include { COVERAGE as COVERAGE_COLL; COVERAGE as COVERAGE_UNCOLL } from './modules/bedtools/coverage/main.nf'
include { MUTECT2_COLL; MUTECT2_UNCOLL; DICT_GEN } from './modules/gatk/mutect2/main.nf'
include { VARDICT as VARDICT_COLL; VARDICT as VARDICT_UNCOLL} from './modules/vardict/variantcall/main.nf'
include { MPILEUP as MPILEUP_COLL; MPILEUP as MPILEUP_UNCOLL } from './modules/samtools/mpileup/main.nf'
include { VARSCAN as VARSCAN_COLL; VARSCAN as VARSCAN_UNCOLL } from './modules/varscan/variantcall/main.nf'
include { ANNOVAR as ANNOVAR_COLL_MUTECT2; ANNOVAR as ANNOVAR_COLL_VARDICT; ANNOVAR as ANNOVAR_COLL_VARSCAN } from './modules/annovar/annotate/main.nf'
include { ANNOVAR as ANNOVAR_UNCOLL_MUTECT2; ANNOVAR as ANNOVAR_UNCOLL_VARDICT; ANNOVAR as ANNOVAR_UNCOLL_VARSCAN } from './modules/annovar/annotate/main.nf'
include { COMBINE_CALLERS as COMBINE_CALLERS_COLL; COMBINE_CALLERS as COMBINE_CALLERS_UNCOLL } from './modules/python/combine/main.nf'
include { ERROR_CORRECTN } from './modules/python/error_correctn/main.nf'


workflow CEBPA_MRD {
	Channel.fromPath(params.input)
		.splitCsv(header:false)
		.map { row ->
			def sample = row[0].trim()
			def r1 = file("${params.sequences}/${sample}_S*_R1_*.fastq.gz", checkIfExists: false)
			def r2 = file("${params.sequences}/${sample}_S*_R2_*.fastq.gz", checkIfExists: false)

			if (!r1 && !r2) {
				r1 = file("${params.sequences}/${sample}*_R1.fastq.gz", checkIfExists: false)
				r2 = file("${params.sequences}/${sample}*_R2.fastq.gz", checkIfExists: false)
			}
			tuple(sample, r1, r2)
		}
		.set { samples_ch }

	main:

		DICT_GEN(genome_file)
		TRIM(samples_ch, truseq_adapters, nextera_adapters)
		ADD_UMI(TRIM.out) 
		MAPBAM(ADD_UMI.out, genome_file, index_files)
		ZIPPERBAM(MAPBAM.out, genome_file, index_files)
		SORT_INDEX(ZIPPERBAM.out)
		GROUPREADSBYUMI(SORT_INDEX.out)
		PLOT_FAMILY_SIZE(GROUPREADSBYUMI.out.family_sizes)
		CALLMOLCONSREADS(GROUPREADSBYUMI.out.grouped_bam_ch, genome_file, index_files)
		MAPBAM_CONS(CALLMOLCONSREADS.out, genome_file, index_files)
		FILTERCONSBAM(MAPBAM_CONS.out, genome_file, index_files)
		ADDGROUPS(FILTERCONSBAM.out)
		SORT_INDEX_CONS(ADDGROUPS.out)

		COVERAGE_COLL(SORT_INDEX_CONS.out, fai_index_file, bedfile, collapsed)
		COVERAGE_UNCOLL(SORT_INDEX.out, fai_index_file, bedfile, uncollapsed)

		MUTECT2_COLL(SORT_INDEX_CONS.out, bedfile, genome_file, DICT_GEN.out, index_files, collapsed)
		MUTECT2_UNCOLL(SORT_INDEX.out, bedfile, genome_file, DICT_GEN.out, index_files, uncollapsed)

		VARDICT_COLL(SORT_INDEX_CONS.out, bedfile, genome_file, DICT_GEN.out, index_files, collapsed)
		VARDICT_UNCOLL(SORT_INDEX.out, bedfile, genome_file, DICT_GEN.out, index_files, uncollapsed)

		MPILEUP_COLL(SORT_INDEX_CONS.out, bedfile, genome_file, DICT_GEN.out, index_files, collapsed)
		VARSCAN_COLL(SORT_INDEX_CONS.out.join(MPILEUP_COLL.out), bedfile, genome_file, DICT_GEN.out, index_files, collapsed)
		MPILEUP_UNCOLL(SORT_INDEX.out, bedfile, genome_file, DICT_GEN.out, index_files, uncollapsed)
		VARSCAN_UNCOLL(SORT_INDEX.out.join(MPILEUP_UNCOLL.out), bedfile, genome_file, DICT_GEN.out, index_files, uncollapsed)

		ANNOVAR_COLL_MUTECT2(MUTECT2_COLL.out, mutect2)
		ANNOVAR_COLL_VARDICT(VARDICT_COLL.out, vardict)
		ANNOVAR_COLL_VARSCAN(VARSCAN_COLL.out, varscan)

		ANNOVAR_UNCOLL_MUTECT2(MUTECT2_UNCOLL.out, mutect2)
		ANNOVAR_UNCOLL_VARDICT(VARDICT_UNCOLL.out, vardict)
		ANNOVAR_UNCOLL_VARSCAN(VARSCAN_UNCOLL.out, varscan)

		COMBINE_CALLERS_COLL(ANNOVAR_COLL_MUTECT2.out.join(ANNOVAR_COLL_VARDICT.out.join(ANNOVAR_COLL_VARSCAN.out.join(COVERAGE_COLL.out))), collapsed)
		COMBINE_CALLERS_UNCOLL(ANNOVAR_UNCOLL_MUTECT2.out.join(ANNOVAR_UNCOLL_VARDICT.out.join(ANNOVAR_UNCOLL_VARSCAN.out.join(COVERAGE_UNCOLL.out))), uncollapsed)

		ERROR_CORRECTN(COMBINE_CALLERS_COLL.out)
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
	println "Completed at: ${workflow.complete}"
	println "Total time taken: ${workflow.duration}"
}
