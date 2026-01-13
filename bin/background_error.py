#!/usr/bin/env python
# This script will calculate the background error as described in Blood. 2018;132(16):1703-1713

import argparse
import pandas as pd
import statistics
import re

def median_cal(input_list):
	# Remove -1 from the list
	return float(statistics.median([ values for values in input_list if values >= 0 ]))

def SNP_INDEL(ref, alt, min_length):
	ref = re.sub("-","", ref, flags = re.IGNORECASE)
	alt = re.sub("-","", alt, flags = re.IGNORECASE)
	# Segregating variants based on length
	if abs(len(alt) - len(ref)) <= min_length:
		return 0
	else:
		return 1

def MRD_conditions(updated_df, background_error, background_std_dev):
	lvaf = "LVAF%"
	vartype = "variant_type"
	limit = background_error + ( 3 * background_std_dev )
	mrd_status = "MRD_STATUS"
	peak_limit = 3 * background_std_dev
	background_column = "BACKGROUND + 3SD"
	
	for index, row in updated_df.iterrows():
		if row[vartype] == 0:
			updated_df.loc[index, background_column] = limit
			# Condition 1 for MRD status
			if row[lvaf] > limit: 
				updated_df.loc[index, mrd_status] = "PASS"
				# Check condition 2 for deciding MRD status
				# lower_base_location = int(row['Start']) - 20
				# upper_base_location = int(row['Start']) + 20
				# peak_vaf_limit = row[lvaf] - peak_limit
				# subset_df = updated_df[
				# 	 (updated_df['Start'] > lower_base_location) & 
				# 	 (updated_df['Start'] < upper_base_location) & 
				# 	 (updated_df[vartype] == 0) & 
				# 	 ( updated_df[lvaf] > peak_vaf_limit)]
				# if subset_df.empty:
				# 	updated_df.loc[index, mrd_status] = "PASS"
				# else:
				# 	updated_df.loc[index, mrd_status] = "FAIL, conditn 2"
			else:
				updated_df.loc[index, mrd_status] = "FAIL"
		else:
			updated_df.loc[index, mrd_status] = "NA"
	return updated_df

def BackgroundRate(args):
	infile = args.input_collapsed_excel
	outfile = args.output_excel
	sheet_to_modify = 'Variants'
	median_VAF = "VAF%_median"
	median_REF = "REF_COUNT_median"
	median_ALT = "ALT_COUNT_median"
	mrd_status = "MRD_STATUS"
	background_column = "BACKGROUND + 3SD"
	lvaf = "LVAF%"
	min_INDEL_length = 2 # The original value as per the paper was 3
	vartype = "variant_type"

	df = pd.read_excel(infile, sheet_name=None)
	merged_df = df[sheet_to_modify].copy()	# Making a copy of the original sheet
	merged_df[median_VAF] = int(-1)
	merged_df[median_REF] = int(-1)
	merged_df[median_ALT] = int(-1)
	merged_df[mrd_status] = int(-1)			# default mrd status
	merged_df[background_column] = float(-1)	# New column for Background
	merged_df[lvaf] = int(-1)
	merged_df[vartype] = int(0)				# default variant type = 0 (SNP)

	for index, row in merged_df.iterrows():
		ref_vals = [row['REF_COUNT_Mutect2'], row['REF_COUNT_VarScan2'], row['REF_COUNT_VarDict']]
		alt_vals = [row['ALT_COUNT_Mutect2'], row['ALT_COUNT_VarScan2'], row['ALT_COUNT_VarDict']]
		vaf_vals = [row['VAF%_Mutect2'], row['VAF%_VarScan2'], row['VAF%_VarDict']]
		ref = row['Ref']
		alt = row['Alt']
		merged_df.loc[index, median_REF] = median_cal(ref_vals)
		merged_df.loc[index, median_ALT] = median_cal(alt_vals)
		median_vaf = median_cal(vaf_vals)
		merged_df.loc[index, median_VAF] = median_vaf
	# merged_df = merged_df[ df[sheet_to_modify].columns.to_list() + [median_REF,median_ALT, median_VAF, background_column]]
	# Sorting the dataframe based on VAF values at positions to obtain the LVAF
	merged_df.sort_values(by=['Chr', 'Start', median_VAF ], inplace=True, ascending=[True, True, False])

	chr_list = []
	position_list = dict()
	for index, row in merged_df.iterrows():
		chr = re.sub("chr","", row['Chr'], flags = re.IGNORECASE)	# Removing chr from the chromosome id
		pos = int(row['Start'])
		ref = row['Ref']
		alt = row['Alt']
		median_vaf = row[median_VAF]
		variant_type = SNP_INDEL(ref, alt, min_INDEL_length)
		merged_df.loc[index, vartype] = variant_type
		if (variant_type == 0):
			if chr in chr_list:
				chr_list.append(chr)
				# Update the LVAF at each position only for INDELS less than equal to min_INDEL_length
				if pos in position_list:
					if median_vaf > position_list[pos]:
						position_list[pos] = median_vaf
					else:
						pass
				else:
					position_list[pos] = median_vaf
			else:
				if position_list:
					average = statistics.mean(position_list.values())
					print ("average of LVAFs for chr",chr, "is", average)
				chr_list = []
				position_list = dict()	# 2 dictionaries one each for SNP+short indels and large indels
				chr_list.append(chr)
				position_list[pos] = median_vaf
			merged_df.loc[index,lvaf] = position_list[pos]
		else:
			merged_df.loc[index,lvaf] = median_vaf

	average = statistics.mean(position_list.values())
	standard_dev = statistics.stdev(position_list.values())
	lower_limit = average - (2.5 * standard_dev)
	upper_limit = average + (2.5 * standard_dev)
	# print ("average of LVAFs at the end of df for chr",chr, "is", average)
	# print ("std dev of LVAFs at the end of df for chr",chr, "is", standard_dev)
	# print ("The limits are", lower_limit, upper_limit, average, standard_dev, min(position_list.values()), max(position_list.values()))

	# In this loop we are assuming that only one chr (chr 19) is being analysed
	chr_list = []
	position_list = dict()
	subset_lvafs = []
	for index, row in merged_df.iterrows():
		chr = re.sub("chr","", row['Chr'], flags = re.IGNORECASE)	# Removing chr from the chromosome id
		pos = int(row['Start'])
		if row[vartype] == 0:				# Extracting only the SNPs+short INDELs
			if pos not in position_list:	# Taking only the first occurence of LVAFs
				position_list[pos] = row[lvaf]
				# Filtering LVAFs which deviate 2.5 times the std dev from the average
				if row[lvaf] > lower_limit and row[lvaf] < upper_limit:
					subset_lvafs.append(row[lvaf])
			else:
				pass
		else:
			pass

	background_error_rate = statistics.mean(subset_lvafs)
	background_std_dev = statistics.stdev(subset_lvafs)

	modified_df = MRD_conditions(merged_df, background_error_rate, background_std_dev)
	unique_cols = [background_column, mrd_status]
	modified_df = modified_df[ df[sheet_to_modify].columns.to_list() + unique_cols ]

	with pd.ExcelWriter(outfile , engine="openpyxl", mode="w") as writer:
		modified_df.to_excel(writer, sheet_name=sheet_to_modify, index=False)
		for sheet_name, df in df.items():
			if sheet_name != sheet_to_modify:	# Not printing the original sheet
				df.to_excel(writer, sheet_name=sheet_name, index=False)

def parse_arguments():
	parser = argparse.ArgumentParser(description="Calculation of background error")
	parser.add_argument('--input_collapsed_excel', help='output of combine_callers.py for collapsed bam')		# 
	parser.add_argument('--output_excel', help='same format as input')					# 
	return parser.parse_args()

def main ():
	args = parse_arguments()
	BackgroundRate(args)

if __name__ == "__main__":
	main()
