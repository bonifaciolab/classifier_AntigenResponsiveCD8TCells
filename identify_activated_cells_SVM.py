## Virag Sharma, June 2019
## A wrapper around the lib-SVM package to idenfity the activated cells given a matrix containing normalised counts of the following genes:
## XCL2,TNFRSF9,XCL1,HSP90AB1,PRDX1,PARK7,CRTAM,HSPA8,FABP5,MIR155HG (see Fuchs et al. https://www.biorxiv.org/content/10.1101/673707v1 )

import sys
import subprocess
import os

def run_command(command, return_val = False):
	output = subprocess.run(command, shell = True, capture_output = True, text = True)
	if output.returncode != 0:
		print("Error running the call {}. Aborting".format(command))
		sys.exit()
		
	if return_val == True:
		return output.stdout.rstrip()
####################################################################
def replace_commas_with_tab(input_file, output_file):
	
	with open(input_file, "r") as FI, open(output_file, "w") as FO:
		for line in FI:
			line = line.replace(',','\t')
			FO.write(line)
####################################################################	
def get_transpose(input_file, output_file):
	## from https://gist.github.com/mikeboers/4997162

	all_data = []
	with open(input_file, "r") as FI:
		for line in FI:
			all_data.append(line.strip().split('\t'))
	# Transpose it.
	transposed = []
	for row_i in range(len(all_data[0])):
		new_row = []
		for col_i in range(len(all_data)):
			new_row.append(all_data[col_i][row_i])
		transposed.append(new_row)
	# and print the transpose to a file:
	with open(output_file, "w") as FO:
		for row in transposed:
			row_formatted = '\t'.join(row)
			row_formatted = row_formatted.replace('"','')			
			FO.write(row_formatted + "\n")
####################################################################
def format_line(line_list):
	''' formats a line so that lib-SVM is happy '''

	label = "-1"  ## randomly assign -1 to all cells
	tmp_new = list()
	tmp_new.append(label)
	ct = 1
	for element in line_list:
		element = "{}:{}".format(ct, element)
		tmp_new.append(element)
		ct += 1

	return " ".join(tmp_new)
####################################################################	
def label_svm_output(svm_output, cell_labels, final_output):

	ct = 0
	with open(svm_output, "r") as FI, open(final_output, "w") as FO:	
		for pred in FI:
			pred = int(pred.rstrip())
			if pred == 1:
				FO.write("{}\t1\n".format(cell_labels[ct]))
			else:
				FO.write("{}\t0\n".format(cell_labels[ct]))
			ct += 1

####################################################################

if __name__ == "__main__":

	input_file = sys.argv[1]
	output_dir = sys.argv[2]
	svm_path   = sys.argv[3]

	if not os.path.isdir(output_dir):
		os.mkdir(output_dir)

	input_1 = "{}/input_tab_separated_1.txt".format(output_dir)
	input_2 = "{}/input_transposed_2.txt".format(output_dir)
	input_3 = "{}/input_formatted_3.txt".format(output_dir)
	input_4 = "{}/input_formatted_scaled_4.txt".format(output_dir)
	final_output = "{}/svm_output.txt".format(output_dir)
	cells_labelled = "{}/cell_labels_predicted.txt".format(output_dir)

	## step 1: replace the commas in the input-line with tabs
	replace_commas_with_tab(input_file, input_1)

	## step 2: get transpose of the input file
	get_transpose(input_1, input_2)

	## get cell labels/barcodes from the transposed matrix:
	cell_labels = dict()

	cell_label_ct = 0
	with open(input_2, "r") as FI, open(input_3, "w") as FO:
		for line in FI:
			line = line.rstrip()
			if "XCL1" in line:
				continue
			
			line_list = line.split('\t')
			barcode = line_list.pop(0)
			cell_labels[cell_label_ct] = barcode
			cell_label_ct += 1

			line_formatted = format_line(line_list)
			FO.write(line_formatted + "\n")

	## step 4: scale the input data, run svm-scale
	svm_scale_call = "{}/svm-scale {} > {}".format(svm_path, input_3, input_4)
	run_command(svm_scale_call)

	## step 5: run svm-predict
	svm_predict_call = "{}/svm-predict {} training_model {}".format(svm_path, input_4, final_output)
	run_command(svm_predict_call)

	## step 6: label the predictions
	label_svm_output(final_output, cell_labels, cells_labelled)
