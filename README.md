# Introduction:

classifier_AntigenResponsiveCD8TCells is a python wrapper around the libSVM package to identify antigen-responsive CD8+ Tcells from a pool of cells using the expression of 10 marker genes. See the accompanying preprint for more details: (https://www.biorxiv.org/content/10.1101/673707v1)

# Dependencies:

The classifier requires the following:
1. An installation of lib-SVM (Please see https://www.csie.ntu.edu.tw/~cjlin/libsvm/ for details)
2. python (version >= 3.7)
3. The following R-packages: DESeq2, Seurat, sctransform

# Installation:
Run the following on your terminal (has been tested on macOS Mojave -version 10.14.6 and Ubuntu 16.04.6 LTS)
git clone https://github.com/bonifaciolab/classifier_AntigenResponsiveCD8TCells/

# Pre-requisities:
The classifier requires a counts-matrix that contains information about the gene expression values of the marker genes
If a plate based method such as SMART-SEQ2 was used, the counts matrix (produced by aligning the raw reads to the reference human genome using an aligner such as STAR and associating counts to genes using featureCounts) is normalised using DESeq2's normTransform function

A generalized workflow using a Mac-OS or a Linux terminal is shown below
```
export input_counts= ## the path to the file containing the counts matrix
export input_counts_normalized=  ## the path to the file containing the normalized counts matrix
## create the counts matrix
Rscript input_for_Classifier_SMARTSEQ2.R $input_counts $input_counts_normalized
```

The input counts matrix could also come from the run of a 10X experiment (using the Cell Ranger software). In that case, the normalized counts matrix is generated using the Seurat package
```
export input_counts= ## the path to the file containing the counts matrix
export input_counts_normalized=  ## the path to the file containing the normalized counts matrix
## create the counts matrix
Rscript input_for_Classifier_10X.R $input_counts $input_counts_normalized
```
### Now you are ready to run the classifier
```
export output_dir= ## the directory where your results will be produced. If the directory does not exist, then it will be created.
export svm_path= ## the path to the directory that contains the svm binaries

python3 identify_activated_cells_SVM.py $input_counts_normalized $output_dir $svm_path
```
The output directory contains different files. The final output is contained in the file called cell_labels_predicted.txt
This file is a tab-separated file where each row has two fields: cell_barcode followed by a 0 or 1.
1 indicates that the cell is predicted to be responsive to the antigen, 0 indicates otherwise.
