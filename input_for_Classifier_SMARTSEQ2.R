args = commandArgs(trailingOnly = TRUE)

input_file      <- args[1]
output_file     <- args[2]

if (length(args) < 2) {
	stop("Aborting. Insufficient number of arguments. The script requires atleast 2 arguments -- the input counts matrix and an output file to which the normalized counts are written")	
}

library("DESeq2")

counts_matrix <- read.table(input_file,header = TRUE, sep = "\t", row.names = 1)

## Arbitrarily split the cells into 2 groups, this is to ensure that a dds object can be created.
group_1_list = rep("group1", ceiling(length(colnames(counts_matrix))/2))
group_2_list = rep("group2", length(colnames(counts_matrix)) - ceiling(length(colnames(counts_matrix))/2))
all_cell_groups = c(group_1_list, group_2_list)
cell_groups   <- factor(all_cell_groups,levels = unique(all_cell_groups))

dds <- DESeqDataSetFromMatrix(counts_matrix,DataFrame(cell_groups), ~cell_groups)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- estimateSizeFactors(dds)

## Do the normalised transformation
nt <- normTransform(dds) # defaults to log2(x+1)

### check if the rownames are Ensembl Identifiers or gene symbols:
genes_list_symbols = c("XCL2", "TNFRSF9", "XCL1", "HSP90AB1", "PRDX1", "PARK7", "CRTAM", "HSPA8", "FABP5", "MIR155HG")

if (grepl("^ENSG",rownames(counts_matrix)[1])) {
	genes_list = c("ENSG00000143185","ENSG00000049249","ENSG00000143184","ENSG00000096384","ENSG00000117450","ENSG00000116288","ENSG00000109943","ENSG00000109971","ENSG00000164687","ENSG00000234883")
	nt_counts <- assay(nt)[genes_list,]	
	nt_counts_selected = nt_counts[genes_list,]
	rownames(nt_counts_selected) = genes_list_symbols
	write.table(nt_counts_selected, file = output_file, sep = ",", quote = FALSE, col.names = NA)
} else {
	nt_counts <- assay(nt)[genes_list,]	
	write.table(nt_counts[genes_list_symbols,], file = output_file, sep = ",", quote = FALSE, col.names = NA)
}
