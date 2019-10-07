args = commandArgs(trailingOnly = TRUE)

input_dir   = args[1]
output_file = args[2]

if (length(args) < 2) {
	stop("Aborting. Insufficient number of arguments. The script requires atleast 2 arguments -- the input directory produced by Cell-Ranger and an output file to which the normalized counts are written")	
}

library(Seurat)

input.data <- Read10X(data.dir = input_dir)
input_object <- CreateSeuratObject(counts = input.data[["Gene Expression"]], project = "my_project", min.cells = 8, min.features = 200)
## these values can be changed depending on the input dataset
mito.features <- grep(pattern = "^MT-", x = rownames(x = input_object), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = input_object, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = input_object, slot = 'counts'))

input_object[['percent.mito']] <- percent.mito
input_object <- subset(x = input_object, subset = nFeature_RNA > 500 & percent.mito < 0.1) ## setting a gate, this could be variable based on the dataset

# run sctransform
input_object <- SCTransform(input_object, vars.to.regress = "percent.mito", assay = "RNA")

# Subsetting information for activation genes for the classifier
subset_marker_genes.matrix <- input_object@assays$SCT@data[c("XCL2", "TNFRSF9", "XCL1", "HSP90AB1", "PRDX1", "PARK7", "CRTAM", "HSPA8", "FABP5", "MIR155HG"),]
write.csv(subset_marker_genes.matrix, file = output_file)
