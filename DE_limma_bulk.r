library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(limma)
library(edgeR)
source('RA_asgard_scripts/voomByGroup.r')
set.seed(2024)

bulk <- LoadH5Seurat("/home/icb/zihe.zheng/projects/asgard/RA_data/rna_obj_clinical_bulk.h5Seurat")

# Get case and control sample names
Control <- unique(bulk@meta.data$sample[bulk@meta.data$disease == 'OA'])
Case <- unique(bulk@meta.data$sample[bulk@meta.data$disease == 'RA'])

print(Control)
print(Case)

# Initialize lists to store results
Gene.list <- list()
C_names <- c()

# Define the minimum number of cells required for analysis
min.cells <- 3  # Adjust based on your requirements

# Take cell type index as an argument
args <- commandArgs(trailingOnly = TRUE)
idx <- as.numeric(args[1])

i <- c("B cell/plasma cell", "T cell", "Stromal cell", "Endothelial cell", "Myeloid cell", "NK")[idx]
# print(unique(bulk@meta.data$cell_type))
print(i)
Idents(bulk) <- "cell_type"
c_cells <- subset(bulk, cell_type == i)

Idents(c_cells) <- "disease"
Samples <- c_cells@meta.data

colnames(Samples)

Controlsample <- row.names(subset(Samples, sample %in% Control))
Casesample <- row.names(subset(Samples, sample %in% Case))

if (length(Controlsample) > min.cells & length(Casesample) > min.cells) {
    new_expr <- c_cells@assays$RNA@data  # Keep in sparse format

    # Check the size of the expression matrix
    print(dim(new_expr))

    new_sample <- data.frame(Samples = c(Casesample, Controlsample), 
                             type = Samples[c(Casesample, Controlsample), "disease"],
                             age = Samples[c(Casesample, Controlsample), "age"],
                             sex = Samples[c(Casesample, Controlsample), "sex"],
                             site = Samples[c(Casesample, Controlsample), "site"]
                            )
    row.names(new_sample) <- new_sample$Samples
    
    # Remove genes with zero variance across all samples
    non_zero_var_genes <- apply(new_expr, 1, var) != 0
    new_expr <- new_expr[non_zero_var_genes, ]

    # Model matrix for limma with adjustments for age, sex and site
    design <- model.matrix(~0 + type + age + sex + site, data = new_sample)

    colnames(design) <- make.names(colnames(design))
    print(head(design))
    
    # voomByGroup transformation
    y_vbg <- voomByGroup(new_expr, design = design, group = new_sample$type, plot = "combine")
    
    # Fit the linear model
    fit <- lmFit(y_vbg, design)
    contr <- makeContrasts(typeRA - typeOA, levels = colnames(coef(fit)))
    tmp <- contrasts.fit(fit, contrasts = contr)
    tmp <- eBayes(tmp)
    
    # Extract top differentially expressed genes
    C_data <- topTable(tmp, sort.by = "P", n = nrow(tmp))
    # print(head(C_data))
    C_data_for_drug <- data.frame(row.names = row.names(C_data), 
                                  score = C_data$t, 
                                  adj.P.Val = C_data$adj.P.Val, 
                                  P.Value = C_data$P.Value)
    Gene.list[[i]] <- C_data_for_drug
    C_names <- c(C_names, i)
}

names(Gene.list) <- C_names
if (i == 'B cell/plasma cell'){i = 'B_plasma_cell'}
filename <- sprintf('/home/icb/zihe.zheng/projects/asgard/RA_data/DE_gene_list_bulk_%s.rds', i)
saveRDS(Gene.list, filename)
