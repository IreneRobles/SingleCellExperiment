module SingleCellExperiment

using CSV
using PyCall, DataFrames
using NoLongerProblems, NoLongerProblems_Pandas, ScikitLearn
import Pandas
using Seaborn, PrettyPlotting
using Statistics, MultipleTesting, HypothesisTests,ProgressMeter

@sk_import linear_model: LogisticRegression


export SingleCellExp
export genes, show_gene_data
export apply_function_cell_vectors, apply_function_gene_vectors
export n_cells, n_genes
export find_gene_index, find_cell_index,find_cell_condition 
export get_gene_distribution, get_cells_with_this_characteristic, get_cells_with_this_characteristics
export get_sample, get_genotype
export total_cells_per_gene_rowDataby
export cpm_transform, ln_cpm_plus1_transform, transform_compressexpression1to0, add_assay
export cpm_transform_assay, ln_cpm_plus1_transform_assay
export number_cells_per_sample, total_count_and_genes_per_cell, total_cells_per_gene
export cpm_transform, add_assay

include("SingleCellExperimentObject.jl")
include("ArrayOperations.jl")
include("ModuleScore.jl")
include("PlottingFunctions.jl")

#include("SingleCellExperimentJulia.jl")

end # module
