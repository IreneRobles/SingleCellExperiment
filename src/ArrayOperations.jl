function scale_between_0_1(array_numbers)
    array_numbers
    MAX = maximum(array_numbers)
    MIN = minimum(array_numbers)
    f(x) = (x - MIN)/(MAX-MIN)
    return [f(i) for i in array_numbers]
end

##########################################################
## FUNCTIONS TO APPLY FUNCTIONS TO GENE OR CELL VECTORS ##
##########################################################

#20thSept2018

function apply_function_cell_vectors(sceexp::SingleCellExp, f; assay = "CPM")
    [f(sceexp.assays[assay][:, ii]) for ii in 1:n_cells(sceexp)]
end

function apply_function_gene_vectors(sceexp::SingleCellExp, f; assay = "CPM")
    [f(sceexp.assays[assay][ii, :]) for i in 1:n_genes(sceexp)]
end


#################################################
## FUNCTIONS TO MANIPULATE ROWDATA AND COLDATA ##
#################################################
function add_genecounts_to_colData(singlecellexperiment::SingleCellExp, gene; assay = "lnCPMplus1")
    
    gene_ind = find_gene_index(gene, singlecellexperiment)
    gene_counts = singlecellexperiment.assays[assay][gene_ind, :]
    singlecellexperiment.colData[!,Symbol(gene)] = gene_counts
    
    return singlecellexperiment
end
    
function total_count_and_genes_per_cell(singlecellexperiment::SingleCellExp)
    coldata = singlecellexperiment.colData
    cols = coldata[!,:CellID]
    # Count the total number of transcripts detected
    totalcounts = [sum(singlecellexperiment.counts[:, i]) for i in 1:size(singlecellexperiment.counts)[2]]
    # Count the number of genes expressed in the cell
    genesdetected = [count(!iszero, singlecellexperiment.counts[:, i]) for i in 1:size(singlecellexperiment.counts)[2]]
    
    coldata[!,:TotalCounts] = totalcounts
    coldata[!,:TotalGenes] = genesdetected
    
    singlecellexperiment.colData = coldata
    
    return singlecellexperiment
end

function total_cells_per_gene(singlecellexperiment::SingleCellExp)
    rowdata = singlecellexperiment.rowData
    # Count the total number of transcripts detected
    totalcounts = [sum(singlecellexperiment.counts[i, :]) for i in 1:size(singlecellexperiment.counts)[1]]
    # Count the number of cells with at least one transcript detected
    cellsdetected = [count(x -> x != 0, singlecellexperiment.counts[i, :]) for i in 1:size(singlecellexperiment.counts)[1]]
    
    rowdata[!,:TotalCounts] = totalcounts
    rowdata[!,:TotalCells] = cellsdetected
    rowdata[!,:FractionCells] = [i/size(singlecellexperiment.counts)[1] for i in cellsdetected]
    
    singlecellexperiment.rowData = rowdata
    
    return singlecellexperiment
end
    
function total_cells_per_gene_by(singlecellexperiment::SingleCellExp, bycolumn)
    singlecellexperiment = SingleCellExpJulia.total_cells_per_gene(singlecellexperiment)
    
    bysplits = unique(singlecellexperiment.colData[bycolumn])
    
    for bysplit in bysplits
        bool = singlecellexperiment.colData[bycolumn].== bysplit
        subset = SingleCellExpJulia.filter_cells(bool, singlecellexperiment)
        subset = SingleCellExpJulia.total_cells_per_gene(subset)
        singlecellexperiment.rowData[Symbol(string(bysplit, "_TotalCounts"))] =  subset.rowData[!,:TotalCounts]
        singlecellexperiment.rowData[Symbol(string(bysplit, "_TotalCells"))] =  subset.rowData[!,:TotalCells]
        singlecellexperiment.rowData[Symbol(string(bysplit, "_FractionCells"))] =  subset.rowData[!,:FractionCells]
    end
    
    return singlecellexperiment
end


    
function total_cells_per_gene_rowDataby(singlecellexperiment::SingleCellExp, bycolumn)
    singlecellexperiment = SingleCellExpJulia.total_cells_per_gene(singlecellexperiment)
    
   bysplits = unique(singlecellexperiment.rowData[bycolumn])
    
    for bysplit in bysplits
        bool = singlecellexperiment.rowData[bycolumn].== bysplit
        subset = SingleCellExpJulia.filter_genes(bool, singlecellexperiment)
        subset = SingleCellExpJulia.total_count_and_genes_per_cell(subset)
        singlecellexperiment.colData[!,Symbol(string(bysplit, "__TotalCounts"))] =  subset.colData[!,:TotalCounts]
        singlecellexperiment.colData[!,Symbol(string(bysplit, "__TotalGenes"))] =  subset.colData[!,:TotalGenes]
        singlecellexperiment.colData[!,Symbol(string(bysplit, "__FractionCounts"))] =  subset.colData[!,:TotalCounts] ./ singlecellexperiment.colData[!,:TotalCounts]
        singlecellexperiment.colData[!,Symbol(string(bysplit, "__FractionGenes"))] =  subset.colData[!,:TotalGenes] ./ singlecellexperiment.colData[!,:TotalGenes]
    end
    
    return singlecellexperiment
end

function get_max_and_min_value_per_gene(singlecellexperiment; 
        assay = "CPM", 
        mincol = Symbol(string(assay, "__Min")), 
        maxcol = Symbol(string(assay, "__Max")))
        
    Assay = singlecellexperiment.assays[assay]

    singlecellexperiment.rowData[!,maxcol] = [maximum(Assay[i, :]) for i in 1:size(Assay)[1]]
    singlecellexperiment.rowData[!,mincol] = [minimum(Assay[i, :]) for i in 1:size(Assay)[1]]
    return singlecellexperiment
    
end  

function percent_mt(sceexp; pat = "mt")
    boolmt = startswith.(sceexp.rownames, pat) 
    mtcounts = sceexp.counts[boolmt, :]
    counts = [sum(sceexp.counts[:,ii]) for ii in 1:size(sceexp.counts)[2]]
    mtcounts = [sum(mtcounts[:,ii]) for ii in 1:size(mtcounts)[2]]
    sceexp.colData[!,"percent.mt"] = mtcounts./counts.*100
    return sceexp
end


###################################
## FUNCTIONS TO TRANSFORM COUNTS ## 
###################################

function cpm_transform(singlecellexperiment)
    # 28th Feb 2018
    # Calculate the size factors (total counts per cell)
    scexp = total_count_and_genes_per_cell(singlecellexperiment)
    cpm = deepcopy(scexp.counts)
    
    siz_fact = scexp.colData[!,:TotalCounts]
    
    for cell in 1:singlecellexperiment.dim[2]
        # Divide counts by size factor and then multiply by 10e6 to obtain cpm (counts per million)
        cpm[:, cell] = [i*10e6/siz_fact[cell] for i in cpm[:, cell]]
    end
    
    singlecellexperiment.assays["CPM"] = cpm
    
    return singlecellexperiment 
end

function ln_cpm_plus1_transform(singlecellexperiment)
    # 1st March 2018
    # Calculate the size factors (total counts per cell)
    scexp = total_count_and_genes_per_cell(singlecellexperiment)
    cpm = deepcopy(scexp.counts)
    
    siz_fact = scexp.colData[!,:TotalCounts]
    
    for cell in 1:singlecellexperiment.dim[2]
        # Divide counts by size factor and then multiply by 10e6 to obtain cpm (counts per million)
        cpm[:, cell] = [log(i*10e6/siz_fact[cell] + 1) for i in cpm[:, cell]]
    end
    
    singlecellexperiment.assays["lnCPMplus1"] = cpm
    
    return singlecellexperiment 
end


    
function cpm_transform_assay(sceexp, assay)
    counts_per_cell = [sum(sceexp.assays[assay][:, ii]) for ii in 1:size(sceexp.assays[assay])[2]]
    genes_per_cell = [count(x -> x != 0 ,sceexp.assays[assay][:, ii]) for ii in 1:size(sceexp.assays[assay])[2]]
    sceexp.colData[!,Symbol("TotalCounts_"*assay)] = counts_per_cell
    sceexp.colData[!,Symbol("TotalGenes_"*assay)] = genes_per_cell
    
    cpm = float(deepcopy(sceexp.assays[assay]))
    siz_fact = sceexp.colData[!,Symbol("TotalCounts_"*assay)]
        
    for cell in 1:sceexp.dim[2]
        # Divide counts by size factor and then multiply by 10e6 to obtain cpm (counts per million)
        cpm[:, cell] = [ii*10e6/siz_fact[cell] for ii in cpm[:, cell]]
    end
    
    sceexp.assays["CPM_"*assay] = cpm
        
    return sceexp
end

function ln_cpm_plus1_transform_assay(sceexp, assay)
    counts_per_cell = [sum(sceexp.assays[assay][:, ii]) for ii in 1:size(sceexp.assays[assay])[2]]
    genes_per_cell = [count(x -> x != 0 ,sceexp.assays[assay][:, ii]) for ii in 1:size(sceexp.assays[assay])[2]]
    sceexp.colData[!,Symbol("TotalCounts_"*assay)] = counts_per_cell
    sceexp.colData[!,Symbol("TotalGenes_"*assay)] = genes_per_cell
    
    cpm = float(deepcopy(sceexp.assays[assay]))
    siz_fact = sceexp.colData[!,Symbol("TotalCounts_"*assay)]
        
    for cell in 1:sceexp.dim[2]
        # Divide counts by size factor and then multiply by 10e6 to obtain cpm (counts per million)
        cpm[:, cell] = [log(i*10e6/siz_fact[cell] + 1) for i in cpm[:, cell]]
    end
    
    sceexp.assays["lnCPMplus1_"*assay] = cpm
        
    return sceexp
end

    
function compress_between_1to0(value, min, max)
    if Float64(min) != Float64(max)
        return  (value-min)/(max-min)
    else 
        return 0.0
    end
end


function transform_compressexpression1to0(singlecellexperiment; assay = "CPM")
    # Compress the expression range between 1 and 0
    assayname = string(assay, "_1to0")
    singlecellexperiment = get_max_and_min_value_per_gene(singlecellexperiment, assay = assay)

    new_assay = deepcopy(singlecellexperiment.assays[assay])
    
    mincol = Symbol(string(assay, "__Min"))
    maxcol = Symbol(string(assay, "__Max"))
    
    mins = singlecellexperiment.rowData[mincol]
    maxs = singlecellexperiment.rowData[maxcol]
    
    for gene in 1:size(new_assay)[1] 
        min = mins[gene]
        max = maxs[gene]
        
         new_assay[gene, :] = [compress_between_1to0(value, min, max) for value in new_assay[gene, :]]
        
    end
    singlecellexperiment.assays[assayname] = new_assay
    
    return singlecellexperiment 
end
    
function calculate_cv_cellexp(sceexp; min_cells_expressing_gene = 0, assay = "lnCPMplus1")
    
    scep = deepcopy(SingleCellExpJulia.select_expressed_genes(sceexp,min_cells_expressing_gene = min_cells_expressing_gene))
    ncell = scep.dim[2]
    #Coefficient of variation only measures on the genes detected in each cell
   sceexp.colData[:CV_gene_expression] = [variation(scep.assays[assay][.!iszero.(scep.assays[assay][:, ii]), ii]) for ii in 1:ncell]
    return sceexp
end


    
function add_assay(sceexp, newassay; assay_name = "new_assay", geneid_col = :GeneID, cellid_col = :CellID)
    # geneid_col is the column in the counts dataframe that contains the gene identifiers
    # cellid_col is the column in the colData dataframe that contains the cell identifiers
    
    newassay = dropmissing(newassay)
    
    # Find out which cells and wich genes appear in both the assay and the sceexp
    cells = intersect(sceexp.colData[cellid_col], names(newassay))
    genes = intersect(sceexp.rowData[!,:GeneID], newassay[geneid_col])
    
    # Reorder cells and genes in the new_assay
    new_assay_gene_ind = [findfirst(x-> x ==ii, newassay[geneid_col]) for ii in genes]
    new_assay_cell_ind = [findfirst(x-> x ==ii, names(newassay)) for ii in cells]
    new_assay = newassay[new_assay_gene_ind, cells]
    
    # Reorder cells and genes in the sce
    
    new_indexes_cells = [findfirst(x -> x == ii, sceexp.colData[:CellID]) for ii in cells]
    new_indexes_genes = [findfirst(x -> x == ii, sceexp.rowData[:GeneID]) for ii in genes]
    
    new_sce = reorder_cells!(sceexp, new_indexes_cells)
    new_sce = reorder_genes!(new_sce, new_indexes_genes)
    
    # Add reordered assay to sceexp
    
    new_sce.assays[assay_name] = Matrix(new_assay)
    
    return new_sce
    
end

