mutable struct SingleCellExp
    dim
    counts
    metadata
    assays
    rownames
    rowData
    colnames
    colData
    reducedDimNames
    spikeNames
end

function Base.show(io::IO, z::SingleCellExp)
    # 1st March 2018
    println("SingleCellExp")
    genes = size(z.counts)[1]
    cells = size(z.colData)[1]
    println("Cells = $cells, Genes = $genes")
    print("assays :")
    for key in keys(z.assays)
        print(" $key")
    end
    println()
    print("colData :")
    for key in names(z.colData)
        print(" $key")
    end
    println()
    print("rowData :")
    for key in names(z.rowData)
        print(" $key")
    end
    println()
    print("reducedDimNames :")
    for key in keys(z.reducedDimNames)
        print(" $key")
    end
end

function SingleCellExp(counts, colData; 
    # geneid_col is the column in the counts dataframe that contains the gene identifiers
    # cellid_col is the column in the colData dataframe that contains the cell identifiers
    geneid_col = :GeneSymbol, cellid_col = :RowName, 
    )
    
    rownames = counts[!, geneid_col]
    colnames = []
    
    try 
        colnames = [Symbol(replace(string(i), "." => "_"))  for i in colData[!,cellid_col]]

        colData[!,:CellID] = colnames

        # This bit reorder the columns in counts to have the same order as colData
    counts = Array{Float64, 2}(counts[:, colnames])
        catch 
        
         colnames = [Symbol(i)  for i in colData[!,cellid_col]]
        colData[!,:CellID] = colnames
        # This bit reorder the columns in counts to have the same order as colData
        counts = Array{Float64, 2}(counts[!, colnames])
    end
    metadata = Dict()
    assays = Dict()
    
    rowData = DataFrames.DataFrame()
    rowData[!,:GeneID] = rownames
    

    colData = colData
    reducedDimNames = Dict()
    spikeNames = Dict()
    
    dim = length(rownames), length(colnames)
    
    SingleCellExp(
    
    dim,
    counts,
    metadata,
    assays,
    rownames,
    rowData,
    colnames,
    colData,
    reducedDimNames,
    spikeNames
    
    )
end


######################################################
## FUNCTIONS TO ACCESS TO SINGLECELLEXPERIMENT DATA ##
######################################################

function genes(singlecellexperiment::SingleCellExp)
    # 2nd March 2018
    return singlecellexperiment.rownames
end

function show_gene_data(sceexp::SingleCellExp, gene)
    # Get gene data
    f(x) = [i == gene for i in x]
    d = sceexp.rowData[f(sceexp.rowData[:GeneID]), :]
    # Get samples
    cols = [string(i) for i in names(d)]
    samples = []
    measures = []
    for col in cols
        s = split(col, "__")
        if length(s) > 1
            push!(samples, s[1])
            push!(measures, s[2])
        end
    end
    
    usamples = unique(samples)
    umeasures = [Symbol(i) for i in unique(measures)]
    
    df = DataFrames.DataFrame()
    df[!,:Sample] = usamples
    
    for umeasure in umeasures
        dat = []
        for usample in usamples
            push!(dat, d[1, Symbol(string(usample, "__", umeasure))])
        end
        df[umeasure] = dat
    end
    df
    
end


###############################################################
## FUNCTIONS TO GET THE NUMBER OF CELLS OR GENES IN THE DATA ##
###############################################################

function n_cells(singlecellexperiment::SingleCellExp)
    # 14th March 2018
    singlecellexperiment.dim[2]
end

function n_genes(singlecellexperiment::SingleCellExp)
    # 14th March 2018
    singlecellexperiment.dim[1]
end

############################################################################
## FUNCTIONS TO FIND GENE AND CELL INDEXES IN SINGLECELLEXP OBJECT ##
############################################################################

function find_gene_index(geneid, singlecelldata)
    # 1st March 2018
    genesid = singlecelldata.rownames
    geneindex = findfirst(x -> x == geneid,genesid)
    return geneindex
end

function find_cell_index(cellid, singlecelldata)
    # 1st March 2018
    cellsid = singlecelldata.colData[!,:CellID]
    cellindex = findfirst(x -> x == cellid,cellsid)
    return cellindex
end

########################################################
## FUNCTIONS TO FIND OUT PROPERTIES OF CELL AND GENES ##
########################################################

function find_cell_condition(cellid, conditioncol, singlecelldata)
    # 1st March 2018
    ind = find_cell_index(cellid, singlecelldata)
    return singlecelldata.colData[ind, conditioncol]
end

####################################################
## FUNCTIONS TO SORT SINGLECELLEXP OBJECT ##
###################################################

export sort_cells!, reorder_cells!, sort_genes!, reorder_genes!

function sort_cells!(singlecellexperiment::SingleCellExp; kwargs...)
    celldatabefore = singlecellexperiment.colData
    
    celldataafter = celldatabefore
    
    try
    celldataafter = sort(celldatabefore, kwargs[:cols], rev = kwargs[:rev])
    catch
    celldataafter = sort(celldatabefore, kwargs[:cols])
    end

    
    new_indexes_cells = [findfirst(x -> x == i, celldatabefore[!,:CellID]) for i in celldataafter[!,:CellID]]
    
    new_sce = deepcopy(singlecellexperiment)
    
    # Update colData and colnames
    new_sce.colData = singlecellexperiment.colData[new_indexes_cells, :]
    new_sce.colnames = singlecellexperiment.colData[new_indexes_cells, :CellID]
    
    # Update counts
    new_sce.counts = new_sce.counts[:, new_indexes_cells]
   
    # Update Assays
    for assay in keys(new_sce.assays)
        new_sce.assays[assay] = new_sce.assays[assay][:, new_indexes_cells]
    end
    
    # Update Dimensionality reductions
    for dim in keys(new_sce.reducedDimNames)
        new_sce.reducedDimNames[dim] = new_sce.reducedDimNames[dim][new_indexes_cells, :]
    end
    
    # Update dims
    new_sce.dim = length(new_sce.rownames), length(new_sce.colnames)
    
    return new_sce
end


function reorder_cells!(singlecellexperiment::SingleCellExp, new_indexes_cells)  
    # new_indexes_cells is an array which contains the previous indexes, in the order that we want them after
    new_sce = deepcopy(singlecellexperiment)
    
    # Update colData and colnames
    new_sce.colData = singlecellexperiment.colData[new_indexes_cells, :]
    new_sce.colnames = singlecellexperiment.colData[new_indexes_cells, :CellID]
    
    # Update counts
    new_sce.counts = new_sce.counts[:, new_indexes_cells]
   
    # Update Assays
    for assay in keys(new_sce.assays)
        new_sce.assays[assay] = new_sce.assays[assay][:, new_indexes_cells]
    end
    
    # Update Dimensionality reductions
    for dim in keys(new_sce.reducedDimNames)
        new_sce.reducedDimNames[dim] = new_sce.reducedDimNames[dim][:, new_indexes_cells]
    end
    
    # Update dims
    new_sce.dim = length(new_sce.rownames), length(new_sce.colnames)
    
    return new_sce
end

function sort_genes!(singlecellexperiment::SingleCellExp; kwargs...)
    
    genedatabefore = singlecellexperiment.rowData
    genedataafter = genedatabefore
    
     try
    genedataafter = sort(genedatabefore, kwargs[:cols], rev = kwargs[:rev])
    catch
    genedataafter = sort(genedatabefore, kwargs[:cols])
    end

    
    new_indexes_genes = [findfirst(x -> x == i, genedatabefore[:CellID]) for i in genedataafter[!,:CellID]]
    
    new_sce = deepcopy(singlecellexperiment)
    
    # Update colData and colnames
    new_sce.rowData = singlecellexperiment.rowData[new_indexes_genes, :]
    new_sce.rownames = singlecellexperiment.rowData[new_indexes_genes, :GeneID]
    
    # Update counts
    new_sce.counts = new_sce.counts[new_indexes_genes, :]
   
    # Update Assays
    for assay in keys(new_sce.assays)
        new_sce.assays[assay] = new_sce.assays[assay][new_indexes_genes, :]
    end
    
    # Update Dimensionality reductions
    for dim in keys(new_sce.reducedDimNames)
        new_sce.reducedDimNames[dim] = new_sce.reducedDimNames[dim][new_indexes_genes, :]
    end
    
    # Update dims
    new_sce.dim = length(new_sce.rownames), length(new_sce.colnames)
    
    return new_sce

    
end


function reorder_genes!(singlecellexperiment::SingleCellExp, new_indexes_genes)
    # new_indexes_genes is an array which contains the previous indexes, in the order that we want them after
    
    new_sce = deepcopy(singlecellexperiment)
    
    # Update colData and colnames
    new_sce.rowData = singlecellexperiment.rowData[new_indexes_genes, :]
    new_sce.rownames = singlecellexperiment.rowData[new_indexes_genes, :GeneID]
    
    # Update counts
    new_sce.counts = new_sce.counts[new_indexes_genes, :]
   
    # Update Assays
    for assay in keys(new_sce.assays)
        new_sce.assays[assay] = new_sce.assays[assay][new_indexes_genes, :]
    end
    
    # Update Dimensionality reductions
    for dim in keys(new_sce.reducedDimNames)
        new_sce.reducedDimNames[dim] = new_sce.reducedDimNames[dim][new_indexes_genes, :]
    end
    
    # Update dims
    new_sce.dim = length(new_sce.rownames), length(new_sce.colnames)
    
    return new_sce

end

#############################
## Save Data to load in R ##
############################
    
function drop_non_unique_geneid(singlecellexperiment)
    isnotduplicate = .!nonunique(singlecellexperiment.rowData, :GeneID)
    new_sce = filter_genes(isnotduplicate, singlecellexperiment)
end


function SingleCellExperiment.save_sce_for_R(singlecellexperiment; dir = pwd())
    mkpath(dir)
    singlecellexperiment = SingleCellExperiment.drop_non_unique_geneid(singlecellexperiment)
    CSV.write(normpath(dir, "rowData.csv"), singlecellexperiment.rowData)
    CSV.write(normpath(dir, "colData.csv"), singlecellexperiment.colData)
    CSV.write(normpath(dir, "counts.csv"), Tables.table(singlecellexperiment.counts), writeheader=false)
    CSV.write(normpath(dir, "geneid.csv"), Tables.table(singlecellexperiment.rownames), writeheader=false)
end



######################################################
## FUNCTIONS TO SUBSET SINGLECELLEXP OBJECT ##
#####################################################

export filter_cells, filter_genes, select_these_genes, select_expressed_genes

function filter_cells(cells_bool, singlecellexperiment::SingleCellExp)
    # 2nd March 2018
    # cells_bool is an Array{Bool} which says which cells should be included and excluded in the same order as it was
    
    new_sce = deepcopy(singlecellexperiment)
    
    coldata = singlecellexperiment.colData
    
    # Create a Boolean of the cells that should be included
    col_bool = cells_bool
    
    # Update colData and colnames
    new_sce.colData = coldata[col_bool, :]
    new_sce.colnames = coldata[col_bool, :CellID]
    
    # Update counts
    new_sce.counts = new_sce.counts[:, col_bool]
   
    # Update Assays
    for assay in keys(new_sce.assays)
        new_sce.assays[assay] = new_sce.assays[assay][:, col_bool]
    end
    
    # Update Dimensionality reductions
    for dim in keys(new_sce.reducedDimNames)
        new_sce.reducedDimNames[dim] = new_sce.reducedDimNames[dim][col_bool, :]
    end
    
    # Update dims
    new_sce.dim = length(new_sce.rownames), length(new_sce.colnames)
    
    return new_sce

end

function filter_genes(genes_bool, singlecellexperiment::SingleCellExp)
    # 2nd March 2018
    # genes_bool is an Array{Bool} which says which genes should be included and excluded in the same order as it is
    #  in the SingleCellExp object
    
    new_sce = deepcopy(singlecellexperiment)
    
    # Create a Boolean of the cells that should be included
    col_bool = genes_bool
    
    # Update colData and colnames
    new_sce.rowData = new_sce.rowData[col_bool, :]
    new_sce.rownames = new_sce.rownames[col_bool]
    
    # Update counts
    new_sce.counts = new_sce.counts[col_bool, :]
   
    # Update Assays
    for assay in keys(new_sce.assays)
        new_sce.assays[assay] = new_sce.assays[assay][col_bool, :]
    end
    
    # Update dims
    new_sce.dim = length(new_sce.rownames), length(new_sce.colnames)
    
    return new_sce

end


function select_these_genes(gene_list, singlecellexperiment)
    # Returns a singlecellexperiment with only the genes that appear in the gene list. 
    # If there is no match between a gene and the singlecellexperiment.rownames
    # it will not return an error, simply it will not return that gene.

    
    # This is slower
    # bool = [in(i, gene_list) for i in genes(singlecellexperiment)]
    
    # This is faster
    # 26th Sept 18
   
    bool = fill(false,singlecellexperiment.dim[1])
    
    gene_inds = Dict(singlecellexperiment.rownames, 1:singlecellexperiment.dim[1])
    
    for gene in gene_list
        try
            index = gene_inds[gene]
            bool[index] = true
        catch
            
        end
    end
    
    new_sce = filter_genes(bool, singlecellexperiment)
    return new_sce
    
end

function select_expressed_genes(sceexp; min_cells_expressing_gene = 114)
    
    sceexp = total_cells_per_gene(sceexp)
    f(x) = x > min_cells_expressing_gene
    bool  = [f(i) for i in sceexp.rowData[!,:TotalCells]] 
    new_sce = filter_genes(bool, sceexp)
    
end
###################################################

function get_gene_distribution(geneid, condition, conditioncol::Symbol, singlecelldata; assay = "lnCPMplus1")
    new_sce = get_cells_with_this_characteristic(condition, conditioncol, singlecelldata)
    gene_index = find_gene_index(geneid, singlecelldata)
    gene_dist = new_sce.assays[assay][gene_index, :]
end


function get_cells_with_this_characteristic(characteristic, characteristiccolumn::Symbol, singlecellexperiment::SingleCellExp)
    # characteristiccolumn must be a column in singlecellexperimment.colData
    # characteristic is what all the cells must have in common
    
    new_sce = deepcopy(singlecellexperiment)
    
    coldata = singlecellexperiment.colData
    
    # Create a Boolean of the cells that should be included
    col_bool = coldata[!,characteristiccolumn] .== characteristic
    
    # Update colData and colnames
    new_sce.colData = coldata[col_bool, :]
    new_sce.colnames = coldata[col_bool, :CellID]
    
    # Update counts
    new_sce.counts = new_sce.counts[:, col_bool]
   
    # Update Assays
    for assay in keys(new_sce.assays)
        new_sce.assays[assay] = new_sce.assays[assay][:, col_bool]
    end
    
    # Update Dimensionality reductions
    # Fixed 17 July 18
    for dim in keys(new_sce.reducedDimNames)
        new_sce.reducedDimNames[dim] = new_sce.reducedDimNames[dim][col_bool, :]
    end
    
    # Update dims
    new_sce.dim = length(new_sce.rownames), length(new_sce.colnames)
    
    return new_sce

end

function get_cells_with_this_characteristics(characteristics, characteristiccolumn::Symbol, singlecellexperiment::SingleCellExp)
    # characteristiccolumn must be a column in singlecellexperimment.colData
    # characteristic is what all the cells must have in common
    
    new_sce = deepcopy(singlecellexperiment)
    
    coldata = singlecellexperiment.colData
    
    # Create a Boolean of the cells that should be included
    col_bool = [in(i, characteristics) for i in coldata[!,characteristiccolumn]]
    
    # Update colData and colnames
    new_sce.colData = coldata[col_bool, :]
    new_sce.colnames = coldata[col_bool, :CellID]
    
    # Update counts
    new_sce.counts = new_sce.counts[:, col_bool]
   
    # Update Assays
    for assay in keys(new_sce.assays)
        new_sce.assays[assay] = new_sce.assays[assay][:, col_bool]
    end
    
    # Update Dimensionality reductions
    # Fixed 17 July 18
    for dim in keys(new_sce.reducedDimNames)
        new_sce.reducedDimNames[dim] = new_sce.reducedDimNames[dim][col_bool, :]
    end
    
    # Update dims
    new_sce.dim = length(new_sce.rownames), length(new_sce.colnames)
    
    return new_sce

end


function get_sample(sampl::Union{String,SubString{String}}, singlecellexperiment::SingleCellExp)
    # 8Jul18
    # sampl::String changed to sampl::Union{String,SubString{String}} 
    get_cells_with_this_characteristic(sampl,:Sample, singlecellexperiment)
end

# Returns all the cells with the same genotype
function get_genotype(genotype::String, singlecellexperiment::SingleCellExp)
    get_cells_with_this_characteristic(genotype, :Genotype, singlecellexperiment)
end



function number_cells_per_sample(df::DataFrames.DataFrame; sergitimepointsample = false)
    # This version of the code works with just a colData provided, o need for a SCE object
    coldata = df
    # Calulate the number of cells per sample
    df_summary = DataFrames.DataFrame()
    df_summary[!,:Sample] = unique(coldata[!,:Sample])
    
    if sergitimepointsample == true
        df_summary[!,:Timepoint] = [split(i, "_")[end] for i in df_summary[!,:Sample]]
        df_summary[!,:Genotype] = [split(i, "_")[1] for i in df_summary[!,:Sample]]
    end
    
    function get_s(str, df)
        coldata = df
        col_bool = [in(i, [str]) for i in coldata[!,:Sample]]
        return coldata[col_bool, :]
    end

    df_summary[!,:Number_Cells] = [size(get_s(string(i), df))[1] for i in df_summary[!,:Sample]]

    df_summary
end
    
function number_cells_per_sample(singlecellexperiment::SingleCellExp; sergitimepointsample = false)
    # This one works for a SCE object
    number_cells_per_sample(singlecellexperiment.colData)
end
    
