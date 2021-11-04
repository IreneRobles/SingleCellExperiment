#########################################
## FUNCTIONS TO CALCULATE MODULE SCORE ## 
#########################################


function fit_mu_std_alpha(singlecellexperiment; splitdataby = :Sample, assay = "lnCPMplus1")
    # Fixed 21st March 2018
    conditions = unique(singlecellexperiment.colData[!,splitdataby])
    genes = singlecellexperiment.rowData[!,:GeneID]
    
    for condition in conditions
        condition_sce = get_cells_with_this_characteristic(condition, splitdataby, singlecellexperiment)
        ncells = condition_sce.dim[2]
        
        mu_name = Symbol(string(condition, "__mu"))
        std_name = Symbol(string(condition, "__std"))
        alpha_name = Symbol(string(condition, "__alpha"))
        average_name = Symbol(string(condition, "__averagepop"))
        var_name = Symbol(string(condition, "__var"))
        
        alphas = []
        mus = []
        stds = []
        average_pops = [] 
        vars = []
           
        for gene in genes
            gene_index = find_gene_index(gene, condition_sce)
            gene_dist = condition_sce.assays[assay][gene_index, :]
            # Find cells with transcripts detected 
            cells_tc_detected = findall(x -> x > 1, gene_dist)
            cells_tc_detected = gene_dist[cells_tc_detected]
            # Calculate alpha, fraction of cells with transcript detected
            alpha = length(cells_tc_detected)/ncells
            # From the cells that have transcripts detected calculate mu and std
            standardeviation = std(cells_tc_detected)
            varian = var(gene_dist)
            
            mu = mean(cells_tc_detected)
            average_pop = mean(gene_dist)
            if isnan(mu) == true
                mu = 0.0
                standardeviation = 0.0
            end
            
            push!(alphas, alpha)
            push!(mus, mu)
            push!(stds, standardeviation)
            push!(average_pops, average_pop)
            push!(vars, varian)
        end
        singlecellexperiment.rowData[!,mu_name] = mus
        singlecellexperiment.rowData[!,std_name] = stds
        singlecellexperiment.rowData[!,alpha_name] = alphas
        singlecellexperiment.rowData[!,average_name] = average_pops
        singlecellexperiment.rowData[!,var_name] = vars
    end
    condition = "All"
    condition_sce = singlecellexperiment
    ncells = condition_sce.dim[2]
        
        mu_name = Symbol(string(condition, "__mu"))
        std_name = Symbol(string(condition, "__std"))
        alpha_name = Symbol(string(condition, "__alpha"))
        average_name = Symbol(string(condition, "__averagepop"))
        var_name = Symbol(string(condition, "__var"))
        
        alphas = []
        mus = []
        stds = []
        average_pops = [] 
        vars = []
           
        for gene in genes
            gene_index = find_gene_index(gene, condition_sce)
            gene_dist = condition_sce.assays[assay][gene_index, :]
            # Find cells with transcripts detected 
            cells_tc_detected = findall(x -> x > 1, gene_dist)
            cells_tc_detected = gene_dist[cells_tc_detected]
            # Calculate alpha, fraction of cells with transcript detected
            alpha = length(cells_tc_detected)/ncells
            # From the cells that have transcripts detected calculate mu and std
            standardeviation = std(cells_tc_detected)
            varian = var(gene_dist)
            
            mu = mean(cells_tc_detected)
            average_pop = mean(gene_dist)
            if isnan(mu) == true
                mu = 0.0
                standardeviation = 0.0
            end
            
            push!(alphas, alpha)
            push!(mus, mu)
            push!(stds, standardeviation)
            push!(average_pops, average_pop)
            push!(vars, varian)
        end
        singlecellexperiment.rowData[!,mu_name] = mus
        singlecellexperiment.rowData[!,std_name] = stds
        singlecellexperiment.rowData[!,alpha_name] = alphas
        singlecellexperiment.rowData[!,average_name] = average_pops
        singlecellexperiment.rowData[!,var_name] = vars
    
    return singlecellexperiment
    
end
    

function alpha_table(singlecellexperiment::SingleCellExp, gene; 
        name = "alpha"
    )
    
    #This is to keep my color code
    pal = Dict(
            "WT" => 1, 
            "RAD21" => 2,
            "CTCF" => 4,
            "RAD21CTCF" => 3
        )
    
     gene_ind = find_gene_index(gene, singlecellexperiment)
    
    gene_data = singlecellexperiment.rowData[gene_ind, :]
    
     cols_interest = columns_containing(gene_data, name)
    
    gene_data = gene_data[cols_interest]
    
    new_dict = Dict()
    
    genotypes = []
    timepoints = []
    ds = []
    
    for col in cols_interest
        split_name = split(string(col), "_")
        
        pop!(split_name)
        pop!(split_name)
        
        genotype = split_name[1]
        push!(genotypes, genotype)
        timepoint = split_name[end]
         push!(timepoints, timepoint)
        
       d =  gene_data[col]
        
        push!(ds, d)
        
     
    end
    
    df = DataFrames.DataFrame()
    df[!,:Genotype] = genotypes
    df[!,:Timepoint] = timepoints
    df[!,:Sample] = [string(df[i, :Genotype], "_", df[i, :Timepoint]) for i in 1:nrow(df)]
    df[!,Symbol(string(name, "_", gene))] = ds
    df
end

function mu_table(singlecellexperiment::SingleCellExp, gene; 
        name = "mu"
    )
    alpha_table(singlecellexperiment, gene, name = name)
end
function averagepop_table(singlecellexperiment::SingleCellExp, gene; 
        name = "averagepop"
    )
    alpha_table(singlecellexperiment, gene, name = name)
end
function std_table(singlecellexperiment::SingleCellExp, gene; 
        name = "std"
    )
    alpha_table(singlecellexperiment, gene, name = name)
end


function fit_single_cell_logistic_regression(singlecelldata; splitdataby = :Sample, assay = "CPM", fitparameter = "mu")
    cells = singlecelldata.colData[!,:CellID]
    b0 = []
    b1 = []
    #p = Progress(length(cells))
    
    for cell in cells
        cell_index = find_cell_index(cell, singlecelldata)
        condition = find_cell_condition(cell, splitdataby, singlecelldata)
        genes__mus_in_condition = Array{Float64, 2}(singlecelldata.rowData[:, [Symbol(string(condition, "__", fitparameter))]])
         expression_cell_bool = Array{Float64,1}([i > 0 for i in singlecelldata.assays[assay][:, cell_index]])
        # fit each cell to a logistic regression
        model_ = LogisticRegression(penalty = "l2").fit(reshape(genes__mus_in_condition, length(genes__mus_in_condition), 1), expression_cell_bool)
        #next!(p)
        push!(b0, model_.intercept_[1])
        push!(b1, model_.coef_[1])
    end
    
    singlecelldata.colData[!,:B0] = b0
    singlecelldata.colData[!,:B1] = b1
    
    return singlecelldata
    
end

function Shalek2014_weight(cellid, geneid, singlecelldata; conditioncol = :Sample, fitparameter = "mu")
    gene_index = find_gene_index(geneid, singlecelldata)
    cell_index = find_cell_index(cellid, singlecelldata)
    condition = find_cell_condition(cellid, conditioncol, singlecelldata)
    b0 = singlecelldata.colData[cell_index, :B0]
    b1 = singlecelldata.colData[cell_index, :B1]
    u = singlecelldata.rowData[gene_index, Symbol(string(condition, "__", fitparameter))]
    return w = 1 / (1+ e^(-(b0 -b1*u)))
end


function transform_induction_values(singlecellexperiment; untreated_pattern = "UT", 
            assay = "lnCPMplus1", 
            comparedtothissample = "",
            fitparameter = "mu")
        
    cells = singlecellexperiment.colnames 
    
    new_array = similar(singlecellexperiment.assays[assay])
    
    for cell in cells
        cell_index = find_cell_index(cell, singlecellexperiment)
        cell_genotype = singlecellexperiment.colData[cell_index, :Genotype]
        cell_sample = singlecellexperiment.colData[cell_index, :Sample]
        
        # This solves a bug on SergiSingleCellData
        if cell_genotype == "RAD21CTCF"
            mu_untreated_col =  Symbol(string(cell_genotype, "_DKO_", untreated_pattern,  "__", fitparameter))
        else
            mu_untreated_col =  Symbol(string(cell_genotype, "_", untreated_pattern,  "__", fitparameter))
        end
        
        if comparedtothissample == "itself"
            mu_untreated_col = Symbol(string(cell_sample,  "__", fitparameter))
        elseif comparedtothissample != "" && untreated_pattern != ""
            mu_untreated_col = Symbol(string(comparedtothissample, "_", untreated_pattern,  "__", fitparameter))
        elseif comparedtothissample != "" && untreated_pattern == ""
            mu_untreated_col = Symbol(string(comparedtothissample,  "__", fitparameter))
        end
        
        cell_expression = singlecellexperiment.assays[assay][:, cell_index]
        genes_mu_untreated_popu = singlecellexperiment.rowData[!,mu_untreated_col]
        
        induction_values = [cell_expression[i]/(genes_mu_untreated_popu[i]+1) for i in 1:length(cell_expression)]
        
        new_array[:, cell_index] = induction_values
    end
    singlecellexperiment.assays["induction_values"] = new_array
    
    return singlecellexperiment
end


function Shalek2014_weight_modulescore(cellid, geneid, singlecelldata; 
        conditioncol = :Sample,
        fitparameter = "mu")
    
    gene_index = find_gene_index(geneid, singlecelldata)
    cell_index = find_cell_index(cellid, singlecelldata)
    condition = find_cell_condition(cellid, conditioncol, singlecelldata)
    b0 = singlecelldata.colData[cell_index, :B0]
    b1 = singlecelldata.colData[cell_index, :B1]
    u = singlecelldata.rowData[gene_index, Symbol(string(condition,  "__", fitparameter))]
    if singlecelldata.counts[gene_index, cell_index] > 0
        w = 1.0
    else
        w = 1 - (1 / (1+ MathConstants.e^(-(b0 -b1*u))))
    end
    
end


function Shalek2014_module_score(gene_list, singlecellexperiment; modulescore_name = :ModuleScore, 
            untreated_pattern = "UT", 
            assay = "lnCPMplus1",
            comparedtothissample = "",
            fitparameter = "mu")
        
    subset_sce = select_these_genes(gene_list, singlecellexperiment)
    cells = singlecellexperiment.colnames
    
    """
    In order to quantitate the activation of the antiviral, peaked, and sustained inflammatory genes within 
    any individual cell, we developed cell-specific module scores. 
    
    First, for all genes in each module, we transformed single-cell expression values into “induction values” 
    by dividing by the average expression of the gene in un-stimulated cells. 
    """
    subset_sce = transform_induction_values(subset_sce, 
            untreated_pattern = untreated_pattern, 
            assay = assay, comparedtothissample = comparedtothissample,
            fitparameter = "mu")
    """
    Second, to ensure that scores were robust to technical differences between single-cell libraries, 
    we assigned a weight of 1 to all induction scores where the gene was detected, 
    and a weight of  for all non-detected genes, where w(C,G) is defined as in Eq. (1). """
    
    modulescores = Array{Float64, 1}()
    
    for cell in cells
        cell_index = find_cell_index(cell, singlecellexperiment)
        
        induction_values = subset_sce.assays["induction_values"][:, cell_index]
        weights = [Shalek2014_weight_modulescore(cell, gene, subset_sce,  conditioncol = :Sample, fitparameter = fitparameter) for gene in genes(subset_sce)]
        
        weighted_induction_values = [induction_values[i]*weights[i] for i in 1:length(weights)]
        weighted_induction_values = [if isnan(i)||isinf(i) 0.0 else i end for i in weighted_induction_values]
        
        """The “module score” of a cell was then calculated as the weighted mean of “induction values” 
        over all module genes within that individual cell. 
        """
        
        tot_weight = sum(weights)
        
        
        modulescore = sum(weighted_induction_values)/tot_weight
        
        push!(modulescores, modulescore)
    end
    
    
    singlecellexperiment.colData[!,modulescore_name] = modulescores
    
    return singlecellexperiment
    
end