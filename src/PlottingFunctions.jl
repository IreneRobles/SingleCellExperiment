#########################
### PLOTTING FUNCTIONS ##
#########################
    
    
function show_me_gene_violin(singlecellexperiment::SingleCellExp, gene; 
        x = "Timepoint", hue = "Genotype", assay = "lnCPMplus1"
    )
    
    #This is to keep my color code
    pal = Dict(
            "WT" => 1, 
            "RAD21" => 2,
            "CTCF" => 4,
            "RAD21CTCF" => 3
        )
    
    sce_show = add_genecounts_to_colData(singlecellexperiment, gene, assay = assay)
    sce_show.colData[:Order] =[ pal[gen] for gen in sce_show.colData[:Genotype]]
    sce_show = sort_cells!(sce_show, cols = :Order)
    
    pd = Pandas.DataFrame(sce_show.colData)
    
    Seaborn.violinplot(data = pd, y = gene, x = x, hue = hue, cut = 0, inner = nothing)
    Seaborn.stripplot(data = pd, y = gene, x = x, hue = hue, dodge = true, jitter = 0.3, alpha = 0.2, color = "gray")
    legend_out_of_plot()
    
    ax = gca()
    handles, labels = ax[:get_legend_handles_labels]()
    n = length(unique(labels))
    legend_help(1:n)
    ylabel(string(gene, " ", assay))
end
    
function show_me_gene_alphas( singlecellexperiment::SingleCellExp, gene; x = "Timepoint", hue = "Genotype", assay = "lnCPMplus1", name = "alpha")
        
    #This is to keep my color code
    pal = Dict(
            "WT" => 1, 
            "RAD21" => 2,
            "CTCF" => 4,
            "RAD21CTCF" => 3
        )
    
    table = alpha_table(singlecellexperiment, gene, name = name)
    table[:Order] = [pal[gen] for gen in table[:Genotype]]
    sort!(table, cols = :Order)
    
    
    
    pd = Pandas.DataFrame(table)
    y = string(name, "_", gene)

    Seaborn.barplot(data = pd, y = y, x = x, hue = hue)
end

function show_me_gene_mus( singlecellexperiment::SingleCellExp, gene; x = "Timepoint", hue = "Genotype", assay = "lnCPMplus1", name = "mu")
  show_me_gene_alphas(singlecellexperiment::SingleCellExp, gene; x = x, hue = hue, assay = assay, name = name)
end
    
function show_me_gene_averages( singlecellexperiment::SingleCellExp, gene; x = "Timepoint", hue = "Genotype", assay = "lnCPMplus1", name = "averagepop")
  show_me_gene_alphas(singlecellexperiment::SingleCellExp, gene; x = x, hue = hue, assay = assay, name = name)
end
    
function show_me_gene_stds( singlecellexperiment::SingleCellExp, gene; x = "Timepoint", hue = "Genotype", assay = "lnCPMplus1", name = "std")
  show_me_gene_alphas(singlecellexperiment::SingleCellExp, gene; x = x, hue = hue, assay = assay, name = name)
end

function show_me_gene(singlecellexperiment::SingleCellExp, gene; assay = "CPM")
    fig = figure(figsize = (20, 10))
    
    subplot(2, 1, 1)
        show_me_gene_violin(singlecellexperiment, gene, assay = assay)
    
    subplot(2, 4, 5)
        show_me_gene_alphas(singlecellexperiment, gene, assay = assay)
        legend_removal()
    
    subplot(2, 4, 6)
        show_me_gene_mus(singlecellexperiment, gene, assay = assay)
        legend_removal()
    
    subplot(2, 4, 7)
        show_me_gene_stds(singlecellexperiment, gene, assay = assay)
        legend_removal()
    
    subplot(2, 4, 8)
        show_me_gene_averages(singlecellexperiment, gene, assay = assay)
        legend_removal()

end

function figure_Shalek2014_antiviralcoremodule_clustermaps(singlecellexperiment::SingleCellExp, genesmodule; assay = "lnCPMplus1", 
        allow_cell_clustering = false,
        allow_gene_clustering = true,
        reordering_of_genes = 1:n_genes(select_these_genes(genesmodule, singlecellexperiment)),
        timepoint = "UT",
        modulename = "",
        figsize = (20, 15), folder = "ClusterMaps_OrderBy_ModuleScore",
        comparedtothissample = "")

        # 14th March 2018
    
    
    singlecellexperiment = Shalek2014_module_score(genesmodule, singlecellexperiment, comparedtothissample = comparedtothissample)
    gene_index_ifnb1 = find_gene_index("Ifnb1", singlecellexperiment)
    singlecellexperiment.colData[:Ifnb1] = singlecellexperiment.assays[assay][gene_index_ifnb1, :]
    
    subset_genes = select_these_genes(genesmodule, singlecellexperiment)
    
    
    
    subset_genes = reorder_genes!(subset_genes, reordering_of_genes[1:end])
    
    # Attach Ifnb1 expression to cell table

    subset_genes = sort_cells!(subset_genes, cols = [:Timepoint, :Genotype, :ModuleScore], rev = false)
    
    if timepoint != "All"
        subs = get_cells_with_this_characteristic(timepoint, :Timepoint, subset_genes)
        else 
        subs = subset_genes
    end
    
    
    subs = sort_cells!(subs, cols = [:Timepoint, :Genotype, :ModuleScore], rev = false)
    cells = [string(i) for i in subs.colData[:CellID]]
    genes = subs.rownames
    pd_df = Pandas.DataFrame(subs.assays[assay], columns = cells, index = genes)
    df_df = DataFrames.DataFrame(pd_df)
    
    
    # Make folder for figures in case it has not been made previously
    if !in(folder, readdir())
        mkdir(folder)
    end
    # Produce figure name   
    figname = string(folder, "/", timepoint, "_clustermap_cellsorderedbymodulescore", modulename ,".png")
    

    CSV.write("data_forPython_assay.csv", df_df)
    CSV.write("data_forPython.csv", subs.colData)
    
    gene_list = py"""
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import csv

# Assay matrix

assay = 'data_forPython_assay.csv'
df1 = pd.read_csv(assay)
df1 = df1.rename(index=str, columns={"x": "GeneSymbol"})
genes = df1.pop('GeneSymbol')
df1 = df1.set_index(genes)


# Data set
url2 = 'data_forPython.csv'
df2 = pd.read_csv(url2)
df2 = df2.set_index(df2.pop("CellID"))
genotypes = df2.pop("Genotype")

lut = dict(zip(genotypes.unique(), ["red","yellow","green", "blue"]))
row_colors = genotypes.map(lut)

cg = sns.clustermap(df1, col_colors=row_colors, row_cluster=$allow_gene_clustering, col_cluster=$allow_cell_clustering, cmap = "plasma", yticklabels = 1, figsize = (50, 50))

plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize = 7)
plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90, fontsize = 7)
ax = plt.gca()
ax.set_xlabel("log(CPM + 1)")


plt.savefig($figname)

if $allow_gene_clustering == True:
    a = cg.dendrogram_row.reordered_ind
"""
    rm("data_forPython_assay.csv")
    rm("data_forPython.csv")    
    # Python uses 0 indexing so we need to reformat to Julia indexis
    if allow_gene_clustering == true
        return gene_list = [i + 1 for i in py"a"]
    end
end


function figure_Shalek2014_modulescoreordered_ifnb1expression(singlecellexperiment::SingleCellExp, genesmodule; assay = "lnCPMplus1", 
        timepoint = "UT",
        figsize = (20, 10), folder = "AntiviralScore_Ifnb1", 
        comparedtothissample = "")
        
        # 14th March 2018

    
    
    singlecellexperiment = Shalek2014_module_score(genesmodule, singlecellexperiment, 
            comparedtothissample = comparedtothissample)
    gene_index_ifnb1 = find_gene_index("Ifnb1", singlecellexperiment)
    singlecellexperiment.colData[:Ifnb1] = singlecellexperiment.assays[assay][gene_index_ifnb1, :]
    
    singlecellexperiment.colData[:ModuleScoreScaled] = scale_between_0_1(singlecellexperiment.colData[:ModuleScore])
    
    
    subset_genes = select_these_genes(genesmodule, singlecellexperiment)

    
    # Attach Ifnb1 expression to cell table

    subset_genes = sort_cells!(subset_genes, cols = [:Timepoint, :Genotype, :ModuleScore], rev = false)
    
    
    subs = get_cells_with_this_characteristic(timepoint, :Timepoint, subset_genes)
    colors = subs.colData[:Genotype]
    colors = [if g == "WT" "blue" elseif g=="RAD21" "orange" elseif g=="CTCF" "red" else "green" end for g in colors]
    
    fig = figure(figsize = figsize) 
    subplot(2, 1, 1)
    title(timepoint)
    bar(0:(nrow(subs.colData)-1), subs.colData[:ModuleScoreScaled], color = colors)
    ylim(0, 1)
    xlim(0, nrow(subs.colData))
    axes = gca()
    axes[:spines]["top"][:set_visible](false) # Hide the top edge of the axis
    axes[:spines]["right"][:set_visible](false)
    xticks([])
    ylabel("ModuleScore")
    subplot(2, 1, 2)
    scatter(0:(nrow(subs.colData)-1), subs.colData[:Ifnb1], color = colors)
    xlim(0, nrow(subs.colData))
    axes = gca()
    axes[:spines]["top"][:set_visible](false) # Hide the top edge of the axis
    axes[:spines]["right"][:set_visible](false)
    ylim(0, 15)
    xticks([])
    ylabel("Ifnb1 expression")
    
    #Make folder for figures in case it has not been made previously
    if !in(folder, readdir())
        mkdir(folder)
    end
    #Produce figure name   
    figname = string(folder, "/", timepoint, "_modulescore_cellsorderedbymodulescore.png")
    savefig(figname)

end
    
#######################
## STATISTICAL TESTS ##
#######################
    
function MultipleTesting.adjust(df::DataFrames.DataFrame, method_adjust_pvalues)
    
    pvalues = collect(Float64, df[:pvalue])
    df[:padj] = adjust(pvalues, method_adjust_pvalues)
    df
    
end
    
function fisher_test_alphas(singlecellexperiment, sample1, sample2)
    # Tests for each pair of genes whether the fraction of expressing cells is different    
    
    ncells = number_cells_per_sample(singlecellexperiment)
    f(x, sam) = [i == sam for i in x]
    # Number of cells in each sample
    nsample1 = ncells[f(ncells[:Sample], sample1), :Number_Cells][1]
    nsample2 = ncells[f(ncells[:Sample], sample2), :Number_Cells][1]
    
    
    colalpha1 = Symbol(string(sample1, "__alpha"))
    colalpha2 = Symbol(string(sample2, "__alpha"))
    
    genes = singlecellexperiment.rowData[:, [:GeneID, colalpha1, colalpha2]]
    
    pvals = []
    confints =  []
    diffs = []
    
    p = Progress(nrow(genes), 0.5)   # minimum update interval: 1 second
    
    for gene in 1:nrow(genes)
        
        x1 = Int(round(genes[gene, colalpha1]*nsample1, 0))
        x2 = Int(round(genes[gene, colalpha2]*nsample2, 0))
        
        n1 = nsample1
        n2 = nsample2
        
        fisher_test = FisherExactTest(x1, x2, n1, n2)
        
        pval = pvalue(fisher_test)
        push!(pvals, pval)
        conf = confint(fisher_test)
        push!(confints, conf)
        diff = genes[gene, colalpha1] - genes[gene, colalpha2]
        push!(diffs, diff)
        next!(p)
    end
    
    genes[:pvalue] = pvals
    genes[:confidence_interval] = confints
    genes[:difference__alpha] = diffs
    
    return genes
    
end
    
export GeneClustermap
    
function GeneClustermap(singlecellexperiment::SingleCellExp, genesmodule, modulescoresymbol::Symbol; assay = "lnCPMplus1", 
        allow_cell_clustering = false,
        allow_gene_clustering = true,
        reordering_of_genes = 1:n_genes(select_these_genes(genesmodule, singlecellexperiment)),
        timepoint = "UT",
        modulename = "",
        figsize = (20, 15), folder = "ClusterMaps_OrderBy_ModuleScore",
        comparedtothissample = "")

        # 14th March 2018
    
    
   # singlecellexperiment = Shalek2014_module_score(genesmodule, singlecellexperiment, comparedtothissample = comparedtothissample)
    
    subset_genes = select_these_genes(genesmodule, singlecellexperiment)
    
    
    
    subset_genes = reorder_genes!(subset_genes, reordering_of_genes[1:end])
    
    # Attach Ifnb1 expression to cell table

    subset_genes = sort_cells!(subset_genes, cols = [:Timepoint, :Genotype, modulescoresymbol], rev = false)
    
    if timepoint != "All"
        subs = get_cells_with_this_characteristic(timepoint, :Timepoint, subset_genes)
        else 
        subs = subset_genes
    end
    
    
    subs = sort_cells!(subs, cols = [:Timepoint, :Genotype, modulescoresymbol], rev = false)
    cells = [string(i) for i in subs.colData[:CellID]]
    genes = subs.rownames
    pd_df = Pandas.DataFrame(subs.assays[assay], columns = cells, index = genes)
    df_df = DataFrames.DataFrame(pd_df)
    
    
    # Make folder for figures in case it has not been made previously
    if !in(folder, readdir())
        mkdir(folder)
    end
    # Produce figure name   
    figname = string(folder, "/", timepoint, "_clustermap_cellsorderedbymodulescore", modulename ,".png")
    

    CSV.write("data_forPython_assay.csv", df_df)
    CSV.write("data_forPython.csv", subs.colData)
    
    gene_list = py"""
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import csv

# Assay matrix

assay = 'data_forPython_assay.csv'
df1 = pd.read_csv(assay)
df1 = df1.rename(index=str, columns={"x": "GeneSymbol"})
genes = df1.pop('GeneSymbol')
df1 = df1.set_index(genes)


# Data set
url2 = 'data_forPython.csv'
df2 = pd.read_csv(url2)
df2 = df2.set_index(df2.pop("CellID"))
genotypes = df2.pop("Genotype")

lut = dict(zip(["CTCF", "RAD21", "RAD21CTCF", "WT"], ["red","yellow","green", "blue"]))
row_colors = genotypes.map(lut)

cg = sns.clustermap(df1, col_colors=row_colors, row_cluster=$allow_gene_clustering, col_cluster=$allow_cell_clustering, cmap = "plasma", yticklabels = 1, figsize = (50, 50))

plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize = 7)
plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90, fontsize = 7)
ax = plt.gca()
ax.set_xlabel("log(CPM + 1)")


plt.savefig($figname)

if $allow_gene_clustering == True:
    a = cg.dendrogram_row.reordered_ind
"""
    rm("data_forPython_assay.csv")
    rm("data_forPython.csv")    
    # Python uses 0 indexing so we need to reformat to Julia indexis
    if allow_gene_clustering == true
        return gene_list = [i + 1 for i in py"a"]
    end
end
    
export GeneClustermap_noorder
    
function GeneClustermap_noorder(singlecellexperiment::SingleCellExp, genesmodule::Array; assay = "lnCPMplus1", 
        allow_cell_clustering = false,
        allow_gene_clustering = true,
        reordering_of_genes = 1:n_genes(select_these_genes(genesmodule, singlecellexperiment)),
        timepoint = "UT",
        modulename = "",
        figsize = (20, 15), folder = "ClusterMaps_noOrderBy_ModuleScore",
        comparedtothissample = "")

        # 14th March 2018
    
    
   # singlecellexperiment = Shalek2014_module_score(genesmodule, singlecellexperiment, comparedtothissample = comparedtothissample)
    
    subset_genes = select_these_genes(genesmodule, singlecellexperiment)
    
    
    
    subset_genes = reorder_genes!(subset_genes, reordering_of_genes[1:end])
    
    # Attach Ifnb1 expression to cell table

    
    if timepoint != "All"
        subs = get_cells_with_this_characteristic(timepoint, :Timepoint, subset_genes)
        else 
        subs = subset_genes
    end
    
    
    subs = sort_cells!(subs, cols = [:Timepoint, :Genotype], rev = false)
    cells = [string(i) for i in subs.colData[:CellID]]
    genes = subs.rownames
    pd_df = Pandas.DataFrame(subs.assays[assay], columns = cells, index = genes)
    df_df = DataFrames.DataFrame(pd_df)
    
    
    # Make folder for figures in case it has not been made previously
    if !in(folder, readdir())
        mkdir(folder)
    end
    # Produce figure name   
    figname = string(folder, "/", timepoint, "_clustermap_cells", modulename ,".png")
    

    CSV.write("data_forPython_assay.csv", df_df)
    CSV.write("data_forPython.csv", subs.colData)
    
    gene_list = py"""
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import csv

# Assay matrix

assay = 'data_forPython_assay.csv'
df1 = pd.read_csv(assay)
df1 = df1.rename(index=str, columns={"x": "GeneSymbol"})
genes = df1.pop('GeneSymbol')
df1 = df1.set_index(genes)


# Data set
url2 = 'data_forPython.csv'
df2 = pd.read_csv(url2)
df2 = df2.set_index(df2.pop("CellID"))
genotypes = df2.pop("Genotypes")

lut = dict(zip(["CTCF", "RAD21", "RAD21CTCF", "WT"], ["red","yellow","green", "blue"]))
row_colors = genotypes.map(lut)

cg = sns.clustermap(df1, col_colors=row_colors, row_cluster=$allow_gene_clustering, col_cluster=$allow_cell_clustering, cmap = "plasma", yticklabels = 1, figsize = (50, 50))

plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize = 7)
plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90, fontsize = 7)
ax = plt.gca()
ax.set_xlabel("log(CPM + 1)")


plt.savefig($figname)

if $allow_gene_clustering == True:
    a = cg.dendrogram_row.reordered_ind
"""
    rm("data_forPython_assay.csv")
    rm("data_forPython.csv")    
    # Python uses 0 indexing so we need to reformat to Julia indexis
    if allow_gene_clustering == true
        return gene_list = [i + 1 for i in py"a"]
    end
end

function genes_per_count(scedata; timepoint = "ALL", genotypes = ["ALL"], genes = "", p = "")
    if genes == ""
        coldata = deepcopy(scedata.colData)
    else
        counts = scedata.colData[:TotalCounts]
         sceexp = select_these_genes(genes,scedata)
        sceexp = SingleCellExpJulia.total_count_and_genes_per_cell(sceexp)
         sceexp.colData[:TotalCounts] = counts
        coldata = deepcopy(sceexp.colData)
    end
    coldata[:Diversity_Genes_per_Counts] = coldata[:TotalGenes]./coldata[:TotalCounts]
    if timepoint != "ALL"
        coldata = coldata[coldata[:Timepoint].== timepoint, :]
    end
    if genotypes != ["ALL"]
        f(x) = [in(ii, genotypes) for ii in x]
         coldata = coldata[f(coldata[:Genotype]), :]
    end
    pd = Pandas.DataFrame(coldata)
    
    if p != ""
    
        Seaborn.boxplot(x = "Timepoint", y = "Diversity_Genes_per_Counts", hue = "Genotype",data = pd, showfliers = false, palette = p)
        Seaborn.stripplot(x = "Timepoint", y = "Diversity_Genes_per_Counts", hue = "Genotype",data = pd, dodge = true, jitter = 0.3, palette = p, s = 3)
    else
        Seaborn.boxplot(x = "Timepoint", y = "Diversity_Genes_per_Counts", hue = "Genotype",data = pd, showfliers = false, )
        Seaborn.stripplot(x = "Timepoint", y = "Diversity_Genes_per_Counts", hue = "Genotype",data = pd, dodge = true, jitter = 0.3, s = 3)
    
        
        
    end
        
        pretty_axes2()
        ylabel("Number of genes detected per cell count")
    
        
end


function plot_expressionranges_singlecell(scedata; timepoint = "ALL", genotypes = ["ALL"], min_cells_expressing_gene = 30, assay = "CPM")
    coldata = deepcopy(scedata)
    if timepoint != "ALL"
        coldata = SingleCellExpJulia.get_cells_with_this_characteristics([timepoint], :Timepoint, coldata)
    end
    if genotypes != ["ALL"]
        coldata = SingleCellExpJulia.get_cells_with_this_characteristics(genotypes, :Genotype, coldata)
    end
    
    coldata = calculate_cv_cellexp(coldata, assay = assay).colData
    
    pd = Pandas.DataFrame(coldata)
    
    Seaborn.boxplot(x = "Timepoint", y = "CV_gene_expression", hue = "Genotype",data = pd, showfliers = false)
    Seaborn.stripplot(x = "Timepoint", y = "CV_gene_expression", hue = "Genotype",data = pd, dodge = true, jitter = 0.3)
    pretty_axes2()
    
end