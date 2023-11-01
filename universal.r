####################################################################################################################
# Applying the MoSBI package on proteomics data (https://bioconductor.org/packages/release/bioc/html/mosbi.html)
# Author: André Michael Bembennek
# Contact: andre.bembennek@studium.uni-hamburg.de
# Date: 06.06.2023
# new features: GO enrichment via gprofiler2
####################################################################################################################


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")


if (!require("mosbi", quietly = TRUE))
    BiocManager::install("mosbi")
#################################################### DISCLAIMER ####################################################
# All passages that may need to be adjusted are marked with a (*).
# Please adjust them according to your data.
# If you have any questions, please contact me.

#################################################### USAGE ####################################################
# Rscript universal_mosbi4genevention2.r savepath path2data path2meta algorithms(optional) min_size(optional) protein_mapping(optional) timepoint(optional) 

# savepath: path to save the results

# path2data: path to the data (.csv, sep=",")

# path2meta: path to the meta data (.csv, sep="\t")

# algorithms: "base", "extra" or "all"
# "base"     --> Fabia, isa2, BCplaid, QUBIC
# "extra"   --> akmbi, bimax, cc, quest, spectral, xmotifs
# "all"     --> all algorithms

# min_size: minimum size of a community (optional, default=3)

# timepoint: timepoint to be analyzed (optional, default=FALSE)

# path2mapping: path to the mapping between Protein.ID and e.g. genes (optional, default=FALSE)





#install.packages("gprofiler2")

############# load libraries #############
library(mosbi)
library(igraph)
library(gprofiler2)


############# read in arguments #############
args <- commandArgs(trailingOnly = TRUE)
#print("teeeeeeeesttttt--------------testtesttesttesttesttesttesttesttest")
if (length(args) < 3) {

    stop("Usage: Rscript universal.r savepath path2data path2meta algorithms(optional: 'base', 'extra' or 'all') min_size(optional) protein_mapping(optional) timepoint(optional)", call.=FALSE)

}else{
    savepath <- args[1]
    data <- args[2]
    meta <- args[3]
    algos <- "base"
    min_size <- 3
    timepoint <- FALSE
    path2mapping <- FALSE
}
if (length(args) == 4){
    algos <- args[4]
}
if (length(args) == 5){
    algos <- args[4]
    min_size <- args[5]  
}
if (length(args) == 6){
    algos <- args[4]
    min_size <- args[5]
    path2mapping <- args[6]
}
if (length(args) == 7){
    algos <- args[4]
    min_size <- args[5]
    path2mapping <- args[6]
    timepoint <- args[7]
    
}

############# read in data and meta data (*) #############

#(*) May need to be adjusted. The test data have the gene names in the 1st column and the protein IDs in the 2nd column.
data <- read.csv(data, header = TRUE, sep = ",", row.names = 1,check.names = FALSE)
datamatrix <- as.matrix(data)
# Saving the connection between gene and protein

if (path2mapping != FALSE){
    protein_mapping <- read.csv(path2mapping, header = TRUE, sep = ",")
}

#(*) colnames of meta: Column, Batch, Animal, Timepoint, Group
#(*) tab separated
meta <- read.csv(meta, header = TRUE, sep = "\t")                                   

# (optional) filtering for timepoint
if (timepoint != FALSE){
    time <- meta[meta$Timepoint == timepoint,]
    time <- time$Column
    datamatrix <- datamatrix[,time]
}



########################## run bicluster algorithms #######################
#' apply bicluster algorithms on datamatrix
#' 
#' @param algos "base", "extra" or "all"
#' "base" --> Fabia, isa2, BCplaid, QUBIC
#' "extra" --> akmbi, bimax, cc, quest, spectral, xmotifs
#' "all" --> all algorithms
#' @return List of bicluster objects
algo <- function(algos, own_data_matrix){
if (algos == "base"){
# Fabia
fb <- mosbi::run_fabia(own_data_matrix)

# isa2
BCisa <- mosbi::run_isa(own_data_matrix)

# BCplaid
BCplaid <- mosbi::run_plaid(own_data_matrix)

# QUBIC
BCqubic <- mosbi::run_qubic(own_data_matrix)

## merge base
all_bics <- c(fb, BCisa, BCplaid, BCqubic)
}

if (algos == "extra"){
## additional
akmbi <- mosbi::run_akmbiclust(own_data_matrix)
bimax <- mosbi::run_bimax(own_data_matrix)
cc <- mosbi::run_cc(own_data_matrix)
quest <- mosbi::run_quest(own_data_matrix)
spectral <- mosbi::run_spectral(own_data_matrix)
xmotifs <- mosbi::run_xmotifs(own_data_matrix)

## merge additional
all_bics <- c(akmbi,bimax,cc,quest,spectral,xmotifs)
}

if (algos == "all"){
    # Fabia
    fb <- mosbi::run_fabia(own_data_matrix)

    # isa2
    BCisa <- mosbi::run_isa(own_data_matrix)

    # BCplaid
    BCplaid <- mosbi::run_plaid(own_data_matrix)

    # QUBIC
    BCqubic <- mosbi::run_qubic(own_data_matrix)

    ## additional
    akmbi <- mosbi::run_akmbiclust(own_data_matrix)
    bimax <- mosbi::run_bimax(own_data_matrix)
    cc <- mosbi::run_cc(own_data_matrix)
    quest <- mosbi::run_quest(own_data_matrix)
    spectral <- mosbi::run_spectral(own_data_matrix)
    # unibic <- mosbi::run_unibic(own_data_matrix)
    xmotifs <- mosbi::run_xmotifs(own_data_matrix)

    all_bics <- c(fb, BCisa, BCplaid, BCqubic, akmbi,bimax,cc,quest,spectral,xmotifs)
}
return(all_bics)
}
algos <- "all"
all_bics <- algo(algos, datamatrix)



########################## compute bicluster network #######################
bic_net <- mosbi::bicluster_network(all_bics, # List of biclusters
    datamatrix, # Data matrix
    n_randomizations = 5,
    # Number of randomizations for the
    # error model
    MARGIN = "both",
    # Use datapoints for metric evaluation
    metric = 4, # Fowlkes–Mallows index
    # For information about the metrics,
    # visit the "Similarity metrics
    # evaluation" vignette
    n_steps = 1000,
    # At how many steps should
    # the cut-of is evaluated
    plot_edge_dist = TRUE
    # Plot the evaluation of cut-off estimation
)


min_size <- 3
coms <- mosbi::get_louvain_communities(bic_net,
    min_size = min_size,
    bics = all_bics
)

########################## Functions to convert bicluster objects to igraph objects(*) ##########################

#' convert bicluster object to igraph object
biclust2igraph <- function(biclustObj){
    # get igraph object from adjacency matrix
    GraphML <- graph.adjacency(biclustObj@adjacency_matrix, mode = "undirected", weighted = TRUE)

    # delete edges with weights < threshold
    # from mosbi: "Estimated threshold for the bicluster similarity adjacency matrix. All values lower than that in the matrix should be discarded."
    GraphML <- delete.edges(GraphML, which(E(GraphML)$weight < biclustObj@threshold))
    
    # delete edges with weights == 1
    GraphML <- delete.edges(GraphML, which(E(GraphML)$weight ==1))

    return(GraphML)
}

#' get Batch-Info for each Sample via meta data
loadBatch <- function(column_names){
    batches <- c()
    for(col in column_names){
        index <- which(meta$Column == col)
        batch <- meta[index,"Batch"]                                            #(*)
        batches <- c(batches,batch)
    }
    batches <- paste0(batches, collapse = ",")
    return(batches)
}

#' get Time-Info for each Sample via meta data
loadTime <- function(column_names){
    timepoints <- c()
    for(col in column_names){
        index <- which(meta$Column == col)
        time <- meta[index,"Timepoint"]                                         #(*)
        timepoints <- c(timepoints,time)
    }
    timepoints <- paste0(timepoints, collapse = ",")
    return(timepoints)
}

#' get Animal-Info for each Sample via meta data
loadAnimal <- function(column_names){
    animals <- c()
    for(col in column_names){
        index <- which(meta$Column == col)
        animal <- meta[index,"Animal"]                                          #(*)
        animals <- c(animals,animal)
    }
    animals <- paste0(animals, collapse = ",")
    return(animals)
}

#' get Group-Info for each Sample via meta data
loadGroup <- function(column_names){
    groups <- c()
    for(col in column_names){
        index <- which(meta$Column == col)
        group <- meta[index,"Group"]                                            #(*)
        groups <- c(groups,group)
    }
    groups <- paste0(groups, collapse = ",")
    return(groups)
}

#' get Genename for each Protein via protein_mapping
loadGenenames <- function(row_names){
    genenames <- c()
    for(row in row_names){
        index <- which(protein_mapping$Protein.ID == row)
        genename <- protein_mapping[index,"Mapping"]                        #(*)      
        genenames <- c(genenames,genename)
    }
    genenames <- paste0(genenames, collapse = ",")
    return(genenames)
}



find_com <- function(bicluster){
    for(i in 1:length(coms)){
        comy <- coms[[i]]
        if(bicluster %in% rownames(comy@adjacency_matrix)){
            return(i)
        }
    }
    return(NA) #if bicluster is not in any community
}
















#' load all meta information into igraph object
loadColor <- function(igraphObj){
    
    # for each node in igraph object
    for(i in 1:vcount(igraphObj)){

        # get node name (e.g. "Bicluster1"")
        nodename <- V(igraphObj)[i]$name

        # extract number from node name (e.g. "1")
        my_match <- regexpr("\\d+", nodename, perl=TRUE)
        number <- as.numeric(substr(nodename, my_match, my_match + attr(my_match, "match.length") - 1))
        
        # get bicluster object from list of biclusters (number as index)
        bicluster <- all_bics[[number]]
        
        # get meta information from bicluster object
        row_name <- paste(bicluster@rowname, collapse = ",")
        V(igraphObj)[i]$rowname <- row_name

        col_name <- paste(bicluster@colname, collapse = ",")
        V(igraphObj)[i]$colname <- col_name

        rows <- paste(bicluster@row, collapse = ",")
        V(igraphObj)[i]$row <- rows

        columns <- paste(bicluster@column, collapse = ",")
        V(igraphObj)[i]$column <- columns
        
        algo <- bicluster@algorithm
        V(igraphObj)[i]$algorithm <- algo


        # get meta information from meta data
        V(igraphObj)[i]$batch <- loadBatch(bicluster@colname)
        V(igraphObj)[i]$timepoint <- loadTime(bicluster@colname)
        V(igraphObj)[i]$animal <- loadAnimal(bicluster@colname)
        V(igraphObj)[i]$group <- loadGroup(bicluster@colname)
        if (path2mapping != FALSE){
            V(igraphObj)[i]$mapping <- loadGenenames(bicluster@rowname) 
        }
        V(igraphObj)[i]$community <- find_com(nodename)
        
    }
    return(igraphObj)
}



########################## save bicluster network as graphml ##########################

# convert  bicluster network to igraph object
GraphML_bic_net <- biclust2igraph(bic_net)

# load all meta information into igraph object
GraphML_bic_net <- loadColor(GraphML_bic_net)

# save igraph object as graphml
path <- paste0(savepath,"/bicluster_network2.graphml")
write_graph(GraphML_bic_net, path, format = "graphml")



########################## get and save communities #######################

# convert all communities to igraph objects
# load all meta information into igraph object
# save igraph objects as graphml

connect_strings <- function(lst) {
  for (i in 1:length(lst)) {
    if (length(lst[[i]]) > 1) {
      lst[[i]] <- paste(lst[[i]], collapse = ", ")
    }
  }
  return(lst)
}


find_query<-function(query,input_list){
    number <- str_extract(query, "\\d+")
    number <- as.numeric(number)
    return(input_list[[number]])
}


for(i in 1:length(coms)){

    ordner <- paste0(savepath,"/community",i)
    dir.create(ordner)
    com <- biclust2igraph(coms[[i]])
    com <- loadColor(com)
    com_path <- paste0(ordner,"/community",i,".graphml")
    write_graph(com, file = com_path, format = "graphml")


    if(path2mapping != FALSE){
    gene_list <- V(com)$mapping
    
    vector <- unlist(strsplit(gene_list, ",|;"))
    vector <- gsub("[()]", "", vector)
    vector <- gsub("N/A", "", vector)
    
        
    
    results <- gost(vector, organism = "rnorvegicus")
    export_df <- results$result

    export_df$parents <- unlist(connect_strings(export_df$parents))
   

    csv_path <- paste0(ordner,"/enrichment_table.csv")
    write.csv(export_df, file = csv_path, row.names = FALSE)
    }
 }


