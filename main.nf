params.output = "./output/"
params.input = "./example/"
params.meta_file = "columnDescription.txt"
params.count_file = "data_withoutGenes.csv"
params.algrorithm = "base" //"extra" or "all"
params.min_size = 3
params.timepoint = true 
params.protein_mapping = true 

meta_file = Channel.fromPath("${params.input}/${params.meta_file}")
count_file = Channel.fromPath("${params.input}/${params.count_file}")
script_file = Channel.fromPath("universal.r")


process mosbi {
    container "bioandre/mosbi_container" // use docker conatainer
    memory "8 GB"
    publishDir params.output, mode: "copy"

    input:
    path script_file
    path count_file
    path meta_file
    //Rscript $script_file $count_file $meta_file ./

    

    //output:
    //path "community*"
    //path "GraphML_*"
    //path "Rplots.pdf"

    file('current_lol2.txt').withWriter { writer ->
    writer.println "Script file: ${script_file}"
}


    """
    pwd > current_directory.txt
    
    """
}


workflow {
  mosbi(script_file, count_file, meta_file)
}
