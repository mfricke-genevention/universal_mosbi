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

process mosbi {
    container "bioandre/mosbi_container" // use docker conatainer
    memory "8 GB"
    publishDir params.output, mode: "copy"

    input:
    path count_file
    path meta_file

    output:
    path "community*"
    path "GraphML_*"
    path "Rplots.pdf"

    """
    Rscript universal.r $count_file $meta_file ./ $params.algrorithm $params.min_size $timepoint
    """
}


workflow {
  mosbi(count_file, meta_file)
}
