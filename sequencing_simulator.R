# reads simulator for nanopore sequencing
nanopore_simulator = function(genomefile, reads, readlength, err, read_output,align_output, error_output){
    raw_data = readLines(genomefile) # load data
    genome = paste(raw_data[-1], collapse = "") # ignore the description of data in first line and transfer data into one long string
    circlegenome = paste(genome,substr(genome, 1, readlength), sep='') # add last few fake base to imitate the circle genome
    length = nchar(genome) # get the genome length
    read_start = sample(length, reads, replace = T) # select read start position randomly
    perfect_simulate = sapply(read_start, function(x){substr(circlegenome, x, x+readlength-1)}) # the "perfect" reads without error
    perfect_reads = paste(perfect_simulate, collapse = "") # transfer perfect reads into one long string
    error_position = sample(nchar(perfect_reads),floor(nchar(perfect_reads)*err)) # select error position randomly
    error_base = sample(3, floor(nchar(perfect_reads)*err), replace = T) # select error base randomly where we place ACGTACG in order so that the number represents the count from original base to error base.
    replaced_base = sapply(error_position, function(x){substr(perfect_reads, x, x)})
    replaced_number = as.numeric(gsub("A", "1", gsub("C", "2", gsub("G","3", gsub("T", "4", replaced_base)))))
    simulate = perfect_reads # add error on the perfect reads
    replaced_base = c() # record the replaced base
    for(j in 1:length(error_position)){
        substr(simulate, error_position[j], error_position[j])= as.character(error_base[j] + replaced_number[j])
    }
    simulate_merge = gsub("1", "A", gsub("2", "C", gsub("3","G", gsub("4", "T", gsub("5", "A", gsub("6", "C", gsub("7", "G", simulate)))))))# transfer error label into base letter
    error = sapply(error_position, function(x){substr(simulate_merge, x, x)}) # record all error base
    split.start = seq(1, nchar(simulate_merge), by = readlength)
    simulate_result = sapply(split.start, function(x){substr(simulate_merge, x, x+readlength-1)}) # split one long reads with error into the length we want

    # build an error table to summarize the error information
    read_number = error_position%/%readlength +1
    position = error_position%%readlength
    read_number[position[position==0]]=read_number[position[position==0]]-1
    position[position==0]=readlength
    error_table = cbind(read_number, position, replaced_base, error)
    colnames(error_table)=c("error_read_ID", "error_position", "correct_base", "error_base")
    # build an alignment table 
    alignment = cbind(read_start, read_start+readlength-1)
    colnames(alignment) = c("start", "end")
    # save the output file
    write.table(simulate_result, file = read_output, sep="\t", col.names = F)
    write.table(alignment, file = align_output, sep="\t", col.name = T) 
    write.table(error_table, file = error_output, sep="\t", col.name = T) 
    print(paste("results have saved as", read_output, align_output, ", and", error_output, sep = ' '))
}
