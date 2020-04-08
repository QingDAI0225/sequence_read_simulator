# reads simulator for nanopore sequencing
nanopore_simulator = function(genomefile, coverage, readlength, err, read_output,align_output, error_output){
    raw_data = readLines(genomefile) # load data
    genome = paste(raw_data[-1], collapse = "") # ignore the description of data in first line and transfer data into one long string
    length = nchar(genome) # get the genome length
    read_start = sample(length, coverage*length/readlength, replace = T) # select read start position randomly
    perfect_simulate = c() # the "perfect" reads without error
    for (i in 1:length(read_start)){
        perfect_simulate[i] = substr(genome, read_start[i] , read_start[i] + readlength -1)
    }
    perfect_reads = paste(perfect_simulate, collapse = "") # transfer perfect reads into one long string
    error_position = sample(nchar(perfect_reads),floor(nchar(perfect_reads)*err)) # select error position randomly
    error_base = sample(4, floor(nchar(perfect_reads)*err), replace = T) # select error base randomly where 1 represents A, 2 represents C, 3 represents G and 4 represents T.
    simulate = perfect_reads # add error on the perfect reads
    replaced_base = c() # record the replaced base
    for(j in 1:length(error_position)){
        replaced_base[j] = substr(simulate, error_position[j], error_position[j])
        substr(simulate, error_position[j], error_position[j])=as.character(error_base[j])
    }
    simulate_merge = gsub("1", "A", gsub("2", "C", gsub("3","G", gsub("4", "T", simulate))))# transfer error label into base letter
    error = gsub("1", "A", gsub("2", "C", gsub("3","G", gsub("4", "T", error_base)))) # record all error base
    simulate_result = c() # split one long reads with error into the length we want
    for (k in 1:(nchar(simulate_merge)/readlength)){
        simulate_result[k] = substr(simulate_merge, (k-1)*readlength+1, k*readlength)
    } 
    # build an error table to summarize the error information
    read_number = error_position%/%readlength
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
    print(paste(c("results have saved as", read_output, "and", align_output)))
}
