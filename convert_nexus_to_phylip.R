library(ape)
library(tools)

# the directories
setwd("~/Desktop/Berkeley/Phylogenetic_analysis/BPP/Phased_run_files")
source_dir <- "Phased_data_excluded"
out_dir    <- "Phased_data_BPP"

# create the output folder
dir.create(out_dir, recursive = TRUE)

# list all nexus files
nexus_files <- list.files(source_dir, full.names = TRUE)
nexus_files <- nexus_files[file_ext(nexus_files) %in% c("nexus", "nex")]

# read all the nexus files
alns <- vector("list", length(nexus_files))
for(i in 1:length(alns)) {

    # get the file
    this_file <- nexus_files[i]

    cat("Processing file ", this_file, "\n")

    # read the alignment
    this_aln <- try(read.nexus.data(this_file))

    if ( class(this_aln) != "try-error" ) {

        # get the new filename
        new_file <- gsub(source_dir, out_dir, this_file)
        new_file <- gsub(".nexus", ".fasta", new_file)
        new_file <- gsub(".nex", ".fasta", new_file)

        # drop excluded characters
        aln_lines <- readLines(this_file)
        exclude_line <- grep("EXSET", aln_lines, value = TRUE)
        exclude_line <- strsplit(exclude_line, "=")[[1]][2]
        exclude_line <- strsplit(exclude_line, " ")[[1]]
        exclude_line <- gsub(";", "", exclude_line)
        exclude_line <- exclude_line[exclude_line != ""]
        if ( length(exclude_line) > 0 ) {
            excludes <- unlist(lapply(exclude_line, function(x) {
                tmp <- as.numeric(strsplit(x,"-")[[1]])
                tmp <- seq(head(tmp, 1), tail(tmp, 1))
                return(tmp)
            }))
            
            for(j in 1:length(this_aln)) {
                this_aln[[j]] <- this_aln[[j]][-excludes]
            }
        }
        
        # reformat
        this_aln <- as.DNAbin(this_aln)

        # write the alignment
        write.dna(this_aln, file = new_file, format = "fasta", nbcol = -1, colsep="")

    } else {
        cat("Problem with aln ", this_file, "\n")
    }

}

