library(ape)
library(tools)
library(stringr)

# point to the delimitation file and set parameters of the model

#ancestral
delimitation_file <- "~/Desktop/Berkeley/Phylogenetic_analysis/SVDQ/Freund_population_ancestral.csv"
outgroup_population <- c("Isoetes_bolanderii", "Isoetes_howellii")
num_quartets <- 25000000
num_bootstrap_reps <- 3000
num_quartets_bootstrap <- 100000

#sites
# delimitation_file <- "~/Desktop/Berkeley/Phylogenetic_analysis/SVDQ/Freund_population_sites.csv"
# outgroup_population <- c("Isoetes_bolanderii", "Isoetes_howellii")
# num_quartets <- 50000000
# num_bootstrap_reps <- 1000
# num_quartets_bootstrap <- 570000

#Individuals 
# delimitation_file <- "~/Desktop/Berkeley/Phylogenetic_analysis/SVDQ/Freund_population_individuals.csv"
# outgroup_population <- c("Isoetes_bolanderii", "Isoetes_howellii")
# num_quartets <- 50000000
# num_bootstrap_reps <- 1000
# num_quartets_bootstrap <- 570000
# 

# point to the dataset directory
dir <- "~/Desktop/Berkeley/Phylogenetic_analysis/BPP/Phased_run_files/Phased_data_nex/"

# analysis settings
model <- "lset nst=6 rmatrix=estimate basefreq=empirical rates=gamma ncat=4 shape=estimate pinvar=estimate;"

num_tree_runs <- 1
num_bootstrap_runs <- 1

# create the new directory
new_dir <- file_path_sans_ext(delimitation_file)
dir.create(new_dir, showWarnings = FALSE)

# read the delimitation file
delimitation <- read.csv(delimitation_file)

# list all the files there
if ( any(grepl("nexus", list.files(dir))) ) {
  format <- "nexus"
  phased_alignments <- list.files(dir, pattern = "phased.nexus", full.names = TRUE)
} else {
  format <- "fasta"
  phased_alignments <- list.files(dir, pattern = "phased.fasta", full.names = TRUE)  
}

# create the alignment
alignments <- vector("list", length(phased_alignments))
for(i in 1:length(phased_alignments)) {
  
  cat(i,"\n")
  
  # read the alignment
  if ( format == "nexus" ) {
    
    # read the nexus file
    aln <- try(read.nexus.data(phased_alignments[i]))
    if (class(aln) == "try-error") {
      warning("Could not read alignment ", phased_alignments[i], ". Make sure it is a nexus file")
      next
    }
    
    # drop excluded characters
    aln_lines <- readLines(phased_alignments[i])
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
      
      for(j in 1:length(aln)) {
        aln[[j]] <- aln[[j]][-excludes]
      }
    }
    
    aln <- as.DNAbin(aln)
    
  } else {
    aln <- read.FASTA(phased_alignments[i])  
  }
  
  # get rid of spaces in sequence name
  names(aln) <- sapply(strsplit(names(aln), " "), head, n=1)
  
  # drop the gene name
  aln <- as.matrix(aln)
  old_names <- rownames(aln)
  new_names <- sapply(strsplit(old_names, "-"), head, n=1)
  new_names <- paste0(new_names, ifelse(grepl("_h1", old_names), ".A", ".B"))
  rownames(aln) <- new_names
  
  # record the alignment
  alignments[[i]] <- aln
  
}

# remove nulls
alignments <- alignments[sapply(alignments, is.null) == FALSE]

# concatenate the alignments
bind_aln <- function(...) cbind.DNAbin(..., fill.with.gaps = TRUE)
concatenated_aln <- do.call(bind_aln, alignments)

# relabel the taxa
old_names <- rownames(concatenated_aln)
new_names <- paste0("^", old_names)
rownames(concatenated_aln) <- new_names

# get the groups
groups  <- unique(delimitation$Assignment.to.species)
ngroups <- length(groups)
tips_per_group <- vector("list", ngroups)
names(tips_per_group) <- groups
for(i in 1:ngroups) {

  # get the group
  this_group <- groups[i]
  
  # find all sequences in this group
  these_sequences <- delimitation[delimitation$Assignment.to.species == this_group,]$Specimen
  these_sequences <- as.vector(outer(c("A","B"), these_sequences, function(x, y) paste0("^",y, ".", x) ))
  
  # store
  tips_per_group[[i]] <- these_sequences
    
}

# get all the samples
all_samples <- unlist(tips_per_group)

# write the alignment
write.nexus.data(concatenated_aln,
                 file = paste0(new_dir, "/seq.nex"),
                 interleaved = FALSE)

# write the SVDQ nexus file
svdq_file <- paste0(new_dir, "/svdq.nex")
cat("#NEXUS\n\n", file = svdq_file, sep = "")

cat("begin paup;
  cd *;
  set hashcomment=y;
end;\n\n", file = svdq_file, sep = "", append = TRUE)

cat("execute seq.nex;\n\n", file = svdq_file, sep = "", append = TRUE)

cat("begin sets;\n", file = svdq_file, sep = "", append = TRUE)
cat("taxpartition species =\n", file = svdq_file, sep = "", append = TRUE)
for(i in 1:length(tips_per_group)) {
  this_group    <- groups[i]
  these_samples <- tips_per_group[[i]]
  these_samples   <- paste0(these_samples, collapse = " ")
  cat("\t", this_group, " : ", these_samples, file = svdq_file, sep = "", append = TRUE)
  if ( i == length(tips_per_group) ) {
    cat(";\n\n", file = svdq_file, sep = "", append = TRUE)
  } else {
    cat(",\n", file = svdq_file, sep = "", append = TRUE)
  }
}

cat("charpartition lociset =\n", file = svdq_file, sep = "", append = TRUE)
current_position <- 0
for(i in 1:length(alignments)) {
  
  # get the locus
  this_aln <- alignments[[i]]
  
  # get the positions
  start <- 1 + current_position
  end   <- ncol(this_aln) + current_position
  
  # add the info
  cat(i, ": ", start, "-", end, file = svdq_file, sep = "", append = TRUE)
  if ( i == length(alignments) ) {
    cat(";\n", file = svdq_file, sep = "", append = TRUE)
  } else {
    cat(",\n", file = svdq_file, sep = "", append = TRUE)
  }
  
  # increment the positions
  current_position <- end
  
}
cat("end;\n\n", file = svdq_file, sep = "", append = TRUE)

# specify the outgroup
outgroup <- unlist(tips_per_group[outgroup_population])
cat("outgroup ", paste0(outgroup, collapse=" "), ";\n", file = svdq_file, sep = "", append = TRUE)

# specify the model
cat(model, "\n\n", file = svdq_file, sep = "", append = TRUE)

# replicates
for(i in 1:num_tree_runs) {
  cat("svdq taxpartition=species evalQuartets=all showScores=no replace=yes;\n", file = svdq_file, sep = "", append = TRUE)
  cat("savetrees brlens file=MyTree", i, ".tre replace;\n\n", file = svdq_file, sep = "", append = TRUE)
}

for(i in 1:num_bootstrap_runs) {
  cat(paste0("svdq taxpartition=species nquartets=", num_quartets_bootstrap, " showScores=no bootstrap nreps=", num_bootstrap_reps, " treeFile=mybootstraptrees",i,".tre replace=yes;\n"), file = svdq_file, sep = "", append = TRUE)
  cat("savetrees brlens file=MyBSTree", i, ".tre replace;\n\n", file = svdq_file, sep = "", append = TRUE)
}

cat("quit;\n", file = svdq_file, sep = "", append = TRUE)

