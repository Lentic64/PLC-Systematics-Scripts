library(ape)
library(stringr)

# point to the dataset directory
dir <- "Phased data/"

# specify the new directory
new_dir <- "Test"
dir.create(new_dir)

# analysis settings
model <- "lset nst=6 rmatrix=estimate basefreq=empirical rates=gamma ncat=4 shape=estimate pinvar=estimate;"
bootstrap_reps <- 2000
n_quartets <- 300000
n_quartest_bootstrap <- 100000

# list all the files there
phased_alignments <- list.files(dir, pattern = "phased.fasta", full.names = TRUE)

# randomly sample loci
num_loci_to_sample <- 150
if ( num_loci_to_sample < length(phased_alignments) ) {
  phased_alignments <- phased_alignments[sort(sample.int(length(phased_alignments), size=num_loci_to_sample))]
}

cat(phased_alignments, sep="\n", file=paste0(new_dir, "/sampled_alignments.txt"))

# in case you want to include only some populations
# populations_to_include <- c("DT13378",
#                             "FF254",
#                             "FF280",
#                             "FF283",
#                             "FF259",
#                             "FF261",
#                             "FF285",
#                             "FF286")

# All populations
populations_to_include <- c("CR5205",
                            "CR5256",
                            "DT13378",
                            "F088",
                            "FF169",
                            "FF254",
                            "FF257",
                            "FF258",
                            "FF259",
                            "FF261",
                            "FF262",
                            "FF263",
                            "FF266",
                            "FF279",
                            "FF280",
                            "FF282",
                            "FF283",
                            "FF284",
                            "FF285",
                            "FF286",
                            "FF290",
                            "KW380")


# relabel prefixes
prefix_map <- c("CR5205"  = "CR5205_CoR8",
                "CR5256"  = "CR5256_SN1",
                "DT13378" = "DT13378_Min",
                "F088"    = "FF088_Out",
                "FF169"   = "FF169_SN2",
                "FF254"   = "FF254_HL",
                "FF257"   = "FF257_SN3",
                "FF258"   = "FF258_SN4",
                "FF259"   = "FF259_SN5",
                "FF261"   = "FF261_SN6",
                "FF262"   = "FF262_SC1",
                "FF263"   = "FF263_CoR1",
                "FF266"   = "FF266_CoR2",
                "FF279"   = "FF279_SC2",
                "FF280"   = "FF280_CoR3",
                "FF282"   = "FF282_CoR4",
                "FF283"   = "FF283_CoR5",
                "FF284"   = "FF284_CoR6",
                "FF285"   = "FF285_SC3",
                "FF286"   = "FF286_SC4",
                "FF290"   = "FF290_CoR7",
                "KW380"   = "KW380_SN7")

# create the outgroup
outgroup_population <- "Out"

# create the alignment
alignments <- vector("list", length(phased_alignments))
for(i in 1:length(phased_alignments)) {
  
  # read the alignment
  aln <- read.FASTA(phased_alignments[i])
  
  # get rid of spaces in sequence name
  names(aln) <- sapply(strsplit(names(aln), " "), head, n=1)
  
  # relabel things
  for(j in 1:length(prefix_map)) {
    names(aln) <- gsub(names(prefix_map[j]), prefix_map[j], names(aln))
  }
  
  # remove populations, if necessary
  if ( exists("populations_to_include") ) {
    
    # check each sequence whether it is in the populations to be included
    includes <- do.call(rbind, lapply(populations_to_include, function(x) {
      grepl(x, names(aln))
    }))
    
    # find the name to keep
    names_to_include <- names(aln)[colSums(includes) > 0]
    
    # throw out everything else
    aln <- aln[names_to_include]
    
  }
  
  # drop the gene name
  aln <- as.matrix(aln)
  old_names <- rownames(aln)
  new_names <- sapply(strsplit(old_names, "-"), head, n=1)
  new_names <- paste0(new_names, ifelse(grepl("_h1", old_names), ".A", ".B"))
  rownames(aln) <- new_names
  
  # TODO: drop exclude block
  
  # record the alignment
  alignments[[i]] <- aln
  
}

# concatenate the alignments
bind_aln <- function(...) cbind.DNAbin(..., fill.with.gaps = TRUE)
concatenated_aln <- do.call(bind_aln, alignments)

# relabel the taxa
old_names <- rownames(concatenated_aln)
new_names <- paste0("^", old_names)
rownames(concatenated_aln) <- new_names

# randomly sample some species
num_samples_to_sample <- 4
populations <- prefix_map[populations_to_include]
samples_to_include <- vector("list", length(populations))
unique_samples <- unique(sapply(strsplit(new_names, ".", fixed=TRUE), head, n = 1))

for(i in 1:length(populations)) {
  
  # get this population
  this_population <- populations[i]
  
  # get all the samples in this population
  individuals_to_include <- grep(this_population, unique_samples, value = TRUE)
  if ( length(individuals_to_include) > num_samples_to_sample ) {
    individuals_to_include <- sample(individuals_to_include, size = num_samples_to_sample) 
  }
  
  # get the phased labels
  individuals_to_include <- c(paste0(individuals_to_include, ".A"), paste0(individuals_to_include, ".B"))
  
  # store the names of the individuals
  samples_to_include[[i]] <- individuals_to_include
  
}

# check the outgroup
renamed_populations <- sapply(strsplit(populations, "_"), tail, n = 1)
if (outgroup_population %in% renamed_populations == FALSE) {
  stop("Error defining outgroup. Please make sure your defined outgroup is in the analysis.")
}

# drop samples from the alignment
all_samples      <- unlist(samples_to_include)
concatenated_aln <- concatenated_aln[rownames(concatenated_aln) %in% all_samples,]

# get the names for each sample
sample_names <- unique(sapply(strsplit(all_samples, "\\."), head, n = 1))
sample_names <- gsub("^", "", sample_names, fixed = TRUE)

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
for(i in 1:length(sample_names)) {
  this_population <- sample_names[i]
  these_samples   <- grep(this_population, all_samples, value = TRUE)
  these_samples   <- paste0(these_samples, collapse = " ")
  cat("\t", this_population, " : ", these_samples, file = svdq_file, sep = "", append = TRUE)
  if ( i == length(sample_names) ) {
    cat(";\n\n", file = svdq_file, sep = "", append = TRUE)
  } else {
    cat(",\n", file = svdq_file, sep = "", append = TRUE)
  }
}

cat("charpartition lociset =\n", file = svdq_file, sep = "", append = TRUE)
current_position <- 0
for(i in 1:length(phased_alignments)) {
  
  # get the locus
  this_aln <- alignments[[i]]
  
  # get the positions
  start <- 1 + current_position
  end   <- ncol(this_aln) + current_position

  # add the info
  cat(i, ": ", start, "-", end, file = svdq_file, sep = "", append = TRUE)
  if ( i == length(phased_alignments) ) {
    cat(";\n", file = svdq_file, sep = "", append = TRUE)
  } else {
    cat(",\n", file = svdq_file, sep = "", append = TRUE)
  }
  
  # increment the positions
  current_position <- end
      
}
cat("end;\n\n", file = svdq_file, sep = "", append = TRUE)

# specify the outgroup
outgroup <- samples_to_include[renamed_populations == outgroup_population][[1]]
cat("outgroup ", paste0(outgroup, collapse=" "), ";\n", file = svdq_file, sep = "", append = TRUE)

# specify the model
cat(model, "\n\n", file = svdq_file, sep = "", append = TRUE)

# replicates
num_tree_reps <- 4
for(i in 1:num_tree_reps) {
  cat("svdq taxpartition=species nquartets=", n_quartets," showScores=no replace=yes;\n", file = svdq_file, sep = "", append = TRUE)
  cat("savetrees brlens file=MyTree", i, ".tre replace;\n\n", file = svdq_file, sep = "", append = TRUE)
}

num_boot_reps <- 2
for(i in 1:num_boot_reps) {
  cat(paste0("svdq taxpartition=species nquartets=", n_quartest_bootstrap, " showScores=no bootstrap nreps=", bootstrap_reps, " treeFile=mybootstraptrees",i,".tre replace=yes;\n"), file = svdq_file, sep = "", append = TRUE)
  cat("savetrees brlens file=MyBSTree", i, ".tre replace;\n\n", file = svdq_file, sep = "", append = TRUE)
}

cat("quit;\n", file = svdq_file, sep = "", append = TRUE)








