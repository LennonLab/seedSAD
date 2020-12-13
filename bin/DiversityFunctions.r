################################################################################
#                                                                              #
# Diversity Functions Source Code: Functions for microbial community analyses  #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Mario Muscarella                                                 #
#                                                                              #
# Created: 2012/08/22                                                          #
#                                                                              #
# Last update: 2014/07/09                                                      #
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code provides numerous functions to be used in the analysis of   #
#        16S rRNA sequence data post mothur anlaysis                           #
#        The ideas in the code are the products of a project between           #
#        Mario Muscarella and Tracy Teal                                       #
#                                                                              #
# Issues: Slow performance reading in OTU tables                               #
#                                                                              #
# Recent Changes:                                                              #
#         1. Added Coverage Function                                           #
#         2. Minor Updates                                                     #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1. Design functions to work with shared files in memory              #
#         2. Add warnings                                                      #
#                                                                              #
################################################################################

# Needed packages
require("vegan")||install.packages("vegan");require("vegan")

# Read OTU File
read.otu <- function(shared = " ",
                     cutoff = "0.03"){
  cols <- count.fields(shared)
  classes <- c(rep("character", 2), rep("numeric", (cols[1] - 2)))
  matrix <- read.delim(shared, header=T, as.is=T, colClasses=classes, comment.char="")
  matrix.cutoff <- subset(matrix, matrix$label == cutoff)
  matrix.out <- as.matrix(matrix.cutoff[1:dim(matrix.cutoff)[1],
    4:(mean(cols))])
  row.names(matrix.out) <- matrix.cutoff$Group
  return(matrix.out)
  }

# Count All Groups
count.groups <- function(otu.matrix = " "){
  counts <- rowSums(otu.matrix)
  return(counts)
  }

# Subsampling wrapper
sub.sample <- function(otu.matrix  = " ",
                       sample.size = "min(count.groups(test))"){
  counts <- count.groups(otu.matrix)
  if (sample.size == " "){
    sample.size = counts}
  statement <- counts > sample.size  # Add warning message
  otu.matrix <- subset(otu.matrix, rowSums(otu.matrix)>sample.size)
  x <- rrarefy(otu.matrix, sample.size)
  return(x)
  }

# Example function for calculating species richness
richness.iter <- function(input  = " ",
                          cutoff = " ",
                          size   = " ",
                          iters  = " ",
                          shared = "TRUE"){
  if(shared == TRUE){
    otu.matrix <- read.otu(input, cutoff)
    }else{
      otu.matrix <- input
      }
  counts <- count.groups(otu.matrix)
  statement <- counts > size  # Add warning message
  if(iters > 1){
    otu.matrix <- subset(otu.matrix, rowSums(otu.matrix)>size)
    rich.matrix <- matrix(NA, dim(otu.matrix)[1], iters)
    rownames(rich.matrix) <- rownames(otu.matrix)
    for (i in 1:iters){
      temp.matrix <- sub.sample(otu.matrix, size)
      rich.matrix[,i] <- rowSums((temp.matrix>0)*1)
      }
  }else{
    rich.matrix <- rowSums((input>0)*1)
    }
  return(rich.matrix)
  }

# Example function for calculating species evenness
# (Pielou's Evenness = J' = H'/H'max = H'/lnS, H' = Shannon)
evenness.iter <- function(input  = " ",
                          cutoff = " ",
                          size   = " ",
                          iters  = " ",
                          shared = "TRUE"){
  if(shared == TRUE){
    otu.matrix <- read.otu(input, cutoff)
    }else{
      otu.matrix <- input
      }
  counts <- count.groups(otu.matrix)
  statement <- counts > size  # Add warning message
  if(iters > 1){
    otu.matrix <- subset(otu.matrix, rowSums(otu.matrix)>size)
    even.matrix <- matrix(NA, dim(otu.matrix)[1], iters)
    rownames(even.matrix) <- rownames(otu.matrix)
    for (i in 1:iters){
      temp.matrix <- sub.sample(otu.matrix, size)
      even.matrix[,i] <- diversity(temp.matrix, "shannon")/
        log(rowSums((temp.matrix>0)*1))
      }
   }else{
     even.matrix <- diversity(input, "shannon")/log(rowSums((input>0)*1))
      }
  return(even.matrix)
  }

# Example function for calculating diversity
diversity.iter <- function(input  = " ",
                           index  = "shannon",
                           cutoff = " ",
                           size   = " ",
                           iters  = " ",
                           shared = "TRUE"){
  if(shared == TRUE){
    otu.matrix <- read.otu(input, cutoff)
    }else{
      otu.matrix <- input
      }
  counts <- count.groups(otu.matrix)
  statement <- counts > size  # Add warning message
  if(iters > 1){
    otu.matrix <- subset(otu.matrix, rowSums(otu.matrix)>size)
    div.matrix <- matrix(NA, dim(otu.matrix)[1], iters)
    rownames(div.matrix) <- rownames(otu.matrix)
    for (i in 1:iters){
      temp.matrix <- sub.sample(otu.matrix, size)
      div.matrix[,i] <- diversity(temp.matrix, index)
      }
  }else{
    div.matrix <- diversity(input, index)
    }
  return(div.matrix)
  }

coverage <- function(input= " ", cutoff = " ", size = " ", shared = "TRUE"){
  if(shared == TRUE){
    otu.matrix <- read.otu(input, cutoff)
    } else {
      otu.matrix <- input
      }
  counts <- count.groups(otu.matrix)
  statement <- counts > size  # Add warning message
  otu.matrix <- subset(otu.matrix, rowSums(otu.matrix)>size)
  cov <- matrix(NA, dim(otu.matrix)[1], 1)
  rownames(cov) <- rownames(otu.matrix)
  colnames(cov) <- "Coverage"
  temp.matrix <- sub.sample(otu.matrix, size)
  for (i in 1:dim(temp.matrix)[1]){
    cov[i,] <- 1 - ((length(which(temp.matrix[1,] == 1))) /
        (length(which(temp.matrix[i,] > 0))))
    }
  return(cov)
  }

  # Function for Calculating Species Richness w/ Resampling
richness.iter <- function(otu.matrix  = " ",
                          sample.size = " ",
                          iters       = " "){
  rich.matrix <- matrix(NA, dim(otu.matrix)[1], iters)
  rownames(rich.matrix) <- rownames(otu.matrix)
  for (i in 1:iters){
    temp.matrix <- rrarefy(otu.matrix, sample.size)
    rich.matrix[,i] <- rowSums((temp.matrix>0)*1)
  }
  rich.mean <- apply(rich.matrix, 1, mean)
  rich.se <- apply(rich.matrix, 1, sem)
  return(data.frame(S = rich.mean, S_se = rich.se))
}

# Function for Simpsons Evenness w/ Resampling
evenness.iter <- function(otu.matrix  = " ",
                          sample.size = " ",
                          iters       = " "){
  simp_even <- function(SAD = " "){
  SAD <- subset(SAD, SAD > 0)
  S <- length(SAD); N <- sum(SAD); X <- rep(NA, S)
  for (i in 1:S){
    X[i] <- (SAD[i]*(SAD[i] - 1)) / (N * (N - 1))
    }
  D <- sum(X); e_d <- (1/D)/S
  return(e_d)
  }
  even.matrix <- matrix(NA, dim(otu.matrix)[1], iters)
  rownames(even.matrix) <- rownames(otu.matrix)
  for (i in 1:iters){
    temp.matrix <- rrarefy(otu.matrix, sample.size)
    even.matrix[,i] <- apply(temp.matrix, 1, simp_even)
    }
  even.mean <- apply(even.matrix, 1, mean)
  even.se <- apply(even.matrix, 1, sem)
  return(data.frame(E = even.mean, E_se = even.se))
}
