################################################################################
#                                                                              #
#	Evolution Canyon Project: Microbial Community PL-SDA                         #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Jay Lennon and Mario Muscarella                                  #
#                                                                              #
#	Last update: 2014/07/01                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/evolution-canyon")
source("./bin/DiversityFunctions.r")
source("./bin/ec.pcoa.r")
source("./bin/DiversityFunctions.r")
require("vegan")
require("mixOmics")
require("som")
require("DiscriMiner")

# Data Files
shared     = "./mothur/EC.bac.final.shared"
design     = "./data/design.txt"

#1 -- Import data (site by OTU matrix) and deal w/ "problem" samples
ec_data <- t(read.otu(shared, "0.03")) # takes a long time (10 mins)
design <- read.delim(design, header=T, row.names=1)

  # Note: owing to amplification issues, we only submitted 76 of 80 samples.
  # The four samples NOT submitted were C-1E-R, EC-2G-R, EC-2J-R, EC-6I-D
  # So, removed their pairs for PLS-DA: EC-1E-D, EC-2G-D, EC-2J-D,EC-6I-R
  # 4 samples have had crappy quality in past: EC-2A-D, EC-2A-R, EC-2C-R, EC-2D-R
  # Automatically exlcuded, but may want to look again
  # Remove their pairs EC-2C-D, EC-2D-R
  # Results in a total of 66 samples (14 problematics)
  # Use following code to identify columns to remove

samps <- (colnames(ec_data))
which(samps == "EC-2A-D")

ec_data_red <- ec_data[,-c(9,20:21,24:27,32,37,74)]
# removes problematic samples and their pairs

samps_red <- colnames(ec_data_red) # recover samples from reduced dataset
design_red <- design[samps_red,] # design matrix w/o problem samples & pairs

#2 -- Options for transformation, relativation, and normalization
Xt <- t(ec_data_red) # transpose sample x OTU matrix; necssary for 'multilelve'
Xlogt <- decostand(Xt,method="log")
Xpa <- (Xt > 0)*1 # Calculate Presense Absence
Xrel <- ec_data_red
  for(i in 1:ncol(ec_data_red)){
    Xrel[,i] = ec_data_red[,i]/sum(ec_data_red[,i])}
Xrelt <- t(Xrel)
Xlog <- decostand(ec_data_red,method="log")
Xrellog <- Xlog
  for(i in 1:ncol(Xlog)){
    Xrellog[,i] = Xlog[,i]/sum(Xlog[,i])}
Xrellogt<-t(Xrellog)
X <- normalize(Xrelt,byrow=TRUE) # normalizing by row (mean = 0, var = 1) with som package, not mixOmics

  # Some comments on normalizatio: In mixOmics "default is PLS centers data by substracting mean of each column
  # (variables) and scaling with the stdev. Result is that each column has a mean zero and a variance of 1
  # personal communication with Kim-Anh Le Cao (mixOmics developer). However, it appears that 'multilevel'
  # procedure doesn't run w/o normalization

#3 -- PLS-DA using DiscriMiner
ECplsDA <- plsDA(X, design_red$paired_within_slope, comps=3)



  slope <- design_red$slope # factor 1
  molecule <- design_red$molecule # factor 2
  slope.molecule <- data.frame(cbind(as.character(slope),as.character(molecule))) # Y matrix with factor 1 and 2
  slope.molecule.concat <- do.call(paste, c(slope.molecule[c("X1", "X2")], sep = "")) # create unique treat ID vector
  paired <- design_red$paired_across_slope # identifies paired samples (i.e., DNA-RNA)
  EC_multilevel <- multilevel(X, cond = slope.molecule, sample = paired, ncomp = 3, method = 'splsda')
  # change up "X" in multilevel to reflect transformation, etc. above

#4 -- Contribution (%) variance of factors (Y) to PLS-DA axes
  # mostly lifed from http://perso.math.univ-toulouse.fr/mixomics/faq/numerical-outputs/
  # not sure what "Rd" means. U are the "variates"...
  Y <- slope # maybe this should be done one-at-a-time (e.g., slope and molecule)
  Rd.YvsU = cor(as.numeric(as.factor(Y)),EC_multilevel$variates$X) # turns categorical vars into quantitative
  Rd.YvsU = apply(Rd.YvsU^2, 2, sum)
  Rd.Y = cbind(Rd.YvsU, cumsum(Rd.YvsU))
  colnames(Rd.Y) = c("Proportion", "Cumulative")
  Rd.Y # percent of variance explained by each component

#4 -- Variation Explained by axis
  # Calculate distance between samples (Bray Curtis or Euclidean?)
  X.dist  <- vegdist(t(ec_data_red),method="euclidean")
  # Calculate distance between samples in reduced (ordination) space
  plsda.1 <- dist(EC_multilevel$variates$X[,1],method="euclidean")
  plsda.2 <- dist(EC_multilevel$variates$X[,2],method="euclidean")
  plsda.3 <- dist(EC_multilevel$variates$X[,3],method="euclidean")
  # Calculate variation explained
  var1 <- cor(X.dist, plsda.1)
  var2 <- cor(X.dist, plsda.2)
  var3 <- cor(X.dist, plsda.3)

#5 -- PLOTTING
  points <- EC_multilevel$variates$X
  par(mar=c(5,5,1,1), oma=c(1,1,1,1)+0.1 )
  plot(points[,1], points[,2], xlab="PLSDA Axis 1 ", ylab="PLSDA Axis 2",
  xlim=c(min(points[,1])+min(points[,1])*0.1,max(points[,1])+max(points[,1])*0.1),
  ylim=c(min(points[,2])+min(points[,2])*0.1,max(points[,2])+max(points[,2])*0.1),
  pch=16, cex=2.0, type="n",xaxt="n", yaxt="n", cex.lab=1.5, cex.axis=1.2)
  axis(side=1, las=1)
  axis(side=2, las=1)
  abline(h=0, lty="dotted")
  abline(v=0, lty="dotted")
  mol.shape <- rep(NA, dim(points)[1])
    for (i in 1:length(mol.shape)){
      if (molecule[i] == "DNA"){mol.shape[i] = 21}
      else {mol.shape[i] = 22}
      }
  slope.color <- rep(NA, dim(points)[1])
    for (i in 1:length(slope)){
      if (slope[i] == levels(slope)[1]) {slope.color[i] = "brown2"}
      else {slope.color[i] = "green3"}
      }
  points(points[,1], points[,2], pch=mol.shape, cex=2.0, col="black", bg=slope.color, lwd=2)
  ordiellipse(cbind(points[,1], points[,2]), slope.molecule.concat, kind="sd", conf=0.95,
    lwd=2, lty=2, draw = "lines", col = "black", label=TRUE)

#6 -- Other stuff: some note and tries based on following website:
  # http://perso.math.univ-toulouse.fr/mixomics/methods/spls-da/
  # calculate the coefficients of the linear combinations
  pred <- predict(EC_multilevel, X[1:2, ])
  pred$B.hat
  # calculate R2 and Q2 values for sPLS-DA?
  # Q2 = "Q2 is the square of the correlation between the actual and predicted response"
  # warning = this is memory intenstive, but didn't crash
  val <- valid(EC_multilevel, criterion = c("R2", "Q2"))








