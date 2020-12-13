################################################################################
#                                                                              #
# Evolution Canyon Project: Microbial Community PL-SDA                         #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Jay Lennon and Mario Muscarella                                  #
#                                                                              #
# Last update: 2014/07/21                                                      #
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code creates the function ec.plsda for analysis of Evolution     #
#        Canyon simulation and emperical data                                  #
#        The function also creates function ec.plsda.plot for plotting         #
#                                                                              #
#                                                                              #
# Recent Changes:                                                              #
#         1. Calculate Variance Explained                                      #
#         2. Confirm Normalizations                                            #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1. ID influential Taxa                                               #
#         2. Load Taxonomy Data                                                #
#         3. Consider Rare Taxa Removal                                        #
#         4. Confirm Problamatic Taxa (EC-2A-D, EC-2A-R, EC-2C-R, EC-2D-R)     #
#         5. Compare with PERMANOVA                                            #
#         6. Paired PERMANOVA                                                  #
#                                                                              #
################################################################################

ec.plsda <- function(shared     = " ",
                     cutoff     = "0.03",
                     design     = " "){

  source("./bin/DiversityFunctions.r")  
  require("vegan")||install.packages("vegan");require("vegan")
  require("mixOmics")||install.packages("mixOmics");require("mixOmics")
  
  # Import Data (Site by OTU Matrix)
  ec_data <- t(read.otu(shared, cutoff))  # rows = sp, col = site
  design  <- read.delim(design, header=T, row.names=1)
  
  # Remove OTUs with <2 individuals
  ec_data <- ec_data[which(rowSums(ec_data) >= 2),] # removes rows (sp) with <2
  
  # Data Transformation: Relativize, Transpose, Log10-Transform
  Xrel     <- ec_data
	              for(i in 1:ncol(ec_data)){
                  Xrel[,i]=ec_data[,i]/sum(ec_data[,i])
                } 
  Xrelt    <- t(Xrel)
  Xlogrelt <- decostand(Xrelt,method="log")[,which(colSums(Xrelt) !=0)]
  
  # Create Factors for Multilevel Model
  slope     <- design$slope       # Factor 1
  molecule  <- design$molecule    # Factor 2 
  site      <- design$site
  station   <- design$station
  paired    <- design$paired
  
  # Define Y Matrix: Factors 1 and 2
  slope.molecule <- data.frame(cbind(as.character(slope),
                                     as.character(molecule))) 

  #4 -- Run multilevel models
  EC_multilevel <- multilevel(Xlogrelt, 
                               cond   = slope.molecule,
                               sample = paired, 
                               ncomp  = 3, 
                               method = 'splsda')
                               
  ec.plsda.out <- list("EC_multilevel"   = EC_multilevel, 
                       "Orig_Log_Dists"  = Xlogrelt) 
  
  return(ec.plsda.out)  
  }



  ec.plsda.plot <- function(plsda.in        = " ",
                            plot.title   = "test",
                            output.dir   = "./plots/"){

  source("./bin/DiversityFunctions.r")  
  require("vegan")||install.packages("vegan");require("vegan")
  require("mixOmics")||install.packages("mixOmics");require("mixOmics")
  
  # Plot PLSDA ordination  
  # Define plot margines and outer margin area
  par(mfrow=c(1,1), 
      mar=c(5,5,1,1), 
      oma=c(1,1,1,1)+0.1 ) 
  layout(rbind(1, 2), 
         height=c(7, 1))
  
  # Define plotting points
  plsda <- get(plsda.in)
  points  <- plsda$EC_multilevel$variates$X

  
  # Variation Explained by axis
  # Calculate distance between samples (Bray Curtis or Euclidean?)
  X.dist  <- vegdist(plsda$Orig_Log_Dists,method="euclidean")
  
  # Calculate distance between samples in reduced (ordination) space
  plsda.1 <- dist(plsda$EC_multilevel$variates$X[,1],method="euclidean")
  plsda.2 <- dist(plsda$EC_multilevel$variates$X[,2],method="euclidean")
  plsda.3 <- dist(plsda$EC_multilevel$variates$X[,3],method="euclidean")
  
  # Calculate variation explained
  var1    <- round(cor(X.dist, plsda.1)*100,2) 
  var2    <- round(cor(X.dist, plsda.2)*100,2) 
  var3    <- round(cor(X.dist, plsda.3)*100,2) 
  
  # Define Plot Limits
  lim1 <- c(min(points[,1])+min(points[,1])*0.1,
            max(points[,1])+max(points[,1])*0.1)
  lim2 <- c(min(points[,2])+min(points[,2])*0.1,
            max(points[,2])+max(points[,2])*0.1) 

  # Base Plot
  plot(points[,1], 
       points[,2], 
       xlim = lim1,
       ylim = lim2,
       pch  = 16, 
       cex  = 2.0, 
       type = "n",
       xaxt = "n", 
       yaxt = "n", 
       cex.lab  = 1.5, 
       cex.axis = 1.2,
       xlab=paste("PLS-DA Axis 1 (",var1, "%)", sep=""), 
       ylab=paste("PLS-DA Axis 2 (",var2, "%)", sep=""),)
        
  axis(side=1, las=1) # add x-axis ticks and labels  
  axis(side=2, las=1) # add y-axis ticks and labels   
  abline(h=0, lty="dotted") # add horizontal dashed line at 0 
  abline(v=0, lty="dotted") # add vertical dashed line at 0 
  mol.shape <- rep(NA, dim(points)[1])
                 for (i in 1:length(mol.shape)){
                   if (plsda$EC_multilevel$name.time[i] == "DNA"){
                     mol.shape[i] = 21}
                   else {mol.shape[i] = 22}
                 } # identifies symbol shape based on molecule
  slope.color <- rep(NA, dim(points)[1])
                   for (i in 1:length(slope.color)){
                     if (plsda$EC_multilevel$name.condition[i] == "North") {
                       slope.color[i] = "brown"}
                     else {slope.color[i] = "green3"}
                   } 
  points(points[,1], 
         points[,2], 
         pch=mol.shape, 
         cex=2.0, 
         col="black", 
         bg=slope.color, 
         lwd=2)

  slope.molecule <- data.frame(cbind(
                    as.character(plsda$EC_multilevel$name.condition),
                    as.character(plsda$EC_multilevel$name.time)))
                    
  slope.molecule.c <- do.call(paste, c(slope.molecule[c("X1", "X2")], sep = ""))
     
  ordiellipse(cbind(points[,1], points[,2]), 
              slope.molecule.c,
              kind="sd",
              conf=0.95, 
              lwd=2, lty=2, 
              draw = "lines", 
              col = "black", 
              label=TRUE, 
              cex=1.2)

  box(lwd=2)
  par(mar=c(0, 3, 0, 0))
  plot.new()
  legend("center", c(
    paste("All; ",levels(slope.molecule[,1])[1]," Slope", sep=""), 
    paste("All; ",levels(slope.molecule[,1])[2]," Slope", sep=""), 
    paste("Active; ",levels(slope.molecule[,1])[1]," Slope", sep=""),
    paste("Active; ",levels(slope.molecule[,1])[2]," Slope", sep="")), 
    pt.lwd=2, col="black", pt.bg=c("brown", "green3", "brown", 
    "green3"), pch=c(21,21,22,22), bty='n', ncol=2, cex=1.2, pt.cex=2)

  }
  
  
  

  
