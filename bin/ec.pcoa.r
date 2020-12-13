################################################################################
#                                                                              #
#	Evolution Canyon Project: Microbial Community PcoA and PERMANOVA           #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                               #
#                                                                              #
#	Last update: 2014/03/11                                                    #
#                                                                              #
################################################################################

ec.pcoa <- function(shared = " ", design = " ", plot.title = "test"){

  source("../bin/DiversityFunctions.r")
  require(vegan)

  # Import Site by OTU Matrix
  ec_data <- t(read.otu(shared, "0.03"))
  design <- read.delim(design, header=T, row.names=1)

  # Note: owing to amplification issues, we only sequenced 76 of 80 samples.
  # The four samples not included are C-1E-R, EC-2G-R, EC-2J-R, EC-6I-D

  # Remove problematic samples
  ec_data_red <- ec_data[,-c(20, 21, 25, 27)] #EC_2A_D, EC_2A_R, EC_2C_R, EC_2D_R

  # Remove 'alona' OTUs
  ec_data_red <- ec_data_red[rowSums(ec_data_red) > (0.001 * sum(colSums(ec_data_red))),]

  # Remove Zero Sum OTUs
  ec_data_red <- ec_data_red[!(rowSums(abs(ec_data_red)) ==0),] # zero occurrences after problem-sample removal

  # Check on fancy code to remove OTU that are rare across samples

  # Calculate Presense Absence
  dataPA <- (ec_data > 0)*1

  # Calculating Relative Abundance
  dataREL <- ec_data
  for(i in 1:ncol(ec_data)){
    dataREL[,i] = ec_data[,i]/sum(ec_data[,i])
    }

  # Create Distance Matrix with bray (deafault), manhattan, euclidean, canberra, bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup, binomial, or chao. Most should be part of vegan, but possilbly 'labdsv' or 'BiodiversityR' packages
  samplePA.dist <- vegdist(t(dataPA),method="bray")
  sampleREL.dist <- vegdist(t(dataREL),method="bray")

  # Principal Coordinates Analysis
  ec_pcoa <- cmdscale(sampleREL.dist,k=3,eig=TRUE,add=FALSE)
    # Classical (Metric) Multidimensional Scaling; returns PCoA coordinates
    # eig=TRUE returns eigenvalues; k = # of dimensions to calculate

  # Percent Variance Explained Using PCoA (Axis 1,2,3)
  explainvar1 <- round(ec_pcoa$eig[1]/sum(ec_pcoa$eig)*100,2)
  explainvar2 <- round(ec_pcoa$eig[2]/sum(ec_pcoa$eig)*100,2)
  explainvar3 <- round(ec_pcoa$eig[3]/sum(ec_pcoa$eig)*100,2)

  pcoap <- merge(as.data.frame(ec_pcoa$points),design,by=0,all.x=T)[,-1]
  rownames(pcoap) <- rownames(ec_pcoa$points)
  pcoap$slope.molecule <- paste(pcoap$slope, pcoap$molecule, sep="")

  # Plot Parameters
  par(mfrow=c(1,1), mar=c(5,5,1,1))
  layout(rbind(1, 2), height=c(7, 1))
  x.dim <- c(min(pcoap$V2)+min(pcoap$V2)*0.2,max(pcoap$V2)+max(pcoap$V2)*0.2)
  y.dim <- c(min(pcoap$V1)+min(pcoap$V1)*0.2,max(pcoap$V1)+max(pcoap$V1)*0.2)

  # Initiate Plot
  plot(pcoap$V2, pcoap$V1, xlab=paste("PCoA Axis 2 (",explainvar2, "%)", sep="")
    , ylab=paste("PCoA Axis 1 (",explainvar1, "%)", sep=""),
    xlim=x.dim,ylim= y.dim, pch=16, cex=2.0, type="n",xaxt="n",
    yaxt="n", cex.lab=1.5, cex.axis=1.2)
  axis(side=1, las=1)
  axis(side=2, las=1)
  abline(h=0, lty="dotted")
  abline(v=0, lty="dotted")
  mol.shape <- rep(NA, dim(pcoap)[1])
    for (i in 1:length(mol.shape)){
      if (pcoap$molecule[i] == "DNA"){mol.shape[i] = 21}
      else {mol.shape[i] = 22}
      }
  slope.color <- rep(NA, dim(pcoap)[1])
    for (i in 1:length(slope.color)){
      if (pcoap$slope[i] == levels(pcoap$slope)[1]) {slope.color[i] = "brown"}
      else {slope.color[i] = "green3"}
      }
  points(pcoap$V2, pcoap$V1, pch=mol.shape, cex=2.0, col="black", bg=slope.color, lwd=2)
  ordiellipse(cbind(pcoap$V2, pcoap$V1), pcoap$slope.molecule, kind="sd", conf=0.95,
    lwd=2, lty=3, draw = "lines", col = "black", label=TRUE)
  # legend("topleft", c(paste("All; ",levels(pcoap$slope)[1]," Slope", sep=""),
  #   paste("All; ",levels(pcoap$slope)[2]," Slope", sep=""),
  #   paste("Active; ",levels(pcoap$slope)[1]," Slope", sep=""),
  #   paste("Active; ",levels(pcoap$slope)[2]," Slope", sep="")),
  #   pt.lwd=2, col="black", pt.bg=c("brown3", "green3", "brown3",
  #   "green3"), pch=c(21,21,22,22), bty='o', box.lty=0, bg="white", cex=1.5)
  box(lwd=2)
  par(mar=c(0, 3, 0, 0))
  plot.new()
  legend("center", c(paste("All; ",levels(pcoap$slope)[1]," Slope", sep=""),
    paste("All; ",levels(pcoap$slope)[2]," Slope", sep=""),
    paste("Active; ",levels(pcoap$slope)[1]," Slope", sep=""),
    paste("Active; ",levels(pcoap$slope)[2]," Slope", sep="")),
    pt.lwd=2, col="black", pt.bg=c("brown", "green3", "brown",
    "green3"), pch=c(21,21,22,22), bty='n', ncol=2, cex=1.5, pt.cex=2)






  dev.copy2pdf(file=paste("../figures/",plot.title,".pdf",sep=""))
  dev.copy(png, file=paste("../figures/",plot.title,".png",sep=""), width=72*(7*4),
    height=72*(8*4), res=72*4)
  dev.off()

  # Adonis (PERMANOVA)
  # Adonis runs a PERMANOVA (Created by Marti J. Anderson)
  # this is very similar to ANOVA but for multivariate data
  # You can make very complex experimental designs with it
  # The default distance measure is bray-curtis, but other measures
  # (Chao, Jaccard, Euclidean) can be used when specified
#  Adonis <- adonis(sampleREL.dist ~ design$molecule*design$slope, method="bray",
#    permutations=1000)
#    return(Adonis)
  }



ec.pcoa.fun <- function(shared = " ", cutoff = " ", design = " "){

  source("../bin/DiversityFunctions.r")
  require("vegan")||install.packages("vegan");require("vegan")


  # Import Data (Site by OTU Matrix)
  ec_data <- t(read.otu(shared, cutoff))  # rows = sp, col = site
  design  <- read.delim(design, header=T, row.names=1)

  # Remove OTUs with <2 individuals
  ec_data <- ec_data[which(rowSums(ec_data) >= 2),]

  # Data Transformation: Relativize, Transpose, Log10-Transform
  Xrel     <- ec_data
	              for(i in 1:ncol(ec_data)){
                  Xrel[,i]=ec_data[,i]/sum(ec_data[,i])
                }
  Xrelt    <- t(Xrel)     # rows = site , col = sp
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

  # Calculate Presense Absence
  dataPA <- (ec_data > 0)*1

  # Calculating Relative Abundance
  dataREL <- ec_data
  for(i in 1:ncol(ec_data)){
    dataREL[,i] = ec_data[,i]/sum(ec_data[,i])
    }

  dataREL.log <- decostand(dataREL,method="log")

  # Create Distance Matrix with bray (deafault), manhattan, euclidean,
  #  canberra, bray, kulczynski, jaccard, gower, altGower, morisita, horn,
  #  mountford, raup, binomial, or chao. Most should be part of vegan,
  #  but possilbly 'labdsv' or 'BiodiversityR' packages
  samplePA.dist <- vegdist(t(dataPA),method="bray")
  sampleREL.dist <- vegdist(Xlogrelt,method="bray")

  # Principal Coordinates Analysis
  ec_pcoa <- cmdscale(sampleREL.dist,k=3,eig=TRUE,add=FALSE)
    # Classical (Metric) Multidimensional Scaling; returns PCoA coordinates
    # eig=TRUE returns eigenvalues; k = # of dimensions to calculate

  # Percent Variance Explained Using PCoA (Axis 1,2,3)
  explainvar1 <- round(ec_pcoa$eig[1]/sum(ec_pcoa$eig)*100,2)
  explainvar2 <- round(ec_pcoa$eig[2]/sum(ec_pcoa$eig)*100,2)
  explainvar3 <- round(ec_pcoa$eig[3]/sum(ec_pcoa$eig)*100,2)

  pcoap <- merge(as.data.frame(ec_pcoa$points),design,by=0,all.x=T)[,-1]
  rownames(pcoap) <- rownames(ec_pcoa$points)

  pcoa.out <- list("pcoap" = pcoap,
                   "explainvar1" = explainvar1,
                   "explainvar2" = explainvar2,
                   "explainvar3" = explainvar3)

  return(pcoa.out)

  }

ec.pcoa.plot <- function(pcoa.in    = "",
                         plot.title = "",
                         output.dir = ""){

  source("../bin/DiversityFunctions.r")
  require("vegan")||install.packages("vegan");require("vegan")

  pcoap <- pcoa.in$pcoap
  explainvar1 <- pcoa.in$explainvar1
  explainvar2 <- pcoa.in$explainvar2
  explainvar3 <- pcoa.in$explainvar3

  # Plot Parameters
  par(mfrow=c(1,1),
      mar=c(5,5,1,1),
      oma=c(1,1,1,1)+0.1 )
  layout(rbind(1, 2),
         height=c(7, 1))

  x.dim <- c(min(pcoap$V1)+min(pcoap$V1)*0.2,max(pcoap$V1)+max(pcoap$V1)*0.2)
  y.dim <- c(min(pcoap$V2)+min(pcoap$V2)*0.2,max(pcoap$V2)+max(pcoap$V2)*0.2)
  mol.shape <- rep(NA, dim(pcoap)[1])
    for (i in 1:length(mol.shape)){
      if (pcoap$molecule[i] == "DNA"){mol.shape[i] = 21}
      else {mol.shape[i] = 22}
      }
  slope.color <- rep(NA, dim(pcoap)[1])
    for (i in 1:length(slope.color)){
      if (pcoap$slope[i] == levels(pcoap$slope)[1]) {slope.color[i] = "brown"}
      else {slope.color[i] = "green3"}
      }

  # Initiate Plot
  plot(pcoap$V1, pcoap$V2,
       xlab=paste("PCoA Axis 1 (",explainvar1, "%)", sep=""),
       ylab=paste("PCoA Axis 2 (",explainvar2, "%)", sep=""),
       xlim=x.dim,ylim= y.dim, pch=16, cex=2.0, type="n",xaxt="n",
       yaxt="n", cex.lab=1.5, cex.axis=1.2)
  points(pcoap$V1, pcoap$V2,
         pch=mol.shape, cex=2.0, col="black", bg=slope.color, lwd=2)


  slope.molecule <- data.frame(cbind(
                    as.character(pcoap$slope),
                    as.character(pcoap$molecule)))

  slope.molecule.c <- do.call(paste, c(slope.molecule[c("X1", "X2")], sep = ""))

  ordiellipse(cbind(pcoap$V1, pcoap$V2),
              slope.molecule.c,
              kind="sd",
              conf=0.95,
              lwd=2, lty=2,
              draw = "lines",
              col = "black",
              label=TRUE,
              cex=1.2)

  axis(side=1, las=1, at=c(seq(-0.2,2,by=0.1)))
  axis(side=2, las=1, at=c(seq(-0.2,2,by=0.1)))
  abline(h=0, lty="dotted")
  abline(v=0, lty="dotted")
  box(lwd=2)

  # Legnd
  par(mar=c(0, 3, 0, 0))
  plot.new()
  legend("center", c(paste("All; ",levels(pcoap$slope)[1]," Slope", sep=""),
    paste("All; ",levels(pcoap$slope)[2]," Slope", sep=""),
    paste("Active; ",levels(pcoap$slope)[1]," Slope", sep=""),
    paste("Active; ",levels(pcoap$slope)[2]," Slope", sep="")),
    pt.lwd=2, col="black", pt.bg=c("brown", "green3", "brown",
    "green3"), pch=c(21,21,22,22), bty='n', ncol=2, cex=1.2, pt.cex=2)

#  dev.copy2pdf(file=paste("../plots/",plot.title,".pdf",sep=""))
#  dev.copy(png, file=paste("../plots/",plot.title,".png",sep=""), width=72*(7*4),
#    height=72*(8*4), res=72*4)
#  dev.off()


  }

