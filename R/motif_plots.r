# script to filter for specific motifs to process with centipede

library(CENTIPEDE)


#' extract one motif from an object generated with centipede_data
#' @param list object generated with centipede_data 
#' @param index indices of motif occurences in list
#' 
extractMotif <- function(list,index) {
  mat <- list[[1]][index,,drop = FALSE]
  regions <- list[[2]][index,]
  
  list(mat=mat,regions=regions)
}

#' function to fit a motif with centipede and directly plot the profile and posterior distribution
#' @param motif A string of a specific motif
#' @param tfs an object 
#' @param score fimo score cutoff
#' @param p.value p.value cutoff
#' @return plot of profile and posterior distribution of a specific motif
#' 
#' @details  tfs is an object generated with centipede_data(). It contains regions and a matrix where each row is a candidate motif match, 
# and each column is a position relative to the center of the motif match. 
# The values in this matrix will be number of read start- sites that occur at that position.
#  
fitAndPlot <- function(motif,tfs,score=0,p.value=0) {

  # extract a single motif 
  motif.sites <- which(tfs$regions$pattern.name == motif)
  motif.tfs <- extractMotif(tfs,motif.sites)
  if ( (score != 0) || (p.value != 0) ) {
    
    if (p.value != 0) { 
      motif.p.value <- which(-log10(motif.tfs$regions$p.value) >= p.value)
      motif.tfs <- extractMotif(motif.tfs,motif.p.value)
    }
    if (score != 0) {
      motif.score <- which(motif.tfs$regions$score >= score)
      motif.tfs <- extractMotif(motif.tfs,motif.score)
    }
  }
  
  # estimate the motif length
  motif.length <- median(sapply(motif.tfs$regions$matched.sequence,FUN = function(x) nchar(as.character(x))))
  # fit the model
  motif.fit <- fitCentipede(Xlist = list(DNase = motif.tfs$mat), 
                                  Y = as.matrix(data.frame(Intercept = rep(1, nrow(motif.tfs$mat)) )))
  
  #PLOTS IMAGE OF CUTSITES RANKED BY CENTIPEDE POSTERIORS
  #imageCutSites(motif.tfs$mat[order(motif.fit$PostPr),][c(1:100, (dim(motif.tfs$mat)[1]-100):(dim(motif.tfs$mat)[1])),])
  #imageCutSitesCombined(motif.tfs$mat[order(motif.fit$PostPr),][c(1:100, (dim(motif.tfs$mat)[1]-100):(dim(motif.tfs$mat)[1])),])
  
  #PLOT ESTIMATED FOOTPRINT
  #imageCutSitesCombined(motif.tfs$mat[order(motif.fit$PostPr),][motif.fit$PostPr >= 0.9,])
  plotProfile(motif.fit$LambdaParList[[1]], Mlen = motif.length)
  hist(motif.fit$PostPr,main = motif,  xlab = "Posterior Probability")
  legend("top",legend = c(paste("Total Sites:", length(motif.fit$PostPr)),paste("Sites >= 0.99:", sum(motif.fit$PostPr>=0.99))))

}



# load the regions and a matrix where each row is a candidate motif match, 
# and each column is a position relative to the center of the motif match. 
# The values in this matrix will be number of read start- sites that occur at that position.
#tfs <- readRDS("tfs_Psoas_Muscle.rds")



plot_all_motifs <- function(cen_dat,folder="motif_plots",score=0,p.value=0) {
  tfs <- readRDS(cen_dat)
  # extract all unique pattern names
  pattern.list <- unique(tfs$regions$pattern.name)
  dir.create(folder)
  pattern.removed <- c()
  
  for (motif in pattern.list ) {
    
    pdf(paste0(folder,"/",motif,".pdf"),paper = "a4")
    par(mfrow=c(2,1))
    
    res <- try( fitAndPlot(motif,tfs,p.value = p.value, score = score),silent = TRUE )
    if(inherits(res, "try-error"))
    {
      #error handling code, maybe just skip this iteration using
      dev.off()
      file.remove(paste0(folder,"/",motif,".pdf"))
      pattern.removed = append(pattern.removed,motif)
      next
    }
    #rest of iteration for case of no error
    
    dev.off()
    #par(old.par)
    
  }
  
  cat("######################\n ",length(pattern.removed)," removed Motifs:\n")
  for ( i in pattern.removed) print(i)
}