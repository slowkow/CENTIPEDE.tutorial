#' Fit a motif with CENTIPEDE and plot the profile and posterior distribution.
#'
#' @param motif A string of a specific motif
#' @param cen an object
#' @param score fimo score cutoff
#' @param p.value p.value cutoff
#' @return plot of profile and posterior distribution of a specific motif
#' @details
#'    cen is a list returned by \code{\link{centipede_data}}. It has a dataframe
#'    called \code{regions} with one row for each putative motif match.
#'    It also has a matrix called \code{mat} with one row for each match,
#'    and each column is a nucleotide position relative to the center of the
#'    motif match. The values in this matrix are the count of read start sites
#'    at that position.
fitAndPlot <- function(motif, cen, score = NA, pvalue = NA) {
  # Match the name of the motif.
  idx <- cen$regions$pattern.name == motif
  if (!is.na(pvalue)) {
    idx <- idx & cen$regions$p.value <= pvalue
  }
  if (!is.na(score)) {
    idx <- idx & cen$regions$score >= score
  }
  motif.tfs <- list(
    mat = cen$mat[idx, , drop = FALSE],
    regions = cen$regions[idx, ]
  )

  # Estimate the motif length.
  motif.length <- median(sapply(
    X = motif.tfs$regions$matched.sequence,
    FUN = function(x) nchar(as.character(x))
  ))

  # Fit the model.
  motif.fit <- CENTIPEDE::fitCentipede(
    Xlist = list(DNase = motif.tfs$mat),
    Y = as.matrix(data.frame(
      Intercept = rep(1, nrow(motif.tfs$mat))
    ))
  )

  # Plot image of cutsites ranked by CENTIPEDE posteriors.
  #imageCutSites(motif.tfs$mat[order(motif.fit$PostPr),][c(1:100, (dim(motif.tfs$mat)[1]-100):(dim(motif.tfs$mat)[1])),])
  #imageCutSitesCombined(motif.tfs$mat[order(motif.fit$PostPr),][c(1:100, (dim(motif.tfs$mat)[1]-100):(dim(motif.tfs$mat)[1])),])

  # Plot estimated footprint.
  #imageCutSitesCombined(motif.tfs$mat[order(motif.fit$PostPr),][motif.fit$PostPr >= 0.9,])
  CENTIPEDE::plotProfile(
    Profile = motif.fit$LambdaParList[[1]],
    Mlen = motif.length
  )
  hist(x = motif.fit$PostPr, main = motif,  xlab = "Posterior Probability")
  legend(
    "top",
    legend = c(
      paste("Total Sites:", length(motif.fit$PostPr)),
      paste("Sites >= 0.99:", sum(motif.fit$PostPr >= 0.99))
    )
  )
}

plot_all_motifs <- function(cen_dat, folder = "motif_plots", score = 0, p.value = 0) {
  tfs <- readRDS(cen_dat)

  motifs <- unique(tfs$regions$pattern.name)
  dir.create(folder)
  pattern.removed <- c()

  for (motif in motifs) {

    filename <- file.path(folder, sprintf("%s.pdf", motif))

    pdf(filename)
    par(mfrow = c(2,1))

    res <- try({
      fitAndPlot(motif, tfs, p.value = p.value, score = score)
    }, silent = TRUE)

    if (inherits(res, "try-error")) {
      dev.off()
      file.remove(filename)
      pattern.removed <- append(pattern.removed, motif)
      next
    }

    dev.off()
  }

  cat(
    "######################\n ", length(pattern.removed), " removed Motifs:\n"
  )
  for (i in pattern.removed) {
    print(i)
  }
}
