#'
#' This function generates a biplot for the results of the PCovR analysis.
#'
#' @param X A matrix containing the data for the X variables.
#' @param Y A matrix containing the data for the Y variables.
#' @param PcovR A list containing the results of the PCovR analysis.
#' @return A list containing the biplot information.
#' @export
#'
Biplot.PCovR <- function(X, Y, PcovR) {
  # Carga de librerías necesarias
  if (!requireNamespace("MultBiplotR", quietly = TRUE)) {
    install.packages("MultBiplotR")
  }

  I = dim(X)[1]
  J = dim(X)[2]
  K = dim(Y)[2]

  Xprep <- switch(PcovR$prepX,
                  stand = nrm2(scale(X, center = TRUE, scale = FALSE)) * sqrt(dim(X)[1]),
                  cent = scale(X, center = TRUE, scale = FALSE))

  Yprep <- switch(PcovR$prepY,
                  stand = nrm2(scale(Y, center = TRUE, scale = FALSE)) * sqrt(dim(X)[1]),
                  cent = scale(Y, center = TRUE, scale = FALSE))

  S <- dim(PcovR$Te)[2]

  Biplot <- list()
  Biplot$Title <- " PCovR - Biplot"
  Biplot$Type <- "PCovR"

  if(results$prepX == "stand") {
    Biplot$Initial_Transformation <- 5
  } else {
    Biplot$Initial_Transformation <- 4
  }
  Biplot$ncols <- dim(X)[2]
  Biplot$nrows <- dim(X)[1]
  Biplot$Dimension <- S
  Biplot$alpha <- 0
  Biplot$Means <- apply(X, 2, mean)
  Biplot$Medians <- apply(X, 2, median)
  Biplot$Deviations <- apply(X, 2, sd)
  Biplot$Minima <- apply(X, 2, min)
  Biplot$Maxima <- apply(X, 2, max)
  Biplot$P25 <- apply(X, 2, quantile)[2, ]
  Biplot$P75 <- apply(X, 2, quantile)[4, ]

  a <- PcovR$Te
  b <- t(PcovR$Px)
  sca <- sum(a^2) / dim(X)[1]
  scb <- sum(b^2) / dim(X)[2]
  scf <- sqrt(sqrt(scb / sca))
  a <- a * scf
  b <- b / scf

  Biplot$RowCoordinates <- a
  Biplot$ColCoordinates <- b

  Cont <- CalculateContributions(Xprep, as.matrix(PcovR$Te), as.matrix(t(PcovR$Px)))
  Biplot$EigenValues <- Cont$Fit
  Biplot$Inertia <- Cont$Fit * 100
  Biplot$CumInertia <- cumsum(Biplot$Inertia)
  Biplot$RowContributions <- Cont$RowContributions
  Biplot$ColContributions <- Cont$ColContributions
  Biplot$Structure <- Cont$Structure
  class(Biplot) <- "ContinuousBiplot"

  YBiplot <- list()
  YBiplot$Title <- " PCovR - Biplot (Y)"
  YBiplot$Type <- "PCovR"
  if(results$prepY == "stand") {
    YBiplot$Initial_Transformation <- 5
  } else {
    YBiplot$Initial_Transformation <- 4
  }
  YBiplot$ncols <- dim(Y)[2]
  YBiplot$Means <- apply(Y, 2, mean)
  YBiplot$Medians <- apply(Y, 2, median)
  YBiplot$Deviations <- apply(Y, 2, sd)
  YBiplot$Minima <- apply(Y, 2, min)
  YBiplot$Maxima <- apply(Y, 2, max)
  YBiplot$P25 <- apply(Y, 2, quantile)[2, ]
  YBiplot$P75 <- apply(Y, 2, quantile)[4, ]
  YBiplot$b0 <- rep(0, dim(Y)[2])
  YBiplot$ColCoordinates <- t(PcovR$Py) / scf
  class(YBiplot) <- "ContSupVarsBiplot"
  Biplot$ContSupVarsBiplot <- YBiplot
  return(Biplot)
}
