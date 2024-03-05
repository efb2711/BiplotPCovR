#' PCovR Regression Biplot
#'
#' This function generates a regression biplot to examine the relationship
#' between the response and explanatory variables.
#'
#' @param PcovR A list containing the results of the PCovR analysis.
#' @return A list containing the biplot information.
#' @export
#'
Biplot.PCovR.Weights <- function(PcovR) {
  # Carga de librerías necesarias
  if (!requireNamespace("MultBiplotR", quietly = TRUE)) {
    install.packages("MultBiplotR")
  }
  Z = as.matrix(PcovR$W) %*% as.matrix(PcovR$Py)
  Biplot = list()
  Biplot$Title = " PCovR - Regression"
  Biplot$Type = "PCovR - Regression"
  Biplot$Initial_Transformation <- 4
  J = dim(PcovR$W)[1]
  K = dim(PcovR$Py)[2]
  S = dim(PcovR$W)[2]
  Biplot$ncols = K
  Biplot$nrows = J
  Biplot$Dimension = S
  Biplot$alpha = 0
  Biplot$Means = apply(Z, 2, mean)
  Biplot$Medians = apply(Z, 2, median)
  Biplot$Deviations = apply(Z, 2, sd)
  Biplot$Minima = apply(Z, 2, min)
  Biplot$Maxima = apply(Z, 2, max)
  Biplot$P25 = apply(Z, 2, quantile)[2, ]
  Biplot$P75 = apply(Z, 2, quantile)[4, ]
  a = PcovR$W
  b = t(PcovR$Py)
  sca = sum(a^2)
  scb = sum(b^2)
  sca = sca/J
  scb = scb/K
  scf = sqrt(sqrt(scb/sca))
  a = a * scf
  b = b/scf
  Biplot$RowCoordinates = a
  Biplot$ColCoordinates = b
  Cont = CalculateContributions(Z, as.matrix(PcovR$W),
                                as.matrix(t(PcovR$Py)))
  Biplot$EigenValues = Cont$Fit
  Biplot$Inertia = Cont$Fit * 100
  Biplot$CumInertia = cumsum(Biplot$Inertia)
  Biplot$RowContributions = Cont$RowContributions
  Biplot$ColContributions = Cont$ColContributions
  Biplot$Structure = Cont$Structure
  class(Biplot) = "ContinuousBiplot"
  return(Biplot)
}
