
#' QQ plot
#'
#' The function creates a QQ plot from observed and simulated data.
#'
#' @param x A data frame holding observations and the corresponding
#'   predictions.
#' @param colObs Name of the column in \code{x} holding the observations.
#' @param colsSim A vector of names identifying the columns in \code{x} that
#'   contain the predictions.
#' @param alphas A vector of values between 0 and 1 defining the quantils to
#'   be compared.
#' @param log If \code{TRUE}, the data are log-transformed before being analyzed.
#' @param log_add A small value to be added to the data before the
#'   log-transformation is done. This is required of the data contain zero values.
#'   Ignored if \code{log=FALSE}.
#' @param xlims,ylims Limits for the x-axis (observation data) and y-axis (simulated data axis). The default will show data
#'   between -3 and 3 standard deviations away from the mean.
#' @param xlabel,ylabel Labels for the x-axis (observation axis) and y-axis (simulated data axis).
#' @param xnumbers,ynumbers Should numerical annotations be printed on the axis? Default is \code{TRUE}.
#' @param title Title to be shown above the chart.
#' @param title.cex Controls font size of title.
#' @param showLegend Should the legend be visible?
#'
#' @return \code{NULL}
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @examples
#' d= data.frame(obs=runif(1000), sim1=runif(1000), sim2=runif(1000))
#' qqplt(x=d, colsSim=c("sim1","sim2"), title="Test data")

qqplt= function(
  x,
  colObs="obs",
  colsSim=c("sim"),
  alphas= c(c(0.01,0.05),seq(0.1, 0.9, 0.1), c(0.95, 0.99)),
  log=FALSE,
  log_add=1.e-03,
  xlims=-3:3,
  ylims=-3:3,
  xlabel="Obs.",
  ylabel="Sim.",
  xnumbers=TRUE,
  ynumbers=TRUE,
  title="",
  title.cex=1,
  showLegend=TRUE
) {
  # Check args
  stopifnot(is.data.frame(x))
  stopifnot(colObs %in% names(x))
  stopifnot(all(colsSim %in% names(x)))
  if (colObs %in% colsSim) stop("Bad selection of column name(s).")
  # Prepare plot
  plot(xlims, ylims, type="n", xlab=xlabel, ylab=ylabel, xaxt="n", yaxt="n")
  axis(side=1, labels=xnumbers)
  axis(side=2, labels=ynumbers)
  abline(0,1, lty=3)
  mtext(side=3, title, cex=title.cex)
  if (length(colsSim) <= 3) {
    pchs= c(1,3,4)
  } else {
    pchs= 1:length(colsSim)
  }
  # Tranform observations
  if (log) {
    obs= log(x[,colObs] + log_add)
  } else {
    obs= x[,colObs]
  }
  obs= (obs - mean(obs)) / sd(obs)
  # Loop through simulated data
  for (i in 1:length(colsSim)) {
    # Tranform sim. data
    if (log) {
      sim= log(x[,colsSim[i]] + log_add)
    } else {
      sim= x[,colsSim[i]]
    }
    sim= (sim - mean(sim)) / sd(sim)
    # Add to plot
    points(quantile(obs,alphas), quantile(sim,alphas), pch=pchs[i])
  }
  if (showLegend)
    legend("bottomright",bty="n",pch=pchs[1:length(colsSim)],legend=colsSim)
  return(invisible(NULL))
}

