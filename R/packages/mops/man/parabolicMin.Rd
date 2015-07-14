\name{parabolicMin}
\alias{parabolicMin}
\title{Location of the minimum of a parabola fitted through 3 points}
\usage{
parabolicMin(x.left, y.left, x.center, y.center, x.right, y.right)
}
\arguments{
  \item{x.left}{x-value of leftmost point}

  \item{y.left}{y at x.left}

  \item{x.center}{x-value of center point. Must be \eqn{>
  x.left}.}

  \item{y.center}{y at x.center. It should be \eqn{y.center
  <= min(y.left, y.right)} for the fitted parabola to have
  a minimum in the range defined by \code{x.left} and
  \code{x.right}.}

  \item{x.right}{x-value of rightmost point. Must be \eqn{>
  x.center}.}

  \item{y.right}{y at x.right}
}
\value{
A list with two components \item{status}{An integer value.
A value of 0 indicates success. Other values indicate a
problem (1: collinearity, 2: local maximum).}
\item{x.min}{The x-value where the parabola has its minimum
if \eqn{status == 0}. Otherwise, the x-value corresponding
to the minimum of \code{y.left}, \code{y.center}, and
\code{y.right}.}
}
\description{
Given 3 points (x.left, y.left), (x.center, y.center),
(x.right, y.right), the function returns the x-value where
a parabola ax^2+bx+c fitted through these 3 points has a
minimum. This is one of the core functions used in
optimization methods based on successive parabolic
interpolation (like \code{\link{optimize}}, for example).
}
\examples{
# Determine minimum of fitted parabola
x= c(-2, 0.5, 5)
y= c(2, 1, 10)
tmp= parabolicMin(x[1], y[1], x[2], y[2], x[3], y[3])
print(tmp)
stopifnot(tmp$status == 0)
# Visual proof (by fitting the coefficients of the parabola)
parabola= function(x, p) { p[1]*x^2 + p[2]*x + p[3] }
sse= function(obs, mod) { sum((obs-mod)^2) }
objFun= function(p, x, y_obs) { sse(obs=y_obs, mod=parabola(x, p)) }
iniGuess= c(1,0,0)
opt= optim(par=iniGuess, fn=objFun, x=x, y_obs=y, method="Nelder-Mead")
stopifnot(opt$convergence == 0)
xplot= seq(from=min(x),to=max(x),length.out=100)
plot(xplot, parabola(xplot, opt$par), type="l")  # The fitted parabola
points(x, y)                                     # The 3 points
abline(v=tmp$x.min, lty=3)                       # Location of the minimum
}
\author{
David Kneis \email{david.kneis@uni-potsdam.de}
}
\references{
Adapted from C code at \url{www.mymathlib.com/optimization/
nonlinear/one_dim/parabolic_interpolation.html}
}

