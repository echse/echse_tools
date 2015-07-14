
#' Improved step-wise line search algorithm
#'
#' This is a simple \emph{deterministic} algorithm to fit the parameters of a
#' model. It was desinged for cases where 'good' initial estimates of the
#' parameters are available. The algorithm is supposed to be useful in the
#' optimization of computationally expensive, multi-parameter models where
#' classical methods may be too fragile and stochatic methods may consume too
#' much computer time.
#' See the below-mentioned reference for details on the original algorithm. This
#' version includes a number of improvements over the original version.
#'
#' @param p A numeric vector holding the initial values of the parameters.
#'   Elements must be named.
#' @param f The objective function to be minimized. If \code{f} returns a vector,
#'   minimization is done with respect to the first element of the return value. See details below.
#' @param ... Arguments other than \code{p} to be passed to \code{f}. See details below.
#' @param steps A vector with the same length as \code{p} defining absolute step sizes for
#'   variation for the individual parameters. In each iteration, a particular parameter with
#'   index \code{i} is varied between \code{p[i]+steps[i]} and \code{p[i]-steps[i]}
#'   in order to test for a reduction (or rise) of \code{f}.
#' @param suspend The number of iterations for which an element of \code{p} is
#'   excluded from further optimization after being identified as insensitive. The
#'   default is \code{suspend=0}, i.e. optimization is attempted in every iteration
#'   irrespective of a sensitivity. Higher values may speed up the optimization
#'   if (locally) non-sensitive parameters are present.
#' @param lower Lower limits for the parameters. Must be of the same length as \code{p}.
#' @param upper Upper limits for the parameters. Must be of the same length as \code{p}.
#' @param maxiter The maximum number of iterations.
#'
#' @return A list with the following components
#'   \item{opt.p}{The 'best' parameter vector found with the method.}
#'   \item{opt.f}{The value of \code{f} corresponding to the 'best' parameter set.}
#'   \item{opt.iter}{The iteration index corresponding to \code{opt.f}. The index 0
#'     corresponds to the evaluation of \code{f} for the initial values of \code{p}, i.e.
#'     this is not considered as an iteration.}
#'   \item{opt.eval}{Index of the function evaluation corresponding to \code{opt.f}.
#'     The index 1 corresponds to the evaluation of \code{f} for the initial values of \code{p}.}
#'   \item{niter}{The total number of iterations used.}
#'   \item{neval}{The total number of evaluations of \code{f}.}
#'   \item{convergence}{A logical value. Set to \code{FALSE} if the algorithm
#'     failed to detect a \emph{local} minimum (or maximum) of \code{f} within no more than
#'     \code{maxiter} iterations. Otherwise \code{TRUE}.}
#'   \item{trace_p}{A data frame with as many columns as there are elements in
#'     \code{p}. The first row lists the initial values and the subsequent rows
#'     hold the new parameter values after each iteration. The total number
#'     of rows is \code{niter}+1.}
#'   \item{trace_f}{A data frame with as many columns as there are elements in
#'     the return value of \code{f}. The first row lists
#'     the value(s) of \code{f} evaluated for the initial parameter values.
#'     The subsequent rows hold the 'optimum' values of \code{f} obtained in each
#'     iteration (optimality refers to the 1st element of \code{f} only; see below).
#'     The total number of rows is \code{niter}+1.}
#'
#' @note The objective function is called as \code{f(p,...)}. It must return
#'   either a scalar result or a numeric vector, preferably named. In the vector case,
#'   minimization is carried out with respect to the first element only but values of
#'   the additional elements are reported in \code{trace_f}.
#'   Errors occuring within \code{f} may be signaled by a non-finite return value
#'   (first element, if it is a vector). In this case, an appropriate error message
#'   is generated which includes the recently tested parameter values.
#'   Note that proper rounding of the result of \code{f}, e.g. using \code{signif}
#'   may help to avoid problems with minor local minima (or maxima).
#'
#'   The success of this method is sensitive to the chosen step sizes and
#'   the initial estimates for the parameters. See the referenced paper for
#'   potential uses and pitfalls.
#'
#'   Note that the result may also be sensitive to the order of the elements
#'   in vector \code{p}. It might be useful to put sensitive parameters first and
#'   to put less sensitive ones at positions with higher element indices. There
#'   is no established recommendation, however.
#'
#'   In the case of known parameter interactions (or constraints applying to parameter
#'   combinations) it may be useful to introduce artificial parameters.
#'   If, for example, two parameters \eqn{a} and \eqn{b} are known to be corellated
#'   in a linear way, one could substitute \eqn{b} by a new parameter
#'   \emph{c} being defined as \eqn{b= a * c}.
#'
#'   It is recommended to check the
#'   \code{convergence} component of the returned list before using the
#'   estimated parameters. Note that the parameter estimates may be 'bad' even
#'   though convergence is signalled. This may happen, in particular, for
#'   objective functions with many local minima (maxima) and in the case of
#'   badly chosen initial values and/or step sizes.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @references The original algorithm was presented by V. Kuzmin et al. (2008), doi:10.1016/j.jhydrol.2008.02.001
#'
#' @examples
#' # Estimate parameters of a linear model
#' model= function(params, prtors) prtors * params["a"] + params["b"]
#' err= function(x,y) { c(rmse= sqrt(mean((x - y)^2)), bias=mean(x-y)) }
#' objFun= function(params, model, prtors, obsrv) err(model(params,prtors), obsrv)
#' params= c(a=1, b=0)  # True parameters
#' prtors= 1:1000       # Predictors (x-values)
#' # Generate noisy data (y-values) being compatible with the model
#' obsrv= model(params=params, prtors=prtors) +
#'   model(params=params, prtors=prtors) * rnorm(length(prtors), mean=0, sd=1)
#' # Optimize initial guess
#' params_ini= c(a=2, b=50)
#' opt= isls(p=params_ini, f=objFun, model=model, prtors=prtors, obsrv=obsrv,
#'   steps=c(0.01,1), maxiter=1000)
#' # Check result
#' if (!opt$convergence) stop("No convergence within max. number of iterations.")
#' layout(matrix(1:4,ncol=2,byrow=TRUE))
#' plot(prtors,obsrv, col="grey")
#' lines(prtors, model(params=params,prtors=prtors), col="blue")
#' lines(prtors, model(params=params_ini,prtors=prtors), col="red", lty=2)
#' lines(prtors, model(params=opt$opt.p,prtors=prtors), col="red")
#' legend("topleft",bty="n",lty=c(1,2,1),col=c("blue","red","red"),
#'   legend=c("True parameters","Initial guess","SLS estimate"))
#' plot(1:(opt$niter+1),opt$trace_f[,1],xlab="Iteration",ylab=names(opt$trace_f)[1])
#' plot(1:(opt$niter+1),opt$trace_p$a,xlab="Iteration",ylab="Parameter a")
#' plot(1:(opt$niter+1),opt$trace_p$b,xlab="Iteration",ylab="Parameter b")

isls= function(
  p,
  f,
  ...,
  steps,
  suspend=0,
  lower= rep(-Inf,length(p)),
  upper= rep(Inf,length(p)),
  maxiter=10*length(p)^2
) {
  # Helper functions
  minmax= function(x, mini, maxi) { min(max(x,mini),maxi) }
  # Check args
  stopifnot(is.numeric(p), length(p)>=1, !is.null(names(p)), all(names(p)!=""))
  stopifnot(is.function(f))
  stopifnot(is.numeric(steps), length(steps)==length(p))
  stopifnot(is.numeric(suspend), length(suspend)==1, suspend==round(suspend))
  stopifnot(is.numeric(lower), length(lower)==length(p))
  stopifnot(is.numeric(upper), length(upper)==length(p))
  stopifnot(all(lower < upper))
  stopifnot(is.numeric(maxiter), length(maxiter)==1)
  # Initializations
  miniFun= f(p,...)
  neval= 1
  if (!is.finite(miniFun[1])) {
    stop("Objective function returned non-finite result for initial parameter values.")
  }
  miniFun_old= Inf
  niter=0
  converged=TRUE
  stepSigns= rep(1, length(p))    # Initially all parameter values are increased
  suspendUntil= rep(1, length(p)) # Initially, all parameters are modified
  trace_p= data.frame(as.list(p))
  trace_f= data.frame(as.list(miniFun))
  opt.iter=0
  opt.eval=1
  # Iterations to minimize f
  while (miniFun[1] < miniFun_old[1]) {
    niter= niter + 1
    if (niter > maxiter) {
      converged=FALSE
      break
    }
    miniFun_old= miniFun
    # Loop over parameters
    for (i in 1:length(p)) {
      if (niter >= suspendUntil[i]) {
        p_new= p
        sensitive= FALSE
        triedDirections=0
        while (triedDirections < 2) {
          triedDirections= triedDirections + 1
          # Modify parameter and evaluate objective function
          p_new[i]= minmax(p[i] + stepSigns[i]*steps[i], mini=lower[i], maxi=upper[i])
          miniFun_new= f(p_new,...)
          neval= neval + 1
          if (!is.finite(miniFun_new[1])) {
            stop(paste("Objective function returned non-finite result in",
              " iteration ",niter,". The tested parameter values were: ",
              paste(names(p),p_new,sep="=",collapse=", "),".",sep=""))
          }
          # Value of f decreased: Accept new estimate & continue with next parameter
          if (miniFun_new[1] < miniFun[1]) {
            p= p_new
            miniFun= miniFun_new
            sensitive= TRUE
            opt.iter=niter
            opt.eval=neval
            break
          # Value of f increased or remained constant: Change search direction
          } else {
            stepSigns[i]= stepSigns[i] * (-1)
            if (miniFun_new[1] > miniFun[1])
              sensitive= TRUE
          }
        }
        if (!sensitive) {
          suspendUntil[i]= niter + 1 + suspend
        }
      }
    } # End of loop over parameters in n^th iteration
    trace_p= rbind(trace_p,as.list(p))
    trace_f= rbind(trace_f,as.list(miniFun))
  } # End of iteration
  return(list(opt.p=p, opt.f=miniFun[1],
    opt.iter=opt.iter, opt.eval=opt.eval,
    niter=niter, neval=neval, convergence=converged,
    trace_p=trace_p, trace_f=trace_f))
}


#' Location of the minimum of a parabola fitted through 3 points 
#'
#' Given 3 points (x.left, y.left), (x.center, y.center), (x.right, y.right),
#' the function returns the x-value where a parabola ax^2+bx+c fitted through
#' these 3 points has a minimum. This is one of the core functions used in
#' optimization methods based on successive parabolic interpolation (like
#' \code{\link{optimize}}, for example).
#'
#' @param x.left x-value of leftmost point
#' @param y.left y at x.left
#' @param x.center x-value of center point. Must be \eqn{> x.left}. 
#' @param y.center y at x.center. It should be
#'   \eqn{y.center <= min(y.left, y.right)} for the fitted parabola to have
#'   a minimum in the range defined by \code{x.left} and \code{x.right}.
#'   
#' @param x.right x-value of rightmost point. Must be \eqn{> x.center}.
#' @param y.right y at x.right
#'
#' @return A list with two components
#'   \item{status}{An integer value. A value of 0 indicates success. Other
#'      values indicate a problem (1: collinearity, 2: local maximum).}
#'   \item{x.min}{The x-value where the parabola has its minimum if
#'     \eqn{status == 0}. Otherwise, the x-value corresponding to the minimum
#'     of \code{y.left}, \code{y.center}, and \code{y.right}.}
#' 
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @references Adapted from C code at \url{www.mymathlib.com/optimization/
#'   nonlinear/one_dim/parabolic_interpolation.html}
#'
#' @examples
#' # Determine minimum of fitted parabola
#' x= c(-2, 0.5, 5)
#' y= c(2, 1, 10)
#' tmp= parabolicMin(x[1], y[1], x[2], y[2], x[3], y[3])
#' print(tmp)
#' stopifnot(tmp$status == 0)
#' # Visual proof (by fitting the coefficients of the parabola)
#' parabola= function(x, p) { p[1]*x^2 + p[2]*x + p[3] }
#' sse= function(obs, mod) { sum((obs-mod)^2) }
#' objFun= function(p, x, y_obs) { sse(obs=y_obs, mod=parabola(x, p)) }
#' iniGuess= c(1,0,0)
#' opt= optim(par=iniGuess, fn=objFun, x=x, y_obs=y, method="Nelder-Mead")
#' stopifnot(opt$convergence == 0)
#' xplot= seq(from=min(x),to=max(x),length.out=100)
#' plot(xplot, parabola(xplot, opt$par), type="l")  # The fitted parabola
#' points(x, y)                                     # The 3 points
#' abline(v=tmp$x.min, lty=3)                       # Location of the minimum

parabolicMin= function(x.left, y.left, x.center, y.center, x.right, y.right) {
  stopifnot(x.left < x.center, x.center < x.right)
  d1 = (x.right - x.center) * (y.left - y.center)
  d2 = (x.center - x.left) * (y.right - y.center)
  denominator = d1 + d2
  numerator = (x.right + x.center) * d1 + (x.center + x.left) * d2
  if (denominator == 0.0) {       # the points are collinear
    return(list(status=1,
      x.min= c(x.left,x.center,x.right)[which.min(c(y.left,y.center,y.right))]))
  } else if (denominator < 0.0) { # local maximum
    return(list(status=2,
      x.min= c(x.left,x.center,x.right)[which.min(c(y.left,y.center,y.right))]))
  } else {                        # the desired 'normal' case
    return(list(status= 0, x.min=(0.5 * numerator / denominator)))
  }
}

################################################################################

#' Approximate 1D minimization with limited number of function evaluations
#'
#' This function tries to minimize the value of a 1-dimensional function using
#' a limited number of function evaluations. It is intended for operational
#' parameter updating of computationally expensive simulation models.
#'
#' @param f Objective function to be minimized with interface \code{f(x,...)}.
#'   The 1st argument of this function must be a scalar numeric model parameter.
#'   The function must return a scalar numeric result (RMSE, for example).
#' @param base Base value of the parameter with respect to which the objective
#'   function is minimized. It is used as the initial value, thus \code{base}
#'   should be chosen close to a suspected minimum of the objective function.
#' @param range1 Vector of 2 elements defining the search range for the 1st try.
#'   The minimum of the objective function is searched in range
#'   \eqn{range1[1] ... range1[2]}. Must fulfill the condition
#'   \eqn{range1[1] < base < range1[2]}. If no minimum is found, the search
#'   range is extended (see argument \code{range2}).
#' @param range2 Vector of 2 elements defining the search range for the 2nd try.
#'   The minimum of the objective function is searched in range
#'   \eqn{range2[1] ... range2[2]} if the 1st try (see argument \code{range1})
#'   failed. Must fulfill
#'   \eqn{range2[1] < range1[1] < base < range1[2] < range2[2]}.
#' @param ... Additional arguments to be passed to the objective function. 
#'
#' @return A list with two components.
#'   \item{min.x}{The parameter value at the approximate minimum.}
#'   \item{min.fx}{The value of the function at \code{min.x}.}
#'
#' @note The algorithm tries to approximate a minimum of \code{f} based on a
#'   single application of parabolic interpolation. If the fitted parabola has
#'   no minimum, the parameter corresponding to the smallest \emph{tested} value
#'   of \code{f} is returned.
#'
#' @seealso Use \code{\link{optimize}} or another appropriate method for
#'   minimization problems that require an accurate solution.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @examples
#' fun= function(x, ...) {
#'   print(paste("Evaluated at",x,"using additional args:",...))
#'   return((x-0.5)^2)
#' }
#' res= fivePointMin(f=fun, base=0, range1=c(-1,1), range2=c(-2,2), "None")
#' print(res)

fivePointMin= function(f, base=1, range1=c(0.5,2), range2=c(0.1,10), ...) {
  # Check args
  stopifnot(length(range1)==2, range1[1] < base, range1[2] > base)
  stopifnot(length(range2)==2, range2[1] < range1[1], range2[2] > range1[2])
  # Try the 3 basic parameters
  tab= data.frame(p= c(range1[1], base, range1[2]), e= rep(NA,3))
  for (i in 1:nrow(tab)) {
    tab$e[i]= f(tab$p[i], ...)
  }
  # Monotonous function with minimum at left boundary
  if ((tab$e[1] < tab$e[2]) && (tab$e[2] < tab$e[3])) {
    tab[3,]= tab[2,]
    tab[2,]= tab[1,]
    tab$p[1]= range2[1]
    tab$e[1]= f(tab$p[1], ...)   # New run with range extended to left
    if (tab$e[1] > tab$e[2]) {  # Determine minimum in extended range
      p_opt= parabolicMin(x.left=tab$p[1], y.left=tab$e[1], x.center=tab$p[2],
        y.center=tab$e[2], x.right=tab$p[3], y.right=tab$e[3])$x.min
      e= f(p_opt, ...)
    } else {# Minimum not in extended range --> use left boundary of ext. range
      p_opt= tab$p[1]
      e=     tab$e[1]
    }
  # Monotonous function with minimum at right boundary
  } else if ((tab$e[1] > tab$e[2]) && (tab$e[2] > tab$e[3])) {
    tab[1,]= tab[2,]
    tab[2,]= tab[3,]
    tab$p[3]= range2[2]
    tab$e[3]= f(tab$p[3], ...)  # New run with range extended to right
    if (tab$e[3] > tab$e[2]) {  # Determine minimum in extended range
      p_opt= parabolicMin(x.left=tab$p[1], y.left=tab$e[1], x.center=tab$p[2],
        y.center=tab$e[2], x.right=tab$p[3], y.right=tab$e[3])$x.min
      e= f(p_opt, ...)
    } else {# Minimum not in extended range --> use right boundary of ext. range
      p_opt= tab$p[3]
      e=     tab$e[3]
    }
  # Maximum at center or constant error (insensitive model)
  } else if (((tab$e[1] < tab$e[2]) && (tab$e[2] > tab$e[3])) ||
             ((tab$e[1] == tab$e[2]) && (tab$e[2] == tab$e[3]))) {
    tab= tab[sort.list(tab$e),] # Put smallest error in 1st row
    tab$p[2]= range2[2]
    tab$e[2]= f(tab$p[2], ...)   # Overwrite pos 2 with range extended to right
    tab$p[3]= range2[1]
    tab$e[3]= f(tab$p[3], ...)   # Overwrite pos 3 with range extended to left
    p_opt= tab$p[which.min(tab$e)] # Just select factor with smallest error
    e= f(p_opt, ...)
  # The 'good' case of a minimum at center (at least 1 boundary value is higher)
  } else {
    p_opt= parabolicMin(x.left=tab$p[1], y.left=tab$e[1], x.center=tab$p[2],
      y.center=tab$e[2], x.right=tab$p[3], y.right=tab$e[3])$x.min
    e= f(p_opt, ...)
  }
  return(list(min.x=p_opt, min.fx=e))
}


