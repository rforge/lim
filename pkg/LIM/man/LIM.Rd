\name{LIM-package}
\alias{LIM-package}
\alias{LIM}
\docType{package}
\title{
  Linear Inverse Model examples and solution methods
}
\description{
  functions that read and solve linear inverse problems (food web problems,
  linear programming problems, flux balance analysis).
  
  These problems find solutions to linear or quadratic functions:
  min or max (f(x)), where \eqn{f(x) = ||Ax-b||^2} or \eqn{f(x) = sum(ai*xi)}

  subject to equality constraints \eqn{Ex=f} and inequality
  constraints \eqn{Gx>=h}.

  Uses package limSolve.
}
\details{
  \tabular{ll}{
    Package: \tab LIM\cr
    Type: \tab Package\cr
    Version: \tab 1.4\cr
    Date: \tab 2009-07-01\cr
    License: \tab GNU Public License 2 or above\cr
  }

  The model problem is formulated in text files in a way that is natural and
  comprehensible. Functions in \code{LIM} then converts this input into the
  required linear equality and inequality conditions, which can be solved
  either by least squares or by linear programming techniques. By
  letting an algorithm formulate the mathematics, it is
  simple to reformulate the model in case a parameter value changes,
  or a component is added or removed.

  Three different types of problems are supported: flow networks,
  reaction networks and other (operations research) problems. The
  first two cases are based on mass balances of the components.

  The package includes many examples
}
\author{
  Karline Soetaert (Maintainer),
  
  Dick van Oevelen
}
%\references{
%% to be filled
%}
\seealso{
 \code{\link{Read}}, \code{\link{Setup}} for reading files and creating the model

 \code{\link{Ldei}}, \code{\link{Lsei}},\code{\link{Linp}},
 \code{\link{Flowmatrix}}, \code{\link{Plotranges}}, \code{\link{Variables}},
 \code{\link{Varranges}}, \code{\link{Xranges}}, \code{\link{Xsample}}.
}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(Lsei)
example(LIMRigaSpring)
example(Ldei)
example(Xsample)
example(Varranges)

## run demos
demo("LIMexamples")

## open the directory with R sourcecode examples
browseURL(paste(system.file(package="LIM"), "/examples/Foodweb", sep=""))
browseURL(paste(system.file(package="LIM"), "/examples/LinearProg", sep=""))
browseURL(paste(system.file(package="LIM"), "/examples/Reactions", sep=""))

## show package vignette with tutorial about how to create input files
vignette("LIM")
browseURL(paste(system.file(package="LIM"), "/doc", sep=""))
}
}

\keyword{ package }