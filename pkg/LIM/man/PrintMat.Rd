\name{PrintMat}
\alias{PrintMat}
\title{Writes linear inverse matrices to the screen}
\description{Prints the linear inverse problem:

Writes to the screen:
\item inverse matrices and vectors \bf{A}, \bf{b} of the equalities \deqn{Ax=b}
\item inverse matrices and vectors \bf{G}, \bf{h} of the inequalities \deqn{Gx>=h}
\item if present, also writes the cost/profit function}
\usage{
PrintMat(lim)
}
\arguments{
  \item{lim }{linear inverse model list, as returned by \code{Setup}}
}
\value{
returns nothing.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{

   PrintMat(LIMBlending)
           }
\keyword{IO}
\keyword{utilities}

