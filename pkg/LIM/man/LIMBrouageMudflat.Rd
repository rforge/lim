\name{LIMBrouageMudflat}
\docType{data}
\alias{LIMBrouageMudflat}
\title{Linear inverse model specification for the Intertidal mudflat food web on the Atlantic coast of France}
\description{Linear inverse model specification for the Intertidal mudflat food web on the Atlantic coast of France
\cr as in Leguerrier et al., 2003. 

The foodweb comprises 16 functional compartments and 3 external compartments,
connected with 95 flows.  Units of the flows are g C/m2/year

The linear inverse model LIMBrouageMudflat is generated from the file BrouageMudflat.input 
which can be found in subdirectory \file{web} of the packages directory

In this subdirectory you will find many foodweb example input files
\cr These files can be read using \code{\link{Read}} and their output processed by 
\code{\link{Setup}} which will produce a linear inverse problem specification similar to LIMBrouageMudflat}
\usage{data(LIMBrouageMudflat)}
\format{a list of matrices, vectors, names and values that specify the linear inverse model problem.
\cr see the return value of \code{\link{Setup}} for more information about this list}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>
   Dick van Oevelen<d.vanoevelen@nioo.knaw.nl>}
\examples{Brouage <- Flowmatrix(LIMBrouageMudflat)
plotweb(Brouage,main="Brouage mudflat food web",sub="gC/m2/yr")
Plotranges(LIMBrouageMudflat,lab.cex=0.7,
               main="Brouage mudflat, Flowranges")
Plotranges(LIMBrouageMudflat,type="V",lab.cex=0.7,
               main="Brouage mudflat, Variable ranges")}
\references{Leguerrier, D., Niquil, N., Boileau, N., Rzeznik, J., Sauriau, P.G., Le Moine, O., Bacher, C., 2003. 
Numerical analysis of the food web of an intertidal mudflat ecosystem on the Atlantic coast of France. Marine Ecology Progress Series 246, 17-37.}
\seealso{\code{\link{LIMTakapoto}}, \code{\link{LIMRigaSummer}} and many others}
\keyword{datasets}
