\name{LIMRigaAutumn}
\docType{data}
\alias{LIMRigaAutumn}
\title{Linear inverse model specification for the Gulf of Riga *autumn* planktonic food web}
\description{Linear inverse model specification for the Gulf of Riga planktonic food web in *autumn*
\cr as in Donali et al. (1999). 

The Gulf of Riga is a highly eutrophic system in the Baltic Sea.
\cr The foodweb comprises 7 functional compartments and two external compartments,
connected with 26 flows.  Units of the flows are mg C/m3/day

The linear inverse model LIMRigaAutumn is generated from the file RigaAutumn.input 
which can be found in subdirectory \file{web} of the packages directory

In this subdirectory you will find many foodweb example input files
\cr These files can be read using \code{\link{Read}} and their output processed by 
\code{\link{Setup}} which will produce a linear inverse problem specification similar to LIMRigaAutumn}
\usage{data(LIMRigaAutumn)}
\format{a list of matrices, vectors, names and values that specify the linear inverse model problem.
\cr see the return value of \code{\link{Setup}} for more information about this list}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>
   Dick van Oevelen<d.vanoevelen@nioo.knaw.nl>}
\examples{rigaAutumn <- Flowmatrix(LIMRigaAutumn)
plotweb(rigaAutumn,main="Gulf of Riga planktonic food web, autumn",
   sub="mgC/m3/day")
Plotranges(LIMRigaAutumn,lab.cex=0.7,
   main="Gulf of Riga planktonic food web, autumn, Flowranges")
Plotranges(LIMRigaAutumn,type="V",lab.cex=0.7,
   main="Gulf of Riga planktonic food web, autumn, Variable ranges")}
\references{Donali, E., Olli, K., Heiskanen, A.S., Andersen, T., 1999. Carbon flow 
patterns in the planktonic food web of the Gulf of Riga, the Baltic Sea: 
a reconstruction by the inverse method. Journal of Marine Systems 23, 251..268.}
\seealso{\code{\link{LIMTakapoto}}, \code{\link{LIMRigaSummer}}, \code{\link{LIMRigaSpring}} and many others}
\keyword{datasets}

