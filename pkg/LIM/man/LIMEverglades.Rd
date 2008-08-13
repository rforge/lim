\name{LIMEverglades}
\docType{data}
\alias{LIMEverglades}
\title{Linear inverse model specification for the herpetological food web of the Everglades}
\description{Linear inverse model specification for the herpetological wet prairie example from the everglades. 
\cr as described in Diffendorfer et al., 2001
 

The everglades are a freshwater wetland in Florida, USA.
\cr The model contains 9 functional compartments and 3 external compartments,
connected with 402 flows.  Units of the flows are gram wet weight / Ha / year

The linear inverse model LIMEverglades is generated from the file Everglades.input 
which can be found in subdirectory \file{web} of the packages directory

In this subdirectory you will find many foodweb example input files
\cr These files can be read using \code{\link{Read}} and their output processed by 
\code{\link{Setup}} which will produce a linear inverse problem specification similar to LIMEverglades}
\usage{data(LIMEverglades)}
\format{a list of matrices, vectors, names and values that specify the linear inverse model problem.
\cr see the return value of \code{\link{Setup}} for more information about this list}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>
   Dick van Oevelen<d.vanoevelen@nioo.knaw.nl>}
\examples{Flows <- Lsei(LIMEverglades)
Flowranges <- Xranges(LIMEverglades)
Everglades <- Flowmatrix(LIMEverglades)
plotweb(Everglades,main="Everglades Herpetological Wet Prairie model",
        sub="g WW/Ha/Yr",lab.size=0.8)}
\references{Diffendorfer, J.E., Richards, P.M., Dalrymple, G.H., DeAngelis, D.L., 2001. 
Applying Linear Programming to estimate fluxes in ecosystems or food webs: 
an example from the herpetological assemblage of the freshwater Everglades. Ecol. Model. 144, 99-120.}
\seealso{\code{\link{LIMTakapoto}}, \code{\link{LIMRigaSummer}} and many others}
\keyword{datasets}

