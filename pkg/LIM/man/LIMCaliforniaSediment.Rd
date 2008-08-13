\name{LIMCaliforniaSediment}
\docType{data}
\alias{LIMCaliforniaSediment}
\title{Linear inverse model specification for the Santa Monica Basin sediment food web}
\description{Linear inverse model specification for the Santa Monica Basin (California) sediment food web
\cr as in Eldridge and Jackson (1993). 

The Santa Monica Basin is a hypoxic-anoxic basin located near California.
\cr The model contains both chemical and biological species. 
\cr The foodweb comprises 7 functional compartments and five external compartments,
connected with 32 flows.  Units of the flows are mg C/m2/day

The linear inverse model LIMCaliforniaSediment is generated from the file 
\cr \file{CaliforniaSediment.input} which can be found in subdirectory \file{web} of the packages directory

In this subdirectory you will find many foodweb example input files
\cr These files can be read using \code{\link{Read}} and their output processed by 
\code{\link{Setup}} which will produce a linear inverse problem specification similar to LIMCaliforniaSediment}
\usage{data(LIMCaliforniaSediment)}
\format{a list of matrices, vectors, names and values that specify the linear inverse model problem.
\cr see the return value of \code{\link{Setup}} for more information about this list}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>
   Dick van Oevelen<d.vanoevelen@nioo.knaw.nl> and Henrik Andersson}
\examples{CaliforniaSediment <- Flowmatrix(LIMCaliforniaSediment)
plotweb(CaliforniaSediment,main="Santa Monica Basin Benthic web",
        sub="mgN/m2/day",lab.size=0.8)
xr <- Xranges(LIMCaliforniaSediment)
min<- xr[,1]
max<- xr[,2]
xv <- Lsei(LIMCaliforniaSediment)$X
i1 <- 1:(nrow(xr)/2)
i2 <- (nrow(xr)/2+1):nrow(xr)
Plotranges(min[i1],max[i1],xv[i1],lab.cex=0.7,
        main="Santa Monica Basin Benthic web, Flowranges - part1")
Plotranges(min[i2],max[i2],xv[i2],lab.cex=0.7,
        main="Santa Monica Basin Benthic web, Flowranges - part2")}
\references{Eldridge, P.M., Jackson, G.A., 1993. Benthic trophic dynamics in California coastal basin and
continental slope communities inferred using inverse analysis. Marine Ecology Progress Series 99, 115-135.}
\seealso{\code{\link{LIMTakapoto}}, \code{\link{LIMRigaSummer}} and many others}
\keyword{datasets}

