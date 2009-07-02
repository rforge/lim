\name{LIMCoralRockall}
\docType{data}
\alias{LIMCoralRockall}
\title{
  Linear inverse model specification for the Deep-water Coral food web
  at Rockall Bank
}
\description{
  Linear inverse model specification for the deep-water coral ecosystem
  at Rockall Bank, North-East Atlantic. See van Oevelen et al. (2009)
}

\usage{
data(LIMCoralRockall)
}

\format{
  a list of matrices, vectors, names and values that specify the linear
  inverse model problem.

  see the return value of \code{\link{Setup}} for more information about
  this list

  A more complete description of this structures is in vignette("LIM")
}

\author{
  Dick van Oevelen<d.vanoevelen@nioo.knaw.nl>

  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\examples{
Coral <- Flowmatrix(LIMCoralRockall)

plotweb(Coral,main="Deep Water Coral Foodweb, Rockall Bank",
        sub="mmolC/m2/day",lab.size=0.8)


\dontrun{
xr <- LIMCoralRockall$NUnknowns
i1 <- 1:(xr/2)
i2 <- (xr/2+1):xr

pm <- par(mfrow=c(1,1))
Simplest <- Ldei(LIMCoralRockall)$X
Ranges <- Xranges(LIMCoralRockall)
Plotranges(Ranges[i1,1],Ranges[i1,2],Simplest[i1],lab.cex=0.7,
        main="Deep Water Coral - ranges part 1")

Plotranges(Ranges[i2,1],Ranges[i2,2],Simplest[i2],lab.cex=0.7,
        main="Deep Water Coral - ranges part 2")
par(mfrow=pm)
}
}

\references{
  van Oevelen, Dick, Gerard Duineveld, Marc Lavaleye, Furu Mienis,
  Karline Soetaert, and Carlo H. R. Heip, 2009.

  The cold-water coral community as hotspot of carbon cycling on
  continental margins: A food web analysis from Rockall Bank
  (northeast Atlantic). Limnology and Oceanography, in press.

}

\seealso{
  \code{\link{LIMTakapoto}}, \code{\link{LIMRigaSummer}} and many others
}

\keyword{datasets}

