\name{NetIndices-package}
\alias{NetIndices-package}
\alias{NetIndices}
\docType{package}
\title{
  Estimates network indices, including trophic structure of foodwebs
}
\description{
  Given a network (e.g. a food web), estimates several network indices.

  These include:
  \itemize{
    \item ascendency network indices,
    \item direct and indirect dependencies,
    \item effective measures,
    \item environ network indices,
    \item general network indices,
    \item pathway analysis,
    \item network uncertainty indices and constraint efficiencies
    \item the trophic level and omnivory indices of food webs.
  }
}

\details{
  \tabular{ll}{
    Package: \tab NetIndices\cr
    Type: \tab Package\cr
    Version: \tab 1.4.1\cr
    Date: \tab 2010-11-01\cr
    License: \tab GNU Public License 2 or above\cr
  }
}

\author{
  Karline Soetaert (Maintainer), Netherlands Institute of Ecology
  
  Julius Kipyegon Kones, University of Nairobi
}

\references{
  Kones, J.K., Soetaert, K., van Oevelen, D. and J.Owino (2009).
  Are network indices robust indicators of food web functioning? a
  Monte Carlo approach.  Ecological Modelling, 220, 370-382.
  
}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(AscInd)
example(TrophInd)
example(Takapoto)

## open the directory with script used to write the Kones et al. (2009) article
browseURL(paste(system.file(package="NetIndices"), "/EcologicalModelling", sep=""))

## open the directory with documents
browseURL(paste(system.file(package="NetIndices"), "/doc", sep=""))

## the vignette
vignette("NetIndices")
}
}

\keyword{ package }