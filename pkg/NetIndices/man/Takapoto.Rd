\name{Takapoto}
\docType{data}
\alias{Takapoto}
\title{Takapoto atoll planktonic food web}
\description{
  Carbon flux matrix of the Takapoto atoll planktonic food web

  as reconstructed by inverse modelling by Niquil et al. (1998).

  The Takapoto Atoll lagoon is located in the French Polynesia of the South
  Pacific

  The food web comprises 7 functional compartments:
  \itemize{
    \item Phytoplankton
    \item Bacteria
    \item Protozoa
    \item Microzooplankton
    \item Mesozooplankton
    \item Detritus
    \item Dissolved organic carbon (DOC)
  }
  one external source:
  \itemize{
    \item CO2
  }
  and three external sinks:
  \itemize{
    \item CO2
    \item Sedimentation
    \item Grazing
  }
  These compartments are connected with 32 flows.

  Units of the flows are mg C/m2/day
}

\usage{Takapoto}

\format{
  matrix with flow values, where element ij denotes flow from compartment i to j

  rownames and columnames are the components.
}

\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\examples{UncInd(Takapoto, Import = "CO2",
       Export = c("CO2", "Sedimentation", "Grazing"))
}

\references{
  Niquil, N., Jackson, G.A., Legendre, L., Delesalle, B., 1998. Inverse
  model analysis of the planktonic food web of Takapoto Atoll
  (French Polynesia). Marine Ecology Progress Series 165, 17..29.
}

\keyword{datasets}

