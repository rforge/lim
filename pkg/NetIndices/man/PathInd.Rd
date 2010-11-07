\name{PathInd}
\alias{PathInd}
\title{
  Pathway analysis
}
\description{
  Calculates the direct and indirect pathways in a network, i.e.
  the total system cycled throughflow, Finn's cycling index and average
  pathlength,...
  
  Based on Finn(1980) (and not Finn (1976))
}
\usage{
PathInd(Flow = NULL, Tij = t(Flow), Import = NULL, Export = NULL)
}

\arguments{
  \item{Flow }{network matrix with Flow[i,j] the flow from i (row) to j
    (column); component positions in rows and columns must be the same;
    if present, rownames or columnnames denote the compartment names.
  }
  \item{Tij }{network matrix where connectance is from column j to row i;
    component positions in rows and columns must be the same ;
    if present, rownames or columnnames denote the compartment names.
  }
  \item{Import }{vector with either the *indices* or the *names* of
    \emph{external} compartmens from where flow enters the network;
    the indices point to the \bold{column} positions in \code{Tij}
    (and the row positions in Flow).
  }
  \item{Export }{vector with either the *indices* or the *names* of
    \emph{external} compartmens to where flow leaves the network; the
    indices point to the \bold{row} positions in \code{Tij} (and the
    column positions in Flow).
  }
}
\value{
   A list with the following items:
  \item{TSTC }{total system cycled throughflow.
  }
  \item{TSTS }{non-cycled throughflow.
  }
  \item{FCI }{Finn's cycling index (1980).
  }
  \item{FCIb }{revised Finn's cycling index, sensu Allesina and Ulanowicz, 2004.
  }
  \item{APL }{average pathlength, also known as Network Aggradation (Sum of
    APLc and APLs in Latham 2006).
  }
}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>,
  Julius Kipyegon Kones<jkones@uonbi.ac.ke>
}
\examples{# The takapoto atoll network
PathInd(Takapoto, Import = "CO2",
        Export = c("CO2", "Sedimentation", "Grazing"))

# Conesprings is the example set 1a from Latham 2006.
as.data.frame(
  PathInd(Tij = Conesprings, Import = "Inflows",
          Export = c("Export", "Dissipation"))
)
}

\references{
  Finn JT. 1980. Flow analysis of models of the Hubbard Brook ecosystem.
  Ecology 61: 562-571.

  Patten BC, Higashi M. 1984. Modified cycling index for ecological
  applications. Ecological Modelling 25: 69-83.

  Patten BC, Bosserman RW, Finn JT, Cale WG. 1976. Propagation of cause in
  ecosystems. Patten BC, editor. Systems Analysis and Simulation in Ecology,
  vol. 4. Academic Press, New York. p457-579.

  Allesina and Ulanowicz, 2004. Cycling in ecological netowrks: Finn's index
  revisited. Computational Biology and Chemistry 28, 227-233.

  Kones, J.K., Soetaert, K., van Oevelen, D. and J.Owino (2009).
  Are network indices robust indicators of food web functioning? a
  Monte Carlo approach.  Ecological Modelling, 220, 370-382.
}

\details{
  The mathematical formulation of these indices can be found in the
  package vignette - \emph{vignette("NetIndices")}.


  The PDF can be found in the subdirectory \file{doc} of the NetIndices package.
}
\keyword{utilities}

