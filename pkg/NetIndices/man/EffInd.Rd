\name{EffInd}
\alias{EffInd}
\title{
  Effective measures (or roles) suite: weighted measures for networks
}
\description{
  Calculates effective connectivity, effective flows, effective nodes and
  effective roles of a network.
}
\usage{
EffInd(Flow = NULL, Tij = t(Flow), 
       Import = NULL, Export = NULL)}
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
    the indices point to the \bold{column} positions in \code{Tij} (and
    the row positions in Flow).
  }
  \item{Export }{vector with either the *indices* or the *names* of
    \emph{external} compartmens to where flow leaves the network; the
    indices point to the \bold{row} positions in \code{Tij} (and the
    column positions in Flow).
  }
}
\value{
  a list with the following items:
  \item{CZ }{Effective connectance}
  \item{FZ }{Effective Flows}
  \item{NZ }{Effective nodes}
  \item{RZ }{Effective roles}
}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>,
  Julius Kipyegon Kones<jkones@uonbi.ac.ke>
}

\examples{# The takapoto atoll network
EffInd(Takapoto, Import = "CO2",
       Export = c("CO2", "Sedimentation", "Grazing"))

# Conesprings is the example set 1a from Latham 2006.
as.data.frame(
  EffInd(Tij = Conesprings, Import = "Inflows",
         Export = c("Export", "Dissipation"))
)
}

\references{
  Latham LG. 2006. Network flow analysis algorithms. Ecological
  Modelling 192: 586-600.

  Zorach and Ulanowicz, 2003. Quantifying the complexity of flow
  networks: how many roles are there?. Complexity 8,68-76.

  Kones, J.K., Soetaert, K., van Oevelen, D. and J.Owino (2009).
  Are network indices robust indicators of food web functioning? a
  Monte Carlo approach.  Ecological Modelling, 220, 370-382.
}

\details{
  The mathematical formulation of these indices can be found in the package
  vignette - \emph{vignette("NetIndices")}.


  The PDF can be found in the subdirectory \file{doc} of the NetIndices
  package.
}
\keyword{utilities}

