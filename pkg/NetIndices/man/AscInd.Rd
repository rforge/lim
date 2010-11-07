\name{AscInd}
\alias{AscInd}
\title{
  Ascendency network indices
}

\description{
  Calculates measures of system growth and development: Ascendency,
  Overhead and Capacity for several (sub)networks.
}
\usage{
AscInd(Flow = NULL, Tij = t(Flow), Import = NULL, 
       Export = NULL, Dissipation = NULL)
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
    \emph{external} compartmens from where flow enters the network; the
    indices point to the \bold{column} positions in \code{Tij} (and the
    row positions in Flow).
  }
  \item{Export }{vector with either the *indices* or the *names* of
    \emph{external} compartmens to where flow leaves the network; the
    indices point to the \bold{row} positions in \code{Tij} (and the
    column positions in Flow).
  }
  \item{Dissipation }{vector with either the *indices* or the *names* to
    \emph{external} compartments that dissipate flows (e.g. respiration); the
    indices point to the \bold{row} positions in \code{Tij} (and the column
    positions in Flow).
  }
}
\value{
   A matrix with ascendency values (columns) for several subnetworks (rows).
   
The subnetworks (rows) are:

  \item{total}{ network}
  \item{internal}{ network (excluding flows from and to external)}
  \item{import}{ flows}
  \item{export}{ flows; this includes the usuable and unusable flows
    (i.e. +dissipation)}
  \item{dissipation}{ flows}

The ascendency indices comprise (columns:)
  \item{asc}{ the ascendency of the network, a measure of growth and
    development.
  }
  \item{overh}{ the overhead of the network.
  }
  \item{cap}{ the development capacity of the network, an upper bound
    on ascendency.
  }
  \item{ACratio}{ the ratio of ascendency and capacity.
  }
}

\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>,
  Julius Kipyegon Kones<jkones@uonbi.ac.ke>
}
\examples{# The takapoto atoll network
AscInd(Takapoto, Import = "CO2",
       Export = c("CO2", "Sedimentation", "Grazing"),
       Dissipation = "CO2")
       
# Conesprings is the example set 1a from Latham 2006.
as.data.frame(
    AscInd(Tij = Conesprings, Import = "Inflows",
    Export = c("Export","Dissipation"), 
    Dissipation = "Dissipation")
    )
}
\references{
  Latham LG. 2006. Network flow analysis algorithms.
  Ecological Modelling 192: 586-600.

  Ulanowicz RE. 2000. Ascendency: a measure of ecosystem performacne.
  Jorgensen SE, Muller F, editors. Handbook of Ecosystem Theories and
  Management. Lewis Publishers, Boca Raton, p303-315.

  Ulanowicz RE, Norden JS. 1990. Symmetrical overhead in flow networks.
  International Journal of System Science 21: 429-437.

  Kones, J.K., Soetaert, K., van Oevelen, D. and J.Owino (2009).
  Are network indices robust indicators of food web functioning? a
  Monte Carlo approach.  Ecological Modelling, 220, 370-382.
}

\details{
  The mathematical formulation of these indices can be found in the package
  vignette - \emph{vignette("NetIndices")}.

  The PDF can be found in the subdirectory \file{doc} of the NetIndices package.
}
\keyword{utilities}

