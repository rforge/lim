\name{GenInd}
\alias{GenInd}
\title{
  General network indices.
}
\description{
  Calculates general network indices such as system throughputs, link
  density, connectance,... of a network.
}

\usage{
GenInd(Flow = NULL, Tij = t(Flow), 
       Import = NULL, Export = NULL, tol = 0)}
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
    \emph{external} compartmens to where flow leaves the network; the indices
    point to the \bold{row} positions in \code{Tij} (and the column positions
    in Flow).
  }
  \item{tol }{flows that are smaller or equal to tol are assumed to be absent.
  }
}

\value{
  A list that contains:
  \item{N }{number of compartments, excluding the externals.
  }
  \item{T.. }{total System Throughput.
  }
  \item{TST }{total System Throughflow.
  }
  \item{Lint }{number of Internal links.
  }
  \item{Ltot }{total number of links.
  }
  \item{LD }{link Density.
  }
  \item{C }{connectance (internal).
  }
  \item{Tijbar }{average Link Weight.
  }
  \item{TSTbar }{average Compartment Throughflow .
  }
  \item{Cbar }{compartmentalization, [0,1], the degree of connectedness of
    subsystems within a network.
  }
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>,
  Julius Kipyegon Kones<jkones@uonbi.ac.ke>
}

\examples{# The takapoto atoll network
(GI<- GenInd(Takapoto, Import = "CO2",
             Export = c("CO2", "Sedimentation", "Grazing")))
as.data.frame(GI)  

# a simple system
as.data.frame(GenInd(diag(5)))

# Conesprings is the example set 1a from Latham 2006.
as.data.frame(
  GenInd(Tij = Conesprings, Import = "Inflows",
         Export = c("Export", "Dissipation"))
)
}

\references{

  Latham LG. 2006. Network flow analysis algorithms. Ecological Modelling
  192: 586-600.

  Hirata H, Ulanowicz RE. 1984. Informational theoretical analysis of
  ecological networks. International journal of systems science
  15 (3): 261-270

  Pimm SL, Lawton JH. 1980. Are food webs divided into compartments?
  Journal of Animal Ecology 49: 879-898.

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

