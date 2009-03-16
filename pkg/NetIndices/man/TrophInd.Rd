\name{TrophInd}
\alias{TrophInd}
\title{The trophic level and omnivory index}
\description{
  Calculates the trophic level and omnivory index of each component of
  a food web.
}

\usage{
TrophInd(Flow=NULL, Tij=t(Flow), 
  Import=NULL, Export=NULL, Dead=NULL)}

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
    \emph{external} compartmens to where flow leaves the network; the indices
    point to the \bold{row} positions in \code{Tij} (and the column positions
    in Flow).
  }
  \item{Dead }{vector with either the *indices* or the *names* of dead matter;
    the indices point to \bold{row} positions in \code{Tij}; the trophic level
    of these components is assumed to be 1.
  }
}
\value{
  a 2-columned data.frame with, for each compartment of the network the
  following:
  \item{TL }{the trophic level of a compartment, defined as 1 + the weighted
    average of the trophic levels of its food items.
  }
  \item{OI }{the omnivory index, the variety in the trophic levels of a
    consumer's food.
  }

}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>,
  Julius Kipyegon Kones<jkones@uonbi.ac.ke>}
\examples{# The takapoto atoll network
# First trophic level without assuming that TL of detritus and DOC is 1
TrophInd(Takapoto,Import ="CO2",
        Export=c("CO2","Sedimentation","Grazing"))

# Now imposing TL=1 for detritus and DOC 
TrophInd(Takapoto,Import ="CO2",
        Export=c("CO2","Sedimentation","Grazing"),
        Dead=c("Detritus","DOC"))}
\references{
  Christensen V, Pauly D. 1992. ECOPATH II - a software for balancing
  steady-state ecosystem models and calculating network characteristics.
  Ecological Modelling 61: 169-185.

  Lindeman RL. 1942. The trophic dynamic aspect of ecology.
  Ecology 23: 399-418.
}

\details{
  Primary producers, defined as those compartments that do not receive matter
  from another internal compartment, will be assigned a trophic level of 1.

  In many trophic level calculations, it is also assumed that TL of detritus,
  dissolved organic matter and other inert material (i.e. that does not feed)
  is also = 1.

  If this is desired, these compartments have to be designated as "Dead"
  (i.e. \code{Dead} should contain an index to row positions in \code{Tij}
  of these compartments. 

  If not specified as "Dead", these compartments will have a TL > 1 and
  consequently the TL of other compartments will be higher too.

  The mathematical formulation of these indices can be found in the package
  vignette - \emph{vignette("NetIndices")}.


  The PDF can be found in the subdirectory \file{doc} of the NetIndices package.
}

\keyword{utilities}

