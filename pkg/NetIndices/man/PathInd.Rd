\name{PathInd}
\alias{PathInd}
\title{Pathway analysis}
\description{Calculates the direct and indirect pathways in a network, i.e.
  the total system cycled throughflow, Finn's cycling index and average pathlength,...
  
  Base on Finn(1980) (and not Finn (1976))}
\usage{PathInd(Flow=NULL, Tij=t(Flow), Import=NULL, Export=NULL)}
\arguments{
  \item{Flow }{Network matrix with Flow[i,j] the flow from i (row) to j (column); component positions in rows and columns must be the same; 
    if present, rownames or columnnames denote the compartment names}
  \item{Tij }{Network matrix where connectance is from column j to row i; component positions in rows and columns must be the same ;
    if present, rownames or columnnames denote the compartment names}
  \item{Import }{Vector with either the *indices* or the *names* of \emph{external} compartmens from where flow enters the network; the indices point to the \bold{column} positions in \code{Tij} (and the row positions in Flow)}
  \item{Export }{Vector with either the *indices* or the *names* of \emph{external} compartmens to where flow leaves the network; the indices point to the \bold{row} positions in \code{Tij} (and the column positions in Flow)}
}
\value{
   A list with the following items:
  \item{TSTC }{Total system cycled throughflow}
  \item{TSTS }{Non-cycled throughflow}
  \item{FCI }{Finn's cycling index (1980)}
  \item{FCIb }{Revised Finn's cycling index, sensu Allesina and Ulanowicz, 2004.}
  \item{APL }{Average pathlength, also known as Network Aggradation (Sum of APLc and APLs in Latham 2006)}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>,
   Julius Kipyegon Kones<jkones@uonbi.ac.ke>}
\examples{# The takapoto atoll network
PathInd(Takapoto,Import ="CO2",
        Export=c("CO2","Sedimentation","Grazing"))
        # Conesprings is the example set 1a from Latham 2006.
as.data.frame(PathInd(Tij=Conesprings,Import="Inflows",
    Export=c("Export","Dissipation")))}
\references{Finn JT. 1980. Flow analysis of models of the Hubbard Brook ecosystem. Ecology 61: 562-571.
  \cr Patten BC, Higashi M. 1984. Modified cycling index for ecological applications. Ecological Modelling 25: 69-83.
  \cr Patten BC, Bosserman RW, Finn JT, Cale WG. 1976. Propagation of cause in ecosystems. Patten BC, editor. Systems Analysis and Simulation in Ecology, vol. 4. Academic Press, New York. p457-579.
  \cr Allesina and Ulanowicz, 2004. Cycling in ecological netowrks: Finn's index revisited. Computational Biology and Chemistry 28, 227-233.
  }
\details{The mathematical formulation of these indices can be found in the package vignette - \emph{vignette("NetIndices")}.


The PDF can be found in the subdirectory \file{doc} of the NetIndices package}
\keyword{utilities}

