\name{EnvInd}
\alias{EnvInd}
\title{Environ network indices}
\description{calculates
  the indices of homogenization, synergism index, dominance of indirect effects,... of a network}
\usage{EnvInd(Flow=NULL, Tij=t(Flow), Import=NULL, 
         Export=NULL, full=FALSE)}
\arguments{
  \item{Flow }{Network matrix with Flow[i,j] the flow from i (row) to j (column); component positions in rows and columns must be the same; 
    if present, rownames or columnnames denote the compartment names}
  \item{Tij }{Network matrix where connectance is from column j to row i; component positions in rows and columns must be the same ;
    if present, rownames or columnnames denote the compartment names}
  \item{Import }{Vector with either the *indices* or the *names* of \emph{external} compartmens from where flow enters the network; the indices point to the \bold{column} positions in \code{Tij} (and the row positions in Flow)}
  \item{Export }{Vector with either the *indices* or the *names* of \emph{external} compartmens to where flow leaves the network; the indices point to the \bold{row} positions in \code{Tij} (and the column positions in Flow)}
  \item{full }{if TRUE, also returns matrices}
}
\value{
  A list with the following items: 
  \item{NAG }{Network aggradation = average path length}
  \item{HP }{Homogenization index}
  \item{BC }{Synergism}
  \item{ID }{Dominance of Indirect effects}
  \item{MN }{Mean of non-dimensional flow-matrix (N)}
  \item{MG }{Mean of direct flow-matrix (G)}
  \item{CVN }{Coefficient of variation of non-dimensional flow-matrix (N)}
  \item{CVG }{Coefficient of variation of direct flow-matrix (G)}
  \item{U }{Only if Full == TRUE: The Utility non-dimensional matrix}
  \item{N1 }{Only if Full == TRUE: The Integral non-dimensional Flow Matrix}
  \item{G }{Only if Full == TRUE: The Normalized direct flow (or transitive closure matrix) matrix}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>,
   Julius Kipyegon Kones<jkones@uonbi.ac.ke>}
\examples{# The takapoto atoll network
EnvInd(Takapoto,Import ="CO2",
         Export=c("CO2","Sedimentation","Grazing"))

as.data.frame(EnvInd(Tij=Conesprings,Import="Inflows",
    Export=c("Export","Dissipation")))
EnvInd(Tij=Conesprings,Import="Inflows",
    Export=c("Export","Dissipation"),full=TRUE)}
\references{Patten BC, Barber MC, Richardson TH. 1982. Path analysis of a reservoir ecosystem model. 
  \cr Fath BD, Patten BC. 1999. Review of the foundations of network environ analysis. Ecosystems 2: 167-179.
  \cr Fath BD, Patten BC. 1999. Quantifying resource homogenization using network flow analysis. Ecological Modelling 123: 193-205.
  \cr Patten BC, Higashi M. 1984. Modified cycling index for ecological applications. Ecological Modelling 25: 69-83.
  \cr Higashi M, Patten BC. 1989. Dominance of indirect causality in ecosystems. The American Naturalist 133: 288-302.}
\details{The mathematical formulation of these indices can be found in the package vignette - \emph{vignette("NetIndices")}.


The PDF can be found in the subdirectory \file{doc} of the NetIndices package}
\keyword{utilities}

