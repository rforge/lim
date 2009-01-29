\name{GenInd}
\alias{GenInd}
\title{General network indices}
\description{Calculates general network indices such as system throughputs, link density, connectance,... of a network}
\usage{GenInd(Flow=NULL, Tij=t(Flow), 
  Import=NULL, Export=NULL, tol=0)}
\arguments{
  \item{Flow }{Network matrix with Flow[i,j] the flow from i (row) to j (column); component positions in rows and columns must be the same; 
    if present, rownames or columnnames denote the compartment names}
  \item{Tij }{Network matrix where connectance is from column j to row i; component positions in rows and columns must be the same ;
    if present, rownames or columnnames denote the compartment names}
  \item{Import }{Vector with either the *indices* or the *names* of \emph{external} compartmens from where flow enters the network; the indices point to the \bold{column} positions in \code{Tij} (and the row positions in Flow)}
  \item{Export }{Vector with either the *indices* or the *names* of \emph{external} compartmens to where flow leaves the network; the indices point to the \bold{row} positions in \code{Tij} (and the column positions in Flow)}
  \item{tol }{Flows that are smaller or equal to tol are assumed to be absent}
}
\value{
  A list that contains:
  \item{N }{Number of compartments, excluding the externals }
  \item{T.. }{Total System Throughput }
  \item{TST }{Total System Throughflow }
  \item{Lint }{Number of Internal links}
  \item{Ltot }{Total number of links}
  \item{LD }{Link Density }
  \item{C }{Connectance (internal)}
  \item{Tijbar }{Average Link Weight }
  \item{TSTbar }{Average Compartment Throughflow }
  \item{Cbar }{Compartmentalization, [0,1], the degree of connectedness of subsystems within a network }
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>,
   Julius Kipyegon Kones<jkones@uonbi.ac.ke>}
\examples{# The takapoto atoll network
(GI<- GenInd(Takapoto,Import ="CO2",
             Export=c("CO2","Sedimentation","Grazing")))
as.data.frame(GI)  
as.data.frame(GenInd(diag(5)))

# Conesprings is the example set 1a from Latham 2006.
as.data.frame(GenInd(Tij=Conesprings,Import="Inflows",
       Export=c("Export","Dissipation")))}
\references{Latham LG. 2006. Network flow analysis algorithms. Ecological Modelling 192: 586-600.
  \cr Hirata H, Ulanowicz RE. 1984. Informational theoretical analysis of ecological networks. The American Naturalist 133: 288-302.
  \cr Pimm SL, Lawton JH. 1980. Are food webs divided into compartments? Journal of Animal Ecology 49: 879-898.}
\details{The mathematical formulation of these indices can be found in the package vignette - \emph{vignette("NetIndices")}.


The PDF can be found in the subdirectory \file{doc} of the NetIndices package}
\keyword{utilities}

