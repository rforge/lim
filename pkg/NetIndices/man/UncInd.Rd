\name{UncInd}
\alias{UncInd}
\title{Network uncertainty indices and constraint effiencies}
\description{calculates the statistical, conditional and realised uncertainty, the average mutual information index,
and the network uncertainty, network constraint and constraint efficiency,...}
\usage{UncInd(Flow=NULL, Tij=t(Flow), 
  Import=NULL, Export=NULL)}
\arguments{
  \item{Flow }{Network matrix with Flow[i,j] the flow from i (row) to j (column); component positions in rows and columns must be the same; 
    if present, rownames or columnnames denote the compartment names}
  \item{Tij }{Network matrix where connectance is from column j to row i; component positions in rows and columns must be the same ;
    if present, rownames or columnnames denote the compartment names}
  \item{Import }{Vector with either the *indices* or the *names* of \emph{external} compartmens from where flow enters the network; the indices point to the \bold{column} positions in \code{Tij} (and the row positions in Flow)}
  \item{Export }{Vector with either the *indices* or the *names* of \emph{external} compartmens to where flow leaves the network; the indices point to the \bold{row} positions in \code{Tij} (and the column positions in Flow)}
}
\value{
  a list with the following items:
  \item{AMI }{The average mutual information; as a system matures to form a web-like pattern, the AMI drops}
  \item{HR }{The statistical uncertainty, upper bound on AMI, a measure of diversity}
  \item{DR }{The conditional uncertainty index, the difference between AMI and HR, a measure of stability}
  \item{RU }{The realised uncertainty index, ratio of AMI and HR}
  \item{Hmax }{Maximum uncertainty}
  \item{Hc }{Constraint information}
  \item{Hsys }{Network uncertainy}
  \item{CE }{Constraint efficiency}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>,
   Julius Kipyegon Kones<jkones@uonbi.ac.ke>}
\examples{# The takapoto atoll network
UncInd(Takapoto,Import ="CO2",
       Export=c("CO2","Sedimentation","Grazing"))

# Conesprings is the example set 1a from Latham 2006.
as.data.frame(UncInd(Tij=Conesprings,Import="Inflows",
    Export=c("Export","Dissipation")))}
\references{Latham LG. 2006. Network flow analysis algorithms. Ecological Modelling 192: 586-600.
  \cr Ulanowicz RE, Norden JS. 1990. Symmetrical overhead in flow networks. International Journal of System Science 21: 429-437.
  \cr Gallager RG. 1968. Information Theory and Reliable Communication. Wiley, New York.
  \cr Shannon CE. 1948. A mathematical theory of communication. Bell System Technical Journal 27: 379-423. 
  \cr Ulanowicz RE. 1997. Ecology, the ascendent perspective. Allen TFH, Roberts DW, editors. Complexity in Ecological Systems Series. Columbia University Press, New York..
  \cr Latham LG, Scully EP. 2002. Quantifying constraint to assess development in ecological networks. Ecological Modelling 154: 25-44.
  \cr Rutledge RW, Basorre BL, Mulholland RJ. 1976. Ecological stability: an information theory viewpoint. Journal of Theoretical Biology 57: 355-371.}
\details{The mathematical formulation of these indices can be found in the package vignette - \emph{vignette("NetIndices")}.


The PDF can be found in the subdirectory \file{doc} of the NetIndices package}
\keyword{utilities}

