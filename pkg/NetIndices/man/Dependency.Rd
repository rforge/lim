\name{Dependency}
\alias{Dependency}
\title{
  Direct and indirect dependency analysis
}
\description{
  Calculates for each component in a flow network the direct+indirect
  dependency on the other components.
}
\usage{Dependency(Flow=NULL, Tij=t(Flow), 
  Import=NULL, Export=NULL)}
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
    \emph{external} compartmens to where flow leaves the network;
    the indices point to the \bold{row} positions in \code{Tij} (and
    the column positions in Flow).
  }
}
\value{
  A matrix with dependency of component i on component j
}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>,
  Julius Kipyegon Kones<jkones@uonbi.ac.ke>}
\examples{# The takapoto atoll network
Dependency(Takapoto,Import ="CO2",
           Export=c("CO2","Sedimentation","Grazing"))

# making "Bacteria" a "primary food source"  
TAK <- Takapoto
TAK[,"Bacteria"]<-c(0,0,0,0,0,0,1,0)
# first two columns
(DD<-Dependency(TAK,Import =c("CO2","DOC"),
         Export=c("CO2","DOC","Sedimentation","Grazing")))
barplot(t(DD[3:nrow(DD),1:2]),beside=TRUE,ylab="-",
       legend=c("Phytoplankton","Bacteria"),
       main="dependency on (primary) food sources")}
\keyword{utilities}

