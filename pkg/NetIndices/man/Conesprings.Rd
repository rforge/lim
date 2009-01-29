\name{Conesprings}
\docType{data}
\alias{Conesprings}
\title{Cone Spring ecosystem}
\description{Cone spring ecosystem (Tilly, 1968) adapted for input - output
analysis by Williams & Crouthamel (1972) in Szyrmer & Ulanowicz (1987, Fig. 1, p.
129) and Ulanowicz & Norden (1990, Fig. 1, p. 435).

This is example 1a from Latham (2006).

The food web comprises 5 functional compartments:
\itemize{
\item  Plants
\item  Detritus
\item  Bacteria
\item  Detritus feeders
\item  Carnivores
}
and two export compartments
\itemize{
\item  usable export
\item  dissipation
}
and one import compartment

\itemize{
\item  Inflows}
}
\usage{Conesprings}
\format{matrix with Tij values, where element (i,j) denotes flow from compartment j to i
\cr rownames and columnames are the components}
\examples{GenInd(Tij=Conesprings,Import="Inflows",
                 Export=c("Export","Dissipation"))
AscInd(Tij=Conesprings,Import="Inflows",
                 Export=c("Export","Dissipation"),
                 Dissipation="Dissipation")
UncInd(Tij=Conesprings,Import="Inflows",
                 Export=c("Export","Dissipation"))
EffInd(Tij=Conesprings,Import="Inflows",
                 Export=c("Export","Dissipation"))
EnvInd(Tij=Conesprings,Import="Inflows",
           Export=c("Export","Dissipation"),full=TRUE)}
\references{Latham LG. 2006. Network flow analysis algorithms. Ecological Modelling 192: 586-600.

Szyrmer, J., & Ulanowicz, R. E. (1987). Total flows in ecosystems. Ecol. Model. 35, 123..136.

Tilly, L. J. (1968). The structure and dynamics of Cone Spring. Ecol. Monogr. 38, 169..197.

Ulanowicz, R. E., & Norden, J. S. (1990). Symmetrical overhead in flow networks. Int. J.
Systems Sci. 21, 429..437.

Williams, M., & Crouthamel, D. (1972). Systems analysis of Cone Spring. Unpublished
manuscript. University of Georgia, Athens, Georgia.}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}

\keyword{datasets}

