\name{Read}
\alias{Read}
\title{Reads an inverse input file}
\description{Reads an inverse input file and creates the inverse problem as a list, of type "liminput"}
\usage{Read(file, verbose=FALSE, checkLinear=TRUE, remtabs=TRUE)}
\arguments{
  \item{file }{name of inverse input file}
  \item{verbose }{if TRUE: prints warnings and messages to the screen}
  \item{checkLinear }{if FALSE: does not check for linearity}
  \item{remtabs }{remove tabs}
}
\value{
   a list containing :
  \item{file }{name of inverse input file}
  \item{pars }{a data.frame with parameter declarations}
  \item{comp }{a data.frame with compartments (or states, stocks)}
  \item{rate }{a data.frame with rate declarations}
  \item{extern }{a data.frame with external declarations}
  \item{flows }{a data.frame with flow declarations}
  \item{vars }{a data.frame with variable declarations}
  \item{cost }{a data.frame with cost declarations}
  \item{profit }{a data.frame with profit declarations}
  \item{equations }{a data.frame with equality declarations}
  \item{constraints }{a data.frame with constraint declarations}
  \item{reactions }{a data.frame with reaction declarations}
  \item{posreac }{a vector with TRUE if reaction or flow is unidirectional (and is thus positive), FALSE if two-way reaction or flow}
  \item{marker }{a data.frame with marker declarations - see vignette("LIM")}
  \item{parnames }{a vector with parameter names}
  \item{varnames }{a vector with variable names}
  \item{compnames }{a vector with compartment names}
  \item{externnames }{a vector with  names of externals}
  \item{Type }{a string; one of "web" (flows are unknowns), "reaction" (reaction rates unknown) and "simple" (compartments are unknowns}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
  # this input has been created with function Read:
  LIMinputBlending
  \dontrun{
  wd <- getwd()
  setwd(paste(system.file(package="LIM"), "/examples/Foodweb", sep=""))
  Read("RigaAutumn.input")
  setwd(wd)
  }
  }
\details{
The structure of an inverse input file is explained in vignette("LIM") which should be consulted.

In short the inverse input file contains the declaration sections enclosed inbetween two lines starting with a \#\#.
	\cr For instance, the following section declares two components
   \cr \# COMP
   \cr State1
   \cr State2
   \cr \# END COMP

   Only the first 4 characters of the section names are read
   \cr The following sections are allowed:
  \itemize{
   \item Parameters   - \#\# PARAMETERS
   \item Components   - \#\# STOCKS or \#\# DECISION VARIABLES or \#\# STATES or \#\# UNKNOWNS
   \item Externals    - \#\# EXTERNALS
   \item Rates        - \#\# RATES
   \item Flows        - \#\# FLOWS
   \item Variables    - \#\# VARIABLES
   \item Cost         - \#\# COST or \#\# MINIMISE
   \item Profit       - \#\# PROFIT or \#\# MAXIMISE
   \item Equalities   - \#\# EQUALITIES
   \item InEqualities - \#\# INEQUALITIES or \#\# CONSTRAINTS
}
   Any (part of a) line starting with a "!" is considered a comment.

   Input is NOT case sensitive

   The output of this function is used as input in function \code{\link{Setup}} which creates the inverse matrices

   By default, only linear problems can be solved, and the function checks whether the input is linear.
   To toggle off this check, set checkLinear to FALSE.

	 Some input files contain tabs, which are converted to spaces, unless this logical is set to FALSE.}
\seealso{\code{\link{Setup}} the function to create inverse matrices, based on output of \code{Read}}
\keyword{IO}

