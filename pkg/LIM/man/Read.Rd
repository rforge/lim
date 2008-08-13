\name{Read}
\alias{Read}
\title{Reads an inverse input file}
\description{Reads an inverse input file and creates the inverse problem as a list.}
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
  \item{state }{a data.frame with state declarations}
  \item{rate }{a data.frame with rate declarations}
  \item{extern }{a data.frame with external declarations}
  \item{flows }{a data.frame with flow declarations}
  \item{vars }{a data.frame with variable declarations}
  \item{cost }{a data.frame with cost declarations}
  \item{profit }{a data.frame with profit declarations}
  \item{equations }{a data.frame with equality declarations}
  \item{constraints }{a data.frame with constraint declarations}
  \item{parnames }{a vector with parameter names}
  \item{varnames }{a vector with variable names}
  \item{statenames }{a vector with state names}
  \item{externnames }{a vector with external names}
  \item{web }{a logical; if TRUE: the input was a web, and flows are unknowns; if FALSE: states are the unknowns}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
  \dontrun{Read("FILERigaAutumn")}
  }
\details{The inverse input file contains the declaration sections enclosed inbetween two lines starting with a \#\#.
	\cr For instance, the following section declares two state variables
   \cr \# STATE
   \cr State1
   \cr State2
   \cr \# END STATE

   Only the first 4 characters of the section names are read
   \cr The following sections are allowed:
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

   Any (part of a) line starting with a "!" is considered a comment.
   \cr Input is NOT case sensitive

   The output of this function is used as input in function \code{\link{Setup}} which creates the inverse matrices

   By default, only linear problems can be solved, and the function checks whether the input is linear.
   To toggle off this check, set checkLinear to FALSE.

	 Some input files contain tabs, which are converted to spaces, unless this logical is set to FALSE.}
\seealso{\code{\link{Setup}} the function to create inverse matrices, based on output of \code{Read}}
\keyword{IO}

