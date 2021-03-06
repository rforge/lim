\name{Plotranges}
\alias{Plotranges}
\alias{Plotranges.double}
\alias{Plotranges.lim}
\alias{Plotranges.character}
\title{
  Plots the minimum and maximum and central values
}
\description{
  Plots minimum and maximum ranges.

  Takes as input either a lim list, as generated by Setup or a set
  of vectors specifying the minimum, maximum and the central value,
  or a data.frame that contains min, max and central values.
  
}
\usage{
Plotranges(...)
\method{Plotranges}{double}(min, max, value = NULL, labels = NULL, log = "",
  pch = 16, pch.col = "black", line.col = "gray", 
  seg.col = "black", xlim = NULL, main = NULL, 
  xlab = NULL, ylab = NULL, lab.cex = 1.0, mark = NULL,...)

\method{Plotranges}{lim}(lim = NULL, labels = NULL, type = "X", log = "",
  pch = 16, pch.col = "black", line.col = "gray", 
  seg.col = "black", xlim = NULL, main = NULL, 
  xlab = NULL, ylab = NULL, lab.cex = 1.0, index = NULL, ...)

\method{Plotranges}{character}(file, ...)

}
\arguments{
  \item{min }{minimum value.
  }
  \item{max }{maximum value.
  }
  \item{value }{median or mean value.
  }
  \item{lim }{a list that contains the linear inverse model
    specification, as  generated by function \code{setup.limfile}.
  }
  \item{file }{name of the inverse input file.
  }
  \item{labels }{names of each value.
  }
  \item{type}{one of "X" or "V" for plotting of unknowns (X) or variables.
  }
  \item{log }{if = x: logarithmic scale for x-axis.
  }
  \item{pch }{pch symbol used for mean value.
  }
  \item{pch.col }{pch color for mean value.
  }
  \item{line.col }{color for each variable, spanning x-axis.
  }
  \item{seg.col }{color for variable range.
  }
  \item{xlim }{limits on x-axis.
  }
  \item{main }{main title.
  }
  \item{xlab }{x-axis label.
  }
  \item{ylab }{y-axis label.
  }
  \item{lab.cex }{label expansion value.
  }
  \item{index}{list of elements to be plotted, a vector of integers;
    default = all elements.
  }
  \item{mark}{list of elements to be marked with a "*", i.e. when range
    is unbounded.
  }
  \item{... }{arguments passed to R-function "text" when writing labels.
  }
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\value{
  Only when a lim list was inputted. A data frame with
  \item{min }{the minimum.
  }
  \item{max }{the maximum.
  }
  \item{values }{the central value.
  }
}

\examples{
# The Takapoto food web.
# some ranges extend to infinity - they are marked with "*"
Plotranges(LIMTakapoto, lab.cex = 0.7, sub = "*=unbounded",
           xlab = "mgC/m2/d",  
           main = "Takapoto atoll planktonic food web, Flowranges")
# ranges of variables, exclude first variable
Plotranges(LIMTakapoto, type = "V", lab.cex = 0.7, 
           index = 2:23, xlab = "mgC/m2/d",
           main = "Takapoto atoll planktonic food web, Variable ranges")
}

\keyword{hplot}

