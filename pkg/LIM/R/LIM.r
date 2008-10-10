################################################################################
##                                                                            ##
## LINEAR INVERSE MODELLING INPUT   -   LIM                                   ##
##    Karline Soetaert                                                        ##
##                                                                            ##
## -----------------------------                                              ##
## part 1: LIM problem, input as a list                                       ##
## -----------------------------                                              ##
##                                                                            ##
## read.limfile :Reads inverse input text file; returns a list                ##
## setup.limfile:Creates inverse model using list returned by read.limfile    ##
## Lsei.limfile:Reads inputfile + solves the linear inverse model with lsei   ##
## Ldei.limfile :Reads inputfile + solves the linear inverse model with ldei  ##
## Linp.limfile:Reads inputfile + solves the linear programming model         ##
##                                                                            ##
## -----------------------------                                              ##
## part 2: solving LIM models                                                 ##
## -----------------------------                                              ##
##                                                                            ##
## Ldei.lim        : Solves inverse model, input as list, uses ldei           ##
## Lsei.lim        : Solves inverse model, input as list, uses lsei           ##
## Linp.lim        : Solves linear programming model, input as list           ##
## Variables       : Creates a list with the variables, input as list         ##
## Flowmatrix      : Creates a matrix with flow magnitudes                    ##
## Xranges         : Calculates ranges of inverse unknowns, input list        ##
## Varranges       : Calculates ranges of inverse variables, input list       ##
##                                                                            ##
## plotranges      : plots variable and flow ranges                           ##
##                                                                            ##
################################################################################

Lsei <- function(...) UseMethod ("Lsei")
Ldei <- function(...) UseMethod ("Ldei")
Linp <- function(...) UseMethod ("Linp")
Plotranges <- function(...) UseMethod ("Plotranges")

Lsei.character <- function(...) Lsei.limfile(...)
Ldei.character <- function(...) Ldei.limfile(...)
Linp.character <- function(...) Linp.limfile(...)

Lsei.double <- function(...) lsei(...)
Ldei.double <- function(...) ldei(...)
Linp.double <- function(...) linp(...)

################################################################################
# reads an input file and solves the linear inverse model using lsei
################################################################################
Lsei.limfile <- function(file,
                     exact =NULL,     # a vector with the equations to be solved EXACTLY
                     parsimonious = FALSE, # if true: also minimises sum of squared unknowns
                     verbose=TRUE,
                     ...)
{
# read file and create inverse matrices
lim    <- Setup.limfile(file,verbose=verbose)
# solve
Lsei.lim(lim,exact,parsimonious=parsimonious,...)

}

################################################################################
# reads an input file and solves the linear programming inverse model
################################################################################

Linp.limfile <- function(file,verbose=TRUE,...)
{
lim    <- Setup.limfile(file,verbose=verbose)
Linp.lim(lim,...)
}   # Linp.limfile

################################################################################
# reads an input file and solves the linear programming inverse model
################################################################################

Ldei.limfile <- function(file,verbose=TRUE,...)
{
lim    <- Setup.limfile(file,verbose=verbose)
Ldei.lim(lim,...)
}

#########################################################################
## Linp.lim : Solves linear programming model                      ##
#########################################################################

Linp.lim <- function(lim,cost=NULL,ispos=lim$ispos,...)

#------------------------------------------------------------------------
# Solves inverse problem, using linear programming
# inverse is a list that contains the problem
#------------------------------------------------------------------------

{
   Nx <- lim$NUnknowns
   B  <-lim$B
   H  <-lim$H
   A  <- NULL
   G  <- NULL
   if(is.null(cost)) Cost   <- lim$Cost else Cost <- cost
   if(is.null(cost)) Profit <- lim$Profit else Profit <- NULL
   if (! is.null(Cost  )&&is.vector(Cost  )) Cost   <- matrix(nr=1,data=Cost  )
   if (! is.null(Profit)&&is.vector(Profit)) Profit <- matrix(nr=1,data=Profit)
if (ispos)
 { A <- lim$A
   G <-lim$G
  } else
   { if (! is.null(lim$A))A <- cbind(lim$A,-1*lim$A)
     if (! is.null(lim$G))G <- cbind(lim$G, -1*lim$G)
     if (! is.null(Cost))   Cost   <- cbind(Cost, -1 * Cost)
     if (! is.null(Profit)) Profit <- cbind(Profit, -1 * Profit)
  }  
res <- NULL
i1 <- 0
if (!is.null(Cost))
 {
 for ( i in 1:nrow(Cost))
   {resi <- linp (A,B,G,H,Cost[i,],...)
    res$residualNorm <- c(res$residualNorm, resi$residualNorm)
    res$solutionNorm <- c(res$solutionNorm, resi$solutionNorm)
    res$X <- rbind(res$X,matrix(nr=1,resi$X))
    rownames(res$X)[i]<-as.character(lim$costnames[i])
   }
   i1 <- i
} else if (is.null(lim$Profit))

{return(NULL)}
if (! is.null(Profit))
{for ( i in 1:nrow(Profit))
   {resi<-linp (A,B,G,H,-1*Profit[i,],...)
    res$residualNorm <- c(res$residualNorm, resi$residualNorm)
    res$solutionNorm <- c(res$solutionNorm, -resi$solutionNorm)
    res$X <- rbind(res$X,matrix(nr=1,resi$X))
    rownames(res$X)[i1+i]<-as.character(lim$profitnames[i])
    }
}
if (!ispos) res$X <-  res$X[,1:Nx]-res$X[,(Nx+1):(2*Nx)]
res$X <- matrix(nc=Nx,data=res$X)
colnames(res$X) <- lim$Unknowns

return(res)
}
 ########## END OF Linp.lim ##########

#########################################################################
## Ldei.lim : Solves inverse model, least distance programming      ##
#########################################################################

Ldei.lim <- function(lim,...)

#------------------------------------------------------------------------
# Solves inverse problem, using ldei
# inverse is a list that contains the problem
#------------------------------------------------------------------------

{
ld<-ldei (E=lim$A,F=lim$B,G=lim$G,H=lim$H,...)
names(ld$X) <- lim$Unknowns
if (ld$IsError) warning("Problem could not be solved")

return(ld)
} ########## END OF Ldei.lim ##########


#########################################################################
## Lsei.lim  : Solves inverse model, lsei                          ##
#########################################################################

Lsei.lim <- function(lim,                  # the linear inverse matrices, a list
                     exact =NULL,          # a vector with the equations to be solved EXACTLY
                     parsimonious = FALSE, # if true: also minimises sum of squared unknowns
                     ...)                  # to print error messages

#------------------------------------------------------------------------
# Solves inverse problem, using lsei
# lsei= least squares with equality and inequality constraints
# lim is a list that contains the input matrices
# if parsimonious is true: sum of squared flows has to be minimised too
#------------------------------------------------------------------------

{
# 0. Setup problem

Ncomp  <- lim$NComponents
Nx     <- lim$NUnknowns

A<-E <- NULL
B<-F <- NULL
Napp <- nrow(lim$A)

if(is.null(exact))
{
# Equalities: all equations are approximate
 Neq  <- 0
 E      <- lim$A[]
 F      <- lim$B[]
 
} else  {                  # minimise sum of squared flows

# Equalities and approximate equations
if (max(exact)>Napp)stop("error: cannot solve Lsei.lim. Equations to be met exactly do not exist")
Neq  <- length(exact)
E    <- lim$A[exact,]
F    <- lim$B[exact]
# Approximate equations : MIN |AX-B|

Napp   <- Napp-Neq    # Number approximate
if (Napp > 0)
{
 A      <- lim$A[-exact,]
 B      <- lim$B[-exact]
}
}

 if (parsimonious) 
 {
 A <- rbind(A,diag(nrow=Nx,ncol=Nx))
 B <- c(B,rep(0,Nx))
 }
 if (is.null(A))warning("Lsei.lim - there are no approximate equations in lsei!")

# Inequalities G*X>=H
G <- lim$G
H <- lim$H


sol <- lsei(A,B,E,F,G,H,...)
names(sol$X) <- lim$Unknowns
if (sol$IsError) warning("Problem could not be solved - least squares solution returned")
return(sol)

} ########## END OF Lsei.lim ##########

#########################################################################
## Xranges  : Calculates ranges of inverse components (linp problem)   ##
#########################################################################

Xranges <- function(lim,...)

{

res   <- xranges (lim$A,lim$B,
                  lim$G,lim$H,...)
rownames(res) <- lim$Unknowns
return(res)

} ########## END OF Xranges ##########

#########################################################################
## Varranges: Calculates ranges of inverse variables            ##
#########################################################################

Varranges <- function(lim,...)

{
if (lim$NVariables == 0) return(NULL)
res   <- varranges (lim$A,lim$B,
                   lim$G,lim$H,
                   lim$VarA,lim$VarB,...)
rownames(res)<-lim$Variables
return(res)

} ########## END OF Varranges ##########
#########################################################################
## Variables : Creates a list with the web variables           ##
#########################################################################

Variables <- function(lim,res=NULL,...)

{
if (lim$NVariables == 0) return(NULL)

if (is.null(res)) res <- Lsei.lim(lim,parsimonious=TRUE,...)$X

variables <- data.frame(values=lim$VarA%*%res-lim$VarB)
rownames(variables)<-lim$Variables

return(variables)

} ########## END OF Variables ##########

#########################################################################
## Flowmatrix: Creates a matrix with the web flows             ##
#########################################################################

Flowmatrix <- function(lim,web=NULL)

{
flowmatrix <- lim$Flowmatrix
if (is.null(web)) web <- Lsei.lim(lim,parsimonious=TRUE)$X

X          <- as.vector(web)
Xpos       <- pmax(0.,X)
ii         <- which(flowmatrix >0,arr.ind=TRUE)
flowmatrix[ii]<-Xpos[lim$Flowmatrix[ii]]

Xneg       <- -1*pmin(0.,X)
if(sum(Xneg)>0)
{
flowmatrix[ii[,c(2,1)]]<-flowmatrix[ii[,c(2,1)]]+Xneg[lim$Flowmatrix[ii]]
}
return(flowmatrix)

} ########## END OF Flowmatrix ##########

#########################################################################
## plotranges   : plots ranges and names                               ##
#########################################################################

Plotranges.double <- function (min,                # minimum value
                               max,                # maximum value
                               value = NULL,       # median or mean value
                              labels=NULL,        # names of each value
                              log="",             # if = x: logarithmic scale for x-axis
                              pch = 16,           # pch symbol used for mean value
                              pch.col="black",     # pch color for mean value
                              line.col = "gray",   # color for each variable, spanning x-axis
                              seg.col="black",     # color for variable range
                              xlim = NULL,
                              main = NULL,
                              xlab = NULL,
                              ylab = NULL,
                              lab.cex = 1.0,
                              mark = NULL,
                              ...)                   # arguments passed to R-function "text" when writing labels
  {
  
  ##-----------------------------------------------------------------
  ## constructing the data 
  ##-----------------------------------------------------------------
      if (! is.vector(value)) value<-as.vector(unlist(value))
      ranges   <- cbind(min,max,value)
      if (log =="x") {

	    minflow<-min(ranges[ranges!=0])     ## minimum, different from 0
      ranges[ranges==0] <- minflow
      min[min==0]       <- minflow        ## replace 0 with minimum	
      max[max==0]       <- minflow	
      value[value==0]   <- minflow	
            }
    numflows <- length(min)

  ##-----------------------------------------------------------------
    if (is.null(labels)) labels <- names(min)
    if (is.null(labels)) labels <- names(max)
    if (is.null(labels)) labels <- as.character(1:numflows)
    
    labelwidth   <- max(strwidth(labels, "inch")*lab.cex, na.rm = TRUE)
    labelheight  <- strheight("M", "inch")
    plot.new()

  ##-----------------------------------------------------------------
  ## new margins
  ##-----------------------------------------------------------------

    nmar         <- nm <- par("mar")
    nmar[2]      <- nmar[4] + (labelwidth + 0.1)/labelheight
    par(mar = nmar)

    y            <- 1:numflows 
    if (is.null (xlim)) xlim <- range(ranges, na.rm=TRUE)
    ylim <- c(0, numflows + 1)

    plot.window(xlim = xlim, ylim = ylim, log = log)

    loffset <- (labelwidth + 0.1)/labelheight
    labs    <- labels[y]
    mtext(labs, side = 2, line = loffset, at = y, adj = 0, 
           las = 2, cex = par("cex") * lab.cex, ...)
 
    abline  (h = y, lty = "dotted", col = line.col)
    if(!is.null(value)) points  (value, y, pch = pch, col = pch.col)
    segments(min,y,max,y,col=seg.col,lty=1)
    if (! is.null(mark)) {
    text(labels=rep("*",length(mark)),
         rep(xlim[2],length(mark)),mark)
    }
    axis(1)
    box()
    title(main = main, xlab = xlab, ylab = ylab, ...)
    invisible()
    par("mar"=nm)

}  ########## END OF plotranges ########## 


#########################################################################
## plotranges   : plots ranges and names                               ##
#########################################################################

Plotranges.lim <- function (lim=NULL,           # lim
                            labels=NULL,        # names of each value
                            type="X",
                            log="",             # if = x: logarithmic scale for x-axis
                            pch = 16,           # pch symbol used for mean value
                            pch.col="black",     # pch color for mean value
                            line.col = "gray",   # color for each variable, spanning x-axis
                            seg.col="black",     # color for variable range
                            xlim = NULL,
                            main = NULL,
                            xlab = NULL,
                            ylab = NULL,
                            lab.cex = 1.0,
                            index=NULL,         # if not NULL, list of values to be plotted...
                            ...)                # arguments passed to function "plotranges"
  {
  
  ##-----------------------------------------------------------------
  ## constructing the data 
  ##-----------------------------------------------------------------
    if (type == "V")
    {ranges <- Varranges(lim)
     value  <- Variables(lim)

    } else
    {ranges <- Xranges(lim)
     value  <- Lsei.lim(lim,parsimonious=TRUE)$X}
    value <- unlist(value)
    if (is.null(index))
     {index <- 1:nrow(ranges)} else {
     
    ranges <- ranges[index,]
    value <- value[index]            }

    infinity <- which (ranges[,2]==1e30)
    if (length(infinity) >0) ranges[infinity,2]<- NA else infinity <- NULL

    if (is.null(labels)) labels <- rownames(ranges)
    if (is.null(labels)) labels <- as.character(1:length(value))

    Plotranges.double(min=ranges[,1],max=ranges[,2],value =value,
   labels=labels,log=log,pch=pch,pch.col=pch.col,line.col=line.col,
   seg.col=seg.col,xlim=xlim,main=main,xlab=xlab,ylab=ylab,
   lab.cex=lab.cex,mark=infinity,...)

}  ########## END OF plotranges.lim ########## 


#########################################################################
## plotranges   : plots ranges and names                               ##
#########################################################################

Plotranges.character <- function (file,           # lim
                            ...)                # arguments passed to function "plotranges"
  {

  ##-----------------------------------------------------------------
  ## constructing the data
  ##-----------------------------------------------------------------

    lim <- Setup(file)
    Plotranges.lim(lim,...)

}  ########## END OF plotranges.character ##########

#
#########################################################################
## Xranges.lim  : Calculates ranges of inverse components (linp problem)   ##
#########################################################################

Xsample <- function(lim,
                     exact =NULL,     # a vector with the equations to be solved EXACTLY)
                     ...) # other arguments passed to function xsample.

{
# 0. Setup problem

Ncomp  <- lim$NComponents
Nx     <- lim$NUnknowns

A<-E <- NULL
B<-F <- NULL
Napp <- nrow(lim$A)
if(is.null(exact))
{
# Equalities: all equations are approximate
 Neq  <- 0
 E      <- lim$A[]
 F      <- lim$B[]
 
} else  {                  # minimise sum of squared flows

# Equalities and approximate equations
if (max(exact)>Napp)stop("error: cannot solve Lsei.lim. Equations to be met exactly do not exist")
Neq  <- length(exact)
E    <- lim$A[exact,]
F    <- lim$B[exact]
# Approximate equations : MIN |AX-B|

Napp   <- Napp-Neq    # Number approximate
if (Napp > 0)
{
 A      <- lim$A[-exact,]
 B      <- lim$B[-exact]
}
}

# Inequalities G*X>=H
G <- lim$G
H <- lim$H



res   <- xsample (A=A,B=B,E=E,F=F,
                  G=G,H=H,...)$X
colnames(res) <- lim$Unknowns
#}
return(res)

} ########## END OF Xsample ##########

PrintMat<-function(lim)
{
A           <- lim$A
colnames(A) <- lim$Unknowns
rownames(A) <- lim$eqnames
G           <- lim$G
colnames(G) <- lim$Unknowns
rownames(G) <- lim$ineqnames
print("A,B")
print(cbind(A,"  B = "=lim$B))

print("G,H")
print(cbind(G,"  H = "=lim$H))

if (!is.null(lim$Cost))
{
Cost           <-  matrix(nc=lim$NUnknowns,data=lim$Cost)
colnames(Cost) <- lim$Unknowns
rownames(Cost) <- lim$costnames
print("Cost")
print(Cost)
}

if (!is.null(lim$Profit))
{
Profit           <-  matrix(nc=lim$NUnknowns,data=lim$Profit)
colnames(Profit) <- lim$Unknowns
rownames(Profit) <- lim$profitnames
print("Profit")
print(Profit)
}
}


