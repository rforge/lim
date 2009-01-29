
#############################################################################
# Implementation of Network indices, as in 
# Latham, LG II 2006 Ecol. Modeling 192: 586-600
#
# Implemented: Julius Kones      - University Nairobi
#              Karline Soetaert  - Netherlands Institute of Ecology
#
#############################################################################


##################################################################################
#_________________________________________________________________________________
# GENERAL SUITE - COMPARTMENTS, TOTAL SYSTEM THROUGHPUT (T..),
#           TOTAL SYSTEM THROUGHFLOW (TST), LINK DENSITY (LD),
#           NUMBER OF LINKS (L), AV. COMPARTMENT THROUGHFLOW (TSTbar),
#           CONNECTANCE (C), AVERAGE LINK WEIGHT (Tijbar), 
#           COMPARTMENTALIZATION (Cbar)
#_________________________________________________________________________________

#############################################################################

GenInd <- function (Flow = NULL,                      # from-to 
                    Tij  = t(Flow),                   # to-from
                    Import =NULL,                     # flow from external (colNr Tij)
                    Export =NULL,                     # flow to external (colNr Tij)
                    tol=0)                            # flow<=tol is assumed absent 
{                      

#------------------------------------------------------------------------
# Flow is a matrix with Flow[i,j] the flow from i (row) to j (column)
# Tij[i,j] contains flow from j to i
# note: component position in rows and columns must be the same - not checked
#------------------------------------------------------------------------

    N          <- InternalNetwork (Tij,Import,Export) 

    # Rate of change of each compartment
    RateComp   <- N$FlowToC-N$FlowFromC

    # number of columns and rows for total flows (Tij)
    ncTij       <- ncol(Tij)
    nrTij       <- nrow(Tij)
    ncomp       <- ncol(N$Tint)              # the number of compartments (without externals)
    compNames   <- rownames(N$Tint)

 
#____________________________________________________
# NUMBER OF TOTAL AND INTERNAL LINKS AND LINK DENSITY
#____________________________________________________
    intlinks    <- length(which(N$Tint>tol))
    links       <- length(which(Tij >tol))
    LD          <- links/ncomp

#____________________________________________________
# THROUGHPUT AND THROUHFLOW  
#____________________________________________________
    ExportSum   <- sum(N$FlowTo[N$export])
    ImportSum   <- sum(N$FlowFrom[N$import])

    #___________________________________________________________
    # THROUGHFLOW - Rates of change taken into consideration
    # Throughflow based on rows and columns has to be the same

    Throughflow <- sum(N$Tint) + ImportSum - sum(RateComp[RateComp<0])
#   Throughflow <- sum(N$Tint) + ExportSum + sum(RateComp[RateComp>0])    

    #___________________________________________________________
    # THROUGHPUT = Sum of all flows, including externals

    Throughput  <- sum(Tij)

#____________________________________________________
# AVERAGE COMPARTMENT THROUGHFLOW (TSTbar)
#____________________________________________________

    Avthrflow <- Throughflow/ncomp

#____________________________________________________
# CONNECTANCE
#____________________________________________________
    
    Connectance <- intlinks/ncomp/(ncomp-1)

#____________________________________________________
# AVERAGE LINK WEIGHT
#____________________________________________________

    Avlinkweight <- Throughput/links

#____________________________________________________
# COMPARTMENTALIZATION
#____________________________________________________
    linkmat             <- N$Tint
    linkmat[linkmat>0]  <- 1   # 1 if there is a link, 0 elsewhere

# The number of components with which both i and j interact divided by
# the number of species by which either i or j interact
    Cij <- matrix(nrow=ncomp,ncol=ncomp,0)
    for (i in 1:ncomp) {
        int_i <- union(which(linkmat[i,]>0),which(linkmat[,i]>0))
        for (j in 1:ncomp) {
            int_j <- union(which(linkmat[j,]>0),which(linkmat[,j]>0))
            sect <- intersect(int_i,int_j)
            uni  <- union    (int_i,int_j)
            Cij[i,j]  <- length(sect)/length(uni)
                           } 
                       }
    Compart <- (sum(Cij)-ncomp) / ncomp /(ncomp-1)


#    "No. of Compartments (n)","Total System Throughput (T..)",
#    "Total System Throughflow (TST)","Link Density (LD)",
#    "No. of Internal links, (Lint)","Total No. of links (Ltot)",
#    "Av. Compartment Throughflow (TSTbar)","Connectance (internal)(C)",
#    "Av. Link Weight (Tijbar)","Compartmentalization (Cbar)"
        
    list(N=ncomp,T..=Throughput,TST=Throughflow,
      Lint=intlinks,Ltot=links,LD=LD,C=Connectance,
      Tijbar=Avlinkweight,TSTbar=Avthrflow,Cbar=Compart)
       
    }              #end generalIndices

###############################################################################

InternalNetwork <- function (Tij  ,                      # to-from
                             Import,                     # flow from external (colNr Tij)
                             Export)                     # flow to external (colNr Tij)

{                      

#------------------------------------------------------------------------
# Tij[i,j] is a matrix with Tij[i,j]  flow from j to i
# note: component position in rows and columns must be the same - not checked
#------------------------------------------------------------------------

    if (is.character(Import)) 
       import <- which(colnames(Tij)%in%Import) else
       import <- Import
    if (length(import) != length(Import)) stop("Import not recognized")
    if (is.character(Export)) 
       export <- which(rownames(Tij)%in%Export) else
       export <- Export
    if (length(import) != length(Import)) stop("Import not recognized")
    if (length(export) != length(Export)) stop("Export not recognized")

#_________________________________________________________________________________
# CHECK THE INPUT            
#_________________________________________________________________________________

    # Flow or Tij should be inputted
    if (is.null(Tij)) stop ("cannot calculate indices - Flow or Tij should be inputted")

#_________________________________________________________________________________
# NUMBER OF COMPARTMENTS, export, import, internal flows,..
#_________________________________________________________________________________

   # Size of the matrices; without the externals, the matrix has to be square

    ncomp     <- ncol(Tij)-length(import)   # the number of compartments (without externals)
    if (ncomp != nrow(Tij)-length(export)) 
         stop ("cannot calculate indices - internal flow input matrix not square ")
 
#_________________________________________________________________________________
# ARRAYS DECLARATION         
#_________________________________________________________________________________

    # indices to elements of T that are internal 
    iN  <- setdiff(1:nrow(Tij),export)  # internal rows    of Tij 
    jN  <- setdiff(1:ncol(Tij),import)  # internal columns of Tij


    # Total internal flows, externals removed.
    Tint        <- Tij
    if (! is.null(export)) Tint <- Tint[-export,]
    if (! is.null(import)) Tint <- Tint[,-import]     

    # Total flows, including flow to/from externals
    FlowFrom   <- colSums(Tij)
    FlowTo     <- rowSums(Tij)
    FlowFromC  <- FlowFrom[jN]    # just the total from internal compartments
    FlowToC    <- FlowTo  [iN]
 
    return(list(Tint=Tint,iN=iN,jN=jN,
                import=import,export=export,
                FlowFrom=FlowFrom,
                FlowTo  = FlowTo,
                FlowFromC=FlowFromC,
                FlowToC  =FlowToC))

    } # END InternalNetwork


 
##################################################################################

    Ascendfun <- function(Tij,irow,icol)

    {
    FlowFrom   <- colSums(Tij)
    FlowTo     <- rowSums(Tij)

    # Sum of all flows, including externals
    Throughput  <- sum(Tij)
    
    Asc   <- 0  # The ascendency of network
    Overh <- 0  # The overhead   of network
    Cap   <- 0  # The capacity   of network

    for (i in irow)   {
          for (j in icol){
            if(Tij[i,j]>0) {
                Asc   <- Asc   + Tij[i,j]*log2(Tij[i,j]*Throughput/FlowTo[i]/FlowFrom[j])
#               Overh <- Overh - Tij[i,j]*log2(Tij[i,j]*Tij[i,j]  /FlowTo[i]/FlowFrom[j])
                Cap   <- Cap   - Tij[i,j]*log2(Tij[i,j]/Throughput)     
                           }
                         }
                       }
                       Overh <- Cap-Asc
    return (c(Asc,Overh,Cap))
    }


##################################################################################
#_________________________________________________________________________________
#   THE ASCENDENCY SUITE
#_________________________________________________________________________________

AscInd <- function (Flow = NULL,                      # from-to 
                        Tij  = t(Flow),                   # to-from
                        Import =NULL,                     # flow from external (colNr Tij)
                        Export =NULL,                     # flow to external (colNr Tij)
                        Dissipation=NULL)                 # external dissipation flow
{                      

#------------------------------------------------------------------------
# Flow is a matrix with Flow[i,j] the flow from i (row) to j (column)
# Tij[i,j] contains flow from j to i
# note: component position in rows and columns must be the same - not checked
#------------------------------------------------------------------------
    N          <- InternalNetwork (Tij,Import,Export) 

    # number of columns and rows for total flows (Tij)
    ncTij       <- ncol(Tij)
    nrTij       <- nrow(Tij)

    # THROUGHPUT = Sum of all flows, including externals
    Throughput  <- sum(Tij)

###################################################################
# THE ASCENDENCY Function:                                      #
#   Calculates Ascendency, Overhead and Capacity                #
#     for list of rows (irow) and columns (icol) in Tij           #
###################################################################

#_________________________________________________________________________________

    Ascend           <- matrix (nrow=5,ncol=4)
    colnames(Ascend) <- c("Ascendency","Overhead","Capacity","ACratio")
    rownames(Ascend) <- c("Total","Internal","Import","Export","Dissipation")

    # ascendency for the entire network (nrTij: number of rows in Tij)
    Ascend[1,1:3] <- Ascendfun(Tij,1:nrTij,1:ncTij)

    # ascendency for the internal network (iN: indices to internal rows in Tij)
    Ascend[2,1:3] <- Ascendfun(Tij,N$iN,N$jN)

    # ascendency for the import (Import: index to external column in Tij)
    Ascend[3,1:3] <- Ascendfun(Tij,1:nrTij,N$import)

    # ascendency for the export
    Export1       <- setdiff(N$export,Dissipation)       # true export flow    
    Ascend[4,1:3] <- Ascendfun(Tij,Export1,1:ncTij)

    # ascendency for the dissipation
    Ascend[5,1:3] <- Ascendfun(Tij,Dissipation,1:ncTij)

    # Ascendency - Capacity Ratio
    Ascend[,4]    <-  Ascend[,1]/  Ascend[,3]

    return(Ascend)
    }
    
##################################################################################
#_________________________________________________________________________________
# THE AMI SUITE - AMI, STATISTICAL UNCERTAINTY (HR), CONDITIONAL UNCERTAINTY (DR)
# REALIZED UNCERTAINTY (AMI/HR)
#_________________________________________________________________________________
# CONSTRAINT EFFICIENCY SUITE
#_________________________________________________________________________________
##################################################################################

UncInd <- function (Flow = NULL,                      # from-to 
                    Tij  = t(Flow),                   # to-from
                    Import =NULL,                     # flow from external (colNr Tij)
                    Export =NULL)                     # flow to external (colNr Tij)
{                      

#------------------------------------------------------------------------
# Flow is a matrix with Flow[i,j] the flow from i (row) to j (column)
# Tij[i,j] contains flow from j to i
# note: component position in rows and columns must be the same - not checked
#------------------------------------------------------------------------
    N          <- InternalNetwork (Tij,Import,Export) 

    # number of columns and rows for total flows (Tij)
    ncTij       <- ncol(Tij)
    nrTij       <- nrow(Tij)
    ncomp       <- ncol(N$Tint)          # the number of compartments (without externals)
    compNames   <- rownames(N$Tint)      # names of internal compartment
 
     # Sum of all flows, including externals
    Throughput  <- sum(Tij)
 
#____________________________________________________
# AVERAGE MUTUAL INFORMATION
#____________________________________________________
    # ascendency for the entire network (nrTij: number of rows in Tij)
    Ascend <- Ascendfun(Tij,1:nrow(Tij),1:ncol(Tij))

    AMI        <- Ascend[[1]] / Throughput

#____________________________________________________
# STATISTICAL UNCERTAINTY
#____________________________________________________
    Q          <- N$FlowFrom/Throughput
    Q          <- Q[Q>0]
    HR         <- -sum(Q*log2(Q))

#____________________________________________________
# CONDITIONAL UNCERTAINTY INDEX
#____________________________________________________
    DR         <- HR - AMI

#____________________________________________________
# REALIZED UNCERTAINTY INDEX
#____________________________________________________
    RU         <- AMI/HR

#____________________________________________________
# MAXIMUM UNCERTAINTY (Hmax)
#____________________________________________________
    Hmax <- ncomp*log2(nrTij)

#____________________________________________________
# CONSTRAINT INFORMATION (Hc)
#____________________________________________________

    blsum <- 0
    for (i in 1:nrTij)    # all the rows (including externals)
     { for (j in N$jN)    # only internal columns
        { if (Tij[i,j]>0) blsum <- blsum + 
                             (Tij[i,j]/N$FlowFrom[j])*log2(Tij[i,j]/N$FlowFrom[j])
        }}

    Hc <- Hmax + blsum
    names(Hc)<-NULL
#____________________________________________________
# CONSTRAINT EFFICIENCY
#____________________________________________________
    
    CE <- Hc/Hmax
#____________________________________________________
# NETWORK EFFICIENCY
#____________________________________________________
    
    Hsys <- Hmax - Hc

    list(AMI=AMI, HR=HR, DR=DR, RU=RU,Hmax=Hmax,Hc=Hc,
     Hsys=Hsys,CE=CE)

    }
##################################################################################
#_________________________________________________________________________________
# THE EFFECTIVE MEASURES SUITE
#_________________________________________________________________________________
##################################################################################
    
EffInd<- function (Flow = NULL,                      # from-to 
                   Tij  = t(Flow),                   # to-from
                   Import =NULL,                     # flow from external (colNr Tij)
                    Export =NULL)                     # flow to external (colNr Tij)
{                      

#------------------------------------------------------------------------
# Flow is a matrix with Flow[i,j] the flow from i (row) to j (column)
# Tij[i,j] contains flow from j to i
# note: component position in rows and columns must be the same - not checked
#------------------------------------------------------------------------
    N          <- InternalNetwork (Tij,Import,Export) 

     # Sum of all flows, including externals

    Throughput  <- sum(Tij)
    ncomp       <- ncol(N$Tint)          # the number of compartments (without externals)
#_________________________________________________________________________________
 

    CZ <-1    # EFFECTIVE CONNECTANCE CZ
    FZ <-1    # EFFECTIVE FLOWS       FZ
    NZ <-1    # EFFECTIVE NODES       NZ
    RZ <-1    # EFFECTIVE ROLES       RZ

    for (i in 1:ncomp) 
      {if (N$FlowToC[i]>0)
         {for (j in 1:ncomp) 
           {if (N$FlowFromC[j]>0)
            NZ <-NZ*(Throughput^2        /N$FlowFromC[j]/N$FlowToC[i])^( 0.5*N$Tint[i,j]/Throughput)
            } # end j
          } # end if flowtoC
       } # end i
    for (i in 1:ncomp) 
       { for (j in 1:ncomp)
          { if (N$Tint[i,j]>0)
            {             
            CZ <-CZ*(N$Tint[i,j]^2         /N$FlowFromC[j]/N$FlowToC[i])^(-0.5*N$Tint[i,j]/Throughput)
            RZ <-RZ*(N$Tint[i,j]*Throughput/N$FlowFromC[j]/N$FlowToC[i])^(     N$Tint[i,j]/Throughput)
            FZ <-FZ*(N$Tint[i,j]           /Throughput)^(-N$Tint[i,j]/Throughput)
             }
          }
        }   
        names(CZ)<-names(RZ)<-names(NZ)<-names(FZ)<-NULL
    list(CZ=CZ,FZ=FZ,NZ=NZ,RZ=RZ)
   }

##################################################################################
#________________________________________________________________________________
# FINN's SUITE   - Pathway analysis
#________________________________________________________________________________
##################################################################################

PathInd <- function (Flow = NULL,                      # from-to 
                     Tij  = t(Flow),                   # to-from
                     Import =NULL,                     # flow from external (colNr Tij)
                     Export =NULL)                     # flow to external (colNr Tij)
{                      

#------------------------------------------------------------------------
# Flow is a matrix with Flow[i,j] the flow from i (row) to j (column)
# Tij[i,j] contains flow from j to i
# note: component position in rows and columns must be the same - not checked
#------------------------------------------------------------------------
    N          <- InternalNetwork (Tij,Import,Export) 

    # Rate of change of each compartment
    RateComp   <- N$FlowToC-N$FlowFromC

    # number of columns and rows for total flows (Tij)
    ncTij       <- ncol(Tij)
    nrTij       <- nrow(Tij)
    ncomp       <- ncol(N$Tint)              # the number of compartments (without externals)
    compNames   <- rownames(N$Tint)

    ExportSum   <- sum(N$FlowTo[N$export])
    ImportSum   <- sum(N$FlowFrom[N$import])

    #___________________________________________________________
    # THROUGHFLOW - Rates of change taken into consideration
    # Throughflow based on rows and columns has to be the same

    Throughflow <- sum(N$Tint) + ImportSum - sum(RateComp[RateComp<0])
    
    #___________________________________________________________
    # THROUGHPUT = Sum of all flows, including externals

    Throughput  <- sum(Tij)

#____________________________________________________
# THE PATHLENGTH
#____________________________________________________
    Pathlength   <- Throughflow/ (ExportSum + sum(RateComp[RateComp>0]))

#____________________________________________________
# TRANSITIVE CLOSURE MATRIX
#____________________________________________________

    CompThroughflow <- pmax(N$FlowFromC,N$FlowToC)  # total flow, internal compartments    
    Qij <- matrix(nrow=ncomp,ncol=ncomp,0)
  
    for (i in 1:ncomp) Qij[i,]  <- N$Tint[i,]/ CompThroughflow[i]

    diagnl    <- diag(nrow=ncomp,1)   # unity matrix
    IQ        <- diagnl-Qij
    M         <- ginv(IQ)           # The generalised inverse

    # Cycled throughflow 
    diaM      <- diag(M)
    TSTC      <- sum((1-1/diaM)*N$FlowFromC)

    # Noncycled throughflow
    TSTS      <- Throughflow - TSTC

    # Finn's cycling index
    FCI       <- TSTC/Throughflow

    # Finn's cycling index revisited
    FCIb      <- TSTC/Throughput

    Finn<-list(TSTC=TSTC,TSTS=TSTS,FCI=FCI,FCIb=FCIb, APL=Pathlength)
    return(Finn)
    }
  
##################################################################################
#________________________________________________________________________________
#   THE ENVIRONS SUITE
#________________________________________________________________________________
    
##################################################################################

EnvInd <- function (Flow = NULL,                      # from-to 
                        Tij  = t(Flow),                   # to-from
                        Import =NULL,                     # flow from external (colNr Tij)
                        Export =NULL,                     # flow to external (colNr Tij)
                        full=FALSE)                       # if True, also returns matrices
{                      

#------------------------------------------------------------------------
# Flow is a matrix with Flow[i,j] the flow from i (row) to j (column)
# Tij[i,j] contains flow from j to i
# note: component position in rows and columns must be the same - not checked
#------------------------------------------------------------------------
    N          <- InternalNetwork (Tij,Import,Export) 

    # number of columns and rows for total flows (Tij)
    ncTij       <- ncol(Tij)
    nrTij       <- nrow(Tij)
    ncomp       <- ncol(N$Tint)              # the number of compartments (without externals)
    compNames   <- rownames(N$Tint)
    ImportSum   <- sum(N$FlowFrom[N$import])
    RateComp    <- N$FlowToC-N$FlowFromC
#____________________________________________________
# Network Aggradation - ratio of internal disorder to generated disorder
#____________________________________________________
    Throughflow  <- sum(N$Tint) + ImportSum - sum(RateComp[RateComp<0])
    TotalOutflow <- sum(Tij[N$export,N$jN])
    Naggr        <- Throughflow/TotalOutflow

#____________________________________________________
# Homogenization Index - measure of evenness of flows
#____________________________________________________
# G1: transitive closure matrix
    G   <- matrix(nrow=ncomp,ncol=ncomp,0)
    for (j in 1:ncomp) G[,j] <- N$Tint[,j]/N$FlowFromC[j]
# N1: Integral nondimensional matrix    
    N1 <- NULL
    if (any(is.nan(G))) stop("Some of the flows are NANs")
    I   <- diag(nrow=ncomp,1)
    IG  <- I-G
    N1  <- ginv(IG)

    # coefficient of variation of G and N1 
    # note: sd of a matrix takes standard deviation per column - overrule by as.vector
    MG  <- mean(G)
    MN  <- mean(N1) 
    CVG <- sd(as.vector(G))/MG
    CVN <- sd(as.vector(N1))/MN
    
    HP  <- CVG/CVN  # Homogenization index
#____________________________________________________
# Dominance of Indirect effects index- ratio i/d
#____________________________________________________

    D1  <- N$Tint-t(N$Tint)
    ID  <- (sum(N1-I-G))/sum(G)    #delta(ij)=1 diagonal elements!


#____________________________________________________
# Utility Matrix (UP)
#____________________________________________________

    D     <- D1/N$FlowToC
    U     <- ginv(I-D)
# synergism index - was wrong in Latham 2006;
# here is the correct(?) value:
    Gamma <- N$FlowToC*U
    bc    <- sum(Gamma[Gamma>0])/abs(sum(Gamma[Gamma<0]))
#    rownames(U)  <-rownames(Tij[-((nrTij-1):nrTij),])
#    colnames(U) <- rownames(U)

        
##################################################################################
  if (full)
  {
    rownames(U) <- rownames(N1) <-rownames(G) <-compNames    
    colnames(U) <- colnames(N1) <-colnames(G)<- compNames    
    
    Res <- list(NAG=Naggr,HP=HP,BC=bc,ID=ID,CVN=CVN,CVG=CVG, U=U,N1=N1,G=G) 
    
    
    }else  Res<-list(NAG=Naggr,HP=HP,BC=bc,ID=ID,MN=MN,MG=MG,CVN=CVN,CVG=CVG)
   
   return(Res)
    }

##################################################################################
#________________________________________________________________________________
# Internal function: estimates the diet composition
#________________________________________________________________________________

Diet <- function (Tint,                   # Calculates diet composition
                  Dead=NULL,              # index from Dead to Tint
                  iN=1:nrow(Tint))      
##################################################################################
#________________________________________________________________________________
# TROPHIC ANALYSIS
#________________________________________________________________________________

{
#____________________________________________________
# p matrix contains the diet composition of predator i
#____________________________________________________

    IntFlowTo   <- rowSums(Tint)  # Total food ingested

    p           <- Tint
    for (i in 1:ncol(Tint)) p[i,] <- Tint[i,]/IntFlowTo[i]
    p[is.na(p)] <- 0

# take into account dead matter; Dead refers to column/row in Tij
# N$iN maps from Tint to Tij
    if (! is.null(Dead)) p[which(iN%in% Dead),]<-0
    return(p)

} # END Diet   

 
#############################################################################
# Trophic analysis function
#############################################################################

TrophInd <- function (Flow = NULL,                      # from-to 
                      Tij  = t(Flow),                   # to-from
                      Import =NULL,                     # flow from external (colNr Tij)
                      Export =NULL,                     # flow to external (colNr Tij)
                      Dead=NULL)                        # flow to dead matter 
 
 

{
    if (is.character(Dead)) 
       dead <- which(rownames(Tij)%in%Dead) else
       dead <- Dead
    if (length(dead) != length(Dead)) stop("Dead not recognized")

# Check input and calculate internal network
    N          <- InternalNetwork (Tij,Import,Export) 
    p          <- Diet(N$Tint,dead,N$iN)
    ncomp      <- ncol(N$Tint)              # the number of compartments (without externals)

#____________________________________________________
# Trophic level TL  TL(i) = 1 + sumj(pij*TL(j))
#                   TL(i)-sumj(pij*TL(j))= 1 
#____________________________________________________
    A       <- -p
    diag(A) <- diag(A)+ 1
    B       <- rep(1,ncomp)
    TL      <- ginv(A) %*% B                   
#____________________________________________________
# Omnivory index: variance of trophic levels of preys
#____________________________________________________

    OI <- vector (length = ncomp)
    for (i in 1:ncomp) OI[i] <-sum((TL-(TL[i]-1))^2*p[i,])


    return(data.frame(TL, OI,row.names=rownames(N$Tint)))
 }               ## END Trophic
 
 
#############################################################################
# Dependency analysis function - dependence on other sources
#############################################################################

Dependency <- function (Flow = NULL,                      # from-to 
                        Tij  = t(Flow),                   # to-from
                        Import =NULL,                     # flow from external (colNr Tij)
                        Export =NULL)                     # flow to external (colNr Tij)
 
##################################################################################
#________________________________________________________________________________
# DEPENDENCY ANALYSIS
#________________________________________________________________________________

{

# Check input and calculate internal network
    N          <- InternalNetwork (Tij,Import,Export) 
    feed       <- Diet(N$Tint)
    
    D          <- solve(diag(nrow=nrow(feed))-feed)        # (I-feed)^-1
          
    return(D)
 }               ## END Dependency
 
  