#######################################################################

# This package uses curve.R too.

#######################################################################

#                            curvature calculation

#######################################################################

curvature.TSestModel <- function (obj, warn=T, ...)  
 {curvature(genD(obj, ...),  warn=warn)}

#-----------------------------------------------------------------------

# S routines for calulating curvatures a la Bates and Watts.

# Notes:
#   Generating the D matrix can be computationally demanding. There are three
#   (or four) versions of the code for generating this matrix. The S version is slow 
#   and suffers from memory problems due to a bug in the way S
#   (S-PLUS Version 3.0 for SUN4) allocates memory in loops. The C version is
#   faster but suffers (even worse) from the memory problem.  
#   Both the S and the C versions take the name
#   of an S function as an arguement and call it. The compiled version is fast and
#   does not suffer from the memory problem, but works only with ARMA 
#   and KF models.
#    
#-----------------------------------------------------------------------

#######################################################################

#                               D matrix calculation

#######################################################################

genD.TSestModel <- function(obj, ...)
   { invisible(genD( TSmodel(obj), TSdata(obj), ...))}

#genD.TSmodel <- function(obj, data, ...) {NextMethod("genD.TSmodel")}

genD.ARMA <- function(obj, data, d=0.01, eps=1e-4, r=6, warn=F){
# Note: m,n,p have different meanings here than they do in 
#  time-series models! ms,ns,ps are use for time-series meaning
   model <- obj
   n <-length(c(output.data(data))) # this has the same length as the residual
   sampleT <-periods(data) 
   ps <-dim(model$A)[3]
      ms <-dim(model$C)[3]
      ns <-ps  # this could be 1 except z0 is used for TREND and F, Q and R for scratch
     loc   <- match(model$location,       c("A","B","C","t"))
     cloc  <- match(model$const.location, c("A","B","C","t"))
   if(is.null(ms))
     {ms<-0
      C<-array(0,c(1,ns,1)) # can't call compiled with 0 length arrays
      u <- matrix(0,sampleT,1)
      G <- matrix(0,ns,1)
     }
   else
     {C <-model$C
      u <- input.data(data)
      G <- matrix(0,ns,ms)
     } 
   x <-model$parms
   zt <- 0.0001
   #Rbug zt <-unix.time(l(model,data,sampleT,sampleT,result="like"))
   zt <- zt[1] *r*2*(length(x)*(length(x) + 3))/2
   if(zt>30) 
       cat("D matrix calculation very roughly estimated to take about ",
            (10*ceiling(zt/10)),"seconds(without other system activity)...\n")
   h0 <- abs(d*x)+eps*(x==0.0)
   D <- matrix(1e20, n,(length(x)*(length(x) + 3))/2) 
   Daprox <-  matrix(0,n,r) 
   Hdiag  <-  matrix(0,n,length(x))
   Haprox <-  matrix(0,n,r)
   storage.mode(D)     <-"double"
   storage.mode(Daprox)<-"double"
   storage.mode(Haprox)<-"double"
   storage.mode(Hdiag )<-"double"
   D <-.Fortran("gend",
            D=D,
            as.integer(is.ARMA(model)), 
            p=as.integer(length(x)),
            x=as.double(x),
            as.double(h0),
            as.integer(n),    #6
            as.integer((length(x)*(length(x)+ 3))/2), #cols of D
            f0=double(n),      
            r=as.integer(r),
            #                       work space for GEND
            Haprox=Haprox,     #10   
            Hdiag=Hdiag,         
            Daprox=Daprox,        
            x=double(length(x)),
            delta=double(length(x)),
            f1=double(n),
            f2=double(n),
            #                   work space for ARMAp / KFp
         #     cov=matrix(1e20,ps,ps),       
            # pred is f0,f1,f2 passed above
            as.integer(ms),     #  input dimension m  #17
            as.integer(ps),     # output dimension p 
            as.integer(sampleT),   
            as.integer(periods(data)), 
            u=as.double(u), 
            y=as.double(output.data(data)),   
            #   model$parm is passed above as x (it is the parameter for curvature calculation)   
            as.integer(loc),   #as.character(model$location), #23
            as.integer(model$i),
            as.integer(model$j),
            as.integer(length(model$const)),
            const=as.double(model$const),
            as.integer(cloc),  #as.character(model$const.location),
            as.integer(model$const.i ), 
            as.integer(model$const.j),
        #  for ARMA models:
            as.integer(model$l),       #31
            as.integer(model$const.l),
            as.integer( dim(model$A)[1]),#1+order of A
            as.integer( dim(model$B)[1]),#1+order of B
            as.integer( dim(C)[1]),#1+order of C
            A=model$A,   
            B=model$B,  
            C=C, 
        #  for state space models:
            as.integer(ns),  # state dim.  #39
        #    state=as.double(matrix(0,sampleT,ns)),  
        #    track=as.double(array(0,c(sampleT,ns,ns))),  
            z0=as.double(rep(0,ps)), # note: this is TREND for ARMA
            P0=as.double(diag(0,ps)),
            F=as.double(matrix(0,ns,ns)),  
            G=as.double(G),  
            H=as.double(matrix(0,ps,ns)),  
            K=as.double(matrix(0,ns,ps)),  
            Q=as.double(matrix(0,ns,ns)),  
            R=as.double(matrix(0,ps,ps)), 
            gain=as.integer(F),   #48
            DUP=.DSEDUP)[c("D","p","f0", "x", "r")] 
   D$d   <- d
   D$eps <- eps
   # D calculation can be done relative to any point (subtracting data does
   #   not matter) but curvature calculation assumes f0 is really a residual.
   D$f0 <- l(model, data, result="pred") - c(output.data(data))
   invisible(classed(D,"Darray")) #constructor
}

#Rbug this does not seem to work as genD.KF or genD.KF.innov

genD.innov <- function(obj, data, d=0.01, eps=1e-4, r=6, warn=F){
# Note: m,n,p have different meanings here than they do in 
#  time-series models! ms,ns,ps are use for time-series meaning
   model <- obj
   n <-length(c(output.data(data))) # this has the same length as the residual
   sampleT <-periods(data) 
   ns <-dim(model$F)[2]
   ps <-dim(model$H)[1]
   ms <-dim(model$G)[2]
   if(is.null(ms))
     {ms<-0
      C<-array(0,c(1,ns,1)) # can't call compiled with 0 length arrays
      u <- matrix(0,sampleT,1)
      G <- matrix(0,ns,1)
     }
   else
     {C <- array(0,c(1,ns,ms))
      u <- input.data(data)
      G <- matrix(0,ns,ms)
     } 
   x <-model$parms
   zt <- 0.0001
   #Rbug zt <-unix.time(l(model,data,sampleT,sampleT,result="like"))
   zt <- zt[1] *r*2*(length(x)*(length(x) + 3))/2
   if(zt>30) 
       cat("D matrix calculation very roughly estimated to take about ",
           (10*ceiling(zt/10)),"seconds(without other system activity)...\n")
   h0 <- abs(d*x)+eps*(x==0.0)
   D <- matrix(1e20, n,(length(x)*(length(x) + 3))/2) 
   Daprox <-  matrix(0,n,r) 
   Hdiag  <-  matrix(0,n,length(x))
   Haprox <-  matrix(0,n,r)
   storage.mode(D)     <-"double"
   storage.mode(Daprox)<-"double"
   storage.mode(Haprox)<-"double"
   storage.mode(Hdiag )<-"double"
   loc   <- match(model$location,       c("f","G","H","K","Q","R","z","P"))
   cloc  <- match(model$const.location, c("f","G","H","K","Q","R","z","P"))
   D <-.Fortran("gend",
            D=D,
            as.integer(is.ARMA(model)), 
            p=as.integer(length(x)),
            x=as.double(x),
            as.double(h0),
            as.integer(n),    #6
            as.integer((length(x)*(length(x)+ 3))/2), #cols of D
            f0=double(n),      
            r=as.integer(r),
            #                       work space for GEND
            Haprox=Haprox,        
            Hdiag=Hdiag,         
            Daprox=Daprox,        
            x=double(length(x)),
            delta=double(length(x)),
            f1=double(n),                # 15
            f2=double(n),
            #                   work space for ARMAp / KFp       
            # cov=matrix(1e20,ps,ps),       
            # pred is f0,f1,f2 passed above
            as.integer(ms),     #  input dimension m  
            as.integer(ps),     # output dimension p 
            as.integer(sampleT),   
            as.integer(periods(data)), 
            u=as.double(u), 
            y=as.double(output.data(data)),   #22
            #   model$parm is passed above as x (it is the parameter for curvature calculation)   
            as.integer(loc),   #as.character(model$location) bug
            as.integer(model$i),
            as.integer(model$j),
            as.integer(length(model$const)),
            const=as.double(model$const),     #28
            as.integer(cloc),  #as.character(model$const.location),
            as.integer(model$const.i ), 
            as.integer(model$const.j),
        #  for ARMA models:
            as.integer(loc),                   #31
            as.integer(loc),
            as.integer(1),#1+order of A
            as.integer(1),#1+order of B
            as.integer(1),#1+order of C
            A=as.double(array(0,c(1,ps,ps))),   
            B=as.double(array(0,c(1,ps,ps))),  
            C=as.double(C), 
        #  for state space models:
            as.integer(ns),  # state dim.     #39
        #    state=as.double(matrix(0,sampleT,ns)),  
        #    track=as.double(array(0,c(sampleT,ns,ns))),  
            z0=as.double(rep(0,ns)),  
            P0=as.double(diag(0,ps)),
            F=as.double(matrix(0,ns,ns)),  
            G=as.double(G),  
            H=as.double(matrix(0,ps,ns)),  
            K=as.double(matrix(0,ns,ps)),  
            Q=as.double(matrix(0,ns,ns)),  
            R=as.double(matrix(0,ps,ps)),
            gain=as.integer(is.innov.SS(model)), #48
            DUP=.DSEDUP)[c("D","p","f0", "x", "r")]
   D$d   <- d
   D$eps <- eps
   # D calculation can be done relative to any point (subtracting data does
   #   not matter) but curvature calculation assumes f0 is really a residual.
   D$f0 <- l(model, data, result="pred") - c(output.data(data))
   invisible(classed(D, "Darray")) #constructor
}

#######################################################################

#              span (dimension of tangent space) calculation

#######################################################################

span.TSestModel <- function (obj, compiled=.DSECOMPILED, ...)  
 {# calculate the singular values of the tangents
  # the compiled version calculates the whole D matrix (which seems like
  # a waste, but the compiled code is much faster, so ...
  if (compiled)
   {D <- genD(obj, ... )$D[,seq(length(parms(obj))),drop=F]
    if (any(is.na(D))) {
       # really should stop here
       warning("D from compiled genD contains NAs. Setting them to zero!!!")
       D[is.na(D)] <- 0
      }
    if (any(is.nan(D))) {
       warning("D from compiled genD contains NANs. Setting them to zero!!!")
       D[is.nan(D)] <- 0
      }
    return(svd(D)$d)
   }
  else {
     # previously used global.assign but that is not necessary 
     func.residual.TSestModel <- function(parms,Shape,data)
      {c(l(set.arrays(Shape,parms=parms),data,result="pred") - output.data(data))}
     #on.exit(remove("func.residual.TSestModel"))
     span.default(func.residual.TSestModel, parms(obj),
               obj.args=list(Shape=TSmodel(obj), data=TSdata(obj), ...))
     }
 }


#######################################################################

#            calculate Fisher Info (Hessian of the likelihood)

#######################################################################

hessian.TSestModel <- function (obj, ...)  
 {# like returns neg. log likelihood
  # previously used global.assign but that is not necessary 
  func.hessian.TSestModel <- function(parms,Shape,data)
  	{l(set.arrays(Shape,parms=parms),data,result="like")}
  #on.exit(rm(func.hessian.TSestModel))
  hessian.default( func.hessian.TSestModel, parms(obj),
   obj.args=list(Shape=TSmodel(obj), data=TSdata(obj), ...))
 }

# was obj.args=append(list(Shape=TSmodel(emodel), data=TSdata(emodel)), list(...)))
#######################################################################

# Test routines for calulating curvatures moved to the tests directory.

#######################################################################
 

#######################################################################

#                    end

#######################################################################

