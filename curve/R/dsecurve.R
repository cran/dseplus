#######################################################################

#######################################################################

#              Time series model curvature calculation

#######################################################################

# S routines for calulating curvatures a la Bates and Watts.

# Notes:
#   Generating the D matrix can be computationally demanding. 
#   The code here uses a compiled version which is fast and
#   does not suffer from the memory problem, but works only with ARMA 
#   and KF models.
#    
#-----------------------------------------------------------------------

curvature.TSestModel <- function (func, x=coef(func),
     method="Richardson", method.args=list(d=0.01, eps=1e-4, r=6, v=2),
     compiled=TRUE, warn=TRUE,
     Shape=TSmodel(func), data=TSdata(func), ...)  
 {if(compiled) 
   curvature(genD(func, x=x,
                  method=method, method.args=method.args,
		  Shape=Shape, data=data, ...),
             warn=warn)
  else
   {func.residual <- function(coefficients,Shape,data)
                      {c(l(setArrays(Shape,coefficients=coefficients),
		           data,result="pred")
                       - outputData(data))}
    curvature.default(func.residual, coef(func), 
                       method="Richardson", method.args=method.args,
                       warn=warn,
		       Shape=Shape, data=data, ...) 
   }
 }

#######################################################################

#                               D matrix calculation

#######################################################################

 
genD.TSestModel <- function(func, x=coef(func),
     method="Richardson", method.args=list(d=0.01, eps=1e-4, r=6, v=2),
     Shape=TSmodel(func), data=TSdata(func), ...)
   { invisible(genD(TSmodel(func), x=x, method=method, method.args=method.args,
                       Shape=Shape, data=data, ...))}

genD.ARMA <- function(func, x=coef(func),
     method="Richardson", method.args=list(d=0.01, eps=1e-4, r=6, v=2),
     Shape=TSmodel(func), data=TSdata(func), ...){
# Note: m,n,p have different meanings here than they do in 
#  time-series models! ms,ns,ps are use for time-series meaning
   if (method != "Richardson") stop("method not implemented.")
   d <- method.args$d
   eps <- method.args$eps
   r <- method.args$r
   if(2 != method.args$v) stop("current code requires v=2 (default)")
   model <- Shape
   n <-length(c(outputData(data))) # this has the same length as the residual
   sampleT <-periods(data) 
   ps <-dim(model$A)[3]
   ms <-dim(model$C)[3]
   #ns <-ps  # this could be 1 except for TREND 
   TREND <- if (is.null(model$TREND)) matrix(0,periods(data), ps) else model$TREND
   is  <- max(ms,ps)
   loc   <- match(model$location,       c("A","B","C","t"))
   cloc  <- match(model$const.location, c("A","B","C","t"))
   if(is.null(ms))
     {ms<-0
      C<-array(0,c(1,ps,1)) # can't call compiled with 0 length arrays
      u <- matrix(0,sampleT,1)
     }
   else
     {C <-model$C
      u <- inputData(data)
     } 
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
   storage.mode(x)     <-"double"
   storage.mode(h0)     <-"double"
   storage.mode(D)     <-"double"
   storage.mode(Haprox)<-"double"
   storage.mode(Hdiag )<-"double"
   storage.mode(Daprox)<-"double"
   storage.mode(u)<-"double"
   storage.mode(outputData(data))<-"double"
   storage.mode(model$const)<-"double"
   storage.mode(TREND)<-"double"
   D <-.Fortran("genda",
            D=D,
            p=as.integer(length(x)),
            x=x,
            h0,
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
            u=u, 
            y=outputData(data),   
            #   model$parm is passed above as x (it is the parameter for curvature calculation)   
            as.integer(loc),   #as.character(model$location), #23
            as.integer(model$i),
            as.integer(model$j),
            as.integer(length(model$const)),
            const=model$const,
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
            CONSTNT=TREND, 
        #  scratch
            as.integer(is),  
            AA=matrix(double(1),is,is),  
            BB=matrix(double(1),is,is),  
            WW=rep(double(1),is),  
            integer(is*is),         # scratch array IPIV
            DUP=.DSEflags()$DUP,
	    PACKAGE="dse1"
	    )[c("D","p","f0", "x", "r")] 
   D$d   <- d
   D$method      <- method
   D$method.args <- method.args
   # D calculation can be done relative to any point (subtracting data does
   #   not matter) but curvature calculation assumes f0 is really a residual.
   D$f0 <- l(model, data, result="pred") - c(outputData(data))
   class(D) <- "Darray"  #constructor
   invisible(D) 
}

#Rbug this does not seem to work as genD.KF or genD.KF.innov

genD.innov <- function(func, x=coef(func),
     method="Richardson", method.args=list(d=0.01, eps=1e-4, r=6, v=2),
     Shape=TSmodel(func), data=TSdata(func), ...){
# Note: m,n,p have different meanings here than they do in 
#  time-series models! ms,ns,ps are use for time-series meaning
   if (method != "Richardson") stop("method not implemented.")
   d <- method.args$d
   eps <- method.args$eps
   r <- method.args$r
   if(2 != method.args$v) stop("current code requires v=2 (default)")
   model <- Shape
   n <-length(c(outputData(data))) # this has the same length as the residual
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
      u <- inputData(data)
      G <- matrix(0,ns,ms)
     } 
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
   loc   <- match(model$location,       c("f","G","H","K","Q","R","z","P"))
   cloc  <- match(model$const.location, c("f","G","H","K","Q","R","z","P"))
   storage.mode(D)     <-"double"
   storage.mode(x)     <-"double"
   storage.mode(h0)     <-"double"
   storage.mode(Daprox)<-"double"
   storage.mode(Haprox)<-"double"
   storage.mode(Hdiag)<-"double"
   storage.mode(u)<-"double"
   storage.mode(outputData(data))<-"double"
   storage.mode(model$const)<-"double"
   storage.mode(C)<-"double"
   storage.mode(G)<-"double"
   IS <- max(ns,ps)
   D <-.Fortran("gendk",
            D=D,
            p=as.integer(length(x)),
            x=x,
            h0,
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
            u=u, 
            y=outputData(data),   #22
            #   model$parm is passed above as x (it is the parameter for curvature calculation)   
            as.integer(loc),   #as.character(model$location) bug
            as.integer(model$i),
            as.integer(model$j),
            as.integer(length(model$const)),
            const=model$const,     #28
            as.integer(cloc),  #as.character(model$const.location),
            as.integer(model$const.i ), 
            as.integer(model$const.j),
        #  for state space models:
            as.integer(ns),  # state dim.     
        #    state=matrix(double(1),sampleT,ns),  
        #    track=array(double(1),c(sampleT,ns,ns)),  
            z0=double(ns),  
            P0=diag(double(1),ps),
            F=matrix(double(1),ns,ns),  
            G=G,  
            H=matrix(double(1),ps,ns),  
            K=matrix(double(1),ns,ps),  
            Q=matrix(double(1),ns,ns),  
            R=matrix(double(1),ps,ps),
            gain=as.integer(is.innov.SS(model)),
	    as.integer(IS),           # scratch arrays for KF, IS
	    matrix(double(1),IS,IS),  #A
	    matrix(double(1),IS,IS),  # AA
	    matrix(double(1),IS,IS),  # PP
	    matrix(double(1),ns,ns),  # QQ
	    matrix(double(1),ps,ps),  # RR 
	    rep(double(1),IS),  # Z
	    rep(double(1),IS), # ZZ
	    rep(double(1),IS), # WW		   
            integer(IS*IS),         # scratch array IPIV
            DUP=.DSEflags()$DUP,
	    PACKAGE="dse1"
	    )[c("D","p","f0", "x", "r")]
   D$d   <- d
   D$method      <- method
   D$method.args <- method.args
   # D calculation can be done relative to any point (subtracting data does
   #   not matter) but curvature calculation assumes f0 is really a residual.
   D$f0 <- l(model, data, result="pred") - c(outputData(data))
   class(D) <- "Darray" #constructor
   invisible(D)
}

#######################################################################

#              span (dimension of tangent space) calculation

#######################################################################

           
span.TSestModel <- function (func, x=coef(func),
     method="Richardson", method.args=list(d=0.01, eps=1e-4, r=6, v=2),
     show.details=FALSE, compiled=.DSEflags()$COMPILED,
     Shape=TSmodel(func), data=TSdata(func), ...){

  # calculate the singular values of the tangents
  # the compiled version calculates the whole D matrix (which seems like
  # a waste, but the compiled code is much faster, so ...
  if (compiled)
   {D <- genD(func, x, method=method, method.args=method.args,
              Shape, data, ...)$D[,seq(length(coef(func))),drop=FALSE]
    if (any(is.na(D))) {
       # really should stop here
       warning("D from compiled genD contains NAs. Setting them to zero. Result is probably not valid!!!")
       D[is.na(D)] <- 0
      }
    if (any(is.nan(D))) {
       warning("D from compiled genD contains NANs. Setting them to zero. Result is probably not valid!!!")
       D[is.nan(D)] <- 0
      }
    return(svd(D)$d)
   }
  else {
     funcTS <- function(coefficients, Shape,data)
      {c(l(setArrays(Shape, coefficients=coefficients),data,result="pred")
           - outputData(data))}
     span.default(funcTS, x, method=method, method.args=method.args,
              show.details=show.details, Shape, data, ...)
     }
 }


#######################################################################

#            calculate Fisher Info (Hessian of the likelihood)

#######################################################################

hessian.TSestModel <- function (func, x=coef(func), 
     method="Richardson", method.args=list(d=0.01, eps=1e-4, r=6, v=2),
     Shape=TSmodel(func), data=TSdata(func), ...)  
 {# like returns neg. log likelihood
  funcTS <- function(coefficients, Shape, data)
  	{l(setArrays(Shape,coefficients=coefficients), data, result="like")}
  # call hessian.default
  hessian(funcTS, x, method=method, method.args=method.args,
           Shape=Shape, data=data, ...)
 }
 

#######################################################################

#                    end

#######################################################################

