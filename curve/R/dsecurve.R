#######################################################################

# This package uses curve.R too.

#######################################################################

#                            curvature calculation

#######################################################################

curvature.TSestModel <- function (emodel, warn=T, ...)  
 {curvature(genD(emodel, ...),  warn=warn)}

#-----------------------------------------------------------------------

# S routines for calulating curvatures a la Bates and Watts.

# Notes:
#   Generating the D matrix can be computationally demanding. There are three
#   (or four) versions of the code for generating this matrix. The S version is slow 
#   and suffers from memory problems due to a bug in the way the current version of S
#   (S-PLUS Version 3.0 for SUN4) allocates memory in loops. The C version is
#   faster but suffers (even worse) from the memory problem. (A fix is promised 
#   in the next release of S-plus.) Both the S and the C versions take the name
#   of an S function as an arguement and call it. The compiled version is fast and
#   does not suffer from the memory problem, but works only with ARMA 
#   and KF models.
#    
#-----------------------------------------------------------------------

#######################################################################

#                               D matrix calculation

#######################################################################

genD.TSestModel <- function(estModel, ...)
   { invisible(genD( TSmodel(estModel), TSdata(estModel), ...))}

#genD.TSmodel <- function(model, data, ...) {NextMethod("genD.TSmodel")}

genD.ARMA <- function(model, data, d=0.01, eps=1e-4, r=6, warn=F){
# Note: m,n,p have different meanings here than they do in 
#  time-series models! ms,ns,ps are use for time-series meaning
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
            DUP=TRUE)[c("D","p","f0", "x", "r")] # Rbug DUP=TRUE seems necessary
   D$d   <- d
   D$eps <- eps
   # D calculation can be done relative to any point (subtracting data does
   #   not matter) but curvature calculation assumes f0 is really a residual.
   D$f0 <- l(model, data, result="pred") - c(output.data(data))
   invisible(classed(D,"Darray")) #constructor
}

#Rbug this does not seem to work as genD.KF or genD.KF.innov

genD.innov <- function(model, data, d=0.01, eps=1e-4, r=6, warn=F){
# Note: m,n,p have different meanings here than they do in 
#  time-series models! ms,ns,ps are use for time-series meaning
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
            DUP=TRUE)[c("D","p","f0", "x", "r")]
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

span.TSestModel <- function (emodel, compiled=.DSECOMPILED, ...)  
 {# calculate the singular values of the tangents
  # the compiled version calculates the whole D matrix (which seems like
  # a waste, but the compiled code is much faster, so ...
  if (compiled)
   {D <- genD(emodel, ... )$D[,seq(length(parms(emodel))),drop=F]
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
     span.default(func.residual.TSestModel, parms(emodel),
               func.args=list(Shape=TSmodel(emodel), data=TSdata(emodel), ...))
     }
 }


#######################################################################

#            calculate Fisher Info (Hessian of the likelihood)

#######################################################################

hessian.TSestModel <- function (emodel, ...)  
 {# like returns neg. log likelihood
  # previously used global.assign but that is not necessary 
  func.hessian.TSestModel <- function(parms,Shape,data)
  	{l(set.arrays(Shape,parms=parms),data,result="like")}
  #on.exit(rm(func.hessian.TSestModel))
  hessian.default( func.hessian.TSestModel, parms(emodel),
   func.args=list(Shape=TSmodel(emodel), data=TSdata(emodel), ...))
 }

# was func.args=append(list(Shape=TSmodel(emodel), data=TSdata(emodel)), list(...)))
#######################################################################

# Test routines for calulating curvatures are (being) moved to the tests directory.

#######################################################################
 

# following from dse1 tests are much slower

dsecurvature.function.testsC <- function( verbose=T, synopsis=T, fuzz.small=1e-12, 
    fuzz.large=1e-6, show.details=F, show.extra.details=F,
    print.values=T, digits=18)
{# Tests of DSE curvature functions
  if      (is.R()) data("eg1.DSE.data.diff", package="dse1")
  else if (is.S()) source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))

  all.ok <- list(ok=T, error=NA)

test <- function(tst, good, all.ok, flag="", fuzz=1e-14, print.values=F)
  {# initialize all.ok <- list(ok=T, error=NA)
   error <- max(abs(good-tst))
   if (fuzz < error) 
     {all.ok$error <- if(is.na(all.ok$error)) error else max(error, all.ok$error)
      all.ok$ok <- F  # assumes ok is initialized T for this to work 
      cat(flag, "failed! error = ", error,"\n")
      if (print.values) print.test.value(c(tst), digits=18)
     }
   all.ok
  }


  if (synopsis & !verbose) cat("All curvature tests ...")

  if (verbose) cat("DSE curvature tests C...\n")
  if (verbose) cat("DSE curvature test C 1 ...")
  VARmodel <- est.VARX.ls(eg1.DSE.data.diff, re.add.means=F)
  SSmodel  <- l(to.SS(VARmodel),  eg1.DSE.data.diff)
  ARMAmodel<- l(to.ARMA(SSmodel), eg1.DSE.data.diff)

  spanVAR <- span(VARmodel)

#  if (print.values) print.test.value(spanVAR, digits=digits)
  # either S or R values work here with fuzz.large but try with something tighter
  # R 1.0.1:
#  if (is.R())
  tst <- c(
     431.324809764142003,  431.324809764133477,  431.324809764125973,  
     29.9938661246376412,  29.9938661245786164,  29.9938661245442368,  
     14.0487108329584771,   14.048710832937374,  14.0487108329158925, 
      8.41716993883251874,  8.41716993883184195,  8.41716993881048836,   
      6.6727488248550193,  6.67274882470962183,  6.67274882470017783,  
      6.23686447689305545,  6.23686447689176671,  6.23686447680381217,  
      0.278224057276565351,  0.278224057275740233,  0.278224057275401837,  
      0.226296510339683066,  0.226296510339373452,  0.226296510338978879,  
      0.22109512238799589,  0.221095122387542337,  0.221095122387375581,  
      0.217635227413440741,  0.217635227413369131,  0.217635227412733029,  
      0.210139509606726088,  0.210139509606644848,  0.210139509606350583,  
      0.193542574580348015,  0.193542574580293447,  0.193542574579900289,  
      0.133100015736768662,  0.133100015736627497,  0.133100015736183325,  
      0.124513362052187021,  0.124513362052168133,  0.124513362051330567,  
      0.121954031479664826,  0.121954031479572386,  0.121954031479351521,  
      0.113732578253665478,  0.113732578253511143,  0.113732578253048111,  
      0.089473111199309549,  0.0894731111990460237,  0.0894731111990405698,  
      0.0871822427356184065,  0.0871822427354376067,  0.0871822427353293739,  
      0.0777910835236551579,  0.0777910835236260007,  0.0777910835233898562,  
      0.0595266223042596881,  0.0595266223042083611,  0.0595266223040415363,  
      0.0499760073414238712,  0.0499760073413548916,  0.0499760073413023087,  
      0.0485695297799133355,  0.0485695297798213396,  0.048569529779736699,  
      0.0474465370323315233,  0.0474465370321883045,  0.0474465370320783369,  
      0.0443198146313130703,  0.0443198146313004693,  0.0443198146312605984 )

#  ok <- fuzz.small >  max(abs( spanVAR - tst ))
#  all.ok <- all.ok & ok
#  if (verbose)  {if (ok) cat("ok\n")  else cat("failed!\n") }
  all.ok <- test(spanVAR, tst, all.ok, flag="1", fuzz=100*fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test C 2 ...")

  spanSS <- span(SSmodel)
#  if (print.values) print.test.value(spanSS, digits=digits)
#  ok <- ok & fuzz.small >  max(abs( spanSS -  tst))
#  all.ok <- all.ok & ok
#  if (verbose)  {if (ok) cat("ok\n")  else cat("failed!\n") }
  all.ok <- test(spanSS, spanVAR, all.ok, flag="2", fuzz=10*fuzz.large,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test C 3 ...")
  spanARMA <- span(ARMAmodel)
  all.ok <- test(spanARMA, 0, all.ok, flag="3", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test C 4 ...")
  curvatureVAR <- curvature(VARmodel)$stats
  tst <- c(24, 100, 0.05, 1.4877780902309642e-06, 3.2049317816656057e-06, 
           1.917679915995849e-06, 4.1310147999845919e-06, 1, 1           )
  all.ok <- test(curvatureVAR, tst, all.ok, flag="4", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test C 5 ...")
  curvatureSS <- curvature(SSmodel)$stats
  tst <- c(48, 100, 0.05, 255.387386434704, 98.33696083718515,
           322.6416248003061, 124.23321787877502, 1, 1.0000000004520695  )
  all.ok <- test(curvatureSS, tst, all.ok, flag="5", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test C 6 ...")
  curvatureARMA <- curvature(ARMAmodel)$stats
  tst <-c(72, 100, 0.05, 31947799733885313024, 1708051.5249938569,
   42274456907383595008,  2260154.101077503, 1, 1.0000267463352852 )
# above looks suspicious
  all.ok <- test(curvatureARMA, tst, all.ok, flag="6", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test C 7 ...")
  hessianVAR <- hessian(VARmodel)
  all.ok <- test(sum(hessianVAR), 7219.717083137912, all.ok, flag="7", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test C 8 ...")
  hessianSS <- hessian(SSmodel)
  all.ok <- test(sum(hessianSS), 7844.3395239153897, all.ok, flag="8", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (verbose) cat("DSE curvature test C 9 ...")
  hessianARMA <- hessian(ARMAmodel)
  all.ok <- test(sum(hessianARMA), 1789846677.6677122, all.ok, flag="9", fuzz=fuzz.small,
                 print.values=print.values)
  if (verbose) cat(" completed\n")


  if (synopsis) 
    {if (verbose) cat("All curvature tests completed")
     if (all.ok$ok) cat(" OK\n")
     else cat(", some FAILED! all.ok$error = ", all.ok$error,"\n")
    }

  if (all.ok$ok) invisible(T)  else stop("FAILED")
}

#######################################################################

#                    end

#######################################################################

