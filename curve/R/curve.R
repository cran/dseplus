############################################################################

#    functions for gradient calculation

############################################################################

gradNumerical <- function(func,x, eps=1e-12) {
#  very simple (crude) numerical approximation
  f <-func(x)
  df <-1:length(x)
  for (i in 1:length(x)) {
    dx <- x
    dx[i] <- dx[i] +eps 
    df[i] <- (func(dx)-f)/eps
   }
df
}


gradRichardson <- function(func, x, d=0.01, eps=1e-4, r=6, show.details=FALSE)
{
#  *  modified by Paul Gilbert from orginal code by XINGQIAO LIU.

  v <- 2               # reduction factor.
  n <- length(x)       # Integer, number of variables.
  a.mtr <- matrix(1, r, n) 
  b.mtr <- matrix(1, (r - 1), n)

#  first order derivatives are stored in the matrix a.mtr[k,i], 
#        where the indexing variables k for rows(1 to r),  i for columns
#       (1 to n),  r is the number of iterations, and n is the number of
#       variables.

  h <- abs(d*x)+eps*(x==0.0)
  for(k in 1:r)  { # successively reduce h                
     for(i in 1:n)  {
         x1.vct <- x2.vct <- x
         x1.vct[i]  <- x[i] + h[i]
         x2.vct[i]  <- x[i] - h[i]
         if(k == 1) a.mtr[k,i] <- (func(x1.vct) - func(x2.vct))/(2*h[i])
         else{
           if(abs(a.mtr[(k-1),i])>1e-20)
                 # some functions are unstable near 0.0              
                 a.mtr[k,i] <- (func(x1.vct)-func(x2.vct))/(2*h[i])
            else  a.mtr[k, i] <- 0
          }
      }
     h <- h/v     # Reduced h by 1/v.
    }	
   if(show.details)  {
        cat("\n","first order approximations", "\n")		
        print(a.mtr, 12)
    }

#------------------------------------------------------------------------
# 1 Applying Richardson Extrapolation to improve the accuracy of 
#   the first and second order derivatives. The algorithm as follows:
#
#   --  For each column of the 1st and 2nd order derivatives matrix a.mtr,
#       say, A1, A2, ..., Ar, by Richardson Extrapolation, to calculate a
#       new sequence of approximations B1, B2, ..., Br used the formula
#
#          B(i) =( A(i+1)*4^m - A(i) ) / (4^m - 1) ,  i=1,2,...,r-m
#
#             N.B. This formula assumes v=2.
#
#   -- Initially m is taken as 1  and then the process is repeated 
#      restarting with the latest improved values and increasing the 
#      value of m by one each until m equals r-1
#
# 2 Display the improved derivatives for each
#   m from 1 to r-1 if the argument show.details=T.
#
# 3 Return the final improved  derivative vector.
#-------------------------------------------------------------------------

  for(m in 1:(r - 1)) {		
#     for(i in 1:(r - m)) b.mtr[i,]<- (a.mtr[(i+1),]*(4^m)-a.mtr[i,])/(4^m-1)
#     a.mtr<- b.mtr
     a.mtr<- (a.mtr[2:(r+1-m),,drop=FALSE]*(4^m)-a.mtr[1:(r-m),,drop=FALSE])/(4^m-1)
     if(show.details & m!=(r-1) )  {
        cat("\n","Richarson improvement group No. ", m, "\n")		
        print(a.mtr[1:(r-m),,drop=FALSE], 12)
      }
   }
a.mtr
}


############################################################################

#    functions for projecting one space on another

############################################################################

project <- function(c1, c2, signif = 0.05, eps=1e-5, fuzz=1e-14,
    show.details=FALSE, warn=TRUE)
{  p1 <- c1$Dlist$p
   N <- length(c1$Dlist$f0)   # dimesion of sample space
   p2 <- c2$Dlist$p
   if (N != length(c2$Dlist$f0)) stop("sample space dimesions do not agree.")
   if (p2 < p1) stop("submodel should be the first argument.")
   if (max(abs(c1$Dlist$f0 - c2$Dlist$f0)) > fuzz)
      warning("projection not computed at the same sample space point as the original curvature calculation.")

   residual <- c(c1$Dlist$f0)   # c() in case it is a column vector
   s.sqr <- sum(residual^2)/(N-2)  # check for time series??

   cat("p1=",p1," p2=",p2,"N=",N,"\n")
   cat("max. span of acceleration space as determined by number of vectors \nin 2nd derivative array:    ")
   cat("   1st model=", p1 *(1+p1)/2,
       ",2nd model=", p2 *(1+p2)/2,"\n" )   

   d1 <- La.svd(c1$Dlist$D[,1:p1])$d
   d2 <- La.svd(c2$Dlist$D[,1:p2])$d
   dj <- La.svd(cbind(c1$Dlist$D[,1:p1], c2$Dlist$D[,1:p2]))$d
   cat("span of tangent vectors:       1st model=", sum(d1 > eps*d1[1]) )
   cat(", 2nd model=",         sum(d2 > eps*d2[1]))
   cat(", jointly=",       sum(dj > eps*dj[1]), "\n")
 
   d1 <- La.svd(c1$Dlist$D[,-(1:p1)])$d
   d2 <- La.svd(c2$Dlist$D[,-(1:p2)])$d
   dj <- La.svd(cbind(c1$Dlist$D[,-(1:p1)], c2$Dlist$D[,-(1:p2)]))$d
   cat("span of acceleration vectors:  1st model=", sum(d1 > eps*d1[1]) )
   cat(", 2nd model=",         sum(d2 > eps*d2[1]))
   cat(", jointly=",       sum(dj > eps*dj[1]), "\n")
   cat("using eps*100 (sensitivity for significance of svd$d )\n")
   cat("span of acceleration vectors:  1st model=", sum(d1 > eps*100*d1[1]) )
   cat(", 2nd model=",         sum(d2 > eps*100*d2[1]))
   cat(", jointly=",       sum(dj > eps*100*dj[1]), "\n")

   m2 <- p2 * (p2 + 3)/2  # m= p+pp is second dimension of D.
   m2 <- min(m2,dim(c2$Dlist$D)[2]) 
   m1 <- p1 * (p1 + 3)/2 

#  c1 projected into c2

   QRofD2 <- qr(c2$Dlist$D)
   T1inT2 <- qr.qty(QRofD2, c1$Dlist$D[,  1:p1 ] )[  1:p2 ,,drop=FALSE]  
   N1inT2 <- qr.qty(QRofD2, c1$Dlist$D[,-(1:p1)] )[  1:p2 ,,drop=FALSE]  
   N1inN2 <- qr.qty(QRofD2, c1$Dlist$D[,-(1:p1)] )[(p2+1):m2,,drop=FALSE]  

#  c2 projected onto c1

   QRofD1 <- qr(c1$Dlist$D)
   N2inN1      <- qr.qty(QRofD1, c2$Dlist$D[,-(1:p2)] )[(p1+1):m1,,drop=FALSE]  
   T2inN1      <- qr.qty(QRofD1, c2$Dlist$D[,  1:p2 ] )[(p1+1):m1,,drop=FALSE] 
   N2andT2inN1 <- qr.qty(QRofD1, c2$Dlist$D           )[(p1+1):m1,,drop=FALSE]
   N2andT2inc1 <- qr.qty(QRofD1, c2$Dlist$D           )

   N1inN1 <- qr.qty(QRofD1, c1$Dlist$D[,(p1+1):m1] )[(p1+1):m1,,drop=FALSE]  
#   v1 <- svd(N2andT2inT1) should use La.svd
#   v2 <- svd(T1inT1)
browser()

   cur2in1 <- relCurvature(s.sqr,N2andT2inc1[1:p2,1:p1], N2andT2inc1[,(p1+1):m1], 
                      show.extra.details=show.details)
   C2in1 <- list(C.parameter=cur2in1[1:p1,,     ,drop=FALSE],
                 C.intrinsic=cur2in1[(p1+1):m1,,,drop=FALSE])

#    z <- relCurvature(s.sqr,N1inN1[1:p1,1:p1], N2andT2inc1[,(p1+1):m1], 
#                       show.extra.details=show.details)
#    zz <- list(C.parameter=z[1:p1,,     ,drop=FALSE],
#                  C.intrinsic=z[(p1+1):m1,,,drop=FALSE])

# svd(   zz$C.intrinsic[1,,])$d   should use La.svd
# svd( c1$C$C.intrinsic[1,,])$d
# svd( c2$C$C.intrinsic[1,,])$d
# svd(C2in1$C.intrinsic[1,,])$d

# svd( c1$C$C.intrinsic[1,,])$d / svd(C2in1$C.intrinsic[1,,])$d[1:2]

# use , symmetric=?, only.values=TRUE
# La.eigen( c1$C$C.intrinsic[1,,])$values
# La.eigen( c2$C$C.intrinsic[1,,])$values
# La.eigen(C2in1$C.intrinsic[1,,])$values

   effective<-effectiveCurvature(cur1on2,QRofD2, residual, s.sqr,
                      show.details=show.details, warn=warn)

   cstats1on2 <-curvatureStats(cur1on2, N, signif=signif)

   R1 <- qr.qty(QRofD1, c2$Dlist$D )[1:m1,,drop=FALSE]  
   cur1 <- relCurvature(s.sqr,R1[1:p1,1:p1], R1[,(p1 + 1):m1], 
                      show.extra.details=show.details)
   effective1<-effectiveCurvature(cur1,QRofD1, residual, s.sqr,
                      show.details=show.details)

   cstats1 <-curvatureStats(cur1, N, signif=signif)
browser()

   result <- rbind(c1$stats,c(p2, N, signif, cstats1on2,
                   effective$min.axis.ratio, effective$max.axis.ratio),
                   c2$stats)
   dimnames(result)<-list(c("original c1", "projected", "original c2"),
                          c("Parms","Sample","Sign. level","RMS Parameter", 
         "RMS Intrinsic","c*sqrt(F) Parameter","c*sqrt(F) Intrinsic",
         "Min Axis Ratio","Max Axis Ratio")) 
   result
}


#######################################################################

#            calculate Hessian 

#######################################################################


hessian <- function (func, x, func.args=NULL, d=0.01, eps=1e-4, r=6) UseMethod("hessian")

hessian.default <- function(func, x, func.args=NULL, d=0.01, eps=1e-4, r=6)
{  D <- genD.default(func, x, func.args=func.args, d=d, eps=eps, r=r)$D
   H <- diag(0,length(x))
   u <- 0
   for(i in 1:length(x))
     {for(j in 1:i) 
        {u <- u + 1
         H[i,j] <- D[,u]
     }  }
   H <- H +t(H)
   diag(H) <- diag(H)/2
   H
}


#######################################################################

#              span (dimension of tangent space) calculation

#######################################################################

span <- function (func, x, func.args=NULL, d=0.01, eps=1e-4, r=6,
      show.details=FALSE, ...)  UseMethod("span")

span.default <- function(func, x, func.args=NULL, d=0.01, eps=1e-4, r=6,
      show.details=FALSE, ...)
{#  (... further arguments, currently disregarded)
 #   v       reduction factor. This could be a parameter but ...
 #             N.B. must be 2 for richarson formula as used.

  v <- 2
  p <- length(x)     #  number of parameters (= dim of M if not over parameterized.
  n <- length(do.call("func",append(list(x), func.args)))  #dimension of sample space.
  D <- array(1, c(n, p, r)) 
  h <- abs(d*x)+eps*(x==0.0)
  for(k in 1:r)  # successively reduce h 
    {if(show.details) cat("\n k=",k, " p: ")
     for(i in 1:p)  # each parameter  - first deriv.
      {if(show.details) cat(" ",i)
       D[,i,k]<-(do.call("func",append(list(x+(i==(1:p))*h), func.args)) -
                 do.call("func",append(list(x-(i==(1:p))*h), func.args)))/(2*h[i])
       NULL 
      }
     h <- h/v     # Reduced h by 1/v.
    NULL 
    }	
  if(show.details)  
      {cat("\n","first order approximations", "\n")		
       for(i in 1:p){ cat(" parameter " ,i,"\n");print(D[,i,(1:r)], 12); NULL}
      }
  for(m in 1:(r - 1)) 		
     {D<- (D[,,2:(r+1-m)]*(4^m)-D[,,1:(r-m)])/(4^m-1)
      if(show.details & m!=(r-1) )  
        {cat("\n","Richarson improvement group No. ", m, "\n")		
         for(i in 1:p){ cat(" parameter " ,i,"\n"); print(D[,i,(1:(r-m))], 12) }
        }
      NULL
     }
   La.svd(D)$d
}


#######################################################################

#                            curvature calculation

#######################################################################

print.curvature <- function(x, ...)  { print(x$stats, ...) }

curvature <- function (func, ...)  
 # calculate the Bates and Watts intrinsic and parameter effects curvature.
  UseMethod("curvature")



curvature.default <- function(func, x, func.args=NULL, d=0.01, eps=1e-4,r=6,
      signif=0.05, show.details=FALSE, warn=TRUE, ...)
{#  (... further arguments, currently disregarded)
  # func is a function
 curvature(genD.default(func,x, func.args=func.args, d=d, eps=eps,r=r),
     signif=0.05, show.details=FALSE, warn=warn)
}


curvature.Darray <- function(func, signif = 0.05,
   show.extra.details=FALSE, show.details=show.extra.details, warn=TRUE, ...)
{#  (... further arguments, currently disregarded)
# Note func is a list of class Darray
#  Page references are to Bates and Watts(1983), Nonlinear Regression 
#   Analysis and Its Applications.
#   Examples are also from this text. e.g.- model A data set 3.

#   modified substantially by Paul Gilbert (May, 1992 and Dec, 1996)
#    from original code by Xingqiao Liu,   May, 1991.
#
   if (all(func$D==0))
     stop("The array of first and second order derivatives is all zeros.")

   p <- func$p
   n <- dim(func$D)[1]   # dimesion of sample space
   pp <- p*(p+1)/2 # p' max number of independent col in the acceleration array.
   m <- p * (p + 3)/2  # m= p+pp is second dimension of D. It is also the max 
                       # span of the combined tangent and accelation spaces.
   if (m!=dim(func$D)[2]) stop("Dimesion of D matrix is not consistent with p")
   residual <- c(func$f0)   # c() in case it is a column vector
   s.sqr <- sum(residual^2)/(length(residual)-2)
   if(show.extra.details) cat("p=",p," pp=",pp," m=",m, " n=",n,"\n")
 
#  Form QR decomposition of D and produce matrix Q and R.
#

   QRofD <- qr(func$D) 

# Now calculate R, which is D rotated onto an othogonal basis with the 
# tangent space first, then acceleration space. 
# In this basis first rows of R are in the tangent and the next rows are in the
# acceleration space. Following rows would be
# all zeros  and so those rows are truncated from R.

#   Q <- qr.qy(QRofD, diag(1,n)) [,1:m]
#   R <- t(Q) %*% func$D   this is the same as
#  R[1:m,]could be used to truncate but it is possible for m to exceed dim(D)[1]
   R <- qr.qty(QRofD, func$D )  

   if(show.extra.details) 
        {cat( "Matrix D, tangent and acceleration vectors in columns ")
         cat("organized as in Bates and Watts p235, 236\n")
         print(func$D)
         cat("Residual Vector\n"); print(residual)
	 cat("Matrix Q from QR decomposition of D\n")
         print( qr.qy(QRofD, diag(n)))
	}
   if(show.details) 
        {cat("Matrix R (organized as in Bates and Watts R1 p236)\n")
         print(R)
	}

#  The matrix R partitions into p by p matrix R11 and m by p'  matrix R2.
#  R2 is Bates and Watts  R12 
#                         R22  p236

#   R11 <- R[1:p,1:p]
#   R2 <- R[1:m,(p + 1):m] 

#                                     R[1:m,(p + 1):m]
   cur <- relCurvature(s.sqr,R[1:p,1:p], R[,(p + 1):m], 
                      show.extra.details=show.extra.details)

   if (m>dim(cur)[1]) 
     {m <- dim(cur)[1]
      warning("acceration space dimension reduced for smaller sample space.")
     }
   C <- classed(list(C.parameter=cur[1:p,,,drop=FALSE],# curvatureArray constructor
             C.intrinsic=cur[(p+1):m,,,drop=FALSE]),"curvatureArray" )

# effective residual curvature. Bates and Watts p260
# Calculate the scaled RMS curvatures relative to a confidence disk 
# radius and extreme axis ratios. ref. Bates and Watts p254 eqn (7.23).
#  and Bates and Watts J.R. Statist.Soc. B (1980).

   effective<-effectiveCurvature(cur,QRofD, residual, s.sqr, 
                              show.details=show.extra.details, warn=warn)
   cstats <-curvatureStats(cur, n, signif=signif)

   result <- matrix(c(p, n, signif, cstats,
               effective$min.axis.ratio, effective$max.axis.ratio),1,9)
   dimnames(result)<-list(NULL,c("Parms","Sample","Sign. level","RMS Parameter",
         "RMS Intrinsic","c*sqrt(F) Parameter","c*sqrt(F) Intrinsic",
         "Min Axis Ratio","Max Axis Ratio")) 
   classed(list(stats=result, Dlist=func, C=C), "curvature")# curvature.Darray constructor
}

print.curvatureArray <- function(x, ...)
  {cat("Parameter effects relative curvature array:\n")
   for (i in 1:dim(x$C.parameter)[1]) 
       {print(x$C.parameter[i,,], ...); cat("\n"); NULL}
   cat("Intrinsic relative curvature array:\n")
   for (i in 1:dim(x$C.intrinsic)[1])
       {print(x$C.intrinsic[i,,], ...); cat("\n"); NULL}
   invisible(x)
  }



curvatureStats <- function(cur, n, signif=0.05)
{# See documentation for curvature.Darray
 p  <- dim(cur)[2]
 pp <- dim(cur)[1] - p
 cstats <- rep(NA,4)
   #  i from 1 to p for the RMS parameter effects curvature.
   #  i from p+1 to m=p+p' for the RMS intrinsic effects curvature.
   di <- 0
   for(i in 1:p) di <- di + sum(diag(cur[i,,]))^2
   cstats[1] <-  2*sum(cur[1:p,,]^2) + di
   cstats[1] <-  sqrt(cstats[1]/(p*(p+2)))
   di <- 0
   for(i in p+(1:pp)) di <- di + sum(diag(cur[i,,]))^2
   cstats[2] <-  2*sum(cur[p+(1:pp),,]^2) + di
   cstats[2] <-  sqrt(cstats[2]/(p*(p+2)))
  value.f <- qf(1 - signif, p, n - p)	# F(p, n-p; signif)
  cstats[3] <- cstats[1] * sqrt(value.f) # c*sqrt(F) Parameter
  cstats[4] <- cstats[2] * sqrt(value.f) # c*sqrt(F) Intrinsic
  cstats
}


effectiveCurvature <- function(cur, QRofD, residual, s.sqr,
			 show.details=FALSE, warn=TRUE)
{
 p  <- dim(cur)[2]
 pp <- dim(cur)[1] - p
 Qz <- qr.qty(QRofD, residual)[p+(1:pp)]/c(sqrt(p * s.sqr))  # p' vector of rotated residuals
 #pxp effective residual (intrinsic) curvature matrix. Bates and Watts p260
   # Qz is a vector of length pp.  
   B <- matrix(NA,p,p) 
   A.intrin <- cur[p+(1:pp),,]
 #   for (i in 1:p)
 #       for (j in 1:p)
 #         B[i,j] <-  A.intrin[,i,j] %*% Qz 
 # following is equivalent to above
   B <- apply(A.intrin * Qz, c(2,3), sum)
   if(show.details) 
     { cat("Scaled and Q-transformed Residual Vector\n"); print(Qz)
       cat("Effective residual curvature matrix B (Bates and Watts p 260)\n")
       print(B)
     }
   if ( max(abs(B- t(B))) > 1e-7 ) warning("B is not symmetric.") 
   eigv <- La.eigen(B, symmetric=TRUE, only.values=TRUE)$values
   if(max(Im(eigv)) >= 1e-15) 
          warning("Calculated eigenvalues have imaginary part.")

   if(show.details)   {cat("Eigenvalues of B\n"); print(eigv )}

   if(warn && (max(Re(eigv)) >= 1.0))
     {warning( paste(
      "N.B. (I-B) is not positive definite as it should be at a local min! ",
      "Axis ratio will evaluate to NA or NaN. ",
      "Eigenvalues of B ", eigv))
     }

 #Mod in the following gets rid of small Im part, but the eigv should be real
# axis.ratios <- 1/sqrt(1-eigv)
# list(B=B, eigenvalues=eigv,
#      min.axis.ratio=min(axis.ratios),
#      max.axis.ratio=max(axis.ratios)  )
# Bates and Watts seem to use the following but the above should be considered.
# It seems possible to have an eigv < -1
# sqrt gives NaN for neg. numbers but test avoids producing additional warnings.
 list(B=B, eigenvalues=eigv,
      min.axis.ratio= if(1 > min(Mod(eigv))) 1/sqrt(1-min(Mod(eigv))) else NaN,
      max.axis.ratio= if(1 > max(Mod(eigv))) 1/sqrt(1-max(Mod(eigv))) else NaN)
}

relCurvature <- function(s.sqr, R11, R2, show.extra.details=FALSE, 
                     eps=sqrt(.Machine$double.eps))
{
 p  <- dim(R11)[2]
 p2 <- dim(R11)[1]  # usually equal p except for projections
 pp <- dim(R2)[2]
 # m will usually be p+pp, but if the number of observations is small relative 
 # to pp, then Q and R2 will be truncated.
 # For smaller problems it would be quicker to do the calculation just to p+pp
 # but small problems are quick anyway.
 m  <- dim(R2)[1] 
 A  <- array(NA,c(m,p,p))  # for accel. array B&W (7.5)
 C  <- array(NA,c(m,p2,p2))
   
# R11.inv <- solve(R11)
 v <-La.svd(R11)
 if (all(v$d == 0)) stop("R11 is completely degenerate.")
 d <- 1/v$d
 illim <- v$d < (v$d[1] * eps)
 if (any(illim)) 
   {warning("eliminating degenerate subspace for R11.")
    d[illim] <- 0
   }
 R11.inv <- Conj(t(v$vt)) %*% diag(d) %*% t(v$u) # switched from svd to La.svd

 R11.invT <- t(R11.inv)
 if(show.extra.details) 
      { cat("Bates and Watts matrix R11 p236)\n")
        print(R11) 
        cat("Matrix R2 (Bates and Watts  R12 \n")
        cat("                            R22 )  p236)\n")
        print(R2)
        cat("Inverse Of Matrix R11  e.g. Bates and Watts p243\n")
        print(R11.inv)
      }
  # expand R2 to make A
  # there may be a better way to do this without expanding
   k <- 0
   for(i in 1:p) 
     for(j in 1:i) 
       {k <- k + 1
        A[,i,j] <- R2[, k]
        A[,j,i] <- R2[, k]
        NULL
       }
   if(show.extra.details)
     {cat("Acceleration array (Bates and Watts matrix A p237).\n")
      cat("Parameter effects relative acceleration array:\n")
      for (i in 1:p) {print(A[i,,]); cat("\n"); NULL}
      cat("Intrinsic relative acceleration array:\n")
      for (i in (p+1):m) {print(A[i,,]); cat("\n"); NULL}
     }	

# compute C  Batts and Watts defn. (7.16) p242 and examples p244 &  p245
   for(i in 1:m)
      {C[i,,]<-(R11.invT %*% A[i,,] %*% R11.inv) * c(sqrt(p * s.sqr)); NULL}
if(show.extra.details)
  {cat("Relative curvature array (Bates and Watts matrix C p242 and eg p 243- 245).\n")
      cat("Parameter affects curvature array:\n")
      for (i in 1:p) {print(C[i,,]); cat("\n")}
      cat("Intrinsic curvature array:\n")
      for (i in (p+1):m) {print(C[i,,]); cat("\n")}
  }	

C
}


#######################################################################

#                               D matrix calculation

#######################################################################

# Notes:

#   Generating the D matrix can be computationally demanding. 
#   The S version is slow and suffers from memory problems 
#   due to the way S allocates memory in loops (but is not so bad in R).
#   The C version is a bit faster but suffers (even worse) from memory problems.
#   It has been moved to defunct. The ARMA and SS methods are much preferred for
#   those models, but do not work with general models.
#   Both the S and the C versions take the name of an S function 
#   as an arguement (and call it multiple times).  
#######################################################################

genD <- function(func, x, func.args=NULL, d=0.01, eps=1e-4, r=6)UseMethod("genD")

genD.default <- function(func, x, func.args=NULL, d=0.01, eps=1e-4, r=6){
#   modified substantially by Paul Gilbert (May, 1992)
#    from original code by Xingqiao Liu,   May, 1991.
#  r       the number of Richardson improvement iterations.
#  v      reduction factor for Richardson iterations. This could
#      be a parameter but the way the formula is coded it is assumed to be 2.

#  This function is not optimized for S speed, but
#  is organized in the same way it is implemented in C (to facilitate checking
#  the code).

#         - THE FIRST ORDER DERIVATIVE WITH RESPECT TO Xi IS 
#           F'(Xi)={F(X1,...,Xi+d,...,Xn) - F(X1,...,Xi-d,...,Xn)}/(2*d)
#                
#         - THE SECOND ORDER DERIVATIVE WITH RESPECT TO Xi IS 
#           F''(Xi)={F(X1,...,Xi+d,...,Xn) - 2*F(X1,X2,...,Xn)+
#                     F(X1,...,Xi-d,...,Xn)}/(d^2)
#
#         - THE SECOND ORDER DERIVATIVE WITH RESPECT TO Xi, Xj IS 
#           F''(Xi,Xj)={F(X1,...,Xi+d,...,Xj+d,...,Xn)-2*F(X1,X2,...,Xn)+
#                       F(X1,...,Xi-d,...,Xj-d,...,Xn)}/(2*d^2) -
#                      {F''(Xi) + F''(Xj)}/2

#  The size of d is iteratively reduced and the Richardson algorithm is applied to
#    improve the accuracy of the computation.

   v <-2
   zt <- .001     #Rbug unix.time has trouble finding func
   # f0 <- func(x)
   f0 <- do.call("func",append(list(x), func.args))
 # zt <-unix.time(f0 <- func(x))   #  f0 is the value of the function at x.
                   # (Real value or point in sample space ie- residual vector)	
   p <- length(x)  #  number of parameters (theta)
   zt <- zt[1] *r*2*(p*(p + 3))/2
   if(zt>500) 
     cat("D matrix calculation roughly estimated to take about ",
          (10*ceiling(zt/10)),"seconds(without other system activity)...\n")
   h0 <- abs(d*x)+eps*(x==0.0)
   D <- matrix(0, length(f0),(p*(p + 3))/2)
   #length(f0) is the dim of the sample space
   #(p*(p + 3))/2 is the number of columns of matrix D.( first
   #   der. & lower triangle of Hessian)
   Daprox <- matrix(0, length(f0),r) 
   Hdiag  <-  matrix(0,length(f0),p)
   Haprox <-  matrix(0,length(f0),r)
   for(i in 1:p)    # each parameter  - first deriv. & hessian diagonal
        {h <-h0
         for(k in 1:r)  # successively reduce h 
           {#f1 <- func(x+(i==(1:p))*h)
            #f2 <- func(x-(i==(1:p))*h) 
            f1 <- do.call("func",append(list(x+(i==(1:p))*h), func.args))
            f2 <- do.call("func",append(list(x-(i==(1:p))*h), func.args))
            Daprox[,k] <- (f1 - f2)  / (2*h[i])    # F'(i) 
            Haprox[,k] <- (f1-2*f0+f2)/ h[i]^2     # F''(i,i) hessian diagonal
            h <- h/v     # Reduced h by 1/v.
            NULL
           }
         for(m in 1:(r - 1))
            for ( k in 1:(r-m))
              {Daprox[,k]<-(Daprox[,k+1]*(4^m)-Daprox[,k])/(4^m-1)
               Haprox[,k]<-(Haprox[,k+1]*(4^m)-Haprox[,k])/(4^m-1)
               NULL
              }
         D[,i] <- Daprox[,1]
         Hdiag[,i] <- Haprox[,1]
         NULL
        }	
   u <- p

   for(i in 1:p)   # 2nd derivative  - do lower half of hessian only
     {for(j in 1:i) 
        {u <- u + 1
            if (i==j) { D[,u] <- Hdiag[,i]; NULL}
            else 
             {h <-h0
              for(k in 1:r)  # successively reduce h 
                {#f1 <- func(x+(i==(1:p))*h + (j==(1:p))*h)
                 #f2 <- func(x-(i==(1:p))*h - (j==(1:p))*h)
                 f1 <- do.call("func", append(
                        list(x+(i==(1:p))*h + (j==(1:p))*h), func.args))
                 f2 <- do.call("func",append(
                        list(x-(i==(1:p))*h - (j==(1:p))*h), func.args))  
                 Daprox[,k]<- (f1 - 2*f0 + f2 -
                                  Hdiag[,i]*h[i]^2 - 
                                  Hdiag[,j]*h[j]^2)/(2*h[i]*h[j])  # F''(i,j)  
                 h <- h/v     # Reduced h by 1/v.
                }
              for(m in 1:(r - 1))
                 for ( k in 1:(r-m))
                   {Daprox[,k]<-(Daprox[,k+1]*(4^m)-Daprox[,k])/(4^m-1); NULL}
              D[,u] <- Daprox[,1]
              NULL
             }
          }  
       }
invisible(classed(list(D=D,  # Darray constructor (genD.default)
             p=length(x),f0=f0,func=func, x=x, d=d, eps=eps, r=r), "Darray"))
}


#######################################################################

#                    end

#######################################################################

