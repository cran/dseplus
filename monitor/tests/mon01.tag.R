   require("monitor")


tagged.function.tests <- function(verbose=TRUE, synopsis=TRUE, fuzz.small=1e-10)
{# A short set of tests of the tagged class methods. 

 data("eg1.DSE.data.diff", package="dse1")

  if (!is.TSdata(eg1.DSE.data.diff))
     stop("Test data not found. Testing stopped.")
  if (synopsis & !verbose) cat("All tagged class tests ...")
  if (verbose) cat("tagged class test 1 ... ")
#  z <- outputData(eg1.DSE.data.diff)
#  tags(z, "tags") <- array("a", dim(z))
#  dseclass(z) <- "tagged"
  z <- outputData(eg1.DSE.data.diff)
  z <- tagged(z, array("a", dim(z)))
  ok <- is.tagged(z)
  all.ok <- ok
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


  if (verbose) cat("tagged class test 2... ")
#  zz <- z
#  tags(zz) <- array("b", dim(z))
  zz <- tagged(z, array("b", dim(z)))
  ok <- testEqual(z,z) & (!testEqual(z,zz))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tagged class test 3... ")
  zz <- tfwindow(z, start=c(1989,1))
  tags(zz) <- array("b", dim(zz))
  zzz <- tbind(tfwindow(z, start=c(1989,1)),zz)
  ok <- (2*sum(tfwindow(outputData(eg1.DSE.data.diff),
           start=c(1989,1)))) ==  sum(zzz)
  ok <- ok & all("a" == tags(zzz)[,1:3]) &  all("b" == tags(zzz)[,4:6]) 
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tagged class test 4... ")
  zzz <- splice(zz, tfwindow(z, end=c(1990,1)))
  ok <- testEqual.matrix(z,zzz) & (!testEqual(z,zzz))
  zzz <- splice(zz, tfwindow(outputData(eg1.DSE.data.diff),
                           end=c(1990,1)), tag2="x")
  ok <- ok & testEqual.matrix(z,zzz) & (!testEqual(z,zzz))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (synopsis) 
    {if (verbose) cat("All tagged class tests completed")
     if (all.ok) cat(" OK\n") else cat(", some FAILED!\n")
    }

  if (all.ok) invisible(TRUE)  else stop("FAILED")
}



   tagged.function.tests(verbose=TRUE)
