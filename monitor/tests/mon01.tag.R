   require("ts",      warn.conflicts=F)
   require("dse2",    warn.conflicts=F)
   require("dsepadi", warn.conflicts=F)
   require("monitor", warn.conflicts=F)
#   A TS PADI server is necessary for these tests.
#   The next line is only necessary to remove this in an old version which set
#     home in frame 0. (I'll never do that again.)
#   if (is.S()) remove("DSE.HOME", where=0) 
   require("padi",    warn.conflicts=F)
   if (is.S()) {
	# the next 2 lines remove old versions of PADI in the search path
 	invisible(if(0!=length(grep("b*/PADI/.Data",search())))
                         detach(grep("b*/PADI/.Data",search()))  )
	attach(paste(getenv("PADI"),"/.Data", sep=""), pos=2)
	#load.padi(from=".")    # this gets the version in pwd
	load.padi()           # this gets the version indicated by PADI
	# load.padi does the following two dynamic loads 
	#dyn.load.shared("/usr/lib/libnsl.so")     # splus 3.3 on SunOS5
	# If the shared library is not loaded then the next has missing symbols
	#dyn.load(paste(getenv("PADI"),"/lib/splusclnt.o", sep=""))
	search()
   }




tagged.function.tests <- function(verbose=T, synopsis=T, fuzz.small=1e-10)
{# A short set of tests of the tagged class methods. 

  if      (is.R()) data("eg1.DSE.data.diff", package="dse1")
  else if (is.S()) source(paste(DSE.HOME, "/data/eg1.DSE.data.diff.R", sep=""))

  if (!is.TSdata(eg1.DSE.data.diff))
     stop("Test data not found. Testing stopped.")
  if (synopsis & !verbose) cat("All tagged class tests ...")
  if (verbose) cat("tagged class test 1 ... ")
#  z <- output.data(eg1.DSE.data.diff)
#  tags(z, "tags") <- array("a", dim(z))
#  dseclass(z) <- "tagged"
  z <- output.data(eg1.DSE.data.diff)
  z <- tagged(z, array("a", dim(z)))
  ok <- is.tagged(z)
  all.ok <- ok
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }


  if (verbose) cat("tagged class test 2... ")
#  zz <- z
#  tags(zz) <- array("b", dim(z))
  zz <- tagged(z, array("b", dim(z)))
  ok <- test.equal(z,z) & (!test.equal(z,zz))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tagged class test 3... ")
  zz <- tfwindow(z, start=c(1989,1))
  tags(zz) <- array("b", dim(zz))
  zzz <- tbind(tfwindow(z, start=c(1989,1)),zz)
  ok <- (2*sum(tfwindow(output.data(eg1.DSE.data.diff),
           start=c(1989,1)))) ==  sum(zzz)
  ok <- ok & all("a" == tags(zzz)[,1:3]) &  all("b" == tags(zzz)[,4:6]) 
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (verbose) cat("tagged class test 4... ")
  zzz <- splice(zz, tfwindow(z, end=c(1990,1)))
  ok <- test.equal.matrix(z,zzz) & (!test.equal(z,zzz))
  zzz <- splice(zz, tfwindow(output.data(eg1.DSE.data.diff),
                           end=c(1990,1)), tag2="x")
  ok <- ok & test.equal.matrix(z,zzz) & (!test.equal(z,zzz))
  all.ok <- all.ok & ok 
  if (verbose) {if (ok) cat("ok\n") else cat("failed!\n") }

  if (synopsis) 
    {if (verbose) cat("All tagged class tests completed")
     if (all.ok) cat(" OK\n") else cat(", some FAILED!\n")
    }

  if (all.ok) invisible(T)  else stop("FAILED")
}



   tagged.function.tests(verbose=T)
