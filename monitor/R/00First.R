.First.lib <- function(library,section){
if(!require("padi"))
   warning("For the database interface this package requires the padi package.")
if(!require("dsepadi"))
 warning("For the database interface this package requires the dsepadi package.")
  if(!require("dse2")) warning("This package requires the dse2 package.")
invisible()}
