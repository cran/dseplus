###################################################
### chunk number 1: 
###################################################
 options(continue="  ")


###################################################
### chunk number 2: 
###################################################
library("monitor")  


###################################################
### chunk number 3: 
###################################################
x11() # in R 
#motif() or something else in Splus 


###################################################
### chunk number 4: 
###################################################
if(require("dsepadi"))
cbps.manuf.data2.ids <- TSPADIdata2(

  output=list(c("ets", "","i37013","percent.change","cbps.prod."),

  c("ets", "", "i37005","percent.change","manuf.prod.")),

  pad.start=FALSE, pad.end=TRUE ) 


###################################################
### chunk number 5: 
###################################################
if(require("padi") && checkPADIserver("ets"))
  cbps.manuf.data2 <- freeze(cbps.manuf.data2.ids) 


###################################################
### chunk number 6: 
###################################################
manuf.data.ids <- TSPADIdata(

input ="lfsa455", input.transforms="percent.change",

input.names="manuf.emp.",

output="I37005", output.transforms="percent.change",

output.names="manuf.prod.",

server="ets", pad.start=FALSE, pad.end =TRUE )

if(require("padi") && checkPADIserver("ets"))
  manuf.data <- freeze(manuf.data.ids) 


###################################################
### chunk number 7: 
###################################################
if(require("padi") && checkPADIserver("ets"))
  tfplot(manuf.data) 


###################################################
### chunk number 8: 
###################################################
if(require("padi") && checkPADIserver("ets"))
  manuf.data <- tfwindow(manuf.data, start=c(1976,2)) 


###################################################
### chunk number 9: 
###################################################
manuf.data.ids <- modify.TSPADIdata(manuf.data.ids, start=c(1976,2)) 


###################################################
### chunk number 10: 
###################################################
if(require("padi") && checkPADIserver("ets"))
  manuf.data <- freeze(manuf.data.ids) 


###################################################
### chunk number 11: 
###################################################
if(require("padi") && checkPADIserver("ets"))
  tfplot(manuf.data, start.=c(1995,11)) 


###################################################
### chunk number 12: 
###################################################
cbps.manuf.data.ids <- TSPADIdata(

input =c("lfsa462","lfsa455"),
input.transforms="percent.change",

input.names=c("cbps.emp.", "manuf.emp"),

output="i37013", output.transforms="percent.change",

output.names="cbps.prod.",

start=c(1976,2),

server="ets", db="", pad.start=FALSE,
pad.end =TRUE )

if(require("padi") && checkPADIserver("ets"))
  cbps.manuf.data <- freeze(cbps.manuf.data.ids) 


###################################################
### chunk number 13: 
###################################################
cbps.manuf.data3.ids <- TSPADIdata(

input ="lfsa462",

input.transforms="percent.change",input.names="cbps.emp.",

output=c("i37013", "i37005"),

output.transforms=c("percent.change", "percent.change"),

output.names=c("cbps.prod.","manuf.prod."),

start=c(1976,2),

server="ets", db ="",
pad.start=FALSE, pad.end =TRUE)

if(require("padi") && checkPADIserver("ets"))
  cbps.manuf.data3 <- freeze(cbps.manuf.data3.ids) 


###################################################
### chunk number 14: 
###################################################
  par(ask=TRUE)


###################################################
### chunk number 15: 
###################################################
if(require("padi") && checkPADIserver("ets"))
manuf.model <- bft(trim.na(manuf.data)) 


###################################################
### chunk number 16: 
###################################################
if(require("padi") && checkPADIserver("ets"))
manuf.model <- bft(trim.na(manuf.data), max.lag=5) 


###################################################
### chunk number 17: 
###################################################
if(require("padi") && checkPADIserver("ets"))
manuf.model <- bft(trim.na(manuf.data), verbose=FALSE, max.lag=5) 


###################################################
### chunk number 18: 
###################################################
if(require("padi") && checkPADIserver("ets"))
manuf.model 


###################################################
### chunk number 19: 
###################################################
if(require("padi") && checkPADIserver("ets")) {
  tfplot(manuf.model)

  tfplot(manuf.model, start.=c(1990,1))

  tfplot(manuf.model, start.=c(1995,1))
} 


###################################################
### chunk number 20: 
###################################################
if(require("padi") && checkPADIserver("ets")) {
  cbps.manuf.model <- bft(trim.na(cbps.manuf.data),verbose=FALSE)
  tfplot(cbps.manuf.model)
  tfplot(cbps.manuf.model, start.=c(1995,1)) 
}


###################################################
### chunk number 21: 
###################################################
if(require("padi") && checkPADIserver("ets")) {
  z <- forecast(TSmodel(manuf.model), tfwindow(manuf.data, end=c(1995,1)),
  conditioning.inputs=tfwindow(input.data(manuf.data), end=c(1996,12)))
  tfplot(z, start.=c(1995,1)) 
}


###################################################
### chunk number 22: 
###################################################
if(require("padi") && checkPADIserver("ets")) {
  forecasts(z)[[1]] 

  tfwindow(forecasts(z)[[1]], start=c(1996,3)) 
}


###################################################
### chunk number 23: 
###################################################
if(require("padi") && checkPADIserver("ets")) {
  fc <- forecastCov(manuf.model)

  tfplot(fc) 
}


###################################################
### chunk number 24: 
###################################################
if(require("padi") && checkPADIserver("ets")) {
  fc <- forecastCov(manuf.model, zero=TRUE, trend=TRUE)

  tfplot(fc) 
}


###################################################
### chunk number 25: 
###################################################
if(require("padi") && checkPADIserver("ets")) {
  outfc <-out.of.sample.forecastCovEstimatorsWRTdata(trim.na(manuf.data),

   estimation.sample=.5,

estimation.methods = list(bft=list(verbose=FALSE), est.VARX.ls=NULL),

trend=TRUE, zero=TRUE)

tfplot(outfc) 
}


###################################################
### chunk number 26: 
###################################################
if(require("padi") && checkPADIserver("ets")) {
outfc <-out.of.sample.forecastCovEstimatorsWRTdata(

trim.na(cbps.manuf.data3),

estimation.sample=.5,

estimation.methods = list(bft=list(verbose=FALSE), est.VARX.ls=NULL),

trend=TRUE, zero=TRUE)

tfplot(outfc) 
}


###################################################
### chunk number 27: 
###################################################
if(require("padi") && checkPADIserver("ets"))
  new.data <- freeze(manuf.data.ids) 


###################################################
### chunk number 28: 
###################################################
if(require("padi") && checkPADIserver("ets"))
z <- l(TSmodel(manuf.model), trim.na(new.data)) 


###################################################
### chunk number 29: 
###################################################
if(require("padi") && checkPADIserver("ets")){
  z <- l(TSmodel(manuf.model), trim.na(tfwindow(freeze(manuf.data.ids),

    start=c(1976,2))))

tfplot(z)

tfplot(z, start.=c(1995,8)) 
}


###################################################
### chunk number 30: 
###################################################
if(require("padi") && checkPADIserver("ets")){
z <- forecast(TSmodel(manuf.model), tfwindow(trim.na(new.data), end=c(1996,1)),
conditioning.inputs=trim.na(input.data(new.data)))

tfplot(z, start.=c(1995,6)) 
}


###################################################
### chunk number 31: 
###################################################
if(require("padi") && checkPADIserver("ets"))
forecasts(z)[[1]] 


###################################################
### chunk number 32: 
###################################################
if(require("padi") && checkPADIserver("ets"))
tfwindow(forecasts(z)[[1]], start=c(1996,2)) 


