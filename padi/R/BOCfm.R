
################################################################################

# The functions in this file are for specific databases 
# and backwards compatibility at the Bank of Canada.

################################################################################

getfm <-function(dbname,series, starty=0,startp=0, endy=0,endp=0, 
    nobs=0,max.obs=2000, transformations=NULL, pad=FALSE)
{# This function is for backwards compatability at the BOC.
 # If dbname[1] begins with "ets" or "/home/ets" then the server is set 
 #   to ets and dbname is 
 #   not used. This will fail for a mix of ets and non ets dbs.
 # Otherwise a server is started to server the dbname from the local host.

 #  The rather crude mapping of startm=startp and endm=endp is used.
 #  This will work for monthly data and for 0 but may cause problems otherwise.

  if (("ets"       == substring(dbname[1],1,3)) | 
      ("/home/ets" == substring(dbname[1],1,9))  )
    {server <-  "ets"  #"bc"  # "padi"
     data <-getpadi(series,server=server,
            server.process="fame.server", dbname="", 
            starty=starty,startm=startp,startd=1,
            endy=endy,endm=endp,endd=1, nobs=nobs,max.obs=max.obs, 
            transformations=transformations, pad=pad)
    }
  else 
    {server <- Sys.info()[["nodename"]]
     data <-getpadi(series,server=server,
            server.process="fame.server",  dbname=dbname, 
            starty=starty,startm=startp,startd=1,
            endy=endy,endm=endp,endd=1, nobs=nobs,max.obs=max.obs, 
            transformations=transformations, pad=pad)
    }
  data
}

putfm <-function(data, dbname, seriesNames)
{# This function is for backwards compatability at the BOC.
 # A server is started on th default (local) host to server dbnames.
 putpadi(data,
      server.process="fame.server", dbname=dbname, series=seriesNames)
}
