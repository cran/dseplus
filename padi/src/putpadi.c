/*UD
 * SccsID: @(#)putpadi.c	1.5		12/1/95 BOC
 * File:   putpadi.c
 * Author: Hope Pioro
 * Date  : Sept 95
 *
 * Put PADI time series
 *
 *LD
 * There are two functions in this file. A main function is used for
 * stand-alone purposes. The other is a control module for obtaining
 * PADI time series from the PADI server.
 *
 *NO
 * A main function exists which is used for debugging purposes
 * and would act as a standalone version when compiled and linked.
 *
 * Compile and Link: make putpadi
 *
 * For all applications the code must be compiled with -DGMALLOC.  This
 * will call specialized malloc() and free() funtions. In the case of SPLUS, 
 * you must also use -DSPLUS. This will call S_alloc instead. (see makefile)
*/

/*IN*/
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "padi.h"
#include <signal.h>
#include <time.h>
#include <string.h>
#include <sys/param.h>
#include "padiutil.h"

/*CO*/
#ifndef L_cuserid
#define L_cuserid 9
#endif

#define DATLEN 16		/* max. number of char in date literal */
#define WRITE_LEN 100
#define MAXMSGBUF 10000
#define TMOUTDEFAULT 60		/* default timeout */
                    /* timeout of 50 is not long enough for some test cases */

private string_t app_name = "putpadi";


/*UF          putpadi function */
int
putpadi(user,password,server_name,  objnam, dbname,frequency, syear, smon, sday, nobs, datadbl, msg, bufsize, timeout)

 /* INPUT variables */
char  **user;			/* user name */
char  **password;		/*  password (not implemented) */
char  **server_name;		/* server name */
char  **objnam;			/*  object name */
char  **dbname;			/*  database name */
int   *syear,			/*  starting year */
       *smon,			/*  starting month */
       *sday;			/* starting day */
int   *frequency;		/* frequency of time series */
int   *nobs;			/* number of observations to write */
double *datadbl;		/* array to hold precision data  */
char  **msg;			/* message buffer */
int   *bufsize;		/* maxsize of buffer, 0 returned on success */
int   *timeout;		/* max time to wait for reply (seconds) */

/*FD
 * This function requests the PADI SERVER to write this data to a database.
 *
 *LD
 * This function receives a server name, database name, an object name,
 * frequency, starting year,month and day, the number of observations,
 * the data buffer and a message buffer.
 *
 * If the database does not exist, it will return one with an error.
 *
 * If the time series already existed in the database, it will be replaced
 * without warning.
 *
 *RV
 * Returns 0
 *
 *NO
*/

{

    PadiSize_t i;
    PadiServer_t server;
    PadiString_t user_name;
    PadiString_t pass_word;
    PadiNewSeries_t new_series;
    PadiResult_tp result;
    char    log_buf[128];
    size_t  nwrite;
    PadiPrecision_t *vector;
    string_t object_name;
    char    option_start[9];
    PadiString_t database_name;

#ifdef DEBUG
    printf("server,dbname,objnam=%s, %s , %s n", *server_name, *dbname, *objnam);
    printf("syear,smon,sday,nobs= %ld,%ld, %ld, %ld\n",
	   *syear, *smon, *sday, *nobs);
    for (i = 0; i < *nobs; i++)
	printf("%f\n", datadbl[i]);
#endif				/* DEBUG */


    (*msg)[0] = '\0';

    user_name = *user;
    pass_word = *password;
    server.name = *server_name;
    object_name = *objnam;
    database_name = *dbname;


    new_series.user = user_name;
    new_series.password = pass_word;
    new_series.dbname = database_name;
    new_series.series.info.name = object_name;
    new_series.series.info.type = "PRECISION";
    new_series.series.info.class = "SERIES";
    new_series.series.info.access = Padi_EMPTY_STR;
    new_series.series.info.frequency = freqstr(*frequency);
    new_series.series.info.desc = Padi_EMPTY_STR;
    new_series.series.info.start = Padi_EMPTY_STR;
    new_series.series.info.end = Padi_EMPTY_STR;
    new_series.series.info.format = Padi_EMPTY_STR;

    nwrite = (*nobs) ? *nobs : WRITE_LEN;
    new_series.series.range.nobs = nwrite;

   if(digitCount((int) *syear) !=  4)
	*syear = *smon = *sday = 0L;

    sprintf(option_start, "%4ld%2ld%2ld", *syear, *smon, *sday);
    new_series.series.range.start = option_start;
    new_series.series.range.end = Padi_EMPTY_STR;
    new_series.series.range.format = "SimpleSeries";
    new_series.series.range.missing_translation[0] = DEFAULT_NC;
    new_series.series.range.missing_translation[1] = DEFAULT_ND;
    new_series.series.range.missing_translation[2] = DEFAULT_NA;
    new_series.series.range.do_missing = TRUE;

    new_series.series.data.data_len = nwrite * sizeof(PadiPrecision_t);
    if (!(new_series.series.data.data_val = (PadiPrecision_t *) malloc(nwrite * sizeof(PadiPrecision_t))))
    {
	sprintf(log_buf, "%s (%s)", app_name, server.name);
	PadiError(stdout, log_buf, Padi_OUT_OF_MEMORY, Padi_FATAL);
    }
    vector = (PadiPrecision_t *) (new_series.series.data.data_val);
    for (i = 0; i < nwrite; i++)
    {
	vector[i] = datadbl[i];

#ifdef DEBUG
	printf("%f\n", vector[i]);
#endif

    }

    PadiSetTimeOut(*timeout, 0);
    if (!(result = PadiNewSeries(&server, &new_series)))
    {

	sprintf(log_buf, "%s PutPadi(%s)", app_name, object_name);
	PadiError(stdout, log_buf, Padi_OUT_OF_MEMORY, Padi_FATAL);
    }

    if (result->status)
    {

#ifdef DEBUG
	sprintf(log_buf, "%s PutPadi(%s)", app_name, object_name);
	PadiError(stdout, log_buf, result->status, Padi_ERROR);
#endif

	if (strlen(PadiStatus(result->status)) < (size_t) *bufsize)
	    strcpy(*msg, PadiStatus(result->status));
	else
	{
	    strncpy(*msg, PadiStatus(result->status), (size_t) *bufsize);
	    (*msg)[(int) *bufsize] = '\0';
	}
    	*bufsize = (int)strlen(*msg);
    }
    else
    {
    	*bufsize = (int)0;
    }


    PadiFreeResult((PadiResult_tp) result);

#ifdef DEBUG
    printf("********* End function putpadi\n");
#endif				/* DEBUG */

    return 0;
}

/* ---------- main program for testing purposes ---------*/

#ifdef MAIN
/*CO*/
#define DATLEN 16		/* max. number of char in date literal */

/*UF    main function */
main(argc, argv)
int     argc;
char   *argv[];

/*FD
 * This main function is used as a debug version in order to test the
 * putpadi modules.
 *
 *LD
 * This function reads the arguments from the command line,  
 * creates data, then passes the server name, database name, object name, 
 * frequency, the starting year, month and day, the number of 
 * observations, the data pointer and the message buffer
 * to the putpadi function.
 *
 * It calls putpadi(). If there was an error, it prints 
 * the error text returned in the message buffer.
 *
 *RV
 * Returns  0.
 *
 *NO
*/

{

    double  datadbl[10000];	/* precision array for data -ap */
    char   *server = NULL;	/* FAME server name */
    char   *dbname = NULL;	/* FAME database name */
    char   *objnam = NULL;	/* FAME object name */
    int    syear,
            smon,
            sday;		/* starting year and period to write data
				   from */
    int    nobs;		/* number of data points to write */
    int    frequency;		/* frequency of object */
    char    msgbuffer[MAXMSGBUF + 1];	/* precision array for data -ap */
    char   *msg = msgbuffer;
    int    bufsize = MAXMSGBUF;
    int     i;
    int    timeout =  TMOUTDEFAULT;   /* time to wait for reply */
    char   *user = NULL;	/* user name */
    char   *password = NULL;	/* password (not implemented) */
    char    userbuf[L_cuserid+1];
    char resolved_name[MAXPATHLEN + 1];
    char *result_name;


    /* Validate command-line arguments */

     if ((argc < 9)) 
     { 
 	fprintf(stderr, 
 	"Usage: putpadi server object database freqency st_year st_mon st_day nobs [user id]\n"); 
 	exit(0); 
     } 

    /* Get command-line arguments */
    server = gstrdup(argv[1]);
    objnam = gstrdup(argv[2]);
    dbname = gstrdup(argv[3]);
    frequency = atol(argv[4]);
    syear = atol(argv[5]);
    smon = atol(argv[6]);
    sday = atol(argv[7]);
    nobs = atol(argv[8]);

   if(digitCount((int) syear) !=  4)
	syear = smon  = sday  = 0L;

     if ((argc == 10)) 
     { 
        user = gstrdup(argv[9]);
     }
     else
        user = (char *) cuserid(userbuf);

    password = gstrdup("(nil)");


#ifdef DEBUG
    printf("dbname,objnam=%s, %s\n", dbname, objnam);
    printf("frequency,syear,smon,day,nobs= %ld, %ld, %ld, %ld, %ld\n",
	   frequency, syear, smon, sday, nobs);
#endif				/* DEBUG */

/* put data in the data buffer for testing purposes */
    for (i = 0; i < nobs; i++)
    {
	datadbl[i] = i + 2.0;

#ifdef DEBUG
	printf("i,datadbl=%d,%f\n", i, datadbl[i]);
#endif				/* DEBUG */
    }

/*
  Forcing putpadi to use a full pathname would ensure consistency between the
   client and server view of the directory structure, but causes problems with
   using putpadi in tests of the public mode fame server (which only has only
   a relative picture of the databases it is serving). 
    if( *dbname != '/' )
    {
 	   printf("Error Message:  %s\n", "Database Name must be full path");
           return(0);
    }
*/
    putpadi(&user,&password,&server,&objnam,&dbname, &frequency, &syear, &smon,
		&sday, &nobs, datadbl, &msg, &bufsize, &timeout);
    if(bufsize != 0)
    {
	printf("Error message: %s\n", msg );
    }

    gfree(server);
    gfree(dbname);
    gfree(objnam);

    return(0);

}

#endif				/* MAIN */

