/*UD
 * SccsID: @(#)getpadi.c	1.6		12/1/95 BOC
 * File:   getpadi.c
 * Author: Hope Pioro
 * Date  : Sept 95
 *
 * Get PADI time series
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
 * Compile and Link: make getpadi
 *
 * For all applications the code must be compiled with -DGMALLOC.  This
 * will call specialized malloc() and free() funtions. In the case of SPLUS,
 * you must also use -DSPLUS. This will call S_alloc instead. (see makefile)
 *
*/

/*IN*/
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <time.h>
#include <string.h>
#include <sys/param.h>
#include <unistd.h>
#include "padi.h"
#include "padiutil.h"

/*CO*/
#ifndef L_cuserid
#define L_cuserid 9
#endif

#define DATLEN 16		/* max. number of char in date literal */
#define MAXMSGBUF 10000		/* max len of message buff */
#define MAXDATA 10000		/* maximum number of data points */
#define TMOUTDEFAULT 90		/* default timeout */
                    /* timeout of 50 is not long enough for some test cases */


private string_t app_name = "getpadi";



/*UF          getpadi function */
int
getpadi(user,password,server_name, objnam, dbname, syear, smon, sday, eyear, emon, eday, frequency,
	     nobs, datadbl, maxsize, msg, bufsize, timeout)

 /* INPUT variables */
char  **user;			/* user id */
char  **password;		/* password (not implemented) */
char  **server_name;			/* server name */
char  **objnam;			/* object name */
char  **dbname;		/* database name */
int   *syear,			/* starting year, if 0 all data returned */
       *smon,			/* starting month, if 0 all data returned  */
       *sday;			/* starting day, if 0 all data returned  */
int   *eyear,			/* ending year, may be 0 if nobs != 0 */
       *emon,			/* ending month, may be 0 if nobs != 0 */
       *eday;			/* ending day, may be 0 if nobs != 0 */
int   *frequency;		/* frequency of time series */
int   *nobs;			/* number of observations to read; set to -1
				   if any error occurs, if 0 all data return */
int   *maxsize;		/* max. size set for datadbl array in calling
				   program */
char  **msg;			/* message buffer */
int   *bufsize;		/* maxsize of buffer, 0 returned on success */
int   *timeout;		/* max time to wait for reply (seconds) */

 /* OUTPUT variable */
double *datadbl;		/* array to hold precision data */
/*FD
 * This function retrieves time series data from the PADI SERVER
 *
 *LD
 * This function receives a server name, an object name and maximum observations
 * and optionally receives starting and ending date range and
 * the number of observations.
 *
 * It attempts to read the object from the server.  If successful, it returns
 * the data, the frequency, starting and ending dates returned,
 * the number of obervations returned and a message.
 *
 * If the starting date and ending date and number of observations is set to
 * zero, it will return the whole series.  If the starting and ending dates
 * are zero but number of observations is set, the LAST numobs observations
 * are returned.  If the starting date is sent, ending date is zero and numobs
 * is non-zero, the FIRST numobs starting from start date are returned.
 *
 * If the number of observations is greater than maxsize,
 * it only returns maxsize number of observations.
 *
 * If not successful, it returns -1 as the number of observations,
 * the appropriate message from the PADI SERVER and bufsize is set to
 * the length of the message in the buffer.
 *
 *RV
 * Returns 0
 *
 *NO
*/

{
    PadiSize_t j;
    PadiServer_t server;
    PadiString_t database_name;
    PadiString_t user_name;
    PadiString_t pass_word;
    PadiRangeArg_t range_arg;
    PadiSeriesResult_tp series_result;
    PadiNewSeries_t new_series;
    PadiSeries_tp series = &(new_series.series);
    char    log_buf[128];
    PadiPrecision_t *vector;
    string_t object_name;
    char    option_start[9];
    char    option_end[9];

#ifdef DEBUG
    printf("user,password=%s %s \n", *user, *password);
    printf("server,objnam=%s %s %s\n", *server_name, *objnam, *dbname);
    printf("syear,smon,sday,eyear,emon,eday,nobs=%ld, %ld, %ld, %ld,%ld, %ld, %ld\n",
	   *syear, *smon, *sday, *eyear, *emon, *eday, *nobs);
#endif				/* DEBUG */

    (*msg)[0] = '\0';

    user_name = *user;
    pass_word = *password;
    server.name = *server_name;
    object_name = *objnam;
    database_name = *dbname;


	if(digitCount((int) *syear) == 4)
	{
	    sprintf(option_start, "%4d%2d%2d", *syear, *smon, *sday);
	    series->range.start = option_start;
	    range_arg.range.start = option_start;
	}
	else
	{
	    series->range.start = "19931231";
	    range_arg.range.start = 0;
	}
	if(digitCount((int) *eyear) == 4)
	{
	    sprintf(option_end, "%4d%2d%2d", *eyear, *emon, *eday);
	    series->range.end = option_end;
	    range_arg.range.end = option_end;
	}
	else
	{
	    series->range.end = 0;
	    range_arg.range.end = 0;
	}


    range_arg.user = user_name;
    range_arg.password = pass_word;
    range_arg.object_name = object_name;
    range_arg.db_name = database_name;
    range_arg.range.nobs = *nobs;
    range_arg.range.format = 0;
    range_arg.range.missing_translation[0] = DEFAULT_NC;
    range_arg.range.missing_translation[1] = DEFAULT_ND;
    range_arg.range.missing_translation[2] = DEFAULT_NA;
    range_arg.range.do_missing = TRUE;

    PadiSetTimeOut(*timeout, 0);
    if (!(series_result = PadiGetSeries(&server, &range_arg)))
    {
	sprintf(log_buf, "%s GetPadi(%s)", app_name, object_name);
	PadiError(stdout, log_buf, Padi_OUT_OF_MEMORY, Padi_FATAL);
    }

    if (series_result->status)
    {

#ifdef DEBUG
	sprintf(log_buf, "%s GetPadi(%s)", app_name, object_name);
	PadiError(stdout, log_buf, series_result->status, Padi_ERROR);
#endif

	if (strlen(PadiStatus(series_result->status)) < (size_t) * bufsize)
	    strcpy(*msg, PadiStatus(series_result->status));
	else
	{
	    strncpy(*msg, PadiStatus(series_result->status), (size_t) * bufsize);
	    (*msg)[(int) *bufsize] = '\0';
	    *bufsize = (int) strlen(*msg);
	}
    }
    else
    {

#ifdef DEBUG
	printf("\n%s: \n", object_name);

	printf("\t%s\n", series_result->series.info.desc);
	printf("\t%s, %s, %s, %s\n", series_result->series.info.class,
	       series_result->series.info.type,
	       series_result->series.info.access,
	       series_result->series.info.frequency);
	if (*(series_result->series.info.start))
	    printf("\t%s to %s\n", series_result->series.info.start,
		   series_result->series.info.end);
	else
	    printf("\n");

	printf("\tdata from %s to %s:\n",
	       series_result->series.range.start,
	       series_result->series.range.end);
#endif

	*frequency = (int) retfreq(series_result->series.info.frequency);
	*syear     = (int) retyear(series_result->series.range.start);
	*smon 	   = (int) retmonth(series_result->series.range.start);
	*sday 	   = (int) retday(series_result->series.range.start);
	*eyear 	   = (int) retyear(series_result->series.range.end);
	*emon      = (int) retmonth(series_result->series.range.end);
	*eday      = (int) retday(series_result->series.range.end);
	*nobs      = series_result->series.range.nobs;
	vector     = (PadiPrecision_t *) (series_result->series.data.data_val);
	for (j = 0; j < series_result->series.range.nobs &&
	     j < *maxsize; j++)
	{
	    datadbl[j] = (double) vector[j];
	}

#ifdef DEBUG
	for (j = 0; j < series_result->series.range.nobs &&
	     j < *maxsize; j++)
	{
	    printf("\t%d:\t%f\n", j + 1, (double) vector[j]);
	}
#endif

	*bufsize = (int) 0;

    }
    PadiFreeResult((PadiResult_tp) series_result);


#ifdef DEBUG
    printf("End function getpadi\n");
#endif				/* DEBUG */

    return (0);

}


#ifdef MAIN
/*UF    main function */
int main(argc, argv)
int     argc;
char   *argv[];

/*FD
 * This main function is used as a stand alone version of the
 * getpadi modules which can be used on the command-line.
 *
 *LD
 * This function reads the input arguments from the command line,
 * then passes the server name and object name
 * and optionally the starting and ending date range and
 * the number of observations to the getpadi function.
 *
 * It calls getpadi(). If there was an error, it prints
 * the error text returned in the message buffer.
 *
 *RV
 * Returns  0.
 *
 *NO
 *
*/

{
    double  datadbl[MAXDATA];	/* precision array for data -ap */
    char   *server = NULL;	/* FAME database name */
    char   *objnam = NULL;	/* FAME object name */
    char   *database = NULL;	/* FAME database name */
    char   *user = NULL;	/* user name */
    char   *password = NULL;	/* password (not implemented) */
    int    frequency = 0;
    int    syear,
            smon,
            sday;		/* starting year and period to read data from */
    int    eyear,
            emon,
            eday;		/* ending year and period to read data to */
    int    nobs;		/* number of data points to read */
    int    maxsize = MAXDATA;	/* max.  size of array  -ap */
    char    msgbuffer[MAXMSGBUF + 1];	/* precision array for data -ap */
    char   *msg = msgbuffer;
    int    bufsize = MAXMSGBUF;/* will be reset to length returned */
    int     i;
    int    timeout =  TMOUTDEFAULT;   /* time to wait for reply */
    char    userbuf[L_cuserid+1];
    int     j, year;

    /* Validate command-line arguments */

    if ((argc == 1) || ((argc > 4) && (argc < 11)))	/* -ap */
    {
	fprintf(stderr, "%s%s\n",
		"Usage: getpadi server object database st_year st_mon st_day end_year end_mon end_day nobs [user id]\n",
		"where st_year, st_mon, st_day, end_year, end_mon, end_day, nobs and user id are optional");
	exit(0);
    }

    /* Get command-line arguments */

    server = gstrdup(argv[1]);
    objnam = gstrdup(argv[2]);
    if (argc > 3)		/* database supplied */
    {
        database = gstrdup(argv[3]);
    }
    else
    {
        database = gstrdup("(nil)");
    }

    if (argc > 5)		/* optional arguments supplied */
    {
	syear = atol(argv[4]);
	smon = atol(argv[5]);
	sday = atol(argv[6]);
	eyear = atol(argv[7]);
	emon = atol(argv[8]);
	eday = atol(argv[9]);
	nobs = atol(argv[10]);
        password = gstrdup("(nil)");

	if(digitCount((int) syear) != 4)
	  syear = smon = sday = 0L;

	if(digitCount((int) eyear) != 4)
	  eyear = emon = eday = 0L;

    }
    else			/* only required arguments supplied */
    {
	syear = 0;
	smon = 0;
	sday = 0;
	eyear = 0;
	emon = 0;
	eday = 0;
	nobs = 0;
    }

     if ((argc == 12)) 
     { 
        user = gstrdup(argv[11]);
     }
     else
        user = (char *) cuserid(userbuf);

    password = gstrdup("(nil)");

#ifdef DEBUG
    printf("user,password=%s %s \n", user, password);
    printf("server,objnam=%s, %s\n", server, objnam);
    printf("database=%s\n", database);
    printf(
	   "syear,smon,sday,eyear,emon,eday,nobs=%ld, %ld, %ld, %ld,%ld, %ld, %ld\n",
	   syear, smon, sday, eyear, emon, eday, nobs);
#endif				/* DEBUG */


    if( !strcmp(database,"(null)") && *database != '/' )
    {
 	   printf("Error Message:  %s\n", "Database Name must be full path");
           return(0);
    }

    getpadi(&user,&password,&server, &objnam, &database, &syear, &smon, &sday, &eyear, &emon, &eday,
		 &frequency, &nobs,
		 datadbl, &maxsize, &msg, &bufsize, &timeout);
    if (bufsize != 0)
    {
	printf("Error Message:  %s\n", msg);
    }
    else
    {
#ifndef X11ARIMA
	printf("syear %d 	smon %d sday %d\n", syear, smon, sday);
	printf("eyear %d 	emon %d 	eday %d\n", eyear, emon, eday);
	printf("frequency %d\n", frequency);
	printf("nobs	     %d\n", nobs);
	for (i = 0; i < nobs && i < maxsize; i++)
	    printf("    %f\n", datadbl[i]);
#endif				/* not X11ARIMA */
#ifdef X11ARIMA
     /* this is a format used as input to X11ARIMA at STC for Monthly data!! */
        year = syear;
	printf("%-8s%4d", objnam, year);
	j = smon - 1;
	for (i = 1; i <= j; i++) printf("           ");
	for (i = 0; i < nobs && i < maxsize; i++)
	    {j = j + 1;
	     printf("%11d", (int) datadbl[i]);
	     if (j == 6) printf("\n            ");
	     if (j == 12) 
	        {printf("\n");
	         j = 0;
		 if (i < nobs) printf("%-8s%4d", objnam, year);
		 year = year + 1;
		}
	    }
	 printf("\n");
#endif				/* X11ARIMA */
    }

    if (objnam != NULL)
	gfree(objnam);
    if (server != NULL)
	gfree(server);
    if (password != NULL)
	gfree(password);

    exit(0);


}

#endif				/* MAIN */


