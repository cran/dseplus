/*UD
 * SccsID: @(#)getinfo.c	1.6		12/1/95 BOC
 * File:   getinfo.c
 * Author: Hope Pioro
 * Date  : Feb 96
 *
 * Get PADI time series
 *
 *LD
 * There are two functions in this file. A main function is used for
 * stand-alone purposes. The other is a control module for obtaining
 * PADI time series information from the PADI server.
 *
 *NO
 * A main function exists which is used for debugging purposes
 * and would act as a standalone version when compiled and linked.
 *
 * Compile and Link: make getinfo
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

private string_t app_name = "getinfo";



/*UF          getinfo function */
int
getinfo(user,password,server_name, objnam, dbname, msg, bufsize)

 /* INPUT variables */
char  **user;			/* user id */
char  **password;		/* password (not implemented) */
char  **server_name;			/* server name */
char  **objnam;			/* object name */
char  **dbname;		/* database name */
char  **msg;			/* message buffer */
long   *bufsize;		/* maxsize of buffer, 0 returned on success */

/*FD
 * This function retrieves informatin about time series data from the public PADI SERVER
 *
 *LD
 * This function receives a server name, an object name 
 * and  starting and ending date range and frequency and format and
 * the number of observations.
 *
 *
 * If not successful, it returns 
 * the appropriate message from the PADI SERVER .
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
    PadiInfoResult_tp info_result;
    PadiInfoArg_t arg_info;
    PadiInfo_t info;
    char    log_buf[128];
    string_t object_name;

#ifdef DEBUG
    printf("user,password=%s %s \n", *user, *password);
    printf("server,objnam=%s %s %s\n", *server_name, *objnam, *dbname);
#endif				/* DEBUG */

    (*msg)[0] = '\0';

    user_name = *user;
    pass_word = *password;
    server.name = *server_name;
    object_name = *objnam;
    database_name = *dbname;


    arg_info.user = user_name;
    arg_info.password = pass_word;
    arg_info.object_name = object_name;
    arg_info.db_name = database_name;

    if (!(info_result = PadiGetInfo(&server, &arg_info)))
    {
	sprintf(log_buf, "%s GetPadi(%s)", app_name, object_name);
	PadiError(stdout, log_buf, Padi_OUT_OF_MEMORY, Padi_FATAL);
    }

    if (info_result->status)
    {

#ifdef DEBUG
	sprintf(log_buf, "%s GetPadi(%s)", app_name, object_name);
	PadiError(stdout, log_buf, info_result->status, Padi_ERROR);
#endif

	if (strlen(PadiStatus(info_result->status)) < (size_t) * bufsize)
	    strcpy(*msg, PadiStatus(info_result->status));
	else
	{
	    strncpy(*msg, PadiStatus(info_result->status), (size_t) * bufsize);
	    (*msg)[(int) *bufsize] = '\0';
	    *bufsize = (long) strlen(*msg);
	}
    }
    else
    {

	info = info_result->info ;
 	printf("\n%s: \n", info.name); 
	printf("\t%s\n", info.desc);
	printf("\t%s, %s, %s, %s\n", info.class,
	info_result->info.type,
	info_result->info.access,
	info_result->info.frequency);
	printf("\t%s\n", info.format);
	printf("\t%s\n", info.start);
	printf("\t%s\n", info.end);



    }
    PadiFreeResult((PadiResult_tp) info_result);


#ifdef DEBUG
    printf("End function getinfo\n");
#endif				/* DEBUG */

    return (0);

}


#ifdef MAIN
/*UF    main function */
main(argc, argv)
int     argc;
char   *argv[];

/*FD
 * This main function is used as a stand alone version of the
 * getinfo modules which can be used on the command-line.
 *
 *LD
 * This function reads the input arguments from the command line,
 * then passes the server name and object name
 *
 * It calls getinfo(). If there was an error, it prints
 * the error text returned in the message buffer.
 *
 *RV
 * Returns  0.
 *
 *NO
 *
*/

{
    char   *server = NULL;	/* FAME database name */
    char   *objnam = NULL;	/* FAME object name */
    char   *database = NULL;	/* FAME database name */
    char   *user = NULL;	/* user name */
    char   *password = NULL;	/* password (not implemented) */
    char    msgbuffer[MAXMSGBUF + 1];	/* precision array for data -ap */
    char   *msg = msgbuffer;
    long    bufsize = MAXMSGBUF;/* will be reset to length returned */
    char    userbuf[L_cuserid+1];

    /* Validate command-line arguments */

    if (argc < 3)
    {
	fprintf(stderr, "%s\n","Usage: getinfo server object database ");
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


    user = (char *) cuserid(userbuf);

    password = gstrdup("(nil)");

#ifdef DEBUG
    printf("user,password=%s %s \n", user, password);
    printf("server,objnam=%s, %s\n", server, objnam);
    printf("database=%s\n", database);
#endif				/* DEBUG */


    if( !strcmp(database,"(null)") && *database != '/' )
    {
 	   printf("Error Message:  %s\n", "Database Name must be full path");
           return(0);
    }

    getinfo(&user,&password,&server, &objnam, &database, &msg, &bufsize);
    if (bufsize != 0)
    {
	printf("Error Message:  %s\n", msg);
    }

    if (objnam != NULL)
	free(objnam);
    if (server != NULL)
	free(server);
    if (password != NULL)
	free(password);

    exit(0);


}

#endif				/* MAIN */


