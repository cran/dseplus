#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "padi.h"


int digitCount(int value);

char *
freqstr( int f)
{
    if(f == 1)
	return("ANNUAL");
    else if(f == 4)
	return("QUARTERLY");
    else if(f == 12)
	return("MONTHLY");
    else if(f == 52)
	return("WEEKLY(WEDNESDAY)");
    else
	return("");

}

int
retfreq( char *frequency)
{
    if( strncmp(frequency,"ANNUAL",strlen("ANNUAL"))==0 )
	return( 1);
    else if(strncmp(frequency,"QUARTERLY",strlen("QUARTERLY"))==0)
	return(4);
    else if(strncmp(frequency,"MONTHLY",strlen("MONTHLY"))==0)
	return(12);
    else if(strncmp(frequency,"WEEKLY",strlen("WEEKLY"))==0)
	return(52);
    else
        return(0);
}


/* 	Sept. 29, 1997
**   	y2k year validation.  
**	Ensure year has four integer digits.
**
**	RETURN 0 if FAIL
*/ 

int
retyear( char *date)
{
    char buffer[5];
    int len;

    if(!date || *date == NULL || (strlen(date) < 4) )
	return 0;
     strncpy(buffer,date,4);
       buffer[4] = '\0';
    for(len = 0;isdigit((int) buffer[len]); len++);
    return( len == 4? atoi(buffer): 0);


}

int
retmonth( char *date)
{
    char buffer[3];

    if(strlen(date)>=(size_t)6)
    {
    strncpy(buffer,date + 4,2);
    buffer[2] = '\0';
    return(atoi(buffer));
    }
    else
    {
	return(0);
    }

}

int
retday( char *date)
{
    char buffer[3];

    if(strlen(date)>=(size_t)8)
    {
    strncpy(buffer,date + 6,2);
    buffer[2] = '\0';
    return(atoi(buffer));
    }
    else
    {
	return(0);
    }

}

char *
gstrdup(char *s)
{
char *p = NULL;

    if( !s )
	return NULL;
    if (!(p = malloc((strlen(s) + (size_t)1) * sizeof(char))))
	return NULL;
    if( !strcpy( p,s ) )
	return NULL;
    return p;

}

int digitCount(int value)
{
   char string[256];
   int i;
   sprintf(string, "%d", value);
   for(i = 0; string[i]; i++)
     if(!isdigit((int) string[i]))
	break;
   return i;
}



