/*
** Protocol for Application -- Database Interface
**
** Simple Server Functions
**
** Copyright 1995, 1996  Bank of Canada.
**
** The user of this software has the right to use, reproduce and distribute it.
** Bank of Canada makes no warranties with respect to the software or its 
** fitness for any particular purpose. The software is distributed by the Bank
** of Canada solely on an "as is" basis. By using the software, user agrees to 
** accept the entire risk of using this software.
**
*/

#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <pwd.h>
#include "padi.h"

/* data types */
typedef int booln_t;
typedef int date_t;
typedef int freq_t;
typedef int integer_t;
typedef float numeric_t;
typedef double precision_t;


private char *base_dir = NULL;

export PadiStatus_t initialize_1
PARAM1(PadiInitArg_tp, object)
{
    PadiStatus_t status = Padi_SUCCEED;

    base_dir = object->object_dbname;;
    return status;
}


export PadiInfoResult_tp getinfo_1
PARAM2(PadiInfoArg_tp, argp, CLIENT *, rqstp)
{
    PadiInfoResult_tp result;
    return result;
}

export PadiSeriesResult_tp getseries_1
PARAM2(PadiRangeArg_tp, arg, CLIENT *, rqstp)
{
    PadiSeriesResult_tp result;
    return result;
}


export PadiResult_tp putseries_1
PARAM2(PadiSeries_tp, series, CLIENT *, rqstp)
{
    PadiResult_tp result;
    return result;
}


export PadiResult_tp newseries_1
PARAM2(PadiNewSeries_tp, new, CLIENT *, rqstp)
{
    PadiResult_tp result;
    return result;
}


export PadiResult_tp destroy_1
PARAM2(PadiDestroyArg_tp, argp, CLIENT *, rqstp)
{
    PadiResult_tp result;

    return result;
}



export PadiStatus_t terminate_1
PARAM1(PadiTermArg_tp, object)
{
    PadiStatus_t status = Padi_SUCCEED;
    return status;
}

export char *status_1
PARAM1(int, status)
{
    static char    buf[64];
    switch (status)
    {
      case 0:
	return "Success";
    }

    sprintf(buf, "Unknown error");
    return buf;
}

export PadiSeriesResult_tp getlocal_1
PARAM2(PadiRangeArg_tp, local, CLIENT *, rqstp)
{
    PadiString_t user = local->user;
    PadiString_t password = local->password;
    PadiRange_tp in_range = &(local->range);
    PadiString_t object_name = local->object_name;
    PadiString_t db_name = local->db_name;
    PadiRange_tp data_range;
    struct passwd *ps;
    char   *svcuser;
    char    name_buf[Padi_MAX_STR_LEN + 1];	/* fame writes to this
						   string! */
    char    dbpath[Padi_MAX_STR_LEN + 1];

    char    desc[Padi_MAX_STR_LEN + 1];
    char    class[Padi_MAX_STR_LEN + 1];
    char    type[Padi_MAX_STR_LEN + 1];
    char    freq[Padi_MAX_STR_LEN + 1];
    char    access[Padi_MAX_STR_LEN + 1];
    char    format[Padi_MAX_STR_LEN + 1];
    char    tmpbuf[Padi_MAX_STR_LEN + 1];
    int     start,
            end,
            nobs,
            numobs;
    int     in_syear,
            in_sprd,
            in_eyear,
            in_eprd;
    size_t  struct_len,
            data_len,
            buf_len;
    char   *s;
    PadiString_t *strp;
    PadiPrecision_t *vector;
    PadiSeriesResult_tp result;
    PadiInfo_tp info;
    DIR    *dir;
    FILE   *fp;
    int     i,
            j;

    struct_len = (size_t) ALIGN(align_t, sizeof(PadiSeriesResult_t));
    if (!(result = (PadiSeriesResult_tp) malloc(struct_len)))
	return NULL;

    bzero((char *) result, struct_len);
    data_range = &(result->series.range);
    data_range->start = Padi_EMPTY_STR;
    data_range->end = Padi_EMPTY_STR;
    data_range->format = Padi_EMPTY_STR;
    info = &(result->series.info);
    for (strp = (PadiString_t *) info;
	 strp < ((PadiString_t *) info + Padi_NUM_INFO_STR);
	 strp++)
	*strp = Padi_EMPTY_STR;

    ps = getpwuid(geteuid());
    svcuser = ps->pw_name;

#ifdef DEBUG
    printf("svcuser %s  passed  %s \n", svcuser, user);
#endif

    if (strcmp(svcuser, user))
    {
	result->status = Padi_UNATH_USER;
	return result;

    }

    if (!(local->object_name && *(local->object_name)))
    {
	result->status = Padi_MISSING_OBJECT;
	return result;

    }

    if (!(local->db_name && *(local->db_name)))
    {
	result->status = Padi_MISSING_OBJECT;
	return result;
    }

    /* check that dbname is under base_dir */
    if( !strstr( db_name, base_dir))
    {
	result->status = Padi_MISSING_DBNAME;
	return result;
    }

    /* build the full path */
    strcpy(dbpath, db_name);


    /* open the object's directory */
    dir = opendir(dbpath);

    if (dir == NULL)
    {
	result->status = Padi_MISSING_OBJECT;
	return result;
    }
    else
	closedir(dir);

#ifdef DEBUG
    printf("OpenDb : %d\n", result->status);
#endif

    strcat(dbpath, "/");
    strcat(dbpath, object_name);
    if ((fp = fopen(dbpath, "r")) == NULL)
    {

#ifdef DEBUG
	perror("Open");
#endif

	result->status = Padi_MISSING_OBJECT;
	return result;
    }

#ifdef DEBUG
    printf("OpenObject : %d\n", result->status);
#endif

    fgets(name_buf, Padi_MAX_STR_LEN, fp);
    if (strncasecmp(name_buf, local->object_name, (int) strlen(local->object_name)))
    {
	result->status = Padi_UNKNOWN_OBJECT;
	return result;
    }
    fgets(desc, Padi_MAX_STR_LEN, fp);

    fgets(tmpbuf, Padi_MAX_STR_LEN, fp);
    sscanf(tmpbuf, "%s*", type);

    fgets(tmpbuf, Padi_MAX_STR_LEN, fp);
    sscanf(tmpbuf, "%s*", class);

    fgets(tmpbuf, Padi_MAX_STR_LEN, fp);
    sscanf(tmpbuf, "%s*", freq);
    fgets(tmpbuf, Padi_MAX_STR_LEN, fp);
    sscanf(tmpbuf, "%s*", access);
    fgets(tmpbuf, Padi_MAX_STR_LEN, fp);
    sscanf(tmpbuf, "%s*", format);
    fgets(tmpbuf, Padi_MAX_STR_LEN, fp);
    start = atoi(tmpbuf);
    fgets(tmpbuf, Padi_MAX_STR_LEN, fp);
    end = atoi(tmpbuf);
    fgets(tmpbuf, Padi_MAX_STR_LEN, fp);
    nobs = atoi(tmpbuf);

    if( (digitCount(start) != 4)  || ( digitCount(end) != 4) )
    {
	result->status = Padi_UNSUPPORTED_DATE;
	return result;
    }

    if (strcmp(class, "SERIES"))
    {
	result->status = Padi_UNSUPPORTED_CLASS;
	return result;
    }

    if (strcmp(type, "PRECISION"))
    {
	result->status = Padi_UNSUPPORTED_TYPE;
	return result;
    }

    if (strcmp(freq, "ANNUAL"))
    {
	result->status = Padi_UNSUPPORTED_FREQ;
	return result;
    }
    else
    {

	if (*(in_range->start))
	{
	    in_syear = atoi(in_range->start);
	    if (in_syear < start)
	    {
		in_syear = -1;
	    }

	}
	else
	{
	    in_syear = -1;
	}

	if (*(in_range->end))
	{
	    in_eyear = atoi(in_range->end);
	    if (in_eyear > end)
	    {
		in_eyear = -1;
	    }

	}
	else
	{
	    in_eyear = -1;
	}
    }

    if (in_range->nobs)
    {
	numobs = in_range->nobs;
    }
    else
    {
	numobs = -1;
    }

    if (numobs != -1)
    {
	if (in_eyear != -1)
	{
	    in_syear = in_eyear - numobs + 1;
	}
	else if (in_syear != -1)
	{
	    in_eyear = in_syear + numobs - 1;
	}
	else
	{
	    in_eyear = end;
	    in_syear = in_eyear - numobs + 1;
	}
    }
    else
    {
	if (in_syear == -1)
	{
	    in_syear = start;
	}
	if (in_eyear == -1)
	{
	    in_eyear = end;
	}
	numobs = in_eyear - in_syear + 1;
    }
    


    struct_len = (size_t) ALIGN(align_t, sizeof(PadiSeriesResult_t));
    data_len = (size_t) ALIGN(align_t, numobs * sizeof(PadiPrecision_t));
    buf_len = (size_t) ALIGN(align_t, Padi_MAX_STR_LEN + 1);
    if (!(result = (PadiSeriesResult_tp) realloc((char *) result,
		struct_len + data_len + (3 + Padi_NUM_INFO_STR) * buf_len)))
	return NULL;

    s = (char *) result;
    bzero(s, struct_len + data_len + (3 + Padi_NUM_INFO_STR) * buf_len);
    data_range = &(result->series.range);
    data_range->nobs = numobs;
    s += struct_len;
    result->series.data.data_len = numobs;
    result->series.data.data_val = (PadiPrecision_t *) s;
    vector = (PadiPrecision_t *) s;
    s += data_len;
    data_range->start = s;
    s += buf_len;
    data_range->end = s;
    s += buf_len;
    data_range->format = s;
    s += buf_len;
    info = &(result->series.info);
    for (strp = (PadiString_t *) info;
	 strp < ((PadiString_t *) info + Padi_NUM_INFO_STR);
	 strp++, s += buf_len)
	*strp = s;


    strcpy(info->name, local->object_name);
    strcpy(info->desc, desc);
    strcpy(info->class, class);
    strcpy(info->type, type);
    strcpy(info->access, access);
    strcpy(info->frequency, freq);
    strcpy(info->format, format);
    sprintf(info->start, "%d0101", start);
    sprintf(info->end, "%d0101", end);

    sprintf(data_range->start, "%d0101", in_syear);
    sprintf(data_range->end, "%d0101", in_eyear);

    if (numobs > 0)
    {
	/* read data */
	for (i = start, j = 0; i <= end; i++)
	{
	    fgets(tmpbuf, Padi_MAX_STR_LEN, fp);
	    if (i >= in_syear && i <= in_eyear)
	    {
		if( i >= start && i <= end)
		{
		    sscanf(tmpbuf, "%lf*\n", (double *) &vector[j]);
		    j++;
		}
	    }
	}
    }


    /* store the range translation */
    data_range->do_missing = in_range->do_missing;
    memcpy((char *) (data_range->missing_translation),
	   (char *) (in_range->missing_translation),
	   sizeof(in_range->missing_translation));

    /* close primary database */
    fclose(fp);

    return result;
}

export PadiResult_tp putlocal_1
PARAM2(PadiNewSeries_tp, new, CLIENT *, rqstp)
{
    PadiString_t user = new->user;
    PadiString_t password = new->password;
    PadiSeries_tp series = &(new->series);
    PadiString_t dbname = new->dbname;
    PadiInfo_tp info = &(series->info);
    PadiString_t object_name = info->name;
    PadiRange_tp in_range = &(series->range);
    char    dbpath[Padi_MAX_STR_LEN + 1];
    PadiResult_tp result;
    size_t  struct_len;
    precision_t *pvector = (precision_t *) (series->data.data_val);
    char   *svcuser;
    struct passwd *ps;
    DIR    *dir;
    FILE   *fp;
    int     i;


    struct_len = sizeof(*result);
    if (!(result = (PadiResult_tp) malloc(struct_len)))
	return NULL;
    bzero((char *) result, struct_len);

    ps = getpwuid(geteuid());
    svcuser = ps->pw_name;
    if (strcmp(svcuser, user))
    {
	result->status = Padi_UNATH_USER;
	return result;

    }

    if (!(new->dbname && *(new->dbname)))
    {
	result->status = Padi_MISSING_DBNAME;
	return result;

    }

    /* check that dbname is under base_dir */
    if( !strstr( dbname, base_dir))
    {
	result->status = Padi_MISSING_DBNAME;
	return result;
    }

    /* build the full path */
    strcpy(dbpath, dbname);

    /* open the object's directory (if not exist, create it) */
    dir = opendir(dbpath);
    if (dir == NULL)
    {
	if (mkdir(dbpath, S_IRWXU))
	{

#ifdef DEBUG
	    perror("mkdir");
#endif

	    result->status = Padi_MISSING_OBJECT;
	    return result;
	}
    }
    else
	closedir(dir);


    strcat(dbpath, "/");
    strcat(dbpath, object_name);
    if ((fp = fopen(dbpath, "w")) == NULL)
    {

#ifdef DEBUG
	perror("Open");
#endif

	result->status = Padi_MISSING_OBJECT;
	return result;
    }

    if (strcmp(info->class, "SERIES"))
    {
	result->status = Padi_UNSUPPORTED_CLASS;
	return result;
    }

    if (strcmp(info->type, "PRECISION"))
    {
	result->status = Padi_UNSUPPORTED_TYPE;
	return result;
    }

    if (strcmp(info->frequency, "ANNUAL"))
    {
	result->status = Padi_UNSUPPORTED_FREQ;
	return result;
    }

    if(digitCount( atoi(in_range->start) ) != 4)
    {
	result->status = Padi_UNSUPPORTED_DATE;
	return result;
    }
    fprintf(fp, "%s\n", info->name);
    if (*(info->desc))
	fprintf(fp, "%s\n", info->desc);
    else
	fprintf(fp, "NO DESCRIPTION FOR %s\n", info->name);

    fprintf(fp, "%s\n", info->type);
    fprintf(fp, "%s\n", info->class);
    fprintf(fp, "%s\n", info->frequency);
    if (*(info->access))
	fprintf(fp, "%s\n", info->access);
    else
	fprintf(fp, "%s\n", "READ");
    if (*(info->format))
	fprintf(fp, "%s\n", info->format);
    else
	fprintf(fp, "%s\n", "TIME");
    fprintf(fp, "%d\n", atoi(in_range->start));
    fprintf(fp, "%d\n", atoi(in_range->start) + in_range->nobs - 1);
    fprintf(fp, "%d\n", in_range->nobs);

    for (i = 0; i < in_range->nobs; i++)
    {
	fprintf(fp, "%f\n", (double) pvector[i]);
    }

    /* close primary database */
    fclose(fp);

    return result;
}
