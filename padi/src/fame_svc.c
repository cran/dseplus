/* 
** Protocol for Application -- Database Interface
** 
** Fame Server Functions
**
** Copyright 1995, 1996  Bank of Canada.
** The user of this software has the right to use, reproduce and distribute it.
** Bank of Canada makes no warranties with respect to the software or its 
** fitness for any particular purpose. The software is distributed by the Bank
** of Canada solely on an "as is" basis. By using the software, user agrees to 
** accept the entire risk of using this software.
**
*/

#include <hli.h>
#include <padi.h>
#include <errno.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <pwd.h>

/* #define  DEBUG 1 */

/* FAME data types */
typedef int booln_t;
typedef int date_t;
typedef int freq_t;
typedef int integer_t;
typedef float numeric_t;
typedef double precision_t;

typedef struct freqname_t {
  int freq;
  string_t name;
} freqname_t, *freqname_tp;

private freqname_t freqname[] = {
  HDAILY,	"DAILY",
  HBUSNS,	"BUSINESS",
  HWKSUN,	"WEEKLY(SUNDAY)",
  HWKMON,	"WEEKLY(MONDAY)",
  HWKTUE,	"WEEKLY(TUESDAY)",
  HWKWED,	"WEEKLY(WEDNESDAY)",
  HWKTHU,	"WEEKLY(THURSDAY)",
  HWKFRI,	"WEEKLY(FRIDAY)",
  HWKSAT,	"WEEKLY(SATURDAY)",
  HTENDA,	"TENDAY",
  HWASUN,	"BIWEEKLY(ASUNDAY)",
  HWAMON,	"BIWEEKLY(AMONDAY)",
  HWATUE,	"BIWEEKLY(ATUESDAY",
  HWAWED,	"BIWEEKLY(AWEDNESDAY)",
  HWATHU,	"BIWEEKLY(ATHURSDAY)",
  HWAFRI,	"BIWEEKLY(AFRIDAY)",
  HWASAT,	"BIWEEKLY(ASATURDAY)",
  HWBSUN,	"BIWEEKLY(BSUNDAY)",
  HWBMON,	"BIWEEKLY(BMONDAY)",
  HWBTUE,	"BIWEEKLY(BTUESDAY)",
  HWBWED,	"BIWEEKLY(BWEDNESDAY)",
  HWBTHU,	"BIWEEKLY(BTHURSDAY)",
  HWBFRI,	"BIWEEKLY(BFRIDAY)",
  HWBSAT,	"BIWEEKLY(BSATURDAY)",
  HTWICM,	"TWICEMONTHLY",
  HMONTH,	"MONTHLY",
  HBMNOV,	"BIMONTHLY(NOV)",
  HBIMON,	"BIMONTHLY(DEC)",
  HQTOCT,	"QUARTERLY(OCT)",
  HQTNOV,	"QUARTERLY(NOV)",
  HQTDEC,	"QUARTERLY(DEC)",
  HQTDEC,	"QUARTERLY",
  HANJAN,	"ANNUAL(JAN)",
  HANFEB,	"ANNUAL(FEB)",
  HANMAR,	"ANNUAL(MAR)",
  HANAPR,	"ANNUAL(APR)",
  HANMAY,	"ANNUAL(MAY)",
  HANJUN,	"ANNUAL(JUN)",
  HANJUL,	"ANNUAL(JUL)",
  HANAUG,	"ANNUAL(AUG)",
  HANSEP,	"ANNUAL(SEP)",
  HANOCT,	"ANNUAL(OCT)",
  HANNOV,	"ANNUAL(NOV)",
  HANDEC,	"ANNUAL(DEC)",
  HANDEC,	"ANNUAL",
  HSMJUL,	"SEMIANNUAL(JUL)",
  HSMAUG,	"SEMIANNUAL(AUG)",
  HSMSEP,	"SEMIANNUAL(SEP)",
  HSMOCT,	"SEMIANNUAL(OCT)",
  HSMNOV,	"SEMIANNUAL(NOV)",
  HSMDEC,	"SEMIANNUAL(DEC)",
  HSMDEC,	"SEMIANNUAL",
  HAYPP,	"YPP",
  HAPPY,	"PPY",
  HSEC,		"SECONDLY",
  HMIN,		"MINUTELY",
  HHOUR,	"HOURLY",
  HCASEX,	"CASE",
  HUNDFX,	NULL
};

private string_t FreqToStr PARAM1(int, freq) 
{
  freqname_tp p;

  for (p = freqname; p->name; p++)
    if (freq == p->freq) break;

  return (p->name) ? p->name : Padi_EMPTY_STR;
}

private int StrToFreq PARAM1(string_t, s) 
{
  freqname_tp p;

  for (p = freqname; p->name; p++)
    if (!strcmp(s, p->name)) break;

  return p->freq;
}

typedef struct typename_t {
  int type;
  string_t name;
} typename_t, *typename_tp;

private typename_t typename[] = {
HNUMRC,    "NUMERIC"     ,
HNAMEL,    "NAMELIST"    ,
HBOOLN,    "BOOLEAN"     ,
HSTRNG,    "STRING"      ,
HPRECN,    "PRECISION"   ,
HDATE ,    "General DATE",
HRECRD,	   "RECORD"      ,
HUNDFT,    NULL
};

private string_t TypeToStr PARAM1(int, type) 
{
  typename_tp p;

  for (p = typename; p->name; p++)
    if (type == p->type) break;

  /* date types are frequencies */
  return (p->name) ? p->name : FreqToStr(type);
}

private int StrToType PARAM1(string_t, s) 
{
  typename_tp p;

  for (p = typename; p->name; p++)
    if (!strcmp(s, p->name)) break;

  return p->type;
}

typedef struct classname_t {
  int class;
  string_t name;
} classname_t, *classname_tp;

private classname_t classname[] = {
HSERIE,	"SERIES"   ,
HSCALA,	"SCALAR"   ,
HFRMLA,	"FORMULA"  ,
HITEM ,	"ITEM	"    ,
HGLNAM,	"GLNAME"   ,
HGLFOR,	"GLFORMULA",
     0, NULL
};

private string_t ClassToStr PARAM1(int, class) 
{
  classname_tp p;

  for (p = classname; p->name; p++)
    if (class == p->class) break;

  return (p->name) ? p->name : Padi_EMPTY_STR;
}

private int StrToClass PARAM1(string_t, s) 
{
  classname_tp p;

  for (p = classname; p->name; p++)
    if (!strcmp(s, p->name)) break;

  return p->class;
}

typedef struct modename_t {
  int mode;
  string_t name;
} modename_t, *modename_tp;

private modename_t modename[] = {
  HRMODE, "READ",	 
  HUMODE, "UPDATE",	 
  HSMODE, "SHARED",	 
      -1, "-1",
       0, NULL
};

private int StrToMode PARAM1(string_t, s) 
{
  modename_tp p;

  for (p = modename; p->name; p++)
    if (!strcmp(s, p->name)) break;

  return p->mode;
}

#ifdef DEBUG
private string_t ModeToStr PARAM1(int, m)
{
  modename_tp p;

  for (p = modename; p->mode; p++)
    if (m == p->mode) break;

  return p->name;
}
#endif

typedef struct dblist_t {
  char dbpath[256];
  int mode;
  int dbkey;
} dblist_t, *dblist_tp;

private int dblist_initialized = 0;
private unsigned long dblist_size;
private dblist_tp dblist;

private int InitializeDblist PARAMVOID
{
/*
   This function initializes the dblist structure to contain as many elements
   equal to 75% of the soft amount of the file descriptors resource. See the
   limit (csh) or ulimit (not available under SunOS Bourne & Korn shell) commands.

   The 75% amount was chosen arbitrarily. However, it must be noted that the
   size of this structure must be a bit smaller than the soft limit to allow for
   other files/databases to be opened - log files, FAME internal files,
   databases opened via FAME SERVER, etc.
*/

  unsigned long soft_limit;
  int i, ret;
  struct rlimit *rlp;
  dblist_tp p;

  if (dblist_initialized != 0) return 0;

  if (!(rlp = (struct rlimit *) malloc(sizeof(struct rlimit))))
    return Padi_OUT_OF_MEMORY;

  if ((ret=getrlimit(RLIMIT_NOFILE,rlp)) != 0)
    return Padi_FAIL;
  else {
    soft_limit = rlp->rlim_cur;
    dblist_size = (unsigned long) (soft_limit * 3 / 4);

#ifdef DEBUG
    printf("soft_limit is %d;   dblist_size is %d;   sizeof(dblist_t) is %d\n",
            soft_limit,dblist_size,sizeof(dblist_t));
#endif

    /* Allocate dblist structure to dblist_size elements */
    if (!(dblist = (dblist_tp) malloc(sizeof(dblist_t) * dblist_size)))
      return Padi_OUT_OF_MEMORY;

    /* Initialize dblist structure */
    for(i=0, p=dblist; i<dblist_size; i++, p++) {
      strcpy(p->dbpath,"");
      p->mode=-1;
      p->dbkey=-1;
    }    /* for */
    dblist_initialized = 1;
  }    /* if */
  return 0;
}  

private PadiStatus_t OpenDb PARAM3(char *, path, int, mode, int *, dbkey)
{
  dblist_tp p;
  PadiStatus_t status = Padi_SUCCEED;
  int ret;

  /* Initialize dblist structure */
  if (!dblist_initialized)
    if (ret = InitializeDblist()) return ret;

  *dbkey = -1;

#ifdef DEBUG
  printf("OpenDb args: %s  %s\n",path,ModeToStr(mode));
  for (p = dblist; *(p->dbpath); p++) 
    printf("%s %s %d\n",p->dbpath,ModeToStr(p->mode),p->dbkey);
  printf("Number of elements in dblist array is %d\n",dblist_size);
#endif

  /* make sure this is in our list */
  for (p = dblist; *(p->dbpath); p++) 
    if (!strncmp(path, p->dbpath, 255)) break;

  if (!(*(p->dbpath))) {
    if (p - dblist == dblist_size - 1) {
/*    if (p - dblist == NELEMENT(dblist) - 1) { */

#ifdef DEBUG
      printf("Too many databases opened. p - dblist is %d\n",p - dblist);
#endif

      return Padi_TOO_MANY_OPEN;

    }
    strncpy(p->dbpath, path, 255); 
  }
  
  if ((*dbkey = p->dbkey) != -1) {

    /* database is open -- check access mode */
    switch (p->mode) {
    case HRMODE:
      if (mode == HRMODE)
	return Padi_SUCCEED;
	break;

    case HSMODE:
      if (mode == HSMODE || mode == HRMODE)
	return Padi_SUCCEED;
	break;

    case HUMODE:
	return Padi_SUCCEED;
	break;

    }

#ifdef DEBUG
    printf("Close database key %d\n",*dbkey);
#endif

    /* close database so that it can be re-opened with new access mode */
    cfmcldb(&status, *dbkey);
    p->dbkey = -1;
    p->mode = -1;
  }

  if (*path == '*') {
    /* open the work database */ 

#ifdef DEBUG 
    printf("Opened WORK database\n");
#endif
 
    cfmopwk(&status, dbkey);
  } else {
    /* open the fame database */

#ifdef DEBUG 
    printf("Try to open %s in %s mode\n",path,ModeToStr(mode)); 
#endif
 
    cfmopdb(&status, dbkey, path, mode);
    if (status == HRNEXI && mode == HUMODE) /* create database if it doesn't exist */

#ifdef DEBUG 
    { printf("Try to create %s\n",path);
#endif
 
      cfmopdb(&status, dbkey, path, HCMODE);

#ifdef DEBUG 
    }
#endif
 
  }

  if (status != HSUCC)
    *dbkey = -1;
 
#ifdef DEBUG 
    printf("%s now set to %s mode and key %d\n\n",p->dbpath,ModeToStr(mode),*dbkey);
#endif
 
  p->dbkey = *dbkey;
  p->mode = mode;
  return status;
}
 
private PadiStatus_t CloseDb PARAM1(int, dbkey)
{
  dblist_tp p, q;
  PadiStatus_t status = Padi_SUCCEED;

  /* make sure this is in our list */
  for (p = dblist; *(p->dbpath); p++) 
    if (p->dbkey == dbkey) break;

  if (!(*(p->dbpath))) 
    return HBKEY;

  /* close the fame database */
  cfmcldb(&status, dbkey);

  for (q = p + 1; *(p->dbpath); p = q++) {
    strcpy(p->dbpath, q->dbpath);
    p->dbkey = q->dbkey;
  }
 
  return status;
}
 
private int obj_dbkey = -1;
private PadiStatus_t CloseAllDb PARAMVOID
{
  dblist_tp p;
  PadiStatus_t status = Padi_SUCCEED;
  char cmd_buf[1024];

  /* loop through database structure */
  for (p = dblist; *(p->dbpath); p++) {

    if (p->dbkey == -1) continue;
  
#ifdef DEBUG
  printf("CloseAllDb: %s;  dbkey=%d;  mode=%d\n",p->dbpath, p->dbkey, p->mode);
#endif

    /* close database */
    cfmcldb(&status, p->dbkey);
#ifdef DEBUG
  printf("cfmcldb: %d\n",status);
#endif
    p->dbkey = -1;
    p->mode = -1;

  }

  obj_dbkey = -1;


 /* besides closing all databases that open by cfmopdb and cfmopwk */
 /* we need to close all databases with channel names that are opened by FAME 4GL for formula series */ 

#ifdef DEBUG
  printf("CloseAllDb: Closing channel names\n");
  sprintf(cmd_buf, "display @open.db");
  cfmfame(&status, cmd_buf);
#endif

  sprintf(cmd_buf, "close all");
  cfmfame(&status, cmd_buf);

  return status;
}
 
private PadiStatus_t ReleaseDb PARAM2(char *, path, int, mode)
{
  dblist_tp p, q;
  PadiStatus_t status = Padi_SUCCEED;

  /* see if this is in our list */
  for (p = dblist; *(p->dbpath); p++) 
    if (!strncmp(path, p->dbpath, 255)) break;

  if (!(*(p->dbpath)) || p->dbkey == -1)
    return status; /* not open */
 

  /* database is open -- check access mode */
  if (p->mode == HRMODE && mode != HUMODE)
    return Padi_SUCCEED; 

  /* close database so that it can be opened with access mode */
  cfmcldb(&status, p->dbkey);
  p->dbkey = -1;
  p->mode = -1;

  return status;
}

private char *obj_dbname = NULL;
private PadiStatus_t OpenObjectDb PARAM1(int *, dbkey) 
{
  PadiStatus_t status;

  if ((*dbkey = obj_dbkey) != -1)
    return Padi_SUCCEED;

  if (!obj_dbname) 
    return HNINIT;

  if (!(status = OpenDb(obj_dbname, HRMODE, dbkey)))
    obj_dbkey = *dbkey;

  return status;
}
 
private PadiStatus_t IsoToDate PARAM3(int, freq, char *, iso_str, int *, date)
{
  char c;
  int iso_len;
  int status, year;
  int month = 1;
  int day = 1;

  if ( (iso_len = strlen(iso_str)) < 4)
     return (Padi_UNSUPPORTED_DATE);

  c = iso_str[4];
  iso_str[4] = '\0';
  year = atoi(iso_str);
  iso_str[4] = c;
  
  if(digitCount(year) != 4)
     return (Padi_UNSUPPORTED_DATE);
	 
  if (iso_len >= 6) {
    c = iso_str[6];
    iso_str[6] = '\0';
    month = atoi(iso_str + 4);
    iso_str[6] = c;
  }

  
  if (iso_len >= 8) {
    c = iso_str[8];
    iso_str[8] = '\0';
    day = atoi(iso_str + 6);
    iso_str[8] = c;
  }

  cfmddat(&status, freq, date, year, month, day);
  if (status && *date)
    status = Padi_SUCCEED;

  return status;
}

private PadiStatus_t DateToIso PARAM3(int, freq, int, date, char *, iso_str)
{
  int status, year, month, day;

  cfmdatd(&status, freq, date, &year, &month, &day);
  if (status != HSUCC)
    return status;

  sprintf(iso_str, "%4.4d%2.2d%2.2d", year, month, day);

  if( digitCount(year) != 4)
     return (Padi_UNSUPPORTED_DATE);

  return Padi_SUCCEED;
}

export void datetest_1 PARAM1(PadiString_t, datetest)
{
  char iso_str[Padi_DATE_TIME_LEN];
  freqname_tp p;
  int date;
  PadiStatus_t status;

  for (p = freqname; p->name; p++) {
    if (!(status = IsoToDate(p->freq, datetest, &date)))
      status = DateToIso(p->freq, date, iso_str);
    printf("%s at %s:\tdate=%d(%s),\tstatus=%d\n", 
            datetest, p->name, date, 
            (status) ? Padi_EMPTY_STR : iso_str, status);
    
  }
}


typedef struct property_t {
  char desc[Padi_MAX_STR_LEN + 1];
  char class[Padi_MAX_STR_LEN + 1];
  char type[Padi_MAX_STR_LEN + 1];
  char access[Padi_MAX_STR_LEN + 1];
  char wk_frequency[Padi_MAX_STR_LEN + 1];
  char frequency[Padi_MAX_STR_LEN + 1];
  char start[Padi_MAX_STR_LEN + 1];
  char end[Padi_MAX_STR_LEN + 1];
  char format[Padi_MAX_STR_LEN + 1];
  char dbdir[Padi_MAX_STR_LEN + 1];
  char dbname[Padi_MAX_STR_LEN + 1];
  char sysdbpath[Padi_MAX_STR_LEN + 1];
  char fame_name[Padi_MAX_STR_LEN + 1];
  char fame_class[Padi_MAX_STR_LEN + 1];
  char fame_ldpath[Padi_MAX_STR_LEN + 1];
} property_t, *property_tp;

private void RightTrim PARAM1(char *, s)
{
  char *t;

  for (t = s + strlen(s) - 1; t >= s && *t == ' '; t--)
    ;

  *(++t) = '\0';

}

private PadiStatus_t GetProperty PARAM3(PadiString_t, object_name, property_tp, result, PadiString_t, dbname)
{
  PadiStatus_t status;
  int objkey;
  char dbpath[Padi_MAX_STR_LEN + 1];
  char fame_name[Padi_MAX_STR_LEN + 1]; /* fame writes to this string! */
  numeric_t index;
  int type, freq, class, syear, sprd, eyear, eprd, range[8], numobs;
  int ismiss, length;
  int date;
  char new_access[Padi_MAX_STR_LEN + 1];

  *new_access = '\0';

  bzero((char *)result, sizeof(*result));

  if (!(object_name && *object_name)) {
    status = Padi_MISSING_OBJECT;
    return status;

  }

  /* open object database */
  if ((status = OpenObjectDb(&objkey)))
    return status;

  /* get the sysdb path */
  strcpy(fame_name, "obj_sysdbpath");
  cfmrstr(&status, objkey, fame_name, NULL, result->sysdbpath, &ismiss, &length);
  if (status)
    return status;

  /* get the database directory */
  strcpy(fame_name, "obj_dbdir");
  cfmrstr(&status, objkey, fame_name, NULL, result->dbdir, &ismiss, &length);
  if (status)
    return status;

  /* get the object case index */
  strncpy(fame_name, object_name, Padi_MAX_NAME_LEN);
  strcat(fame_name, ".obj");
  cfmrrng(&status, objkey, fame_name, NULL, (float *)(&index), HNTMIS, NULL);

#ifdef DEBUG
	printf("\tobj_sysdbpath = %s\n",result->sysdbpath);
	printf("\tobj_dbdir = %s\n",result->dbdir);
#endif

  if (status == HNOOBJ) {

    /* 
    ** check linked databases 
    */

#ifdef DEBUG
	printf("\t%s does not exist\n",fame_name);
#endif

    cfmosiz(&status, objkey, strcpy(fame_name, "obj_link.dbname"), 
            &class, &type, &freq, &syear, &sprd, &eyear, &eprd);

    if (status == HNOOBJ || eprd < 1) {
      status = HNOOBJ; /* empty list */

#ifdef DEBUG
	printf("\t%s does not exist or end-period is less than 1 (=%d)\n",fame_name, eprd);
#endif

    } else {
      int dbkey, sdb, edb, basis, observ;
      int cyear, cmonth, cday, myear, mmonth, mday;
      char junk[2];

#ifdef DEBUG
	printf("\t%s exists : end-period is %d\n",fame_name, eprd);
#endif

      for (sdb = sprd, edb = eprd; sdb <= edb; sdb++) {

        /* set the case range to the current database */
        sprd = sdb;
        numobs = 1;
        syear = eyear = eprd = -1;
        cfmsrng(&status, HCASEX, &syear, &sprd, &eyear, &eprd, range, &numobs);
        if (status)
          return status;
 
        /* get the database name */
        cfmrstr(&status, objkey, strcpy(fame_name, "obj_link.dbname"), 
                range, result->dbname, &ismiss, &length);
        if (status)
          return status;

        /* get the database access */
        cfmrstr(&status, objkey, 
                strcpy(fame_name, "obj_link.access"), 
                range, result->access, &ismiss, &length);
        if (status)
          return status;

        if (dbname && *dbname && !(*new_access)) 
          if (!strcmp(dbname, result->dbname) && 
               strcmp("READ", result->access))
            strcpy(new_access, result->access);

        /* open the database */
        strcpy(dbpath, result->dbdir);
        strcat(dbpath, "/");
        strcat(dbpath, result->dbname);

#ifdef DEBUG
	printf("\tOpening %s in %s mode\n",dbpath,result->access);
#endif

        if ((status = OpenDb(dbpath, StrToMode(result->access), &dbkey)))
          return status;

        /* look for the object */
        memset(result->desc, ' ', Padi_MAX_STR_LEN);
        result->desc[Padi_MAX_STR_LEN] = '\0';
        junk[0] = ' ';
        junk[1] = '\0';
        cfmwhat(&status, dbkey, strcpy(fame_name, object_name), 
                &class, &type, &freq,
                &basis, &observ, &syear, &sprd, &eyear, &eprd,
                &cyear, &cmonth, &cday, &myear, &mmonth, &mday, 
                result->desc, junk);
        if (status != HSUCC) {
#ifdef DEBUG
	printf("\t%s does not exist in %s\n",fame_name, dbpath);
#endif
          continue; /* not found */
	}

#ifdef DEBUG
	printf("\t%s found in %s\n",fame_name, dbpath);
#endif

        /* object found */
        RightTrim(result->desc);
        strcpy(result->fame_name, fame_name);
        strcpy(result->fame_class, ClassToStr(class));
        if (class == HFRMLA) {
          strcpy(result->access, "READ");
        } else {
          strcpy(result->class, result->fame_class);
          strcpy(result->type, TypeToStr(type));
          strcpy(result->frequency, FreqToStr(freq));
        }

#ifdef DEBUG
	printf("\tAttributes of %s are:\n"
	       "\tDbname = %s\n"
	       "\tDesc = %s\n"
	       "\tClass = %s\n"
	       "\tType = %s\n"
	       "\tFreq = %s\n"
	       "\tAccess = %s\n",
		result->fame_name,result->dbname,result->desc,result->fame_class,result->type,
		result->frequency,result->access );
#endif

        break;

      }
      
      if (sdb > edb) {

#ifdef DEBUG
	printf("\t%s does not exist in any database\n",fame_name);
#endif

        status = HNOOBJ; /* not found */
      }
    }

  } else if (status == HSUCC && (int)index > 0) {

    /*
    ** get info from object table
    */

#ifdef DEBUG
	printf("\t1)%s exists with index = %d. Is an exception.\n",fame_name,(int) index);
#endif

    /* set the case range to this index */
    sprd = (int)index;
    numobs = 1;
    syear = eyear = eprd = -1;
    cfmsrng(&status, HCASEX, &syear, &sprd, &eyear, &eprd, range, &numobs);
    if (status)
      return status;
 
    /* get the description */
    cfmrstr(&status, objkey, strcpy(fame_name, "obj.desc"), range, result->desc, &ismiss, &length);
    if (status)
      return status;

    /* get the access */
    cfmrstr(&status, objkey, strcpy(fame_name, "obj.access"), range, result->access, &ismiss, &length);
    if (status)
      return status;

    /* get the class */
    cfmrstr(&status, objkey, strcpy(fame_name, "obj.class"), range, result->class, &ismiss, &length);
    if (status)
      return status;

    /* get the freq */
    cfmrstr(&status, objkey, strcpy(fame_name, "obj.freq"), range, result->frequency, &ismiss, &length);
    if (status)
      return status;

    /* get the type */
    cfmrstr(&status, objkey, strcpy(fame_name, "obj.type"), range, result->type, &ismiss, &length);
    if (status)
      return status;

    /* get the fame object name */
    cfmrstr(&status, objkey, strcpy(fame_name, "obj.fame_name"), range, result->fame_name, &ismiss, &length);
    if (status)
      return status;

    /* get the fame object class */
    cfmrstr(&status, objkey, strcpy(fame_name, "obj.fame_class"), range, result->fame_class, &ismiss, &length);
    if (status)
      return status;

    /* get the fame function ldpath */
    cfmrstr(&status, objkey, strcpy(fame_name, "obj.fame_ldpath"), range, result->fame_ldpath, &ismiss, &length);
    if (status)
      return status;

    /* get the database name */
    cfmrstr(&status, objkey, strcpy(fame_name, "obj.dbname"), range, result->dbname, &ismiss, &length);
    if (status)
      return status;

#ifdef DEBUG
	printf("\tException Object Attributes (obj.*) of %s are:\n"
	       "\tDesc = %s\n"
	       "\tClass = %s\n"
	       "\tType = %s\n"
	       "\tFreq = %s\n"
	       "\tAccess = %s\n",
	       "\tFame_Class = %s\n",
	       "\tFame_ldpath = %s\n",
	       "\tdbname = %s\n",
		result->fame_name,result->desc,result->class,result->type,
		result->frequency,result->access,result->fame_class,
		result->fame_ldpath,result->dbname);
#endif

    /* build the full path */
    strcpy(dbpath, result->dbdir);
    strcat(dbpath, "/");
    strcat(dbpath, result->dbname);

    if (!strcmp(result->fame_class, "SERIES")) {
      int dbkey;

#ifdef DEBUG
	printf("\tobj.fame_class is SERIES. Open %s in %s mode.\n",dbpath,result->access);
#endif

      /* open the object's database */
      if ((status = OpenDb(dbpath, StrToMode(result->access), &dbkey)))
        return status;

      /* get the start and end */
      cfmosiz(&status, dbkey, strncpy(fame_name, object_name, Padi_MAX_NAME_LEN), &class, &type,
              &freq, &syear, &sprd, &eyear, &eprd);

#ifdef DEBUG
	printf("\t%s: syear=%d ; sprd=%d ; eyear=%d ; eprd=%d ; freq=%d ; class=%d ; type=%d\n",
		fame_name,syear,sprd,eyear,eprd,freq,class,type);
#endif

    } else {

#ifdef DEBUG
	printf("\tobj.fame_class is not SERIES\n");
#endif

      eprd = -1;
      freq = HUNDFX;
    }

  }

  else if (status == HSUCC && (int)index == 0) {

    /*
    ** Exclusions
    */

#ifdef DEBUG
	printf("\t2)%s exists with index = %d. Is an exclusion.\n",fame_name,(int) index);
#endif

      status = HNOOBJ;

  }

  if (status != HSUCC) {
    if (status == HNOOBJ) {
      status = Padi_UNKNOWN_OBJECT;
      strcpy(result->dbname, (dbname) ? dbname : "");
      strcpy(result->access, new_access);
    } 
    return status;

  }

  if (eprd > 0 && freq != HUNDFX) {
    /* convert the start and end periods to Iso dates */
    if (freq == HCASEX) {
      sprintf(result->start, "%d", sprd);
      sprintf(result->end, "%d", eprd);
    } else {
      cfmpdat(&status, freq, &date, syear, sprd);
      if (status)
        return status;

      if ((status = DateToIso(freq, date, result->start)))
        return status;

      cfmpdat(&status, freq, &date, eyear, eprd);
      if (status)
        return status;

      if ((status = DateToIso(freq, date, result->end)))
        return status;

    }

#ifdef DEBUG
	printf("\tConverted start and end periods to Iso dates; start=%s; end=%s\n",
		result->start,result->end);
#endif

  }

  strcpy(result->format, "SimpleSeries");
 
#ifdef DEBUG
	printf("\tReturn from GetProperty. Status is %d\n", status);
#endif

  return status;
}


export PadiInfoResult_tp getinfo_1 PARAM2(PadiInfoArg_tp, argp, CLIENT *, rqstp)
{
  PadiString_t user = argp->user;
  PadiString_t password = argp->password;
  PadiString_t object_name = argp->object_name;
  property_t property;
  PadiInfoResult_tp result;
  PadiInfo_tp info;
  size_t struct_len, buf_len;
  PadiString_t s, *strp;
  

  struct_len = (size_t)ALIGN(align_t, sizeof(PadiInfoResult_t));
  buf_len = (size_t)ALIGN(align_t, Padi_MAX_STR_LEN + 1);
  if (!(result = (PadiInfoResult_tp)malloc(struct_len + 
                                             Padi_NUM_INFO_STR*buf_len)))
    return NULL;

  bzero((char *)result, struct_len + Padi_NUM_INFO_STR*buf_len);
  s = (char *)result + struct_len;
  info = &(result->info);
  for (strp = (PadiString_t *)info; 
       strp < ((PadiString_t *)info + Padi_NUM_INFO_STR); 
       strp++, s += buf_len)
    *strp = s;

  if (!(object_name && *object_name)) {
    result->status = Padi_MISSING_OBJECT;
    return result;

  }

#ifdef DEBUG
	printf("\tgetinfo_1: Calling GetProperty.\n");
#endif

  if ((result->status = GetProperty(object_name, &property, 
		(PadiString_t) NULL))) {
    CloseAllDb();
    return result;
  }

  strcpy(info->name, object_name);
  strcpy(info->desc, property.desc);
  strcpy(info->class, property.class);
  strcpy(info->type, property.type);
  if (!(!strcmp(property.fame_class, "SERIES") || 
        !strcmp(property.fame_class, "SCALAR")))
    strcpy(info->access, "READ");
  else 
    strcpy(info->access, property.access);
  strcpy(info->frequency, property.frequency);
  strcpy(info->start, property.start);
  strcpy(info->end, property.end);
  strcpy(info->format, property.format);

  CloseAllDb();
  return result;
}

export PadiSeriesResult_tp getseries_1 PARAM2(PadiRangeArg_tp, arg, CLIENT *, rqstp)
{
  PadiRange_tp in_range = &(arg->range);
  PadiString_t user = arg->user;
  PadiString_t password = arg->password;
  int dbkey;
  char fame_name[Padi_MAX_STR_LEN + 1]; /* fame writes to this string! */
  char dbpath[Padi_MAX_STR_LEN + 1];
  int type, freq, class, syear, sprd, eyear, eprd, range[8], numobs;
  int in_syear, in_sprd, in_eyear, in_eprd;
  int ismiss, length;
  property_t property;
  PadiSeriesResult_tp result;
  PadiInfo_tp info;
  PadiRange_tp data_range;
  size_t struct_len, data_len, buf_len;
  char *s;
  PadiString_t *strp;
  PadiPrecision_t *vector;

  struct_len = (size_t)ALIGN(align_t, sizeof(PadiSeriesResult_t));
  if (!(result = (PadiSeriesResult_tp)malloc(struct_len))) {
    return NULL;
  }

  bzero((char *)result, struct_len);
  data_range = &(result->series.range);
  data_range->start = Padi_EMPTY_STR;
  data_range->end = Padi_EMPTY_STR;
  data_range->format = Padi_EMPTY_STR;
  info = &(result->series.info);
  for (strp = (PadiString_t *)info; 
       strp < ((PadiString_t *)info + Padi_NUM_INFO_STR); 
       strp++)
    *strp = Padi_EMPTY_STR;

  if (!(arg->object_name && *(arg->object_name))) {
    result->status = Padi_MISSING_OBJECT;
    return result;

  }
 
#ifdef DEBUG
	printf("\tgetseries_1: Calling GetProperty.\n");
#endif

  if ((result->status = GetProperty(arg->object_name, &property,(PadiString_t) NULL))) {
    CloseAllDb();
    return result;
  }

  if (in_range->format && 
      *(in_range->format) && 
      !strstr(property.format, in_range->format)) {
    result->status = Padi_UNSUPPORTED_FORMAT;
    CloseAllDb();
    return result;

  }

#ifdef DEBUG
	printf("\tSupported format\n");
#endif

  /* build the full path */
  strcpy(dbpath, property.dbdir);
  strcat(dbpath, "/");
  strcat(dbpath, property.dbname);

  if (!(strcmp(property.fame_class, "SERIES") && strcmp(property.fame_class, "SCALAR"))) {

    /* 
    ** SERIES or SCALAR
    */

    /* open the object's database */
    if ((result->status = OpenDb(dbpath, StrToMode(property.access), &dbkey))) {
      CloseAllDb();
      return result;
    }

#ifdef DEBUG
	printf("\tObject %s has obj.fame_class set to SERIES or SCALAR =%s\n",
		arg->object_name,property.fame_class);
	printf("\tOpened %s in %s mode.\n",dbpath,property.access);
#endif

  } else {
    char cmd_buf[1024];

    /*
    ** FORMULA or FUNCTION
    */

    /* open a fame session work database */
    if ((result->status = OpenDb("*", HUMODE, &dbkey))) {
      CloseAllDb();
      return result;
    }

    sprintf(cmd_buf, "overwrite on; ignore on;");
    cfmfame(&(result->status), cmd_buf);
    if (result->status) {

      CloseAllDb();
      return result;
    }

    strcpy(property.wk_frequency, arg->working_freq);

    if ( arg->working_freq !=NULL ) {   
       freq = StrToFreq(property.wk_frequency);
    } else {
       freq = StrToFreq(property.frequency);
    } 

    if (freq != HUNDFX) {
       if ( arg->working_freq !=NULL ) { 
          sprintf(cmd_buf, "freq %s", property.wk_frequency);
       } 
       else {
          sprintf(cmd_buf, "freq %s", property.frequency);
       }
       cfmfame(&(result->status), cmd_buf);
       if (result->status) {

        CloseAllDb();
        return result;
       }
   
#ifdef DEBUG
        printf("\tgetseries_1: freq = %s \n", property.frequency);
        printf("\tgetseries_1: FAMEfreq = %d \n", freq);
        printf("\tgetseries_1: work_freq = %s \n", property.wk_frequency);
	printf("\tObject %s has obj.fame_class set to FORMULA or FUNCTION = %s\n",
		arg->object_name,property.fame_class);
	printf("\tOpened WORK and executed: %s\n",cmd_buf);
#endif

    }
#ifdef DEBUG
        printf("\txgetseries_1: freq = %s \n", property.frequency);
        printf("\txgetseries_1: FAMEfreq = %d \n", freq);
        printf("\txgetseries_1: work_freq = %s \n", property.wk_frequency);
#endif

    /* open primary database */
    if ((result->status = ReleaseDb(dbpath, StrToMode(property.access)))) {
      CloseAllDb();
      return result;
    }

 
    sprintf(cmd_buf,"try; open <access %s> file(\"%s\") as %s; end try",    
                      property.access, dbpath, property.dbname);

  
    cfmfame(&(result->status), cmd_buf); 
    if (result->status) {  
      CloseAllDb(); 
     return result; 
    }  

#ifdef DEBUG
	printf("\tOpened %s in %s mode and executed: %s\n",dbpath,property.access,cmd_buf);
#endif

    if (!strcmp(property.fame_class, "FORMULA")) {
      int syskey;

      /* open formula databases for the primary database */

      if ((result->status = OpenDb(property.sysdbpath, HRMODE, &syskey))) {
        CloseAllDb();
        return result;
      }

      strcpy(fame_name, "g_");
      strcat(fame_name, property.dbname);
      strcat(fame_name, "_formula_db_list_str");
      cfmosiz(&(result->status), syskey, fame_name,
              &class, &type, &freq, &syear, &sprd, &eyear, &eprd);

      if (result->status == HNOOBJ || eprd < 1) {
        result->status = HNOOBJ; /* empty list */

#ifdef DEBUG
	printf("\tobj.fame_class is FORMULA: Read %s: Is NON_EXISTANT or EMPTY\n",fame_name);
#endif

      } else {
        int sdb, edb;
        char sdbname[Padi_MAX_STR_LEN + 1];

#ifdef DEBUG
	printf("\tobj.fame_class is FORMULA: Read non-empty %s successfully\n",fame_name);
#endif

        for (sdb = sprd, edb = eprd; sdb <= edb; sdb++) {

          /* set the case range to the current database */
          sprd = sdb;
          numobs = 1;
          syear = eyear = eprd = -1;
          cfmsrng(&(result->status), HCASEX, &syear, &sprd, &eyear,
			 &eprd, range, &numobs);
          if (result->status) {
            CloseAllDb();
            return result;
          }
 
          /* get the database name */
          cfmrstr(&(result->status), syskey, fame_name, 
                  range, sdbname, &ismiss, &length);
          if (result->status) {
            CloseAllDb();
            return result;
          }

          /* build the full path */
          strcpy(dbpath, property.dbdir);
          strcat(dbpath, "/");
          strcat(dbpath, sdbname);

          if ((result->status = ReleaseDb(dbpath, HRMODE))) {
            CloseAllDb();
            return result;
          }

          sprintf(cmd_buf, "try; open <access READ> file(\"%s\") as %s; end try", 
                           dbpath, sdbname);

          cfmfame(&(result->status), cmd_buf);

          if (result->status) { 
           CloseAllDb();  
            return result;  
          }   

#ifdef DEBUG
	printf("\tExecuted: %s\n",cmd_buf);
#endif

        }
      }

      /* copy the formula into a work series (scalar) */
     sprintf(cmd_buf, "work'%s = %s'%s", property.fame_name, property.dbname,
		property.fame_name);

      cfmfame(&(result->status), cmd_buf);
      if (result->status) {
        CloseAllDb();
        return result;
      }

      /* RESET NUMOBS to -1: bugfix hp */
      numobs = -1;

#ifdef DEBUG
      printf("\tFormula: Executed: %s\n",cmd_buf);
      printf("\tFormula: freq: %d\n",freq);
      printf("\tFormula: frequency: %s  wk_frequency %s \n",property.frequency, property.wk_frequency );

#endif


    } else {  /* FUNCTION */
   
      /* make sure the function is loaded */ 
      sprintf(cmd_buf, 
              "if missing(lookup(\"%s\", sl(make(namelist,@functions)))); load file(\"%s\"); end", 
              property.fame_name, property.fame_ldpath);
      cfmfame(&(result->status), cmd_buf);
      if (result->status) {
        CloseAllDb();
        return result;
      }

#ifdef DEBUG
	printf("\tFunction: Executed: %s\n",cmd_buf);
#endif
 
      /* copy the function into a work series */
      sprintf(cmd_buf, "work'%s = %s()", property.fame_name, property.fame_name);
      cfmfame(&(result->status), cmd_buf);
      if (result->status) {
        CloseAllDb();
        return result;
      }

#ifdef DEBUG
	printf("\tExecuted: %s\n",cmd_buf);
#endif

    }

    sprintf(cmd_buf, "close all");
    cfmfame(&(result->status), cmd_buf);
    if (result->status) {
      CloseAllDb();
      return result;
    }

#ifdef DEBUG
	printf("\tExecuted: %s\n",cmd_buf);
#endif

  }

  cfmosiz(&(result->status), dbkey, strcpy(fame_name, property.fame_name), 
          &class, &type, &freq, &syear, &sprd, &eyear, &eprd);

#ifdef DEBUG
    printf("\n\tcfmosiz() call: status is %d\n"
	"\tObject: %s\n"
	"\tClass: %s\n"
	"\tType: %s\n"
	"\tFreq: %s\n"
	"\tSyear: %d\n"
	"\tSprd: %d\n"
	"\tEyear: %d\n"
	"\tEprd: %d\n",
	result->status,fame_name,ClassToStr(class),
	TypeToStr(type),FreqToStr(freq),syear,sprd,eyear,eprd);
    printf("\n\t%s has FAME Type %d\n",fame_name,type);
#endif

  /* only support Numeric and Precision FAME Types */
  switch(type) {
	case HNAMEL: case HBOOLN: case HSTRNG: case HDATE:
		result->status = HBOBJT;
                CloseAllDb();
		return result;
		break;
	case HNUMRC: case HPRECN:
		break;
	default:
		result->status = HBOBJT;
                CloseAllDb();
		return result;
		break;
  }

  if (result->status == HNOOBJ) {
    result->status = Padi_SUCCEED; /* non-existent (i.e. empty) series */
    eprd = 0;

#ifdef DEBUG
	printf("\tobj.fame_name, %s is NON-EXISTANT - Set eprd=0\n",fame_name);
#endif

  }

  if (result->status) {
    CloseAllDb();
    return result;
  }

  /* convert the input range */  

  if (freq == HUNDFX) {
    in_syear = in_eyear = -1;
    in_sprd = in_eprd = -1;
  } else if (freq == HCASEX) {
    in_syear = in_eyear = -1;
    in_sprd = (*(in_range->start)) ? atoi(in_range->start) : -1;
    in_eprd = (*(in_range->end)) ? atoi(in_range->end) : -1;
  } else {
    int date; 

    if (*(in_range->start)) {
      if (!(result->status = IsoToDate(freq, in_range->start, &date)))
        cfmdatp(&(result->status), freq, date, &in_syear, &in_sprd);
      if (result->status) {
        CloseAllDb();
        return result;
      }

    } else {
      in_syear = in_sprd = -1;
    }

    if (*(in_range->end)) {
      if (!(result->status = IsoToDate(freq, in_range->end, &date)))
        cfmdatp(&(result->status), freq, date, &in_eyear, &in_eprd);
      if (result->status) {
        CloseAllDb();
        return result;
      }

    } else {
      in_eyear = in_eprd = -1;
    }
  }

  if (in_range->nobs) {
    numobs = in_range->nobs;
    if (in_eprd != -1)
      in_syear = in_sprd = -1;
  } else {
    numobs = -1;
  }

  if (class == HSCALA) {
    numobs = 1;
  } else if (eprd < 1) {
    numobs = 0; /* empty series */
  } else {
    if (numobs != -1) {
      if (in_sprd == -1 && in_eprd == -1) {
        in_eyear = eyear;
        in_eprd = eprd;
      }
    } else {
      if (in_sprd == -1) {
        in_syear = syear;
        in_sprd = sprd;
      }
      if (in_eprd == -1) {
        in_eyear = eyear;
        in_eprd = eprd;
      }
    }

  
    /* set fame range */
    cfmsrng(&(result->status), freq, &in_syear, &in_sprd, &in_eyear, &in_eprd, range, &numobs);

    if (result->status) {
      CloseAllDb();
      return result;
    }
  
  }


  struct_len = (size_t)ALIGN(align_t, sizeof(PadiSeriesResult_t));
  data_len = (size_t)ALIGN(align_t, numobs*sizeof(PadiPrecision_t));
  buf_len = (size_t)ALIGN(align_t, Padi_MAX_STR_LEN + 1);

#ifdef DEBUG
  printf("sizeof(PadiSeriesResult_t) is %d\n",sizeof(PadiSeriesResult_t));
  printf("ALIGN calls: struct_len is %d;  data_len is %d;  buf_len is %d\n",
	struct_len,data_len,buf_len);
#endif

  if (!(result = (PadiSeriesResult_tp)realloc((char *)result,
                   struct_len + data_len + (3 + Padi_NUM_INFO_STR)*buf_len))) {
    CloseAllDb();
    return NULL;
  }

  s = (char *)result;
  bzero(s, struct_len + data_len + (3 + Padi_NUM_INFO_STR)*buf_len);
  data_range = &(result->series.range);
  data_range->nobs = numobs;
  s += struct_len;
  result->series.data.data_len = numobs;
  result->series.data.data_val = (PadiPrecision_t *) s;
  vector = (PadiPrecision_t *)s;
  s += data_len;
  data_range->start = s; s += buf_len;
  data_range->end = s; s += buf_len;
  data_range->format = s; s += buf_len;
  info = &(result->series.info);
  for (strp = (PadiString_t *)info; 
       strp < ((PadiString_t *)info + Padi_NUM_INFO_STR); 
       strp++, s += buf_len)
    *strp = s;

/* adds found database name to info structure so it can be passed back to Series*/
  strcpy(info->dbname, property.dbname);

  strcpy(info->name, arg->object_name);
  strcpy(info->desc, property.desc);
  strcpy(info->class, ClassToStr(class));
  strcpy(info->type, TypeToStr(type));
  if (!(!strcmp(property.fame_class, "SERIES") || 
        !strcmp(property.fame_class, "SCALAR")))
    strcpy(info->access, "READ");
  else 
    strcpy(info->access, property.access);
  strcpy(info->frequency, FreqToStr(freq));
  strcpy(info->start, property.start);
  strcpy(info->end, property.end);
  strcpy(info->format, property.format);

  if (numobs > 0) {
    /* read data */
    if(type == HPRECN) {
    cfmrrng(&(result->status), dbkey, 
            strcpy(fame_name, property.fame_name), 
            range, (double *)vector, 
            (in_range->do_missing) ? HTMIS : HNTMIS, 
            in_range->missing_translation);

#ifdef DEBUG
	printf("\tRead Precision data of %s\n",fame_name);
#endif

    }
    else {
      float *v;
      int i;
      if (!(v = (float *)malloc(numobs*sizeof(float)))) {
        CloseAllDb();
        return NULL;
      }
      cfmrrng(&(result->status), dbkey, 
            strcpy(fame_name, property.fame_name), 
            range, v, 
            (in_range->do_missing) ? HTMIS : HNTMIS, 
            in_range->missing_translation);
      for(i=0;i<numobs;i++) {
	vector[i]=(double)v[i];
      }
      free((void *)v);

#ifdef DEBUG
	printf("\tRead Numeric data of %s\n",fame_name);
#endif

    }

    if (result->status) {
      CloseAllDb();
      return result;
    }
  
  }

  /* store the range start and end */
  if (freq == HCASEX) {
    sprintf(info->start, "%d", sprd);
    sprintf(info->end, "%d", eprd);
    sprintf(data_range->start, "%d", in_sprd);
    sprintf(data_range->end, "%d", in_eprd);
  } else if (freq != HUNDFX) {
    int date;

    cfmpdat(&(result->status), freq, &date, syear, sprd);
    if (!(result->status))
      if ((result->status = DateToIso(freq, date, info->start))) {
        CloseAllDb();
        return result;
      }

    cfmpdat(&(result->status), freq, &date, eyear, eprd);
    if (!(result->status))
      if ((result->status = DateToIso(freq, date, info->end))) {
        CloseAllDb();
        return result;
      }

    cfmpdat(&(result->status), freq, &date, in_syear, in_sprd);
    if (!(result->status))
      if ((result->status = DateToIso(freq, date, data_range->start))) {
        CloseAllDb();
        return result;
      }

    cfmpdat(&(result->status), freq, &date, in_eyear, in_eprd);
    if (!(result->status))
      if ((result->status = DateToIso(freq, date, data_range->end))) {
        CloseAllDb();
        return result;
      }

  }

  /* store the range translation */
  data_range->do_missing = in_range->do_missing;
  memcpy((char *)(data_range->missing_translation), 
         (char *)(in_range->missing_translation),
         sizeof(in_range->missing_translation));

#ifdef DEBUG
    {
    int j;
    printf("\tgetseries_1: Data read for %s (%d observations)\n",
		fame_name,data_range->nobs);
    printf("\tFrequency is %s;  Start Date is %s;  End Date is %s\n",
		FreqToStr(freq),data_range->start,data_range->end);
    for(j=0; j<data_range->nobs; j++)
       printf("data_val[%d] = %lf\n",j,result->series.data.data_val[j]);
    }
#endif

        /* sprintf(cmd_buf, "freq ANNUAL"); */
         cfmfame(&(result->status), "freq ANNUAL");
         if (result->status) {
            CloseAllDb();
            return result;
         }

  CloseAllDb();


#ifdef DEBUG
	puts("AFTER CloseAllDb\n");
#endif

  return result;
}

export PadiResult_tp putseries_1 PARAM2(PadiSeries_tp, series, CLIENT *, rqstp)
{
  PadiInfo_tp info = &(series->info);
  PadiRange_tp in_range = &(series->range);
  int dbkey;
  char fame_name[Padi_MAX_STR_LEN + 1]; /* fame writes to this string! */
  property_t property;
  char dbpath[Padi_MAX_STR_LEN + 1];
  int type, freq, class, syear, sprd, eyear, eprd, range[8], numobs;
  int in_syear, in_sprd, in_eyear, in_eprd;
  int date;
  PadiResult_tp result;
  size_t struct_len;
  precision_t *pvector = (precision_t *)(series->data.data_val);
  numeric_t *nvector = (numeric_t *)(series->data.data_val);


  struct_len = sizeof(*result);
  if (!(result = (PadiResult_tp)malloc(struct_len)))
    return NULL;
  bzero((char *)result, struct_len);

  if (!(info->name && *(info->name))) {
    result->status = Padi_MISSING_OBJECT;
    return result;
  }

  if ((result->status = GetProperty(info->name, &property, (PadiString_t)NULL))) {
    CloseAllDb();
    return result;
  }

  if (!strcmp(property.access, "READ") ||
      (strcmp(property.fame_class, "SERIES") && 
       strcmp(property.fame_class, "SCALAR"))) {
    result->status = Padi_READ_ACCESS;
    CloseAllDb();
    return result;
  }

  if (!(info->type && *(info->type))) {
    result->status = Padi_MISSING_TYPE;
    CloseAllDb();
    return result;
  }

  if (strcmp(property.type, info->type)) {
    result->status = Padi_TYPE_MISMATCH;
    CloseAllDb();
    return result;
  }

  if (!(in_range->format && *(in_range->format))) {
    result->status = Padi_MISSING_FORMAT;
    CloseAllDb();
    return result;
  }

  if (!strstr(property.format, in_range->format)) {
    result->status = Padi_UNSUPPORTED_FORMAT;
    CloseAllDb();
    return result;
  }

  /* build the full path */
  strcpy(dbpath, property.dbdir);
  strcat(dbpath, "/");
  strcat(dbpath, property.dbname);

  if ((result->status = OpenDb(dbpath, StrToMode(property.access), &dbkey))) {
    CloseAllDb();
    return result;
  }

  cfmosiz(&(result->status), dbkey, strcpy(fame_name, property.fame_name), &class, &type,
          &freq, &syear, &sprd, &eyear, &eprd);

  if (result->status) {
    CloseAllDb();
    return result;
  }

  if (freq != HUNDFX) {
    if (freq == HCASEX) {
      in_syear = in_eyear = -1;
      in_sprd = (*(in_range->start)) ? atoi(in_range->start) : -1;
      in_eprd = (*(in_range->end)) ? atoi(in_range->end) : -1;
    } else {
      if (*(in_range->start)) {
        if (!(result->status = IsoToDate(freq, in_range->start, &date)))
          cfmdatp(&(result->status), freq, date, &in_syear, &in_sprd);
        if (result->status) {
          CloseAllDb();
          return result;
        }
  
      } else {
        in_syear = in_sprd = -1;
      }
  
      if (*(in_range->end)) {
        if (!(result->status = IsoToDate(freq, in_range->end, &date)))
          cfmdatp(&(result->status), freq, date, &in_eyear, &in_eprd);
        if (result->status) {
          CloseAllDb();
          return result;
        }
  
      } else {
        in_eyear = in_eprd = -1;
      }
    }
  
    /* set fame range */
    if (in_range->nobs) {
      numobs = in_range->nobs;
      if (in_eprd != -1) 
        in_syear = in_sprd = -1;
    } else {
      numobs = -1;
    }
  
    cfmsrng(&(result->status), freq, &in_syear, &in_sprd, &in_eyear, &in_eprd, range, &numobs);
    if (result->status) {
      CloseAllDb();
      return result;
    }

  } /* if (freq != HUNDFX) */

  /* write data */
  cfmwrng(&(result->status), dbkey, 
          strcpy(fame_name, property.fame_name), 
          range, (float *)pvector, 
          (in_range->do_missing) ? HTMIS : HNTMIS, 
          in_range->missing_translation);
  if (result->status) {
    CloseAllDb();
    return result;
  }

  CloseAllDb();
  return result;
}

export PadiResult_tp newseries_1 PARAM2(PadiNewSeries_tp, new, CLIENT *, rqstp)
{
  PadiString_t user = new->user;
  PadiString_t password = new->password;
  PadiSeries_tp series = &(new->series);
  PadiInfo_tp info = &(series->info);
  PadiRange_tp in_range = &(series->range);
  int dbkey, freq, class, type;
  char fame_name[Padi_MAX_STR_LEN + 1]; /* fame writes to this string! */
  property_t property;
  char dbpath[Padi_MAX_STR_LEN + 1];
  PadiResult_tp result;
  size_t struct_len;


  struct_len = sizeof(*result);
  if (!(result = (PadiResult_tp)malloc(struct_len)))
    return NULL;
  bzero((char *)result, struct_len);


  if ((result->status = GetProperty(info->name, &property, new->dbname)) 
        != Padi_UNKNOWN_OBJECT) {
    if (result->status == HSUCC)
      result->status = Padi_OBJECT_EXISTS;
    CloseAllDb();
    return result;
  }

  if (!(new->dbname && *(new->dbname))) {
    result->status = Padi_MISSING_DBNAME;
    CloseAllDb();
    return result;
  }

  if (!(*property.access)) {
    result->status = Padi_READ_ACCESS;
    CloseAllDb();
    return result;
  }

  if (!(class = StrToClass(info->class))) {
    result->status = Padi_UNSUPPORTED_CLASS;
    CloseAllDb();
    return result;
  }

  if (!(freq = StrToFreq(info->frequency))) {
    result->status = Padi_UNSUPPORTED_FREQ;
    CloseAllDb();
    return result;
  }

  if (!(type = StrToType(info->type))) {
    result->status = Padi_UNSUPPORTED_TYPE;
    CloseAllDb();
    return result;
  }

  /* build the full path */
  strcpy(dbpath, property.dbdir);
  strcat(dbpath, "/");
  strcat(dbpath, property.dbname);

  if ((result->status = OpenDb(dbpath, StrToMode(property.access), &dbkey))) {
    CloseAllDb();
    return result;
  }

  /* create the object */
  cfmalob(&(result->status), dbkey, strcpy(fame_name, info->name), 
          class, freq, type, HBSDAY, HOBEND,
          in_range->nobs, in_range->nobs*16, (float)2.0);
  if (result->status != HSUCC) {
    CloseAllDb();
    return result;
  }

  free((char *)result); 
  
  return putseries_1(series, rqstp);
}

export PadiResult_tp destroy_1 PARAM2(PadiDestroyArg_tp, argp, CLIENT *, rqstp)
{
  PadiString_t user = argp->user;
  PadiString_t password = argp->password;
  PadiString_t object_name = argp->object_name;
  int dbkey;
  char fame_name[Padi_MAX_STR_LEN + 1]; /* fame writes to this string! */
  property_t property;
  char dbpath[Padi_MAX_STR_LEN + 1];
  PadiResult_tp result;
  size_t struct_len;

  struct_len = sizeof(*result);
  if (!(result = (PadiResult_tp)malloc(struct_len)))
    return NULL;
  bzero((char *)result, struct_len);

  if ((result->status = GetProperty(object_name, &property,(PadiString_t) NULL))) {
    CloseAllDb();
    return result;
  }

  if (!strcmp(property.access, "READ")) {
    result->status = Padi_READ_ACCESS;
    CloseAllDb();
    return result;
  }

  /* build the full path */
  strcpy(dbpath, property.dbdir);
  strcat(dbpath, "/");
  strcat(dbpath, property.dbname);

  /* open object's database */
  if ((result->status = OpenDb(dbpath, StrToMode(property.access), &dbkey))) {
    CloseAllDb();
    return result;
  }

  /* destroy the object */
  cfmdlob(&(result->status), dbkey, strcpy(fame_name, property.fame_name)); 
  CloseAllDb();
  return result;
}

export PadiStatus_t initialize_1 PARAM1(PadiInitArg_tp, object)
{
  PadiStatus_t status = Padi_SUCCEED;
  int dbkey;

  /* initialize fame */
  cfmini(&status);
  if (status != HSUCC)
    return status;

  obj_dbname = object->object_dbname;

  /* open object database */
  status = OpenObjectDb(&dbkey);

  return status;
}

export PadiStatus_t terminate_1 PARAM1(PadiTermArg_tp, object)
{
  PadiStatus_t status = Padi_SUCCEED;
  dblist_tp p;

  /* close all open databases */
  for (p = dblist; *(p->dbpath); p++) {
    if (p->dbkey == -1) continue;

    cfmcldb(&status, p->dbkey);
    if (status != HSUCC)
      return status;
  }

  /* terminate fame */
  cfmfin(&status);
  return status;
}

export char * status_1 PARAM1(int, status)
{
  static char buf[64];
  /*
  ** Return message for a HLI status code
  */
  switch (status) {
  case HSUCC:  return "Success";                                        
  case HINITD: return "HLI has already been initialized";               
  case HNINIT: return "HLI has not been initialized";                   
  case HFIN:   return "HLI has already been finished";                  
  case HBFILE: return "Bad file name";                                  
  case HBMODE: return "Bad or unauthorized file access mode";           
  case HBKEY:  return "Bad data base key";                              
  case HBSRNG: return "Bad starting year or period";                    
  case HBERNG: return "Bad ending year or period";                      
  case HBNRNG: return "Bad number of observations";                     
  case HNOOBJ: return "Object does not exist";                          
  case HBRNG:  return "Bad range";                                      
  case HDUTAR: return "Target object already exists";                   
  case HBOBJT: return "Bad object type";                                
  case HBFREQ: return "Bad frequency";                                  
  case HTRUNC: return "Oldest data has been truncated";                 
  case HNPOST: return "Data base not posted or closed";                 
  case HFUSE:  return "File already in use";                            
  case HNFMDB: return "File not a FAME data base";                      
  case HRNEXI: return "Read or update and file does not exist";         
  case HCEXI:  return "Create and file exists";                         
  case HNRESW: return "Name is reserved or not a legal FAME name";      
  case HBCLAS: return "Bad object class";                               
  case HBOBSV: return "Bad OBSERVED attribute";                         
  case HBBASI: return "Bad BASIS attribute";                            
  case HOEXI:  return "Object already exists";                          
  case HBFMON: return "Bad month";                                      
  case HBFLAB: return "Bad fiscal year label";                          
  case HBMISS: return "Bad missing value type";                         
  case HBINDX: return "Bad value index";                                
  case HNWILD: return "Wildcarding has not been initialized";           
  case HBNCHR: return "Bad number of characters";                       
  case HBGROW: return "Bad growth factor";                              
  case HQUOTA: return "Too many files open or no disk space available"; 
  case HOLDDB: return "Can't update or share old data bases";           
  case HMPOST: return "Data base must be posted";                       
  case HSPCDB: return "Can't write to a special data base";             
  case HBFLAG: return "Bad flag";                                       
  case HPACK:  return "Can't perform operation on packed data base";    
  case HNEMPT: return "Data base is not empty";                         
  case HBATTR: return "Bad attribute name";                             
  case HDUP:   return "A duplicate was ignored";                        
  case HBYEAR: return "Bad year";                                       
  case HBPER:  return "Bad period";                                     
  case HBDAY:  return "Bad day";                                        
  case HBDATE: return "Bad date";                                       
  case HBSEL:  return "Bad date selector";                              
  case HBREL:  return "Bad date relation";                              
  case HBCPU:  return "Unauthorized CPU ID or hardware";                
  case HEXPIR: return "Expired dead date";                              
  case HBPROD: return "Unauthorized product";                           
  case HBUNIT: return "Bad number of units";                            
  case HIFAIL: return "HLI internal failure";                           
  }

  sprintf(buf, "Unknown Fame database error");
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
    int     dbkey;
    char    *svcuser;
    struct passwd *ps;
    char    fame_name[Padi_MAX_STR_LEN + 1];	/* fame writes to this
						   string! */
    char    dbpath[Padi_MAX_STR_LEN + 1];
    int     type,
            freq,
            class,
            syear,
            sprd,
            eyear,
            eprd,
            range[8],
            numobs,
	    status;
    int     in_syear,
            in_sprd,
            in_eyear,
            in_eprd;
    int     ismiss,
            length;
    size_t  struct_len,
            data_len,
            buf_len;
    char   *s;
    PadiString_t *strp;
    PadiPrecision_t *vector;
    property_t property;
    PadiSeriesResult_tp result;
    PadiInfo_tp info;

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
	printf("svcuser %s  passed  %s \n",svcuser,user);
#endif

    if ( strcmp(svcuser,user) )
    {
	printf("svcuser %s  passed  %s \n",svcuser,user);
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


    /* build the full path */
    strcpy(dbpath, db_name);

#ifdef DEBUG
printf("Try to open %s in READ mode\n",dbpath);
#endif

    /* open the object's database */
    cfmopdb( &(result->status), &dbkey,dbpath,HRMODE);
    if (result->status )
	return result;

#ifdef DEBUG
printf("cfmopdb : %d\n",result->status);
#endif

    cfmosiz(&(result->status), dbkey, strcpy(fame_name, object_name),
	    &class, &type, &freq, &syear, &sprd, &eyear, &eprd);
#ifdef DEBUG
printf("cfmosiz : %d\n",result->status);
#endif

    if (result->status)
    {
	        cfmcldb(&status,dbkey);
	return result;
    }

    if (class == HFRMLA)
    {
	char    cmd_buf[1024];

	/*
	   FORMULA or FUNCTION
	*/

	/* open a fame session work database */
	if ((result->status = OpenDb("*", HUMODE, &dbkey)))
	{
	    cfmcldb(&status,dbkey);
	    return result;
	}

	sprintf(cmd_buf, "overwrite on; ignore on;");
	cfmfame(&(result->status), cmd_buf);
	if (result->status)
	{
	        cfmcldb(&status,dbkey);
	    return result;
	}


	/* open primary database */
	if ((result->status = ReleaseDb(dbpath, StrToMode(property.access))))
	    return result;

	sprintf(cmd_buf, "try; open <access %s> file(\"%s\") as %s; end try",
		"READ", db_name, "FRM_DB");
	cfmfame(&(result->status), cmd_buf);
	if (result->status)
	{
	        cfmcldb(&status,dbkey);
	    return result;
	}




	/* copy the formula into a work series (scalar) */
	sprintf(cmd_buf, "work'%s = %s'%s", object_name, "FRM_DB",
		object_name);

	cfmfame(&(result->status), cmd_buf);
	if (result->status)
	{
	        cfmcldb(&status,dbkey);
	    return result;
	}


	/* RESET NUMOBS to -1: bugfix hp */
	numobs = -1;


	sprintf(cmd_buf, "close all");
	cfmfame(&(result->status), cmd_buf);
	if (result->status)
	{
	        cfmcldb(&status,dbkey);
	    return result;
	}

	cfmosiz(&(result->status), dbkey, object_name,
		&class, &type, &freq, &syear, &sprd, &eyear, &eprd);

    }


    if (result->status == HNOOBJ)
    {
	result->status = Padi_SUCCEED;	/* non-existent (i.e. empty) series */
	eprd = 0;

    }

    if (result->status)
    {
	        cfmcldb(&status,dbkey);
	return result;
    }

    /* convert the input range */
#ifdef DEBUG
printf("freq : %d\n",freq);
#endif

    if (freq == HUNDFX)
    {
	in_syear = in_eyear = -1;
	in_sprd = in_eprd = -1;
    }
    else if (freq == HCASEX)
    {
	in_syear = in_eyear = -1;
	in_sprd = (*(in_range->start)) ? atoi(in_range->start) : -1;
	in_eprd = (*(in_range->end)) ? atoi(in_range->end) : -1;
    }
    else
    {
	int     date;

	if (*(in_range->start))
	{
	    if (!(result->status = IsoToDate(freq, in_range->start, &date)))
		cfmdatp(&(result->status), freq, date, &in_syear, &in_sprd);
	    if (result->status)
	    {
	        cfmcldb(&status,dbkey);
		return result;
	    }

	}
	else
	{
	    in_syear = in_sprd = -1;
	}

	if (*(in_range->end))
	{
	    if (!(result->status = IsoToDate(freq, in_range->end, &date)))
		cfmdatp(&(result->status), freq, date, &in_eyear, &in_eprd);
	    if (result->status)
	    {
	        cfmcldb(&status,dbkey);
		return result;
	    }

	}
	else
	{
	    in_eyear = in_eprd = -1;
	}
    }

    if (in_range->nobs)
    {
	numobs = in_range->nobs;
	if (in_eprd != -1)
	    in_syear = in_sprd = -1;
    }
    else
    {
	numobs = -1;
    }

    if (class == HSCALA)
    {
	numobs = 1;
    }
    else if (eprd < 1)
    {
	numobs = 0;		/* empty series */
    }
    else
    {
	if (numobs != -1)
	{
	    if (in_sprd == -1 && in_eprd == -1)
	    {
		in_eyear = eyear;
		in_eprd = eprd;
	    }
	}
	else
	{
	    if (in_sprd == -1)
	    {
		in_syear = syear;
		in_sprd = sprd;
	    }
	    if (in_eprd == -1)
	    {
		in_eyear = eyear;
		in_eprd = eprd;
	    }
	}


	/* set fame range */
	cfmsrng(&(result->status), freq, &in_syear, &in_sprd, &in_eyear,
		&in_eprd, range, &numobs);
#ifdef DEBUG
printf("cfmsrng : %d\n",result->status);
#endif


	if (result->status)
	{
	        cfmcldb(&status,dbkey);
	    return result;
	}

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
    strcpy(info->desc, property.desc);
    strcpy(info->class, ClassToStr(class));
    strcpy(info->type, TypeToStr(type));
    if (!(!strcmp(property.fame_class, "SERIES") ||
	  !strcmp(property.fame_class, "SCALAR")))
	strcpy(info->access, "READ");
    else
	strcpy(info->access, property.access);
    strcpy(info->frequency, FreqToStr(freq));
    strcpy(info->start, property.start);
    strcpy(info->end, property.end);
    strcpy(info->format, property.format);


    if (numobs > 0)
    {
	/* read data */
	if (type == HPRECN)
	{
	    cfmrrng(&(result->status), dbkey,
		    fame_name,
		    range, (double *) vector,
		    (in_range->do_missing) ? HTMIS : HNTMIS,
		    in_range->missing_translation);
	}
	else
	{
	    float  *v;
	    int     i;
	    if (!(v = (float *) malloc(numobs * sizeof(float))))
		return NULL;
	    cfmrrng(&(result->status), dbkey,
		    fame_name,
		    range, v,
		    (in_range->do_missing) ? HTMIS : HNTMIS,
		    in_range->missing_translation);
	    for (i = 0; i < numobs; i++)
	    {
		vector[i] = (double) v[i];
	    }
	    free((char *) v);
	}
#ifdef DEBUG
printf("cfmrrng : %d\n",result->status);
#endif
	if (result->status)
	{
	        cfmcldb(&status,dbkey);
	    return result;
	}

    }

    /* store the range start and end */
    if (freq == HCASEX)
    {
	sprintf(info->start, "%d", sprd);
	sprintf(info->end, "%d", eprd);
	sprintf(data_range->start, "%d", in_sprd);
	sprintf(data_range->end, "%d", in_eprd);
    }
    else if (freq != HUNDFX)
    {
	int     date;

	cfmpdat(&(result->status), freq, &date, syear, sprd);
	if (!(result->status))
	    if ((result->status = DateToIso(freq, date, info->start)))
	    {
	        cfmcldb(&status,dbkey);
		return result;
	    }

	cfmpdat(&(result->status), freq, &date, eyear, eprd);
	if (!(result->status))
	    if ((result->status = DateToIso(freq, date, info->end)))
	    {
	        cfmcldb(&status,dbkey);
		return result;
	    }

	cfmpdat(&(result->status), freq, &date, in_syear, in_sprd);
	if (!(result->status))
	    if ((result->status = DateToIso(freq, date, data_range->start)))
	    {
	        cfmcldb(&status,dbkey);
		return result;
	    }

	cfmpdat(&(result->status), freq, &date, in_eyear, in_eprd);
	if (!(result->status))
	    if ((result->status = DateToIso(freq, date, data_range->end)))
	    {
	        cfmcldb(&status,dbkey);
		return result;
	    }

    }

    /* store the range translation */
    data_range->do_missing = in_range->do_missing;
    memcpy((char *) (data_range->missing_translation),
	   (char *) (in_range->missing_translation),
	   sizeof(in_range->missing_translation));

    /* close primary database */
    cfmcldb(&(result->status),dbkey);
    if (result->status)
    {
#ifdef DEBUG
printf("cfmcldb : %d\n",result->status);
#endif

	return result;
    }

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
    PadiString_t series_name = info->name;
    PadiRange_tp in_range = &(series->range);
    int     dbkey,
            freq,
            class,
            type,
            date,
            numobs,
	    status;
    int     range[8];
    int     in_syear,
            in_sprd,
            in_eyear,
            in_eprd;
    int     ismiss,
            length;
    char    fame_name[Padi_MAX_STR_LEN + 1];	/* fame writes to this
						   string! */
    char    dbpath[Padi_MAX_STR_LEN + 1];
    PadiResult_tp result;
    size_t  struct_len;
    precision_t *pvector = (precision_t *) (series->data.data_val);
    numeric_t *nvector = (numeric_t *) (series->data.data_val);
    char    *svcuser;
    struct passwd *ps;


    struct_len = sizeof(*result);
    if (!(result = (PadiResult_tp) malloc(struct_len)))
	return NULL;
    bzero((char *) result, struct_len);

     ps = getpwuid(geteuid());
     svcuser = ps->pw_name;
    if ( strcmp(svcuser,user))
    {
	result->status = Padi_UNATH_USER;
	return result;

    }

    if (!(new->dbname && *(new->dbname)))
    {
	result->status = Padi_MISSING_DBNAME;
	return result;

    }


    if (!(class = StrToClass(info->class)))
    {
	result->status = Padi_UNSUPPORTED_CLASS;
	return result;

    }

    if (!(freq = StrToFreq(info->frequency)))
    {
	result->status = Padi_UNSUPPORTED_FREQ;
	return result;

    }

    if (!(type = StrToType(info->type)))
    {
	result->status = Padi_UNSUPPORTED_TYPE;
	return result;

    }

    /* build the full path */
    strcpy(dbpath, dbname);

#ifdef DEBUG
printf("database,object %s %s\n",dbpath, series_name);
printf("Try to open %s in SHARED mode\n",dbpath);
#endif

    cfmopdb(&(result->status), &dbkey,dbpath,HSMODE);
    if (result->status)
    {
#ifdef DEBUG
printf("cfmopdb (SHARED): %d\n",result->status);
printf("Try to open %s in CREATE mode\n",dbpath);
#endif
        cfmopdb(&(result->status), &dbkey,dbpath,HCMODE);
	if (result->status)
	    return result;
    }
#ifdef DEBUG
printf("cfmopdb (CREATE) : %d\n",result->status);
#endif

    /* create the object */
    cfmalob(&(result->status), dbkey, strcpy(fame_name, info->name),
	    class, freq, type, HBSDAY, HOBEND,
	    in_range->nobs, in_range->nobs * 16, (float) 2.0);

#ifdef DEBUG
printf("cfmalob : %d\n",result->status);
#endif

    if (result->status == HOEXI)
    {
	cfmdlob(&(result->status), dbkey, fame_name);

	if (result->status != HSUCC)
	{
	    cfmcldb(&status,dbkey);
	    return result;
	}

	cfmalob(&(result->status), dbkey, strcpy(fame_name, info->name),
		class, freq, type, HBSDAY, HOBEND,
		in_range->nobs, in_range->nobs * 16, (float) 2.0);
    }

#ifdef DEBUG
printf("cfmalob : %d\n",result->status);
#endif

    if (result->status != HSUCC)
    {
	cfmcldb(&status,dbkey);
	return result;
    }

#ifdef DEBUG
printf("freq : %d\n",freq);
#endif

    if (freq != HUNDFX)
    {
	if (freq == HCASEX)
	{
	    in_syear = in_eyear = -1;
	    in_sprd = (*(in_range->start)) ? atoi(in_range->start) : -1;
	    in_eprd = (*(in_range->end)) ? atoi(in_range->end) : -1;
	}
	else
	{

	    if (*(in_range->start))
	    {
		if (!(result->status = IsoToDate(freq, in_range->start, &date)))
		    cfmdatp(&(result->status), freq, date, &in_syear, &in_sprd);
		if (result->status)
		{
	    	    cfmcldb(&status,dbkey);
		    return result;

		}

	    }
	    else
	    {
		in_syear = in_sprd = -1;
	    }

	    if (*(in_range->end))
	    {
		if (!(result->status = IsoToDate(freq, in_range->end, &date)))
		    cfmdatp(&(result->status), freq, date, &in_eyear, &in_eprd);
		if (result->status)
		{
	    	    cfmcldb(&status,dbkey);
		    return result;

		}

	    }
	    else
	    {
		in_eyear = in_eprd = -1;
	    }
	}


	/* set fame range */
	if (in_range->nobs)
	{
	    numobs = in_range->nobs;
	    if (in_eprd != -1)
		in_syear = in_sprd = -1;
	}
	else
	{
	    numobs = -1;
	}

	cfmsrng(&(result->status), freq, &in_syear, &in_sprd, &in_eyear, &in_eprd, range, &numobs);

#ifdef DEBUG
printf("cfmsrng : %d\n",result->status);
#endif

	if (result->status)
	{
	    cfmcldb(&status,dbkey);
	    return result;

	}

    }				/* if (freq != HUNDFX) */

    /* write data */
    cfmwrng(&(result->status), dbkey,
	    fame_name,
	    range, (float *) pvector,
	    (in_range->do_missing) ? HTMIS : HNTMIS,
	    in_range->missing_translation);

#ifdef DEBUG
printf("cfmwrng : %d\n",result->status);
#endif

    if (result->status)
    {
	    cfmcldb(&status,dbkey);
	return result;

    }

    cfmcldb(&(result->status),dbkey);
    if (result->status)
	return result;

    return result;
}

