/*  
** Protocol for Application -- Database Interface
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

/* type definitions */
typedef string PadiString_t<>;
typedef unsigned int PadiSize_t;
typedef long PadiFunc_tp;
typedef int PadiStatus_t;
typedef char PadiBoolean_t;
typedef double PadiPrecision_t;

struct PadiResult_t {
  PadiFunc_tp xdr_proc;
  PadiStatus_t status;
};


struct PadiInfo_t {
  PadiString_t name;
  PadiString_t desc;
  PadiString_t class;
  PadiString_t type;
  PadiString_t access;
  PadiString_t frequency;
  PadiString_t start;
  PadiString_t end;
  PadiString_t format;
};

struct PadiInfoResult_t {
  PadiFunc_tp xdr_proc;
  PadiStatus_t status;
  PadiInfo_t info;
};
 
struct PadiRange_t {
  PadiString_t start;
  PadiString_t end;
  PadiSize_t nobs;
  PadiBoolean_t do_missing;
  PadiPrecision_t missing_translation[3];
  PadiString_t format;
};

struct PadiRangeArg_t {
  PadiString_t user;
  PadiString_t password;
  PadiString_t object_name;
  PadiString_t db_name;
  PadiRange_t range;
};

struct PadiDestroyArg_t {
  PadiString_t user;
  PadiString_t password;
  PadiString_t object_name;
  PadiString_t db_name;
};

struct PadiInfoArg_t {
  PadiString_t user;
  PadiString_t password;
  PadiString_t object_name;
  PadiString_t db_name;
};

struct PadiInitArg_t {
  PadiString_t user;
  PadiString_t password;
  PadiString_t object_dbname;
};

struct PadiTermArg_t {
  PadiString_t user;
  PadiString_t password;
};


struct PadiSeries_t {
  PadiInfo_t info;
  PadiRange_t range;
  PadiPrecision_t data<>;
};

struct PadiSeriesResult_t {
  PadiFunc_tp xdr_proc;
  PadiStatus_t status;
  PadiSeries_t series;
};

struct PadiNewSeries_t {
  PadiString_t user;
  PadiString_t password;
  PadiString_t dbname;
  PadiSeries_t series;
};


program PADIPROG {
  version PADIVERS {
     PadiInfoResult_t GETINFO (PadiInfoArg_t) = 1;
     PadiSeriesResult_t GETSERIES (PadiRangeArg_t) = 2;
     PadiResult_t NEWSERIES (PadiNewSeries_t) = 3;
     PadiResult_t DESTROY (PadiDestroyArg_t) = 4;
     PadiSeriesResult_t GETLOCAL (PadiRangeArg_t) = 5;
     PadiResult_t PUTLOCAL (PadiNewSeries_t) = 6;
  } = 1;
} = 0x20004444;

