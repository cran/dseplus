/*  
** Protocol for Application -- Database Interface
**
** Copyright (C) 1995, 1996 Bank of Canada
**
** Any person using this software has the right to use, reproduce and 
** distribute it.
**
** The Bank does not warrant the accuracy of the information contained in 
** the software.  User assumes all liability for any loss or damages, 
** either direct or indirect, arising from the use of the software.
*/

#ifndef PADI_H
#define PADI_H

#ifndef _HPUX_SOURCE
#define PORTMAP
#endif

#ifdef MSDOS
typedef unsigned bool_t;
typedef unsigned u_int;
#define _RPC_TYPES_H "msdos_kludge"
#else
#include <rpc/rpc.h>
#include <rpc/types.h>
#endif

#include "rpcx.h"
#include "define.h"

/* status code */
#define Padi_SUCCEED          0
#define Padi_STATUS           1000
#define Padi_FAIL             Padi_STATUS + 0
#define Padi_OUT_OF_MEMORY    Padi_STATUS + 1
#define Padi_TOO_MANY_OPEN    Padi_STATUS + 2
#define Padi_FILE_OPEN        Padi_STATUS + 5
#define Padi_TCP_CREATE       Padi_STATUS + 6
#define Padi_TCP_REGISTER     Padi_STATUS + 7
#define Padi_FREE_ARG         Padi_STATUS + 8
#define Padi_READ_ACCESS      Padi_STATUS + 9
#define Padi_MISSING_OBJECT   Padi_STATUS + 10
#define Padi_UNKNOWN_OBJECT   Padi_STATUS + 11
#define Padi_OVERWRITE        Padi_STATUS + 12
#define Padi_MISSING_FORMAT   Padi_STATUS + 13
#define Padi_UNSUPPORTED_FORMAT Padi_STATUS + 14
#define Padi_MISSING_TYPE     Padi_STATUS + 15
#define Padi_TYPE_MISMATCH    Padi_STATUS + 16
#define Padi_OBJECT_EXISTS    Padi_STATUS + 17
#define Padi_MISSING_DBNAME   Padi_STATUS + 18
#define Padi_UNSUPPORTED_CLASS  Padi_STATUS + 19
#define Padi_UNSUPPORTED_FREQ   Padi_STATUS + 20 
#define Padi_UNSUPPORTED_TYPE   Padi_STATUS + 21 
#define Padi_UNATH_USER         Padi_STATUS + 22
#define Padi_UNSUPPORTED_DATE	Padi_STATUS + 23 
#define Padi_SVC_STATUS       2000
#define Padi_CLNT_STATUS      3000
#define Padi_NO_CONNECTION    Padi_CLNT_STATUS + 0
#define Padi_CLNT_SET_TIMEOUT Padi_CLNT_STATUS + 333
#define Padi_SYSTEM_STATUS    4000
#define Padi_HOST_STATUS      5000
#define Padi_HOST_NOT_FOUND   Padi_HOST_STATUS + 1
#define Padi_SIGNAL_STATUS    9000

/* severity code */
#define Padi_WARNING   1
#define Padi_ERROR     2
#define Padi_FATAL     3

#define DEFAULT_NC   170141345719746060945050695293894393856.000000
#define DEFAULT_ND   170141670238299719371777478449914970112.000000
#define DEFAULT_NA   170141507979022890158414086871904681984.000000


#define Padi_MAX_NAME_LEN 55 /* Fame has limited name length! */
#define Padi_MAX_STR_LEN 255 
#define Padi_DATE_TIME_LEN 16 
#define Padi_NUM_INFO_STR (sizeof(PadiInfo_t)/sizeof(PadiString_t))

typedef char *PadiData_tp;
typedef PadiResult_t *PadiResult_tp;
typedef PadiInfo_t *PadiInfo_tp;
typedef PadiInfoResult_t *PadiInfoResult_tp;
typedef PadiRange_t *PadiRange_tp;
typedef PadiRangeArg_t *PadiRangeArg_tp;
typedef PadiInitArg_t *PadiInitArg_tp;
typedef PadiTermArg_t *PadiTermArg_tp;
typedef PadiSeries_t *PadiSeries_tp;
typedef PadiSeriesResult_t *PadiSeriesResult_tp;
typedef PadiNewSeries_t *PadiNewSeries_tp;
typedef PadiInfoArg_t *PadiInfoArg_tp;
typedef PadiDestroyArg_t *PadiDestroyArg_tp;
 
typedef struct PadiServer_t {
  PadiString_t name;
  u_short port;
} PadiServer_t, *PadiServer_tp;
 
/* public data */
import PadiBoolean_t PadiVerbose;
import PadiString_t Padi_EMPTY_STR;

/* functions */
import PadiStatus_t PadiInitialize PROTO1(PadiInitArg_tp, object);
import PadiStatus_t PadiTerminate PROTO1(PadiTermArg_tp, object);
import char * PadiStatus PROTO1(int, status);
import void PadiError PROTO4(FILE *, fp, PadiString_t, source, PadiStatus_t, status, int, severity);
import void PadiFreeResult PROTO1(PadiResult_tp, result);
import PadiInfoResult_tp PadiGetInfo PROTO2(PadiServer_tp, server, PadiInfoArg_tp, object);
import PadiSeriesResult_tp PadiGetSeries PROTO2(PadiServer_tp, server, PadiRangeArg_tp, range);
import PadiResult_tp PadiNewSeries PROTO2(PadiServer_tp, server, PadiNewSeries_tp, arg);
import PadiResult_tp PadiDestroy PROTO2(PadiServer_tp, server, PadiDestroyArg_tp, dobject);

#include <errno.h>
 
#define PadiSocketError() (Padi_SYSTEM_STATUS + errno)
 
import void PadiChannelClose PROTOVOID;
import void PadiSetServer PROTO1(PadiServer_tp, server);
 
import PadiStatus_t PadiChannelOpen PROTOVOID;
import void PadiSetTimeOut PROTO2(int, second, int, usecond);
import void PadiSetIdle PROTO1(unsigned, idle);
import PadiResult_tp PadiClntCall PROTO6(PadiServer_tp, server, u_long, procnum,
 xdrproc_t, inproc, PadiData_tp, in, xdrproc_t, outproc, PadiSize_t, result_size
);
 
typedef void (* dispatch_t) PROTO2(struct svc_req *, request, SVCXPRT *, xprt);
import PadiStatus_t PadiServe PROTO2(PadiServer_tp, server, dispatch_t, dispatch
);

#ifdef GMALLOC

import size_t ngmalloced;
import char *gmalloc PROTO1(size_t, n);
import char *grealloc PROTO2(char *, p, size_t, n);
import char *gfree PROTO1(char *, p);
#define malloc(n) gmalloc(n)
#define realloc(p, n) grealloc(p, n)
#define free(p) gfree(p)

#endif /* GMALLOC */

#endif /* PADI_H */
