/*  
** Portmapper Server Protocol Channel
**				
** Copyright 1995, 1996  Bank of Canada.
**
** The user of this software has the right to use, reproduce and distribute it.
** Bank of Canada makes no warranties with respect to the software or its
** fitness for any particular purpose. The software is distributed by the Bank
** of Canada solely on an "as is" basis. By using the software, user agrees to
** accept the entire risk of using this software.
*/

#if !defined(FAME_SVC) && !defined(FS_SVC)
#define PADI_CLIENT
#endif
 
#include "padi.h"

#include "sys/types.h"
#include "time.h"

private struct Channel_t {
  PadiServer_t server;
  struct timeval timeout;
  CLIENT * clnt;
  unsigned idle;
} Channel = {
  {NULL, 0},
  {25, 0},
  NULL,
  0
};

export void PadiChannelClose PARAMVOID
{
  if (Channel.clnt) {
    auth_destroy(Channel.clnt->cl_auth);
    clnt_destroy(Channel.clnt);
    Channel.clnt = NULL;
  }
}

export void PadiSetServer PARAM1(PadiServer_tp, server)
{
  PadiServer_t dummy;
 
  if (!server) {
    server = &dummy;
    server->name = Padi_EMPTY_STR;
    server->port = 0;
  }

  if (!(Channel.server.name && server->name) ||
      strcmp(Channel.server.name, server->name)) {
    PadiChannelClose();
  }
  if (Channel.server.name) {
    free(Channel.server.name);
    Channel.server.name = NULL;
  }
  if (server) {
    if (server->name)
      Channel.server.name = strcpy(malloc(strlen(server->name) + 1), 
                                   server->name);
  }
}
 
#ifdef PADI_CLIENT

export PadiStatus_t PadiChannelOpen PARAMVOID
{
  struct ct_data ;
  PadiStatus_t status; 

  if (Channel.clnt) 
    return Padi_SUCCEED;

  if (!(Channel.clnt = clnt_create(Channel.server.name, PADIPROG, PADIVERS, 
                                   "tcp"))) {
    PadiChannelClose();
    return status = Padi_NO_CONNECTION;

  }

  if (!clnt_control(Channel.clnt, CLSET_TIMEOUT, (char *)&(Channel.timeout))) {
    PadiChannelClose();
    return status = Padi_CLNT_SET_TIMEOUT;

  }

  Channel.clnt->cl_auth = authunix_create_default();
  return status = Padi_SUCCEED;
}

export void PadiSetTimeOut PARAM2(int, second, int, usecond)
/*
** Set the time out period.
**
** second (input)
**   The number of seconds to wait for a server response.
**
** usecond (input)
**   The additional number of micro seconds to wait for a server response.
**
*/
{
  Channel.timeout.tv_sec = second;
  Channel.timeout.tv_usec = usecond;
}
 
export void PadiSetIdle PARAM1(unsigned, idle)
{
  Channel.idle = idle;
}

export PadiResult_tp PadiClntCall PARAM6(PadiServer_tp, server, u_long, procnum, xdrproc_t, inproc, PadiData_tp, in, xdrproc_t, outproc, PadiSize_t, result_size)
/*
** Do a RPC call with result structure created using malloc().
** Call PadiFreeResult() to free this structure.
*/
{
  PadiResult_tp result;
  PadiStatus_t status; 

  PadiSetServer(server);

  if (!(result = (PadiResult_tp)malloc(result_size))) 
    return NULL;

  bzero((char *)result, result_size);
  result->xdr_proc = (PadiFunc_tp)NULL;

  if ((result->status = PadiChannelOpen())) {
    return result;

  } 
  if ((status = clnt_call(Channel.clnt, procnum, inproc, (caddr_t)in,
                          outproc, (caddr_t)result, Channel.timeout)) 
                != RPC_SUCCESS) {
    result->xdr_proc = (PadiFunc_tp)NULL;
    result->status = status + Padi_CLNT_STATUS;
  } else {
    result->xdr_proc = (PadiFunc_tp)outproc;
  }

  if (result->status || !(Channel.idle))
    PadiChannelClose();

  return result;
}

#else /* PADI_CLNT */

export PadiStatus_t PadiServe PARAM2(PadiServer_tp, server, dispatch_t, dispatch)
{
  import FILE * logfp;
  SVCXPRT *transp;

  PadiSetServer(server);
 
  (void) pmap_unset(PADIPROG, PADIVERS);
 
  if(!(transp = svctcp_create(RPC_ANYSOCK, 0, 0))) {
    return Padi_TCP_CREATE;
 
  }
 
  if (!svc_register(transp, PADIPROG, PADIVERS, dispatch, IPPROTO_TCP)) {
    return Padi_TCP_REGISTER;
 
  }
 
  svc_run();

  return Padi_SUCCEED;
}

#endif /* PADI_CLIENT */

