#
# dse/src/pac/padi/src/Makefile

VPATH=../exec

CC=gcc
# consider using rpcgen -N
RPCGEN= rpcgen -C

PKGNAME = padi

#VPATH = @srcdir@
#srcdir = @srcdir@
#top_srcdir = @top_srcdir@

#top_builddir = ../../../..
#subdir = src/library/$(pkg)/src
#R_HOME = $(top_builddir)
R_HOME = $(R)

PRGS = getpadi putpadi objectinfo
FAMEPRGS = compare.pc object.pc target.pc sys.pc

LIBSRC = padi.c padiutil.c channel.c rpcx_xdr.c gmalloc.c
SRC = getpadi.c  putpadi.c  $(LIBSRC)
#SRC = padi.c  channel.c  getpadi.c  gmalloc.c  padiutil.c  putpadi.c  rpcx_xdr.c  
OBJS = $(SRC:.c=.o) 

OTHERSRC = rpcx_clnt.c  rpcx_svc.c  objectinfo.c simple_svc.c
OTHEROBJECTS = $(OTHERSRC:.c=.o) 

PKG_CFLAGS = -ansi -pedantic -I$(PADI_HOME)/include 
#PKG_CPPFLAGS = -DGMALLOC -trigraphs -undef
PKG_CPPFLAGS = -DGMALLOC -trigraphs
#PKG_CPPFLAGS = -DGMALLOC 

ifeq (${shell uname}, SunOS)
  LD_LIBS=-lm -lc -lnsl
else
  LD_LIBS= -L/usr/local/lib 
endif

all: rpcx.h padi.so ../exec/padi_simple_svc ../exec/getpwuid  ../exec/x11arima \
       $(PRGS:%=../exec/%)


$(PKGNAME).so: $(OBJS) rpcx.h
	$(CC) -shared $(LD_LIBS) -o $(PKGNAME).so $(OBJS)
#	$(R_HOME)/bin/R CMD SHLIB -o $(PKGNAME).so $(OBJS)

../exec/padi_simple_svc:  simple_svc.c $(LIBSRC)  
	$(RM) $(notdir $@)
	$(CC) -Wall $(PKG_CFLAGS) $(PKG_CPPFLAGS) -DFS_SVC \
	    simple_svc.c $(LIBSRC)  $(LD_LIBS) -o padi_simple_svc
	mv padi_simple_svc ../exec

#needs work
fame_svc:	fame_svc.c plus fame hli
	echo not yet

$(PRGS:%=../exec/%):../exec/%: %.c $(LIBSRC)
	$(RM) $(notdir $@)
	$(CC) -Wall -DMAIN $(LD_LIBS) $(notdir $@).c $(LIBSRC) -o $(notdir $@)
	mv $(notdir $@) ../exec

../exec/x11arima: getpadi.c $(LIBSRC)
	$(RM) $(notdir $@)
	$(CC) -Wall -DMAIN -DX11ARIMA $(LD_LIBS) getpadi.c $(LIBSRC) -o $(notdir $@)
	mv $(notdir $@) ../exec

../exec/getpwuid:	getpwuid.c
	$(RM) $(notdir $@)
	$(CC) -Wall -DMAIN $(LD_LIBS) getpwuid.c gmalloc.c -o getpwuid
	mv $(notdir $@) ../exec
	
rpcx_clnt.c  rpcx_svc.c  rpcx_xdr.c:	rpcx.h

#rpcx.h: rpcx.x
%.h: %.x
	$(RPCGEN) $<

../exec/$(FAMEPRGS):%.pc: %.pro
	fameCompile $<
	mv $(notdir $@) ../exec

clean:
	@$(RM) $(PRGS) getpwuid padi_simple_svc x11arima
	@$(RM) *.so *.o 
	@$(RM) rpcx_clnt.c  rpcx.h  rpcx_svc.c  rpcx_xdr.c

cvsclean:	clean
	@$(RM) $(PRGS:%=../exec/%)
	@$(RM) ../exec/getpwuid ../exec/padi_simple_svc ../exec/x11arima
	@$(RM) ../exec/*.pro ../exec/*server.log.*

distclean: clean
#	@-rm -f Makefile
maintainer-clean: distclean

## Automagically generated dependencies:
