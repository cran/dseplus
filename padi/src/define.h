/*
** Standard Definitions
*/
#ifndef DEFINE_H
#define DEFINE_H

#if defined(sparc) && !defined(__sparc__)
#define __sparc__
#endif

#ifndef STDIO_H
#ifndef _STDIO_H
#include "stdio.h"
#endif
#endif

#ifndef _DEC
#include "stdlib.h"
#else
#include "values.h"
#endif

#include "limits.h"

#ifndef sun
#define _ANSI_C
#endif

#ifdef _ANSI_C
#include "float.h"
#else
#include "floatingpoint.h"
#define FLT_EPSILON 1.192e-07F
#endif

#include "assert.h"
#include "memory.h"
#include "string.h"
#include "ctype.h"

#include "stdtyp.h"

#ifndef private
#define import extern
#define export
#define private static
#define forward extern
#endif

#ifndef OK_EXIT
#define OK_EXIT 0
#endif

#ifndef ERROR_EXIT
#define ERROR_EXIT -1
#endif

#ifndef YES
#define YES 1
#define NO 0
#endif

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#ifndef FOREVER
#define FOREVER for(;;)
#endif

#ifndef NELEMENT
#define NELEMENT(a) (sizeof(a) / sizeof(*a))
#endif

#ifndef MAX
#define MAX(a, b) ((a < b) ? b : a)
#define MIN(a, b) ((a > b) ? b : a)
#endif

#ifndef ABS
#define ABS(a) (((a) < 0) ? -(a) : (a))
#endif

#ifndef FUZZY_EQ
#define FUZZY_EQ(a, b) (fabs(a - b) < FLT_EPSILON)
#endif

#ifndef ALIGN
/* align pointer p to point to type t, e.g.:
       q = align(double, p) where q >= p  */
#define ALIGN(t, p) (t *)(((long)(p) + (sizeof(t)-1)) & ~(sizeof(t)-1))
#endif

#ifndef PROTO1
#ifdef _ANSI_C
#define PROTOVOID (void)
#define PROTO1(t1, p1) (t1 p1)
#define PROTO2(t1, p1, t2, p2) (t1 p1, t2 p2)
#define PROTO3(t1, p1, t2, p2, t3, p3) (t1 p1, t2 p2, t3 p3)
#define PROTO4(t1, p1, t2, p2, t3, p3, t4, p4) (t1 p1, t2 p2, t3 p3, t4 p4)
#define PROTO5(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5) (t1 p1, t2 p2, t3 p3, t4 p4, t5 p5)
#define PROTO6(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6) (t1 p1, t2 p2, t3 p3, t4 p4, t5 p5, t6 p6)
#define PROTO7(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7) (t1 p1, t2 p2, t3 p3, t4 p4, t5 p5, t6 p6, t7, p7)
#define PROTO8(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8) (t1 p1, t2 p2, t3 p3, t4 p4, t5 p5, t6 p6, t7 p7, t8 p8)
#define PROTO9(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8, t9, p9) (t1 p1, t2 p2, t3 p3, t4 p4, t5 p5, t6 p6, t7 p7, t8 p8, t9 p9)
#define PARAMVOID (void)
#define PARAM1(t1, p1) (t1 p1)
#define PARAM2(t1, p1, t2, p2) (t1 p1, t2 p2)
#define PARAM3(t1, p1, t2, p2, t3, p3) (t1 p1, t2 p2, t3 p3)
#define PARAM4(t1, p1, t2, p2, t3, p3, t4, p4) (t1 p1, t2 p2, t3 p3, t4 p4)
#define PARAM5(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5) (t1 p1, t2 p2, t3 p3, t4 p4, t5 p5)
#define PARAM6(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6) (t1 p1, t2 p2, t3 p3, t4 p4, t5 p5, t6 p6)
#define PARAM7(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7) (t1 p1, t2 p2, t3 p3, t4 p4, t5 p5, t6 p6, t7 p7)
#define PARAM8(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8) (t1 p1, t2 p2, t3 p3, t4 p4, t5 p5, t6 p6, t7 p7, t8 p8)
#define PARAM9(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8, t9, p9) (t1 p1, t2 p2, t3 p3, t4 p4, t5 p5, t6 p6, t7 p7, t8 p8, t9 p9)
#else
#define PROTOVOID ()
#define PROTO1(t1, p1) () 
#define PROTO2(t1, p1, t2, p2) ()
#define PROTO3(t1, p1, t2, p2, t3, p3) ()
#define PROTO4(t1, p1, t2, p2, t3, p3, t4, p4) ()
#define PROTO5(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5) ()
#define PROTO6(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6) ()
#define PROTO7(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7) ()
#define PROTO8(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8) ()
#define PROTO9(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8, t9, p9) ()
#define PARAMVOID ()
#define PARAM1(t1, p1) (p1) t1 p1;
#define PARAM2(t1, p1, t2, p2) (p1, p2) t1 p1; t2 p2;
#define PARAM3(t1, p1, t2, p2, t3, p3) ( p1, p2, p3) t1 p1; t2 p2; t3 p3;
#define PARAM4(t1, p1, t2, p2, t3, p3, t4, p4) (p1, p2, p3, p4) t1 p1; t2 p2; t3 p3; t4 p4;
#define PARAM5(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5) (p1,p2,p3,p4,p5) t1 p1; t2 p2; t3 p3; t4 p4; t5 p5;
#define PARAM6(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6) (p1,p2,p3,p4,p5,p6) t1 p1; t2 p2; t3 p3; t4 p4; t5 p5; t6 p6;
#define PARAM7(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7) (p1,p2,p3,p4,p5,p6,p7) t1 p1; t2 p2; t3 p3; t4 p4; t5 p5; t6 p6; t7 p7;
#define PARAM8(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8) (p1,p2,p3,p4,p5,p6,p7,p8) t1 p1; t2 p2; t3 p3; t4 p4; t5 p5; t6 p6; t7 p7; t8 p8;
#define PARAM9(t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7, t8, p8, t9, p9) (p1,p2,p3,p4,p5,p6,p7,p8,p9) t1 p1; t2 p2; t3 p3; t4 p4; t5 p5; t6 p6; t7 p7; t8 p8; t9 p9;
#endif /* _ANSI_C */
#endif /* PROTO1 */

/* bzero() and bcopy() are BSD only! End run here */
#define bzero(b, length) memset((void *)(b), 0, (size_t)(length))
#define bcopy(from, to, length) memcpy((void *)(to), (void *)(from), (size_t)(length))

#endif /* DEFINE_H */
