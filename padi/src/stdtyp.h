#ifndef STDTYP_H
#define STDTYP_H

/*
** Standard Types for portability
*/

/* boolean types */
typedef char tbool_t;        /* 8 bit */
 
#if !(defined(__rpc_types_h) /* sunos */ || defined(_RPC_TYPES_H) /* solaris */ || defined(_HPUX_SOURCE) || defined(__RPC_TYPES_H__) /* dec */ || defined(__RPC_TYPES_H) /* AIX */)
typedef int bool_t;          /* register size (function return value, etc.) */
#endif

/* signed integer types */
typedef char tiny_t;         /* 8 bit */
#define TINY(n) (char)(n)  /* !!! chars are signed !!! */
/* #define TINY(n) (char)(((n) & 0x80) ? (~0x7f | (n)) : (n))
                              !!! chars are unsigned */
                           /* use short for 16 bit */
                           /* use long for 32 bit */

/* unsigned integer types */
typedef char utiny_t;      /* 8 bit */  /* !!! no unsigned char !!! */
#define UTINY(n) (unsigned)((n) & 0xff) /* !!! no unsigned char !!! */
/* typedef unsigned char utiny_t;       8 bit   !!! unsigned char allowed !!! */
/* #define UTINY(n) (unsigned char)(n)     !!! unsigned char allowed !!! */



/* to hold size of anything */
/* typedef unsigned size_t;  defined in ANSI C now */

/* index type */
/* typedef size_t index_t;  defined as short in posix! */

/* string type (null terminated string) */
typedef char *string_t;

/* pointer to anything */
typedef char *data_tp;

/* type of largest alignment boundary */
typedef double align_t;

#endif /* STDTYP_H */
