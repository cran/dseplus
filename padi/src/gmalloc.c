#include <malloc.h>
#include "define.h"

export size_t ngmalloced = 0;

#define STAMP "GMALLOC!"
#define STAMP_SIZE 8

export char *gmalloc PARAM1(size_t, n)
{
  char *p;

#ifdef SPLUS
  if (!(p = (char *)S_alloc((long)(n + STAMP_SIZE),(int)1)))
#else
  if (!(p = malloc(n + STAMP_SIZE)))
#endif
    return NULL;

  ngmalloced++;
  strcpy(p, STAMP);
  return p + STAMP_SIZE;
}

export char *grealloc PARAM2(char *, p, size_t, n)
{
  char *q = p - STAMP_SIZE;
#ifdef SPLUS
	printf("this will cause Splus to core dump\n");
#endif

  if (!p || strncmp(STAMP, q, STAMP_SIZE) || !ngmalloced) 
    abort();

  if (!(p = realloc(q,n + STAMP_SIZE)))
    return NULL;

  strcpy(p, STAMP);
  return p + STAMP_SIZE;
}

export char *gfree PARAM1(char *, p)
{
  char *q = p - STAMP_SIZE;

#ifdef SPLUS
  return p;
#else
  if (!p || strncmp(STAMP, q, STAMP_SIZE) || !ngmalloced) 
    abort();

  free(q);
  ngmalloced--;
  return (char *) NULL;
#endif

}

