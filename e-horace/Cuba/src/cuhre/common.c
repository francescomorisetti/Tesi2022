/*
	common.c
		includes most of the modules
		this file is part of Cuhre
		last modified 2 Aug 13 11 th
*/


#include "ChiSquare.c"
#include "Rule.c"

static inline bool BadDimension(cThis *t)
{
  if( t->ndim > MAXDIM ) return true;
  /* begin HORACE */
  /* allowing also 1-dim integration! */
  /* original line: return t->ndim < 2;*/
  return t->ndim < 1;
  /* end HORACE */
}

static inline bool BadComponent(cThis *t)
{
  if( t->ncomp > MAXCOMP ) return true;
  return t->ncomp < 1;
}

