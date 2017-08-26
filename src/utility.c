#include <math.h>
#include "utility.h"

void print_vec3(const char * msg, const vec3 * i)
{
    printf(msg);
#ifdef __GSL
    printf(" %e %e %e\n", gsl_vector_get (i, 0),
                          gsl_vector_get (i, 1),
                          gsl_vector_get (i, 2));
#else
    printf(" %e %e %e\n", i->x, i->y, i->z);
#endif
}

void print_mat3(const char * msg, const mat3 * i)
{
    printf(msg);
#ifdef __GSL
    printf(" %e %e %e\n", gsl_matrix_get (i,0,0), gsl_matrix_get (i,0,1), gsl_matrix_get (i,0,2));
    printf(" %e %e %e\n", gsl_matrix_get (i,1,0), gsl_matrix_get (i,1,1), gsl_matrix_get (i,1,2));
    printf(" %e %e %e\n", gsl_matrix_get (i,2,0), gsl_matrix_get (i,2,1), gsl_matrix_get (i,2,2));
#else
    printf(" %e %e %e\n", i->a->x, i->a->y, i->a->z);
    printf(" %e %e %e\n", i->b->x, i->b->y, i->b->z);
    printf(" %e %e %e\n", i->c->x, i->c->y, i->c->z);
#endif
}
