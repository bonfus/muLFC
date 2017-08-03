#include <stdlib.h>

#ifndef VEC3_H
#define VEC3_H

typedef double scalar;
#ifndef __GSL
typedef struct {
	scalar x, y, z;
} vec3;
#else
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

typedef gsl_vector vec3;
#endif

vec3 * new_vec3(scalar, scalar, scalar);
vec3 * new_vec3_zero(void);
void vec3_set(vec3*, scalar, scalar, scalar);
void vec3_getp(vec3*, scalar*);
void vec3_get(vec3*, scalar*, scalar*, scalar*);
int vec3_cpy(vec3 *, vec3 *);
double vec3_dot(vec3 * v, vec3 * u);
int vec3_cross(vec3 * u, vec3 * v, vec3 * t);

#ifndef __GSL
int vec3_add(vec3 *, vec3 *);
int vec3_sub(vec3 *, vec3 *);
int vec3_mul(vec3 *, vec3 *);
int vec3_muls(scalar, vec3 *);
double vec3_norm(vec3 * v);
int vec3_daxpy(scalar alpha, vec3 * , vec3*);
void vec3_free(vec3*);

#else

#define vec3_add gsl_vector_add
#define vec3_sub gsl_vector_sub
#define vec3_mul gsl_vector_mul
#define vec3_muls(X,Y) gsl_vector_scale(Y,X)
#define vec3_norm gsl_blas_dnrm2
#define vec3_daxpy gsl_blas_daxpy
#define vec3_free gsl_vector_free
#define vec3_cpy gsl_vector_memcpy
#endif


#endif
