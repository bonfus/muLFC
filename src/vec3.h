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
void vec3_free(vec3*);
int vec3_cpy(vec3 *, vec3 *);
int vec3_add(vec3 *, vec3 *);
int vec3_sub(vec3 *, vec3 *);
int vec3_mul(vec3 *, vec3 *);
int vec3_muls(scalar, vec3 *);
double vec3_norm(vec3 * v);
double vec3_dot(vec3 * v, vec3 * u);
int vec3_cross(vec3 * u, vec3 * v, vec3 * t);
int vec3_daxpy(scalar alpha, vec3 * , vec3*);
#endif
