#include <math.h>
#include "vec3.h"

#ifndef __GSL


vec3 * new_vec3(scalar x, scalar y, scalar z)
{
	vec3 * v;
        v = malloc(3*sizeof(scalar));
	v->x = x;
	v->y = y;
	v->z = z;
	return v;
}

vec3 * new_vec3_zero() 
{
	vec3 * v;
        v = malloc(3*sizeof(scalar));
	v->x = 0;
	v->y = 0;
	v->z = 0;
	return v;
}

void vec3_free(vec3 * v) { free(v); }

void vec3_set(vec3 * v, scalar x, scalar y, scalar z) 
{
	v->x = x;
	v->y = y;
	v->z = z;
}
void vec3_get(vec3 * v, scalar * x, scalar  * y, scalar * z) 
{
	*x = v->x;
	*y = v->y;
	*z = v->z;
}
void vec3_getp(vec3 * v, scalar * x) 
{
	*x = v->x;
	*(x+1) = v->y;
	*(x+2) = v->z;
}


int vec3_cpy(vec3 * v, vec3 * u) 
{
	v->x = u->x;
	v->y = u->y;
	v->z = u->z;
	return 0;
}

int vec3_add(vec3 * v, vec3 * u) 
{
	v->x = v->x + u->x;
	v->y = v->y + u->y;
	v->z = v->z + u->z;
	return 0;
}

int vec3_sub(vec3 * v, vec3 * u) 
{
	v->x = v->x - u->x;
	v->y = v->y - u->y;
	v->z = v->z - u->z;
	return 0;
}

int vec3_mul(vec3 * v, vec3 * u) 
{
	v->x = v->x * u->x;
	v->y = v->y * u->y;
	v->z = v->z * u->z;
	return 0;
}

int vec3_muls(scalar s, vec3 * v) 
{
	v->x = s * v->x;
	v->y = s * v->y;
	v->z = s * v->z;
	return 0;
}

int vec3_daxpy(scalar s, vec3 * v, vec3* r) 
{
	r->x = s * v->x + r->x;
	r->y = s * v->y + r->y;
	r->z = s * v->z + r->z;
	return 0;
}

double vec3_norm(vec3 * v)
{
    return sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
}

double vec3_dot(vec3 * v, vec3 * u)
{
    return v->x*u->x + v->y*u->y + v->z*u->z;
}

int vec3_cross(vec3 * u, vec3 * v, vec3 * t)
{
	  t->x = u->y * v->z - u->z * v->y;
	  t->y = u->z * v->x - u->x * v->z;
	  t->z = u->x * v->y - u->y * v->x;
	  
          return 0;
}


#else

vec3 * _vec3(scalar x, scalar y, scalar z)
{
	vec3 * v;
  v = gsl_vector_alloc(3);
	gsl_vector_set(v,0,x);
	gsl_vector_set(v,1,y);
	gsl_vector_set(v,2,z);
	return v;
}

vec3 * vec3_zero() 
{
	vec3 * v;
  v = gsl_vector_alloc(3);
  gsl_vector_set_zero(v);
	return v;
}

#define vec3_add gsl_vector_add
#define vec3_sub gsl_vector_sub
#define vec3_mul gsl_vector_mul
#define vec3_muls(X,Y) gsl_vector_scale(Y,X)
#define vec3_norm gsl_blas_dnrm2
#define vec3_daxpy gsl_blas_daxpy

scalar vec3_dot (vec3 * x, vec3 * y) {
        scalar r;
        gsl_blas_ddot(x,y,&r);
        return r;
}

// https://gist.github.com/jmbr/668083
int vec3_cross(vec3 *u, vec3 *v, vec3 *t)
{
        double p1 = gsl_vector_get(u, 1)*gsl_vector_get(v, 2)
                - gsl_vector_get(u, 2)*gsl_vector_get(v, 1);

        double p2 = gsl_vector_get(u, 2)*gsl_vector_get(v, 0)
                - gsl_vector_get(u, 0)*gsl_vector_get(v, 2);

        double p3 = gsl_vector_get(u, 0)*gsl_vector_get(v, 1)
                - gsl_vector_get(u, 1)*gsl_vector_get(v, 0);

        gsl_vector_set(t, 0, p1);
        gsl_vector_set(t, 1, p2);
        gsl_vector_set(t, 2, p3);
}

void vec3_get(vec3 * v, scalar * x, scalar  * y, scalar * z) 
{
	*x = gsl_vector_get(v,0);
	*y = gsl_vector_get(v,1);
	*z = gsl_vector_get(v,2);
}
void vec3_getp(vec3 * v, scalar * x) 
{
	*x = gsl_vector_get(v,0);
	*(x+1) = gsl_vector_get(v,1);
	*(x+2) = gsl_vector_get(v,2);
}

#endif
