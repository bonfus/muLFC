#ifndef MAT3_H
#define MAT3_H

#include "vec3.h"

#ifndef __GSL
typedef struct {
	vec3 *a, *b, *c;
} mat3;
#else
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
typedef gsl_matrix mat3;
#endif



mat3* new_mat3_zero(void);
mat3* new_mat3_identity(void);
mat3* new_mat3_diag(scalar a, scalar b, scalar c);
mat3* new_mat3(scalar ax, scalar ay, scalar az, 
                scalar bx, scalar by, scalar bz,
                scalar cx, scalar cy, scalar cz);

void mat3_set(mat3*, scalar ax, scalar ay, scalar az, 
                scalar bx, scalar by, scalar bz,
                scalar cx, scalar cy, scalar cz);
void mat3_get(mat3*, scalar * ax, scalar * ay, scalar * az, 
                scalar * bx, scalar * by, scalar * bz,
                scalar * cx, scalar * cy, scalar * cz);
void mat3_getp(mat3*, scalar *);

int mat3_aangle(vec3 * v, scalar r, mat3*m);
int mat3_inv(mat3 * i, mat3 * minv);

#ifndef __GSL
void mat3_free(mat3* m);
int mat3_mul(mat3 *, mat3*, mat3*); 
int mat3_add(mat3*,  mat3*);
int mat3_mulv( mat3*,  vec3*, vec3*);
int mat3_vmul( vec3*,  mat3*, vec3*);
int mat3_cpy(mat3*,  mat3*);

#else


#define mat3_free gsl_matrix_free
#define mat3_add gsl_matrix_add
#define mat3_mul(A,B,C) gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, A, B, 0.0, C)
#define mat3_mulv(A,x,y) gsl_blas_dgemv (CblasNoTrans, 1.0, (A), (x), 0.0, (y))
#define mat3_vmul(x,A,y) gsl_blas_dgemv (CblasTrans, 1.0, (A), (x), 0.0, (y))
#define mat3_cpy gsl_matrix_memcpy

#endif

#endif
