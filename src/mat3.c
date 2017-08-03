#include <math.h>
#include "mat3.h"

#ifndef __GSL
mat3 * new_mat3_zero()
{
    mat3 * m;
    m = malloc(sizeof(mat3));
    m->a = new_vec3_zero();
    m->b = new_vec3_zero();
    m->c = new_vec3_zero();

    return m;
}

mat3 * new_mat3_identity()
{
    mat3 * m;
    m = new_mat3_zero();

    m->a->x = 1;
    m->b->y = 1;
    m->c->z = 1;

    return m;
}

mat3 * new_mat3_diag(scalar a, scalar b, scalar c)
{
    mat3 * m;
    m = new_mat3_zero();
    m->a->x = a;
    m->b->y = b;
    m->c->z = c;

    return m;
}

mat3 * new_mat3(scalar ax, scalar ay, scalar az,
                scalar bx, scalar by, scalar bz,
                scalar cx, scalar cy, scalar cz)
{
    mat3 * m;
    m = new_mat3_zero();
    m->a->x = ax;
    m->a->y = ay;
    m->a->z = az;
    m->b->x = bx;
    m->b->y = by;
    m->b->z = bz;
    m->c->x = cx;
    m->c->y = cy;
    m->c->z = cz;
    return m;
}
void mat3_set(mat3* m, scalar ax, scalar ay, scalar az,
              scalar bx, scalar by, scalar bz,
              scalar cx, scalar cy, scalar cz)
{
    m->a->x = ax;
    m->a->y = ay;
    m->a->z = az;
    m->b->x = bx;
    m->b->y = by;
    m->b->z = bz;
    m->c->x = cx;
    m->c->y = cy;
    m->c->z = cz;
}
void mat3_get(mat3* m, scalar * ax, scalar * ay, scalar * az,
              scalar * bx, scalar * by, scalar * bz,
              scalar * cx, scalar * cy, scalar * cz)
{
    *ax=m->a->x;
    *ay=m->a->y;
    *az=m->a->z;
    *bx=m->b->x;
    *by=m->b->y;
    *bz=m->b->z;
    *cx=m->c->x;
    *cy=m->c->y;
    *cz=m->c->z;
}
void mat3_getp(mat3* m, scalar * p)
{
    *(p+0) = m->a->x; *(p+1) = m->a->y; *(p+2) = m->a->z;
    *(p+3) = m->b->x; *(p+4) = m->b->y; *(p+5) = m->b->z;
    *(p+6) = m->c->x; *(p+7) = m->c->y; *(p+8) = m->c->z;
}

void mat3_free(mat3* m) {
    vec3_free(m->a);
    vec3_free(m->b);
    vec3_free(m->c);
    free(m);
}

int mat3_add(mat3 * m, mat3 * n)
{

    m->a->x = m->a->x + n->a->x;
    m->a->y = m->a->y + n->a->y;
    m->a->z = m->a->z + n->a->z;
    m->b->x = m->b->x + n->b->x;
    m->b->y = m->b->y + n->b->y;
    m->b->z = m->b->z + n->b->z;
    m->c->x = m->c->x + n->c->x;
    m->c->y = m->c->y + n->c->y;
    m->c->z = m->c->z + n->c->z;

    return 0;
}

int mat3_mul(mat3 * m, mat3 * n, mat3 * r)
{

    scalar ax = m->a->x * n->a->x + m->a->y * n->b->x + m->a->z * n->c->x;
    scalar ay = m->a->x * n->a->y + m->a->y * n->b->y + m->a->z * n->c->y;
    scalar az = m->a->x * n->a->z + m->a->y * n->b->z + m->a->z * n->c->z;
    scalar bx = m->b->x * n->a->x + m->b->y * n->b->x + m->b->z * n->c->x;
    scalar by = m->b->x * n->a->y + m->b->y * n->b->y + m->b->z * n->c->y;
    scalar bz = m->b->x * n->a->z + m->b->y * n->b->z + m->b->z * n->c->z;
    scalar cx = m->c->x * n->a->x + m->c->y * n->b->x + m->c->z * n->c->x;
    scalar cy = m->c->x * n->a->y + m->c->y * n->b->y + m->c->z * n->c->y;
    scalar cz = m->c->x * n->a->z + m->c->y * n->b->z + m->c->z * n->c->z;
    mat3_set(r, ax,ay,az,bx,by,bz,cx,cy,cz);
    return 0;
}

int mat3_mulv(mat3 * m, vec3 * v, vec3 * r)
{
    vec3 u;
    u.x = m->a->x * v->x +        m->a->y * v->y + m->a->z * v->z;
    u.y = m->b->x * v->x +  m->b->y * v->y + m->b->z * v->z;
    u.z = m->c->x * v->x +  m->c->y * v->y + m->c->z * v->z;
    r->x = u.x; r->y = u.y; r->z = u.z;
    return 0;
}

int mat3_vmul(vec3 * v,mat3 * m, vec3 * r)
{
    vec3 u;
    u.x = m->a->x * v->x +  m->b->x * v->y + m->c->x * v->z;
    u.y = m->a->y * v->x +  m->b->y * v->y + m->c->y * v->z;
    u.z = m->a->z * v->x +  m->b->z * v->y + m->c->z * v->z;
    r->x = u.x; r->y = u.y; r->z = u.z;
    return 0;
}

/* Get the rotation matrix representation a rotation 'r' around axis 'v' */
/* make sure 'v' is a unit vector */
int mat3_aangle(vec3 * v, scalar r, mat3*m)
{
    scalar c, s, C;

    c = cos(r);
    s = -sin(r);
    C = 1. - c;

    m->a->x = 1.; m->a->y = 0.; m->a->z = 0.;
    m->b->y = 1.; m->b->x = 0.; m->b->z = 0.;
    m->c->z = 1.; m->c->x = 0; m->c->y = 0.;

    m->a->x = v->x * v->x * C + c;
    m->a->y = v->y * v->x * C + v->z * s;
    m->a->z = v->z * v->x * C - v->y * s;
    m->b->x = v->x * v->y * C - v->z * s;
    m->b->y = v->y * v->y * C + c;
    m->b->z = v->z * v->y * C + v->x * s;
    m->c->x = v->x * v->z * C + v->y * s;
    m->c->y = v->y * v->z * C - v->x * s;
    m->c->z = v->z * v->z * C + c;
    return 0;
}

int mat3_inv(mat3 * i, mat3 * minv)
{
    /* computes the inverse of a matrix m */
    scalar det = i->a->x * (i->b->y * i->c->z - i->c->y * i->b->z) -
                 i->a->y * (i->b->x * i->c->z - i->b->z * i->c->x) +
                 i->a->z * (i->b->x * i->c->y - i->b->y * i->c->x);

    scalar invdet = 1.0 / det;

    minv->a->x = (i->b->y * i->c->z - i->c->y * i->b->z) * invdet;
    minv->a->y = (i->a->z * i->c->y - i->a->y * i->c->z) * invdet;
    minv->a->z = (i->a->y * i->b->z - i->a->z * i->b->y) * invdet;
    minv->b->x = (i->b->z * i->c->x - i->b->x * i->c->z) * invdet;
    minv->b->y = (i->a->x * i->c->z - i->a->z * i->c->x) * invdet;
    minv->b->z = (i->b->x * i->a->z - i->a->x * i->b->z) * invdet;
    minv->c->x = (i->b->x * i->c->y - i->c->x * i->b->y) * invdet;
    minv->c->y = (i->c->x * i->a->y - i->a->x * i->c->y) * invdet;
    minv->c->z = (i->a->x * i->b->y - i->b->x * i->a->y) * invdet;

    return 0;
}

#else

mat3* new_mat3_zero(){
    mat3 * m = gsl_matrix_calloc(3, 3);
    gsl_matrix_set_zero (m);
    return m;
}
mat3 * new_mat3_diag(scalar a, scalar b, scalar c) {
    mat3 * m = gsl_matrix_calloc (3, 3);
    gsl_matrix_set_zero (m);
    gsl_matrix_set (m, 0, 0, a);
    gsl_matrix_set (m, 1, 1, b);
    gsl_matrix_set (m, 2, 2, c);
    return m;
}
mat3 * new_mat3_identity() {
    mat3 * m = gsl_matrix_calloc (3, 3);
    gsl_matrix_set_identity (m);
    return m;
}
mat3 * new_mat3(scalar ax, scalar ay, scalar az,
                scalar bx, scalar by, scalar bz,
                scalar cx, scalar cy, scalar cz) {
    mat3 * m = gsl_matrix_calloc (3, 3);
    gsl_matrix_set_zero (m);
    gsl_matrix_set (m, 0, 0, ax);
    gsl_matrix_set (m, 0, 1, ay);
    gsl_matrix_set (m, 0, 2, az);
    gsl_matrix_set (m, 1, 0, bx);
    gsl_matrix_set (m, 1, 1, by);
    gsl_matrix_set (m, 1, 2, bz);
    gsl_matrix_set (m, 2, 0, cx);
    gsl_matrix_set (m, 2, 1, cy);
    gsl_matrix_set (m, 2, 2, cz);
    return m;
}

void mat3_set(mat3* m, scalar ax, scalar ay, scalar az,
              scalar bx, scalar by, scalar bz,
              scalar cx, scalar cy, scalar cz)
{
    gsl_matrix_set (m, 0, 0, ax);
    gsl_matrix_set (m, 0, 1, ay);
    gsl_matrix_set (m, 0, 2, az);
    gsl_matrix_set (m, 1, 0, bx);
    gsl_matrix_set (m, 1, 1, by);
    gsl_matrix_set (m, 1, 2, bz);
    gsl_matrix_set (m, 2, 0, cx);
    gsl_matrix_set (m, 2, 1, cy);
    gsl_matrix_set (m, 2, 2, cz);
}

int mat3_aangle(vec3 * v, scalar r, mat3 * m)
{
    scalar c, s, C;

    c = cos(r);
    s = -sin(r);
    C = 1. - c;

    gsl_matrix_set_identity (m);

    gsl_matrix_set(m, 0, 0,  gsl_vector_get(v,0) * gsl_vector_get(v,0) * C + c);
    gsl_matrix_set(m, 0, 1,  gsl_vector_get(v,1) * gsl_vector_get(v,0) * C + gsl_vector_get(v,2) * s);
    gsl_matrix_set(m, 0, 2,  gsl_vector_get(v,2) * gsl_vector_get(v,0) * C - gsl_vector_get(v,1) * s);
    gsl_matrix_set(m, 1, 0,  gsl_vector_get(v,0) * gsl_vector_get(v,1) * C - gsl_vector_get(v,2) * s);
    gsl_matrix_set(m, 1, 1,  gsl_vector_get(v,1) * gsl_vector_get(v,1) * C + c);
    gsl_matrix_set(m, 1, 2,  gsl_vector_get(v,2) * gsl_vector_get(v,1) * C + gsl_vector_get(v,0) * s);
    gsl_matrix_set(m, 2, 0,  gsl_vector_get(v,0) * gsl_vector_get(v,2) * C + gsl_vector_get(v,1) * s);
    gsl_matrix_set(m, 2, 1,  gsl_vector_get(v,1) * gsl_vector_get(v,2) * C - gsl_vector_get(v,0) * s);
    gsl_matrix_set(m, 2, 2,  gsl_vector_get(v,2) * gsl_vector_get(v,2) * C + c);
    return 0;
}
int mat3_inv(mat3 * m, mat3 * minv) {
    gsl_permutation * p;
    int s;

    p = gsl_permutation_alloc(3);
    gsl_linalg_LU_decomp (m, p, &s);
    
    s = gsl_linalg_LU_invert (m, p, minv);
    gsl_permutation_free(p);
    return s;
}

void mat3_get(mat3* m, scalar * ax, scalar * ay, scalar * az,
              scalar * bx, scalar * by, scalar * bz,
              scalar * cx, scalar * cy, scalar * cz)
{
    *ax=gsl_matrix_get(m, 0, 0);
    *ay=gsl_matrix_get(m, 0, 1);
    *az=gsl_matrix_get(m, 0, 2);
    *bx=gsl_matrix_get(m, 1, 0);
    *by=gsl_matrix_get(m, 1, 1);
    *bz=gsl_matrix_get(m, 1, 2);
    *cx=gsl_matrix_get(m, 2, 0);
    *cy=gsl_matrix_get(m, 2, 1);
    *cz=gsl_matrix_get(m, 2, 2);
}
void mat3_getp(mat3* m, scalar * p)
{
    *(p+0) = gsl_matrix_get(m, 0, 0);
    *(p+1) = gsl_matrix_get(m, 0, 1);
    *(p+2) = gsl_matrix_get(m, 0, 2);
    *(p+3) = gsl_matrix_get(m, 1, 0);
    *(p+4) = gsl_matrix_get(m, 1, 1);
    *(p+5) = gsl_matrix_get(m, 1, 2);
    *(p+6) = gsl_matrix_get(m, 2, 0);
    *(p+7) = gsl_matrix_get(m, 2, 1);
    *(p+8) = gsl_matrix_get(m, 2, 2);
}

#endif
