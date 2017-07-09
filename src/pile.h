#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <string.h>

#include <math.h>
#include "vec3.h"

typedef struct {
	unsigned int nElements; /**< Number of elements in the pile. */
	double * ranks;  /**< Scalar value to weight the elements.  Must be positive */
	vec3 ** elements; /**< Pointer to the elements. */
} pile;

void pile_init(pile * p, unsigned int nElements);

void pile_add_element(pile * p, double rank, vec3 * v);

void pile_move_elements_from_position(pile * p, unsigned int pos);

void pile_free(pile * p);
