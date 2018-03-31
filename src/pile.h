#include <math.h>
#include "types.h"

typedef struct {
	unsigned int nElements; /**< Number of elements in the pile. */
	VecX ranks;  /**< Scalar value to weight the elements.  Must be positive */
	MatX elements; /**< Pointer to the elements. */
} pile;

void pile_init(pile& p, unsigned int nElements);

void pile_add_element(pile& p, double rank, const Vec3& v);

void pile_move_elements_from_position(pile& p, unsigned int pos);

void pile_free(pile& p);
