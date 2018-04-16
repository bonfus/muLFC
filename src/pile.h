#include <math.h>
#include "types.h"


class pile {
  public: 
    pile(unsigned int nEls);
    ~pile();
    void add_element(double rank, const Vec3& v);
    void reset();
	  VecX ranks;  /**< Scalar value to weight the elements.  Must be positive */
	  MatX elements; /**< Pointer to the elements. */

  private:
    void move_elements_from_position(unsigned int pos);
  	unsigned int nElements; /**< Number of elements in the pile. */
};
