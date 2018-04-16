/**
 * @file   pile.c
 * @author Pietro Bonfa'
 * @date   2016
 * @brief  Pile of vec3
 *
 *
 * 
 *
 * This file defines a set of functions to deal with a pile struct.
 * A pile is a list of (fixed) n elements with positive definite ranks.
 * Only the n elements with the highest rank r are kept.
 * The addition of an element with rank r' can result
 * in either the removal of another element with lower rank or in an
 * unaltered pile if all the elements in the pile have ranks higher than
 * the one that is under consideration.
 */
#include <iostream>
#include "pile.h"


/**
 * This function initializes the pile introducing nElements zero vectors
 * and setting the rank to -1. Since ranks must be positive, all these 
 * elements will be replaced.
 * 
 */
pile::pile(unsigned int nEls)
{
	nElements = nEls;
	ranks.resize(nElements);
	elements.resize(3,nElements);
  
  
  ranks.setOnes();
  ranks *= -1.0;
  elements.setZero();
}

/**
 * This function resets nElements to zero vectors and their rank to -1.
 * Since ranks must be positive, all these elements will be replaced.
 * 
 */
void pile::reset()
{
  ranks.setOnes();
  ranks *= -1.0;
  elements.setZero();
}

/**
 * This function checks if the element v with rank rank should be added
 * to the pile p.
 * If yes shifts all elements with lower rank by one position, removes
 * the last one, and adds the vector.
 * 
 */
void pile::add_element(double rank, const Vec3& v)
{
	unsigned int i;
	for ( i = 0; i < nElements; i++)
	{
		if (ranks(i) < 0.0) {
			ranks(i) = rank;
			elements.col(i) = v;
			break;
		}
		if (ranks(i) > rank) {
			move_elements_from_position(i);
			ranks(i) = rank;
			elements.col(i) = v;
			break;
		}
	}
}

/**
 * This function moves the element in position pos down by one position.
 * The implicit assumption is that the element at position pos will be
 * overwritten.
 * If nElements = 1 nothing should be done. Otherwise copy all the 
 * elements form nElements-2-i to nElements-1-i
 */
void pile::move_elements_from_position(unsigned int pos)
{
	unsigned int i;
	/* the first -1 is for 0 indexing */
	/* the second -1 is becouse if only the last element must be moved it is thrashed! */
	if (nElements < 2) {
		return;
	}
	for (i = (nElements-1) ; i-- > pos ;)
	{
		ranks(i+1) = ranks(i);
		elements.col(i+1) = elements.col(i);
	}
}


/**
 * 
 * Cleanup allocated memory 
 * 
 */
pile::~pile()
{

}
