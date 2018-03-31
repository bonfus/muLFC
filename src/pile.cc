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

#include "pile.h"


/**
 * This function initializes the pile introducing nElements zero vectors
 * and setting the rank to -1. Since ranks must be positive, all these 
 * elements will be replaced.
 * 
 */
void pile_init(pile& p, unsigned int nElements)
{
	unsigned int i;
	p.nElements = nElements;
	p.ranks.resize(nElements);
	p.elements.resize(3,nElements);
  
  
  p.ranks.setOnes();
  p.ranks *= -1.0;
  p.elements.setZero();

}

/**
 * This function resets nElements to zero vectors and their rank to -1.
 * Since ranks must be positive, all these elements will be replaced.
 * 
 */
void pile_reset(pile& p, unsigned int nElements)
{
	unsigned int i;
  
  p.ranks.setOnes();
  p.ranks *= -1.0;
  p.elements.setZero();
}

/**
 * This function checks if the element v with rank rank should be added
 * to the pile p.
 * If yes shifts all elements with lower rank by one position, removes
 * the last one, and adds the vector.
 * 
 */
void pile_add_element(pile& p, double rank, const Vec3& v)
{
	unsigned int i;
	for ( i = 0; i < p.nElements; i++)
	{
		if (p.ranks(i) < 0.0) {
			p.ranks(i) = rank;
			p.elements.col(i) = v;
			break;
		}
		if (p.ranks(i) > rank) {
			pile_move_elements_from_position(p, i);
			p.ranks(i) = rank;
			p.elements.col(i) = v;
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
void pile_move_elements_from_position(pile& p, unsigned int pos)
{
	unsigned int i;
	/* the first -1 is for 0 indexing */
	/* the second -1 is becouse if only the last element must be moved it is thrashed! */
	if (p.nElements < 2) {
		return;
	}
	for (i = (p.nElements-1) ; i-- > pos ;)
	{
		p.ranks(i+1) = p.ranks(i);
		p.elements.col(i+1) = p.elements.col(i);
	}
}


/**
 * 
 * Cleanup allocated memory 
 * 
 */
void pile_free(pile& p)
{

}
