/*
 * polygon.h
 *
 *  Created on: Mar 8, 2013
 *      Author: hamed
 */

#ifndef POLYGON_H_
#define POLYGON_H_


#include "Line.h"
#include <vector>
/* A temporary class just for the case study. it is more complex
 * Assumptions:
 * 	1) A polygon is consisted of one property line and the region's boundaries.
 * 	2) Polygons are convex.
 */
class Polygon {

public:
	std::vector<Segment*>* segments;

	void print();
	Polygon();
	virtual ~Polygon();

	bool intersect(Line* l, Point &p1, Point &p2);
};

#endif /* POLYGON_H_ */
