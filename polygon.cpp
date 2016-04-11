/*
 * polygon.cpp
 *
 *  Created on: Mar 8, 2013
 *      Author: hamed
 */

#include "polygon.h"

#include <iostream>

Polygon::Polygon() {
	segments = new std::vector<Segment*>();
}

void Polygon::print() {
	for (unsigned int i = 0; i < segments->size(); i++){
		std::cout <<  "(" <<segments->at(i)->p1.X << ", " << segments->at(i)->p1.Y << ") " <<
				", (" <<segments->at(i)->p2.X << ", " << segments->at(i)->p2.Y << ")" << std::endl;
	}
}

Polygon::~Polygon() {
	// TODO Auto-generated destructor stub
}

bool Polygon::intersect(Line* l, Point& p1, Point& p2) {
	int intCnt = 0;
	Point p;

	for (int i = 0; i < this->segments->size() && intCnt < 2; i++){
		if (this->segments->at(i)->intersect(*l, p)){
			if (intCnt == 0) p1 = p;
			if (intCnt == 1) p2 = p;
			intCnt++;
		}
	}

	if (intCnt < 2) return false;
	else return true;
}


