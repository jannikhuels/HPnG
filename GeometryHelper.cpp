/*
 * GeometryHelper.cpp
 *
 *  Created on: Mar 6, 2013
 *      Author: hamed
 */

#include "GeometryHelper.h"

Polygon* GeometryHelper::reform(Polygon* poly, Interval interval) {
	return NULL;
}

RegionState GeometryHelper::getTimeAndDirection(DtrmEvent* dtrmRegion,
		Formula* formula, double& t, Direction& dir) {
	/*
	 * \f$ x \geq c => f_0 - \beta*d + d*t \geq c \f$ in which f_0 is the level of fluid place at the moment of entering this region.
	 * and \beta is the time of entering this region. Therefore:
	 * if d = 0:
	 * 		if f_0 \geq c => entire region
	 * 		o.w. Nothing
	 * 	if d > 0
	 * 		t > (c-f_0)/d + \beta
	 * 	if d < 0
	 * 		t < (c-f_0)/d + \beta
	 */

	/* Checking if we have combination of formulas. if the discrete part holds it holds in all the region. and if not in no where.
	 * IMPORTANT: we assume 'and' between discrete and continuous formulas.
	 */


	AtomContFormula* atomConProp = NULL;
	AtomDisFormula* disFormula = NULL;

	if (formula->getType() == DISCRITE) {
		disFormula = ((AtomDisFormula*) formula);
		if (dtrmRegion->postRegionMarking->tokens[model->places[disFormula->getPlaceIndex()].idInMarking] == disFormula->getN()){
//			std::cout << "in getTimeAndDirection ------- " << std::endl;
//			std::cout << "dtrmRegion->postRegionMarking->tokens[disFormula->getPlaceIndex()]: " << dtrmRegion->postRegionMarking->tokens[disFormula->getPlaceIndex()] << std::endl;
//			std::cout << "disFormula->getPlaceIndex(): " << disFormula->getPlaceIndex() << std::endl;
			return ENTIRE;
		}
	}else if (formula->getType() == COMBI) {
		disFormula = ((ADFormula*) formula)->disFormula;
		atomConProp = ((ADFormula*) formula)->contFormula;
	}else
		atomConProp = ((AtomContFormula*) formula);


	if (disFormula != NULL)
		if (dtrmRegion->postRegionMarking->tokens[model->places[disFormula->getPlaceIndex()].idInMarking] != disFormula->getN())
			return NONE;

	if (atomConProp != NULL){
		double f0 =
				dtrmRegion->postRegionMarking->fluid0[model->places[atomConProp->getPlaceIndex()].idInMarking];
		double d =
				dtrmRegion->postRegionMarking->fluidPlaceDeriv[model->places[atomConProp->getPlaceIndex()].idInMarking];

		if (d == 0) {
			if (atomConProp->getType() == CONT_POS) {
				if (f0 >= atomConProp->getC())
					return ENTIRE;
				else
					return NONE;
			} else if (atomConProp->getType() == CONT_NEG) {
				if (f0 < atomConProp->getC())
					return ENTIRE;
				else
					return NONE;
			}
		} else {
//			std::cout << "atomConProp->getC(): " << atomConProp->getC() << "f0: "
//					<< f0 << "d: " << d << "dtrmRegion->time: " << dtrmRegion->time
//					<< std::endl;

			if (atomConProp->getType() == CONT_POS)
				dir = (d > 0) ? UP : DOWN;
			else if (atomConProp->getType() == CONT_NEG)
				dir = (d > 0) ? DOWN : UP;

			t = (atomConProp->getC() - f0) / d + dtrmRegion->time;

		}

		if ((dir == UP && t > dtrmRegion->nextDtrmEvent->time)
				|| (dir == DOWN && t < dtrmRegion->time))
			return NONE;
		else if ((dir == UP && (t <= dtrmRegion->time))
				|| (dir == DOWN && (t >= dtrmRegion->nextDtrmEvent->time)))
			return ENTIRE;
		else
			return PART;
	}
    return NONE;
}

//void GeometryHelper::drawPolygon(cv::Mat& im, Polygon* poly, const int scale,
//		cv::Scalar color) {
//	for (unsigned int i = 0; i < poly->segments->size(); i++) {
//		cv::Point p1(poly->segments->at(i)->p1.X * scale,
//				poly->segments->at(i)->p1.Y * scale);
//		cv::Point p2(poly->segments->at(i)->p2.X * scale,
//				poly->segments->at(i)->p2.Y * scale);
//
//		cv::line(im, p1, p2, color, 2);
//	}
//
//}
//
//void GeometryHelper::drawVerticalLine(cv::Mat& im, double t, Direction dir,
//		const int scale, cv::Scalar color) {
//	cv::Point p1(t * scale, t * scale);
//	cv::Point p2(model->MaxTime * scale, t * scale);
//	cv::Point pCenter((p1.x + p2.x) / 2, t * scale);
//	cv::Point pUp((p1.x + p2.x) / 2, t * scale + 3);
//	cv::Point pDown((p1.x + p2.x) / 2, t * scale - 3);
//	cv::line(im, p1, p2, color, 2);
//	if (dir == UP)
//		cv::line(im, pCenter, pUp, color, 2);
//	else
//		cv::line(im, pCenter, pDown, color, 2);
//
//}

GeometryHelper::~GeometryHelper() {
	// TODO Auto-generated destructor stub
}

Polygon* GeometryHelper::createPolygon(Region* region, Formula* formula) {
	// TODO: Currently I have assumed that the given formula is an atomic continuous property.

	/* \f$ x = f_0 + f_1s \f$ when entering the region and \f$ t = \alpha s + \beta \f$ is the underlying segment,and \f$ d \f$ is the drift.
	 *
	 * if d != 0:
	 * 		the corresponding line for formula \f$ x \geq c\f$ is \f$ t = (\alpha -f_1/d)s + (\beta + c/d - f_0/d) \f$
	 * 		if \f$ d > 0 => t \geq as+b \f$ and if \f$ d < 0 => t \leq as+b \f$
	 *
	 * if d = 0 and \f$ f_1 \neq 0 \f$:
	 * 		if \f$ f_1 > 0 => s \geq (c-f_0)/f_1 \f$ and if \f$ f_1 <0 => s \leq (c-f_0)/f_1\f$
	 *
	 * 		if f_1 = 0:
	 * 			all the region  if f_0 \geq c
	 * 			nothing o.w.
	 *
	 */


	AtomContFormula* atomConProp = NULL;
	AtomDisFormula* disFormula = NULL;

	if (formula->getType() == DISCRITE) {
		disFormula = ((AtomDisFormula*) formula);
		if (region->marking->tokens[model->places[disFormula->getPlaceIndex()].idInMarking] == disFormula->getN())
			return regionToPoly(region);
	}

	if (formula->getType() == COMBI) {
		disFormula = ((ADFormula*) formula)->disFormula;
		atomConProp = ((ADFormula*) formula)->contFormula;
	}else
		atomConProp = ((AtomContFormula*) formula);


	if (disFormula != NULL)
		if (region->marking->tokens[model->places[disFormula->getPlaceIndex()].idInMarking] != disFormula->getN())
			return NULL;


	if (atomConProp != NULL){
		double f0 =
				region->marking->fluid0[model->places[atomConProp->getPlaceIndex()].idInMarking];
		double f1 =
				region->marking->fluid1[model->places[atomConProp->getPlaceIndex()].idInMarking];
		double d =
				region->marking->fluidPlaceDeriv[model->places[atomConProp->getPlaceIndex()].idInMarking];

//		std::cout << "f0: " << f0 << " f1: " << f1 << " d: " << d << std::endl;
		Polygon* poly;
		Line propLine;

		if (IS_ZERO(d)) {
			if (!IS_ZERO(f1)) {
				propLine.X = (atomConProp->getC() - f0) / f1;
				Direction upOrDown;
				if (atomConProp->getType() == CONT_POS)
					upOrDown = (f1 > 0) ? UP : DOWN;
				else if (atomConProp->getType() == CONT_NEG)
					upOrDown = (f1 < 0) ? UP : DOWN;

				poly = splitRegion(region, &propLine, upOrDown);
			} else if (f0 > atomConProp->getC())
				poly = (atomConProp->getType() == CONT_POS) ?
						regionToPoly(region) : NULL;
			else if (f0 < atomConProp->getC())
				poly = (atomConProp->getType() == CONT_NEG) ?
						regionToPoly(region) : NULL;
			else
				poly = regionToPoly(region);

		} else {
			double a = region->lowerBoundry->a - f1 / d;
			double b = region->lowerBoundry->b + (atomConProp->getC() - f0) / d;
			propLine.a = a;
			propLine.b = b;

			Direction upOrDown;
			if (atomConProp->getType() == CONT_POS)
				upOrDown = (d > 0) ? UP : DOWN;
			else if (atomConProp->getType() == CONT_NEG)
				upOrDown = (d < 0) ? UP : DOWN;

			poly = splitRegion(region, &propLine, upOrDown);
		}

		return poly;
	}
    return NULL;
}

Polygon* GeometryHelper::splitRegion(Region* region, Line* line,
		Direction direction) {
	Polygon* poly = new Polygon();

//	region->print(std::cout);

	for (unsigned int i = 0; i < region->eventSegments->size(); i++)
		addToPolygon(poly, region->eventSegments->at(i)->timeSegment, line,
				direction);

	addToPolygon(poly, region->leftBoundry, line, direction);
	addToPolygon(poly, region->rightBoundry, line, direction);
	addToPolygon(poly, region->lowerBoundry, line, direction);

	Point p1, p2;
	if (region->intersect(*line, p1, p2))
		poly->segments->push_back(new Segment(p1, p2));

	return poly;
}

Polygon* GeometryHelper::regionToPoly(Region* region) {
	Polygon* poly = new Polygon();
	for (unsigned int i = 0; i < region->eventSegments->size(); i++) {
		if (region->eventSegments->at(i)->timeSegment->p1
				!= region->eventSegments->at(i)->timeSegment->p2)
			poly->segments->push_back(
					region->eventSegments->at(i)->timeSegment);
	}

	if (region->leftBoundry->p1 != region->leftBoundry->p2)
		poly->segments->push_back(region->leftBoundry);
	if (region->rightBoundry->p1 != region->rightBoundry->p2)
		poly->segments->push_back(region->rightBoundry);
	if (region->lowerBoundry->p1 != region->lowerBoundry->p2)
		poly->segments->push_back(region->lowerBoundry);

	return poly;

}

void GeometryHelper::addToPolygon(Polygon* poly, Segment* seg, Line* line,
		Direction direction) {
	Point intersection;
	if (seg->intersect(*line, intersection)) {
		if ((direction == UP && line->isUp(seg->p1))
				|| (direction == DOWN && !line->isUp(seg->p1))) {
			if (intersection != seg->p1)
				poly->segments->push_back(new Segment(seg->p1, intersection));
		} else if (intersection != seg->p2)
			poly->segments->push_back(new Segment(seg->p2, intersection));
	} else if ((direction == UP && line->isUp(seg->p1) && line->isUp(seg->p2))
			|| (direction == DOWN && !line->isUp(seg->p1)
					&& !line->isUp(seg->p2))) {
		poly->segments->push_back(seg);
	}
}

Polygon* GeometryHelper::cropPolygon(Polygon* poly, double t, Direction dir) {
	if (poly == NULL || poly->segments->size() == 0)
		return NULL;

	Polygon* newPoly = new Polygon();
	Line* line = new Line(0, t);
	Point p1, p2;
	if (poly->intersect(line, p1, p2)) {
		for (unsigned int i = 0; i < poly->segments->size(); i++)
			addToPolygon(newPoly, poly->segments->at(i), line, dir);
		newPoly->segments->push_back(new Segment(p1, p2));
	} else {
		if ((dir == DOWN && !line->isUp(poly->segments->at(0)->p1))
				|| (dir == UP && line->isUp(poly->segments->at(0)->p1))) {
			for (unsigned int i = 0; i < poly->segments->size(); i++)
				newPoly->segments->push_back(poly->segments->at(i));
		} else
			newPoly = NULL;

	}
	return newPoly;

}

bool GeometryHelper::isPolyContainPt(Polygon* poly, Point* p) {
	for (unsigned int i = 0; i < poly->segments->size(); i++) {
		if (poly->segments->at(i)->p1 == *p || poly->segments->at(i)->p2 == *p)
			return true;
	}
	return false;
}

IntervalSet* GeometryHelper::getIntersectionIntervals(Polygon* poly,
		Segment* segment) {
	Point intPt1, intPt2;
	int numIntersection = 0;
//	std::cout << "---------------------------" << std::endl;
//	poly->print();
	for (unsigned int i = 0; i < poly->segments->size(); i++) {
		Point intPt;
		IntersectionStatus inStatus = poly->segments->at(i)->intersect2(
				*segment, intPt);
		if (inStatus == ALL) {
			Interval i1(segment->p1.X, segment->p2.X);
			Interval i2(poly->segments->at(i)->p1.X,
					poly->segments->at(i)->p2.X);
			Interval res;
			i1.intersect(i2, res);
			IntervalSet* ret = new IntervalSet();
			ret->intervals.push_back(res);
			return ret;
		} else if (inStatus == NO)
			continue;
		if (numIntersection == 0) {
			intPt1 = intPt;
			numIntersection++;
		} else if (numIntersection == 1) {
			intPt2 = intPt;
			numIntersection++;
		} else
			break;
	}
	if (numIntersection < 2)
		return NULL; //No intersection
	else {
		Interval I(intPt1.X, intPt2.X);
		IntervalSet* ret = new IntervalSet();
		ret->intervals.push_back(I);
		return ret;
	}
}

IntervalSet* GeometryHelper::getIntersectionIntervals(Polygon* poly1,
		Polygon* poly2) {
	Point intPt1, intPt2;
	int numIntersection = 0;
	IntervalSet* ret = new IntervalSet();
	for (unsigned int i = 0; i < poly1->segments->size(); i++) {
		for (unsigned int j = 0; j < poly2->segments->size(); j++) {
			Point intPt;
			IntersectionStatus inStatus = poly1->segments->at(i)->intersect2(
					*poly2->segments->at(j), intPt);
			if (inStatus == ALL) {
				Interval i1(poly1->segments->at(i)->p1.X,
						poly1->segments->at(i)->p2.X);
				Interval i2(poly2->segments->at(j)->p1.X,
						poly2->segments->at(j)->p2.X);
				Interval res;
				i1.intersect(i2, res);
				ret->intervals.push_back(res);
				break;
			}
			if (inStatus == POINT && numIntersection == 0) {
				intPt1 = intPt;
				numIntersection++;
			} else if (inStatus == POINT && numIntersection == 1) {
				intPt2 = intPt;
				numIntersection++;
			}
		}
		if (numIntersection == 2) {
			Interval I(intPt1.X, intPt2.X);
			ret->intervals.push_back(I);
		}
	}
//	std::cout << "ret: ";
//	ret->print();
	ret->removeRedundency();
	return ret;
}

