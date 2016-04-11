/*
 * GeometryHelper.h
 *
 *  Created on: Mar 6, 2013
 *      Author: hamed
 */

#ifndef GEOMETRYHELPER_H_
#define GEOMETRYHELPER_H_

//#include "./clipper/clipper.hpp"
#include "Formula.h"
#include "IntervalSet.h"
#include "Region.h"
#include "polygon.h"

//#include <opencv2/opencv.hpp>

enum Direction {UP, DOWN};
enum RegionState {ENTIRE, PART, NONE};


/* IMPORTANT: This class is customized for a simple case!
 * Refer to the header file polygon.h .
 */
class GeometryHelper {
	Model* model;

public:
	//Stochastic part functions
	Polygon* createPolygon(Region* region, Formula* formula);
	Polygon* splitRegion(Region* region, Line* line, Direction direction);
	Polygon* regionToPoly(Region* region);
	/*
	 * Helper function. if the line has an intersection with seg, according to the direction intersection point
	 * and appropriate end point of seg is added to polygon.
	 */
	void addToPolygon(Polygon* poly, Segment* seg, Line* line, Direction direction);
	bool isPolyContainPt(Polygon* poly, Point* p);
	IntervalSet* getIntersectionIntervals(Polygon* poly, Segment* segment);
	IntervalSet* getIntersectionIntervals(Polygon* poly1, Polygon* poly2);
	Polygon* reform(Polygon* poly, Interval interval);

	//Discrete part functions.
	/*
	 * return true if formula holds in entire region. o.w. t and dir represent the validity area.
	 */
	RegionState getTimeAndDirection(DtrmEvent* dtrmRegion, Formula* psi1, double& t, Direction& dir);


	Polygon* cropPolygon(Polygon* poly, double t, Direction dir);

//	void drawPolygon(cv::Mat& im, Polygon* poly, const int scale, cv::Scalar color = cv::Scalar(0, 0, 255));
//	void drawVerticalLine(cv::Mat& im, double t, Direction dir,const int scale, cv::Scalar color =  cv::Scalar(0, 0, 255));

	GeometryHelper(Model* _model): model(_model){};
	virtual ~GeometryHelper();
};

#endif /* GEOMETRYHELPER_H_ */
