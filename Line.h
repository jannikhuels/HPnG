/*
 * Line.h
 *
 *  Created on: 3 nov. 2011
 *      Author: GhasemiehH
 */

#ifndef LINE_H_
#define LINE_H_

#define ZERO_PREC  1.0e-10
//inline bool IS_ZERO(double x) {
//	return x < ZERO_PREC && x > -ZERO_PREC);
//}
#include "DFPN2.h"

#define INF +1.0e100

enum IntersectionStatus {NO, ALL, POINT};

struct Point{
	double X, Y;

	Point(){X = Y = 0;};
	Point(double x, double y) : X(x), Y(y){};
	Point & operator =(const Point &p){X = p.X; Y = p.Y; return *this;};
	bool operator !=(const Point &p){return (!IS_ZERO(X - p.X) || !IS_ZERO(Y - p.Y));};
	bool operator ==(const Point &p){return IS_ZERO(X - p.X) && IS_ZERO(Y - p.Y);};
};

/**
 * Represents typical line equation: \f$ y = ax + b\f$ (our case: \f$t = as + b\f$)
 *
 */
class Line {
public:
	double a;
	double b;
	double X; //in case a = INF, the equation is x=X.

	Line();
	Line(double _a, double _b): a(_a), b(_b) {};
	Line(double _X): X(_X), a(INF){};
	virtual ~Line();

	Line & operator =(const Line &l);
	Line operator -(const Line &l);


	/**
	 * for sorting lines by their slop.
	 */
//	bool operator< (const Line &l);
//	bool operator> (const Line &l);

	double getY(double x);
	double getX(double y);
	bool isUp(Point& p);


	//TODO: This Function is should be modified for the case of a = INF.
	virtual bool intersect(Line &l, Point& intPoint);

};


class Segment : public Line {
public:
	Point p1, p2;

	Segment(Point &p1, Point &p2);

	Segment(double a, double b, Point &p1, Point &p2);
	Segment(double a, double b, double start, double end);
	Segment(Line &line, double start, double end);

	Segment& operator =(const Segment &s);

	bool isUp(Point& p);

	bool intersect(Segment &s, Point &intPoint);
	/*
	 * intPoint has a value if tow segments has an intersection on only one point.
	 */
	IntersectionStatus intersect2(Segment &s, Point &intPoint);
	bool intersect(Line &l, Point &intPoint);

	bool isContainPoint(Point& p);

	void print();
};



#endif /* LINE_H_ */
