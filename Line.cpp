/*
 * Line.cpp
 *
 *  Created on: 3 nov. 2011
 *      Author: GhasemiehH
 */

#include "Line.h"

#include <iostream>

Line::Line() {
	// TODO Auto-generated constructor stub

}

Line::~Line() {

	// TODO Auto-generated destructor stub
}

Line & Line::operator =(const Line &l){
	a = l.a; b = l.b;
	return *this;
};

Line Line::operator -(const Line &l){
	return Line(a - l.a, b - l.b);
}

//bool Line::operator <(Line &l){
//	return this->a < l.a;
//}
//
//bool Line::operator >(Line &l){
//	return this->a > l.a;
//}

double Line::getY(double x){
	if (a == INF)
		return -1;

	return a*x + b;
}


double Line::getX(double y){
	if (a == INF)
		return X;
	if (a == 0)
		return -1;
	else
		return (y - b)/a;

}

bool Line::intersect(Line &l, Point& intPoint){
	if (this->a == l.a) return false;

	intPoint.X = -(b - l.b)/(a - l.a);
	intPoint.Y = a*intPoint.X + b;

	return true;
}

bool Line::isUp(Point& p){
	if (this->a == INF)
		return (p.X <= X);

	return  ( p.Y - this->getY(p.X) > ZERO_PREC);
}
	

Segment::Segment(Point &p1, Point &p2) {
	this->p1 = p1;
	this->p2 = p2;
	if (IS_ZERO(p1.X - p2.X)){
		a = INF;
		b = INF;
		X = p1.X;
	} else{
		this->a = (p1.Y - p2.Y)/(p1.X - p2.X);
		this->b = p1.Y - a*p1.X;
	}
}

Segment::Segment(double a, double b, Point &p1, Point &p2): Line(a, b){
	this->p1 = p1;
	this->p2 = p2;
}

Segment::Segment(Line &line, double start, double end): Line (line.a, line.b) {
	p1.X = start; p2.Y = line.getY(start);
	p2.X = end; p2.Y = line.getY(end);
}

Segment::Segment(double a, double b, double start, double end): Line (a, b) {
	p1.X = start; p1.Y = getY(start);
	p2.X = end; p2.Y = getY(end);
}

Segment& Segment::operator =(const Segment &s){
	p1 = s.p1;
	p2 = s.p2;
	a = s.a;
	b = s.b;

	return *this;
}

bool Segment::isContainPoint(Point& p){
	if (this->a == INF)
		if (IS_ZERO(this->X - p.X))
			if (p.Y >= std::min(this->p1.Y, this->p2.Y) && p.Y <= std::max(this->p1.Y, this->p2.Y))
				return true;
	if (IS_ZERO(this->getY(p.X) - p.Y))
		if (p.X >= std::min(this->p1.X, this->p2.X) && p.X <= std::max(this->p1.X, this->p2.X))
			return true;

	return false;
}

bool Segment::intersect(Segment &s, Point &intPoint){
	if (this->a == s.a){
		//TODO: [IMPORTANT] here we may have an infinite intersection, if the two segments are lying on each others.
//		if (IS_ZERO(this->getY(s.p1.X) - s.p1.Y)){
//			intP
//			return true;
//		}
//		else
			return false;
	}

	double x ,y;
	if (this->a == INF){
		x = this->p1.X;
		y = s.getY(x);
	}else if (s.a == INF){
		x = s.p1.X;
		y = this->getY(x);
	}else{
		x = -(b - s.b)/(a - s.a);
		y = s.getY(x);
	}

	//if (x < p1.X || x > p2.X || x < s.p1.X || x > s.p2.X ||
	//	y < p1.Y || y > p2.Y || y < s.p1.Y || y > s.p2.Y)
	//	return false;

	Point p (x, y);
	if (this->isContainPoint(p) && s.isContainPoint(p)){
		intPoint = p;
		return true;
	}

/******************************************************************************************
//	if (s.p1 == this->p1 || s.p1 == this->p2 || s.p2 == this->p1 || s.p2 == this->p2){
//		intPoint.X = x;
//		intPoint.Y = y;
//		return true;
//	}
//
//	if (((s.isUp(this->p1) && !s.isUp(this->p2)) || (!s.isUp(this->p1) && s.isUp(this->p2))) &&
//		((this->isUp(s.p1) && !this->isUp(s.p2)) || (!this->isUp(s.p1) && this->isUp(s.p2))) ){
//		intPoint.X = x;
//		intPoint.Y = y;
//		return true;
//	}
********************************************************************************************/
	return false;
}

IntersectionStatus Segment::intersect2(Segment& s, Point& intPoint) {
	if (IS_ZERO(this->a - s.a)){
		if (IS_ZERO(this->getY(s.p1.X) - s.p1.Y)){
			return ALL;
		}
		else
			return NO;
	}

	double x ,y;
	if (this->a == INF){
		x = this->p1.X;
		y = s.getY(x);
	}else if (s.a == INF){
		x = s.p1.X;
		y = this->getY(x);
	}else{
		x = -(b - s.b)/(a - s.a);
		y = s.getY(x);
	}

	//if (x < p1.X || x > p2.X || x < s.p1.X || x > s.p2.X ||
	//	y < p1.Y || y > p2.Y || y < s.p1.Y || y > s.p2.Y)
	//	return false;

	if (s.p1 == this->p1 || s.p1 == this->p2 || s.p2 == this->p1 || s.p2 == this->p2){
		intPoint.X = x;
		intPoint.Y = y;
		return POINT;
	}

//	if (IS_ZERO(y - this->getY(x)) || IS_ZERO(y - s.getY(x))){
//		intPoint.X = x;
//		intPoint.Y = y;
//		return POINT;
//	}

	if (((s.isUp(this->p1) && !s.isUp(this->p2)) || (!s.isUp(this->p1) && s.isUp(this->p2))) &&
		((this->isUp(s.p1) && !this->isUp(s.p2)) || (!this->isUp(s.p1) && this->isUp(s.p2))) ){
		intPoint.X = x;
		intPoint.Y = y;
		return POINT;
	}



	return NO;

}

bool Segment::intersect(Line&l, Point &intPoint){
	if (IS_ZERO(this->a - l.a))
		return false;

	double x ,y;
	if (this->a == INF){
		x = this->p1.X;
		y = l.getY(x);
	}else if (l.a == INF){
		x = l.X;
		y = this->getY(x);
	}else{
		x = -(b - l.b)/(a - l.a);
		y = l.getY(x);
	}

	//if (x < p1.X || x > p2.X || x < s.p1.X || x > s.p2.X ||
	//	y < p1.Y || y > p2.Y || y < s.p1.Y || y > s.p2.Y)
	//	return false;
	if ((l.isUp(this->p1) && !l.isUp(this->p2)) || (!l.isUp(this->p1) && l.isUp(this->p2))){
		intPoint.X = x;
		intPoint.Y = y;
		return true;
	}

	return false;
}


bool Segment::isUp(Point& p){

	if (p.Y - this->getY(p.X) == 0) return false;

	if (this->a == INF)
		return (p.X <= this->p1.X);

	return  ( p.Y - this->getY(p.X) > 0);
}

void Segment::print() {
	std::cout << "a:"<< a << "b:" << b <<  " p1:(" << p1.X << "," << p1.Y << ")" << " p2:(" << p2.X << "," << p2.Y << ")" << std::endl;
}



