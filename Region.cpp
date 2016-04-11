/*
 * Region.cpp
 *
 *  Created on: 28 okt. 2011
 *      Author: GhasemiehH
 */

#include <iostream>

#include "Region.h"

Region::Region(std::vector<StochasticEvent*> * eventList, Segment * lowerBoundray) {
	successors = new std::vector<Region*>();
	eventSegments = eventList;
	this->lowerBoundry = lowerBoundray;
	for (int i = 0; i < eventSegments->size(); i++){
		if (IS_ZERO(eventSegments->at(i)->timeSegment->p1.X - lowerBoundray->p1.X))
			this->leftBoundry = new Segment(lowerBoundray->p1, eventSegments->at(i)->timeSegment->p1);
		if (IS_ZERO(eventSegments->at(i)->timeSegment->p2.X - lowerBoundray->p2.X))
			this->rightBoundry = new Segment(lowerBoundray->p2, eventSegments->at(i)->timeSegment->p2);
	}


}

void Region::print(std::ostream &out, Model* model ) {

	out << "Lower Boundry: "; 
	out <<  "(" <<lowerBoundry->p1.X << ", " << lowerBoundry->p1.Y << ")";
	out << "--";
	out << "(" << lowerBoundry->p2.X << ", " << lowerBoundry->p2.Y << ")";
	out << std::endl;

	std:: cout << "*****************************************" << std::endl;
	printFluidLevels(model, this->marking);
	std:: cout << "*****************************************" << std::endl;
	printFluidRates(model, this->marking);
	std:: cout << "*****************************************" << std::endl;

	char* eventTypeMap [] = {"TRANSITION", "PLACE_LOWER_BOUNDRY", "PLACE_UPPER_BOUNDRY", "GAURD_ARC", "MAX_TIME_REACHED", "FIRST_NULL_EVENT"};

	for (int i = 0; i < eventSegments->size(); i++) {
		char* id;
		if (model != 0){
			if (eventSegments->at(i)->eventType == TRANSITION)
				id = model->transitions[eventSegments->at(i)->id].id;
			else if (eventSegments->at(i)->eventType == PLACE_LOWER_BOUNDRY || eventSegments->at(i)->eventType == PLACE_UPPER_BOUNDRY)
				id = model->places[eventSegments->at(i)->id].id;
			else if (eventSegments->at(i)->eventType == GAURD_ARC)
				id = model->arcs[eventSegments->at(i)->id].id;
			else
				id = "max or first";
		}else
			id = "not known";

		out <<  "(" <<eventSegments->at(i)->timeSegment->p1.X << ", "
				<< eventSegments->at(i)->timeSegment->p1.Y << ")";
		out << "--";
		out << "(" << eventSegments->at(i)->timeSegment->p2.X << ", " << eventSegments->at(i)->timeSegment->p2.Y << ")";
		out << "-->event: ["<< eventTypeMap[eventSegments->at(i)->eventType] << "]-id: [" << id << "]";
		out << std::endl;
	}

	out << "Left Boundry: ";
	out <<  "(" <<leftBoundry->p1.X << ", " << leftBoundry->p1.Y << ")";
	out << "--";
	out << "(" << leftBoundry->p2.X << ", " << leftBoundry->p2.Y << ")";
	out << std::endl;


	out << "Right Boundry: ";
	out <<  "(" <<rightBoundry->p1.X << ", " << rightBoundry->p1.Y << ")";
	out << "--";
	out << "(" << rightBoundry->p2.X << ", " << rightBoundry->p2.Y << ")";
	out << std::endl;
}

bool Region::intersect(Segment& s, Point &p1, Point &p2){
	// Our regions are convex so they will have at most 2 intersection points.
	int intCnt = 0;
	Point p;

	for (int i = 0; i < eventSegments->size() && intCnt < 2; i++){
		if (eventSegments->at(i)->timeSegment->intersect(s, p)){
			if (intCnt == 0) p1 = p;
			if (intCnt == 1) p2 = p;
			intCnt++;
		}
	}
	
	if (leftBoundry->intersect(s, p)){
		if (intCnt == 0) p1 = p;
		if (intCnt == 1) p2 = p;
		intCnt++;
	}

	if (rightBoundry->intersect(s, p)){
		if (intCnt == 0) p1 = p;
		if (intCnt == 1) p2 = p;
		intCnt++;
	}

	if (lowerBoundry->intersect(s, p)){
		if (intCnt == 0) p1 = p;
		if (intCnt == 1) p2 = p;
		intCnt++;
	}

	if (intCnt < 2) return false;
	else return true;
}

bool Region::intersect(Line& s, Point &p1, Point &p2){
	// Our regions are convex so they will have at most 2 intersection points.
	int intCnt = 0;
	Point p;

	for (int i = 0; i < eventSegments->size() && intCnt < 2; i++){
		if (eventSegments->at(i)->timeSegment->intersect(s, p)){
			if (intCnt == 0) p1 = p;
			if (intCnt == 1) p2 = p;
			intCnt++;
		}
	}

	if (leftBoundry->intersect(s, p)){
		if (intCnt == 0) p1 = p;
		if (intCnt == 1) p2 = p;
		intCnt++;
	}

	if (rightBoundry->intersect(s, p)){
		if (intCnt == 0) p1 = p;
		if (intCnt == 1) p2 = p;
		intCnt++;
	}

	if (lowerBoundry->intersect(s, p)){
		if (intCnt == 0) p1 = p;
		if (intCnt == 1) p2 = p;
		intCnt++;
	}

	if (intCnt < 2) return false;
	else return true;
}

