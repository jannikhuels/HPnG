/*
 * ModelChecker.cpp
 *
 *  Created on: Mar 1, 2013
 *      Author: hamed
 */

#include "ModelChecker.h"
using namespace std;

IntervalSet* ModelChecker::visitStocRegion(Region* region, Formula* psi1, Formula* psi2, IntervalSet* potentialSatSet, double t, Interval bound) {
	//TODO:[IMPORTANT assumption] t=0 and bound = [0, MAX_TIME]

	if (region == NULL) return NULL;
	//if line t is within the region



	//if line t is NOT within the region
	Polygon *psi1_poly = geometryHelper->createPolygon(region, psi1);
	Polygon *psi2_poly = geometryHelper->createPolygon(region, psi2);
	psi1_poly = geometryHelper->cropPolygon(psi1_poly, bound.end, DOWN);
	psi2_poly = geometryHelper->cropPolygon(psi2_poly, bound.end, DOWN);



//	if (!debugImage.empty()){
////		region->print(cout, model);
////		if (psi1_poly != NULL) psi1_poly->print();
////		if (psi2_poly != NULL) psi2_poly->print();
//		if (psi1_poly != NULL) geometryHelper->drawPolygon(debugImage, psi1_poly, scale, cv::Scalar(0, 255, 0));
//		if (psi2_poly != NULL) geometryHelper->drawPolygon(debugImage, psi2_poly, scale, cv::Scalar(255, 0, 0));
//		cv::Mat flipped;
//		cv::flip(debugImage, flipped, 0);
//		cv::imshow("visitStocRegion: psi1 - psi2", flipped);
//		cv::waitKey(0);
//	}


	//potentialSatSet->print();

	// First interval set.
	if (psi2_poly != NULL){
		for (unsigned int i = 0; i < potentialSatSet->intervals.size(); i++){
			Point p1(potentialSatSet->intervals.at(i).start, region->lowerBoundry->getY(potentialSatSet->intervals.at(i).start));
			Point p2(potentialSatSet->intervals.at(i).end, region->lowerBoundry->getY(potentialSatSet->intervals.at(i).end));
			Segment* seg = new Segment(p1, p2);

			satSet = satSet->unionWith(geometryHelper->getIntersectionIntervals(psi2_poly, seg));
		}
	}

	if (psi1_poly != NULL){
		IntervalSet* psi1_Int = new IntervalSet();

		//Second interval set.
		for (unsigned int i = 0; i < potentialSatSet->intervals.size(); i++){
			Point p1(potentialSatSet->intervals.at(i).start, region->lowerBoundry->getY(potentialSatSet->intervals.at(i).start));
			Point p2(potentialSatSet->intervals.at(i).end, region->lowerBoundry->getY(potentialSatSet->intervals.at(i).end));
			Segment* seg = new Segment(p1, p2);

			psi1_Int = psi1_Int->unionWith(geometryHelper->getIntersectionIntervals(psi1_poly, seg));
			//In general case here the psi1_poly should be reformed.
		}

		if (psi2_poly != NULL){
			IntervalSet* psi12_Int = geometryHelper->getIntersectionIntervals(psi1_poly, psi2_poly);
			if (psi12_Int != NULL){
				satSet = satSet->unionWith(psi1_Int->intersect(psi12_Int));
			}
		}

		//Third Set of intervals.
		for (unsigned int i = 0; i < region->successors->size(); i++){
			IntervalSet* potPsi1_int = geometryHelper->getIntersectionIntervals(psi1_poly, region->successors->at(i)->lowerBoundry);
			if (potPsi1_int != NULL){
				potPsi1_int = potPsi1_int->minus(satSet);
				potPsi1_int = potPsi1_int->intersect(psi1_Int);
				if (!potPsi1_int->isEmpty())
					 satSet = satSet->unionWith(visitStocRegion(region->successors->at(i), psi1, psi2, potPsi1_int, t, bound));
			}
		}
	}

//	satSet->print(cout);
	return satSet;

}

ModelChecker::~ModelChecker() {
	// TODO Auto-generated destructor stub
}

IntervalSet* ModelChecker::visitDtrmRegion(DtrmEvent* dtrmRegion, Formula* psi1, Formula* psi2, double t, Interval bound) {
//	satSet->print();

	if (dtrmRegion->time >= bound.end || dtrmRegion->time >= model->MaxTime || dtrmRegion == NULL)
		return satSet;

	Direction d1 = UP;
	Direction d2 = UP;
	double t1 = INF;
	double t2 = INF;

	RegionState state1 = geometryHelper->getTimeAndDirection(dtrmRegion, psi1, t1, d1);
	RegionState state2 = geometryHelper->getTimeAndDirection(dtrmRegion, psi2, t2, d2);

	std::cout << "dtrmRegion: " << dtrmRegion->time << std::endl;
	std::cout << "psi1: " << t1 << " - state: " << state1 << std::endl;
	std::cout << "psi2: " << t2 << " - state: " << state2 << std::endl;

//	if (!debugImage.empty()){
//		if (state1 != NONE)
//			geometryHelper->drawVerticalLine(debugImage, t1, d1, scale, cv::Scalar(0, 0, 255));
//		if (state2 != NONE)
//			geometryHelper->drawVerticalLine(debugImage, t2, d2, scale, cv::Scalar(255, 0, 0));
//
//		cv::Mat flipped;
//		cv::flip(debugImage, flipped, 0);
//		cv::imshow("visitDtrmRegion: psi1 - psi2", flipped);
//		cv::waitKey(0);
//	}


	if ((state2 != NONE && d2 == DOWN) || state2 == ENTIRE){
		satSet = satSet->unionWith(Interval(dtrmRegion->time, model->MaxTime));
//		std::cout << "psi2 holds -------------------" << std::endl;
//		satSet->print(std::cout);
		return satSet;
	}

	if (state1 != ENTIRE && d1 == UP)
		return NULL;

	if (state2 != NONE && t2 <= t1)
		satSet = satSet->unionWith(Interval(t2, model->MaxTime));
	else if(state1 == ENTIRE && (state2 != NONE && d2 == UP))
		satSet = satSet->unionWith(Interval(t2, model->MaxTime));
	else if(state1 == ENTIRE)
		satSet = satSet->unionWith(visitDtrmRegion(dtrmRegion->nextDtrmEvent, psi1, psi2, t, bound));

//	satSet->print();

	IntervalSet* iSet = new IntervalSet();
	iSet->intervals.push_back(Interval(dtrmRegion->time, (state1 == ENTIRE)? dtrmRegion->nextDtrmEvent->time : t1));
	iSet = iSet->minus(satSet);

//	iSet->print();


	for (unsigned int i =0; i < dtrmRegion->nextRegions->size(); i++)
		satSet = satSet->unionWith(visitStocRegion(dtrmRegion->nextRegions->at(i), psi1, psi2, iSet, t, bound));

	return satSet;

}

IntervalSet* ModelChecker::until(Formula* psi1, Formula* psi2, Interval bound, double time) {
	// TODO: Currently I have assumed that the given bound starts at 0 ie. T_1 = 0.

//	if (satSet != NULL)
//		satSet->clear();
//	else
	satSet = new IntervalSet();

	double t = time; // - std->getTrEnabledTime();
	Line timeLine(0, t);
	Point p1, p2;

//	if (!debugImage.empty()){
//		geometryHelper->drawVerticalLine(debugImage, bound.end, UP, scale, cv::Scalar(0, 128, 128));
//		cv::Mat flipped;
//		cv::flip(debugImage, flipped, 0);
//		cv::imshow("bound", flipped);
//		cv::waitKey(0);
//
//	}


//	satSet->print();
	//iterating over all deterministic events.
	// the last event regards the maximum time reached. so should not be considered.
	for (unsigned int i = 0; i < std->dtrmEventList.size() - 1; i++){
//		std::cout << i << std::endl;
		if (std->dtrmEventList[i]->time <= time && std->dtrmEventList[i]->nextDtrmEvent->time > time){
			satSet = satSet->unionWith(visitDtrmRegion(std->dtrmEventList[i], psi1, psi2, t, bound));
//			satSet->print();
			break;
		}
	}

	//iterating over all stochastic regions.
	for (unsigned int i = 0; i < std->regionList.size(); i++){
		if (std->regionList[i]->intersect(timeLine, p1, p2)){
			IntervalSet* potSet = new IntervalSet();
			potSet->intervals.push_back(Interval(p1.X, p2.X));

			satSet = satSet->unionWith(visitStocRegion(std->regionList[i], psi1, psi2, potSet, t, bound));
		}
	}



	return satSet;
}

double ModelChecker::calcProb(IntervalSet* iSet, double(*cdf)(double), double shift){
	double prob = 0;
	for (unsigned int i = 0; i < iSet->intervals.size(); i++)
		prob += ( cdf(iSet->intervals[i].end - shift) - cdf(iSet->intervals[i].start - shift));
	return prob;

}

bool ModelChecker::propertyXleqCTest(Model* model, Marking* marking, double t0, double t1, double &s1, double &s2, int pIndex, double amount) {
    /**
     * fluid level in a place is (as+b) + (t1s+t0)d. t1s+t0 is time that has passed after entering this region.
     */
    double b = marking->fluid0[model->places[pIndex].idInMarking] + t0*marking->fluidPlaceDeriv[model->places[pIndex].idInMarking];
    double a = marking->fluid1[model->places[pIndex].idInMarking] + t1*marking->fluidPlaceDeriv[model->places[pIndex].idInMarking];
    if (IS_ZERO(a)) {
        if (b <= amount)
            return true;
        else
            return false;
    }
    double p = (amount-b)/a;
    if (p < s2 && p > s1){
        if (a < -ZERO_PREC) s1 = p;
        if (a > +ZERO_PREC) s2 = p;
        return true;
    } else if (p > s2){
        if (a * s2 + b < amount) return true;
    } else if (p < s1){
        if (a * s1 + b < amount) return true;
    }

    return false;
}

IntervalSet* ModelChecker::calcAtomDisISetAtTime(double time, int pIndex, double amount) {

    double s1, s2;
    IntervalSet *iSet = new IntervalSet();
    Point p1, p2;

    /* First iterate over all the regions in the stochastic area. */
    if (time > TimedDiagram::getInstance()->getTrEnabledTime()) {
        for (int i = 0; (unsigned)i < TimedDiagram::getInstance()->regionList.size(); i++) {
            double sFrameTime = time - TimedDiagram::getInstance()->getTrEnabledTime();

            Segment timeSeg(0, sFrameTime , 0, model->MaxTime - TimedDiagram::getInstance()->getTrEnabledTime());
            if (TimedDiagram::getInstance()->regionList[i]->intersect(timeSeg, p1, p2)) {
                s1 = p1.X; s2 = p2.X;
                if (TimedDiagram::getInstance()->regionList[i]->marking->tokens[model->places[pIndex].idInMarking] == amount) {
                    //std::cout << "Interval : ["<< s1 << "," << s2 << "]" << std::endl;
                    Interval I(s1, s2);
                    //std::cout << "Interval : [" << I.end << "," << I.start << "]" << std::endl;
                    iSet->intervals.push_back(I);
                    //iSet->print(std::cout);
                }
            }
        }
    }

    /* Second decide if the deterministic regions hold or not. */
    int cc;
    for (cc = 0; (unsigned)cc < TimedDiagram::getInstance()->dtrmEventList.size() && time > TimedDiagram::getInstance()->dtrmEventList[cc]->time ; cc++);
    s1 = time > TimedDiagram::getInstance()->getTrEnabledTime() ? time : 0;
    s2 = INFINITY;/*model->MaxTime;*/ // Rather infinity right?

    double t;
    if (cc == 0)
        t = time;
    else
        t = time - TimedDiagram::getInstance()->dtrmEventList[cc - 1]->time;
    if (TimedDiagram::getInstance()->dtrmEventList[cc]->preRegionMarking->tokens[model->places[pIndex].idInMarking] == amount) {
        Interval I(time > TimedDiagram::getInstance()->getTrEnabledTime() ? s1 - TimedDiagram::getInstance()->getTrEnabledTime() : s1, time > TimedDiagram::getInstance()->getTrEnabledTime() ? s2 - TimedDiagram::getInstance()->getTrEnabledTime() : s2);
        iSet->intervals.push_back(I);
    }
    return iSet;
}

IntervalSet* ModelChecker::calcAtomContISetAtTime(double time, int pIndex, double amount) {
    //TODO: by sorting region in more genuine way this function could be more effecient.

    double s1, s2;
    IntervalSet *iSet = new IntervalSet();
    Point p1, p2;

    if (time > TimedDiagram::getInstance()->getTrEnabledTime()) {
        //stochastic part : only if at the given time g-transition could have been enabled.
        /**
        * IMPORTANT note: we have computed everything for the stochastic part in the frame with origin at (gTrEnabledTime, gTrEnabledTime)
        * So we should move to that frame.
        */
        double sFrameTime = time - TimedDiagram::getInstance()->getTrEnabledTime();

        Segment timeSeg(0, sFrameTime , 0, model->MaxTime - TimedDiagram::getInstance()->getTrEnabledTime());
        for (int i = 0; (unsigned)i < TimedDiagram::getInstance()->regionList.size(); i++){
            //s1 and s2 are the validity interval for which the probability holds, returned by the function iPropHolds.
            if (TimedDiagram::getInstance()->regionList[i]->intersect(timeSeg, p1, p2)){
                s1 = p1.X; s2 = p2.X;
                double t0 = sFrameTime - TimedDiagram::getInstance()->regionList[i]->lowerBoundry->b /*- regionList[i]->timeBias*/;
                double t1 = - TimedDiagram::getInstance()->regionList[i]->lowerBoundry->a;
                if (this->propertyXleqCTest(model, TimedDiagram::getInstance()->regionList[i]->marking, t0, t1, s1, s2, pIndex, amount)){
                    Interval I(s1, s2);
                    iSet->intervals.push_back(I);
                }
            }
        }
    }

    //determinestic part after g-transition firing
    int cc;
    for (cc = 0; (unsigned)cc < TimedDiagram::getInstance()->dtrmEventList.size() && time > TimedDiagram::getInstance()->dtrmEventList[cc]->time ; cc++);
    s1 = time > TimedDiagram::getInstance()->getTrEnabledTime() ? time : 0;
    s2 = INFINITY;/*model->MaxTime;*/ // Rather infinity right?

    double t;
    if (cc == 0)
        t = time;
    else
        t = time - TimedDiagram::getInstance()->dtrmEventList[cc - 1]->time;
    if (this->propertyXleqCTest(model, TimedDiagram::getInstance()->dtrmEventList[cc]->preRegionMarking, t , 0, s1, s2, pIndex, amount)) {
        Interval I(time > TimedDiagram::getInstance()->getTrEnabledTime() ? s1 - TimedDiagram::getInstance()->getTrEnabledTime() : s1, time > TimedDiagram::getInstance()->getTrEnabledTime() ? s2 - TimedDiagram::getInstance()->getTrEnabledTime() : s2);
        iSet->intervals.push_back(I);
    }
    return iSet;
}

bool ModelChecker::iSetTT(IntervalSet *&res) {
    res = new IntervalSet();
    Interval I(0, INFINITY);
    res->intervals.push_back(I);
    //std::cout << "IntervalSet result TT : ";
    //res->print(std::cout);
    //std::cout << std::endl;
    return true;
}

bool ModelChecker::iSetNeg(IntervalSet *&res, IntervalSet* iset1) {
    //std::cout << "IntervalSet 1 : ";
    //iset1->print(std::cout);
    //std::cout << std::endl;
    if(!iSetTT(res)) return false;
    res = res->minus(iset1);
    //std::cout << "IntervalSet result NEG : ";
    //res->print(std::cout);
    //std::cout << std::endl;
    return true;
}

bool ModelChecker::iSetAnd(IntervalSet *&res, IntervalSet* iset1, IntervalSet* iset2) {
    //std::cout << "IntervalSet 1 : ";
    //iset1->print(std::cout);
    //std::cout << std::endl;
    //std::cout << "IntervalSet 2 : ";
    //iset2->print(std::cout);
    //std::cout << std::endl;
    res = iset1->intersect(iset2);
    //std::cout << "IntervalSet result AND : ";
    //res->print(std::cout);
    //std::cout << std::endl;
    return true;
}

bool ModelChecker::iSetAtomDis(IntervalSet *&res, AtomDisFormula* psi1) {
    res = this->calcAtomDisISetAtTime(this->ttc, psi1->getPlaceIndex(), psi1->getN());
    return true;
}

bool ModelChecker::iSetAtomCont(IntervalSet *&res, AtomContFormula* psi1) {
    res = this->calcAtomContISetAtTime(this->ttc, psi1->getPlaceIndex(), psi1->getC());
    return true;
}
