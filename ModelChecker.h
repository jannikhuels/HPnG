/*
 * ModelChecker.h
 *
 *  Created on: Mar 1, 2013
 *      Author: hamed
 */

#include <vector>

#include "Region.h"
#include "TimedDiagram.h"
#include "GeometryHelper.h"

#ifndef MODELCHECKER_H_
#define MODELCHECKER_H_





class ModelChecker {

private:
	Model* model;
	TimedDiagram* std;
	IntervalSet* satSet;
	GeometryHelper* geometryHelper;

	 

	    /*
	         * Helper functions as in DFPN2, although this one is more generalised independant.
	    */
	    bool propertyXleqCTest(Model* model, Marking* marking, double t0, double t1, double &s1, double &s2, int pIndex, double amount);


	IntervalSet* visitStocRegion(Region* region, Formula* psi1, Formula* psi2, IntervalSet* potentialSatSet, double t, Interval bound);
	IntervalSet* visitDtrmRegion(DtrmEvent* dtrmRegion, Formula* psi1, Formula* psi2, double t, Interval bound);

public:
//	cv::Mat debugImage;
	int scale;
	double ttc;

	ModelChecker(Model* model, TimedDiagram* std) {this->model = model; this->std = std; geometryHelper = new GeometryHelper(model);};
	virtual ~ModelChecker();

	IntervalSet* until(Formula* psi1, Formula* psi2, Interval bound, double time);

	double calcProb(IntervalSet* iSet, double(*cdf)(double), double shift);

	/*
	     * Functions that will calculate the intervalsets by traversal.
	     */
	    bool iSetTT(IntervalSet *&res);
	    bool iSetNeg(IntervalSet *&res, IntervalSet* iset1);
	    bool iSetAnd(IntervalSet *&res, IntervalSet* iset1, IntervalSet* iset2);
	    bool iSetAtomDis(IntervalSet *&res, AtomDisFormula* psi1);
	    bool iSetAtomCont(IntervalSet *&res, AtomContFormula* psi1);
    
    /**
     * @brief calcAtomDisISetAtAtTime calculates the intersection at a specific time to check.
     * @param time The time for which probability calculation is being done.
     * @param pIndex The index number of the discrete place.
     * @param amount The amount of to compare with.
     * @return The intervalset that the given property holds for a given time to check.
     */
    IntervalSet* calcAtomDisISetAtTime(double time, int pIndex, double amount);
    
    /**
     * @brief calcAtomContISetAtAtTime calculates the intersection at a specific time to check.
     * @param time The time for which probability calculation is being done.
     * @param pIndex The index number of the continuous place.
     * @param amount The amount of to compare with.
     * @return The intervalset that the given property holds for a given time to check.
     */
    IntervalSet* calcAtomContISetAtTime(double time, int pIndex, double amount);

};

#endif /* MODELCHECKER_H_ */
