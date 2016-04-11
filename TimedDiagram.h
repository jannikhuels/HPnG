/*
 * TimedDiagram.h
 *
 *  Created on: 28 okt. 2011
 *      Author: GhasemiehH
 */

#ifndef TIMEDDIAGRAM_H_
#define TIMEDDIAGRAM_H_

#include <vector>
#include <math.h>
#include <algorithm>

//#include <opencv2/opencv.hpp>

#include "Event.h"
#include "DFPN2.h"
#include "Region.h"
#include "Line.h"






class TimedDiagram {
public:
	static TimedDiagram* getInstance();

	void setModel(Model* model);

	void clear();

	void generateDiagram(Marking* initialMarking);

	/**
	 * @param potentialEvents List of possible next events.
	 * @param uSegment underlying segment.
	 * @param start Beginning of working interval.
	 * @param end End of working interval.
	 * @param nextEvents List computed next events, depending on the value of \f$ s\f$. \\
	 * Have in mind this list is sorted by the value of \f$ s\f$. This means that the first element of the list is the one happening by lower value of \f$ s\f$
	 *
	 * \brief This function takes the minimum of input event lines set, in the t-s frame. It returns the sorted resulting segments.
	 * TODO: consider multiple events at the same time.
	 */
	void computeNextEvents(std::vector<StochasticEvent*> * potentialEvents, Segment *uSegment, std::vector<StochasticEvent*> * nextEvents);

	void SegmentizeDtrmRegion(Marking* marking);

	/**
	 * @param marking Current marking of the system.
	 * @param eventLine This is the equation of the lower boundary of region we wish to segmentize. We expect this event to contain its parent region.\\
	 * Have in mind that this is a segment so it has the validity interval in itself.
	 *
	 * Each (potential) region to be segmentized at first consists of starting and ending for 's' (like the previous intervals),
	 * and an event line, which specifies the equation for \f$t\f$ (the time after which the event has occurred) depending on \f$s\f$.
	 * This function recursively segmentize such regions, by finding linear equations of all the possible next events.
	 * Note 1 : It is assumed that this function is called only for the upper part of st-diagram (when the general transition is fired.)
	 * Note 2 : Remember that initial marking for this function to start should be such that general transition is fired.
	 */
	void segmentizeStochasticRegion(Marking* marking, StochasticEvent *eventLine,  double timeBias);

	/**
	 * @param time The time for which probability calculation is being done.
	 * @param sPdfInt Pointer to the integral pdf of s.
	 * @param isPropHolds Pointer to the function which determines whether a spesific property for a marking holds. This 
	 *		  function returns the interval that the property holds.
	 * @return The probabilty that the given property holds for a given time.
	 */
	double calProbAtTime(double time, double (*sPdfInt)(double), bool (*isPropHolds)(Model*, Marking*, double t0, double t1, double&, double&));

//	void saveDiagram(std::string filename);

	double getTrEnabledTime() const {
		return gTrEnabledTime;
	}

	int getNumberOfRegions(){return regionList.size() + dtrmEventList.size();};

	Model* model;

	/**
	 * This vector contains list of region after firing of general transition, i.e. at the top of ts-line.
	 */
	std::vector<Region*> regionList;

	/**
	 * Before firing of general transition all events could be determined deterministicly. So they can be maintained as a
	 * sequence points.
	 */
	std::vector<DtrmEvent*> dtrmEventList;

//	cv::Mat debugImage;
	int scale;

private:
	static TimedDiagram* instance;
	unsigned int currentTime;

	/**
	 * keeps the first time g-transition could fire.
	 */
	double gTrEnabledTime;


	void minLines(std::vector<StochasticEvent*> * potentialEvents, Segment* uSegment, std::vector<StochasticEvent*> * nextEvents);
	void createAddRegions(std::vector<StochasticEvent*> * eventList, StochasticEvent * preEvent, Marking* marking);

	void computeNextEventCreatePoly(std::vector<StochasticEvent*> * potentialEvents, StochasticEvent * preEvent, Marking* marking, std::vector<StochasticEvent*> *nextEvents);

	/**
	 * Auxilary drawing functions
	 */

//	void drawSegmet(cv::Mat & image, Segment& seg, const int scale);

	TimedDiagram();

};

#endif /* TIMEDDIAGRAM_H_ */
