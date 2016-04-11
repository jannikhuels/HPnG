#ifndef EVENT_H_
#define EVENT_H_

#include "DFPN2.h"
#include "Line.h"
#include "Region.h"

enum EventType { TRANSITION, PLACE_LOWER_BOUNDRY, PLACE_UPPER_BOUNDRY, GAURD_ARC, MAX_TIME_REACHED, FIRST_NULL_EVENT};



class Region;
class DtrmEvent;

/**
 * Stochastic Event. sd
 */
struct StochasticEvent{
	/**
	 * Describes the equation and the validity s-interval of next potential event, dependent on \f$ s \f$.
	 */
	Segment* timeSegment;

	EventType eventType;

	/**
	 * Transition or place ID;
	 */
	int id;

	/**
	 * marking of region before this event happening.
	 */
	Marking* preRegionMarking;

	/*
	 * If the pre-event for this stochastic event were in stochastic area. this variable has a value.
	 */
	Region* preRegion;

	/*
	 * If the pre-event for this stochastic event were in deterministic area. this variable has a value.
	 */
	DtrmEvent* preDtrmEvent;

	/* TODO[IMPORTANT]: Only one of the above variables is allowed to have value. I know this is NASTY.
	 * the solution is to have a parent region class and then stochastic and deterministic regions as its child class.
	 * and then STD will have the list of region class. this is nicer for model checking. But so far I keep as it is for time being.
	 */

	/**
	 * marking of region before this event happening.
	 */
	Marking* postRegionMarking;

	static bool greaterSlopeFirst(StochasticEvent* e1, StochasticEvent* e2){return e1->timeSegment->a > e2->timeSegment->a;};

	StochasticEvent(Segment *line, EventType type) : timeSegment(line), eventType(type), preRegion(0), preDtrmEvent(0){};
	StochasticEvent(StochasticEvent* se){
		timeSegment = se->timeSegment;
		eventType = se->eventType;
		id = se->id;
		preRegionMarking = se->preRegionMarking;
		postRegionMarking = se->postRegionMarking;
		preRegion = se->preRegion;
		preDtrmEvent = se->preDtrmEvent;
	};

	StochasticEvent(){};

};

/**
 * Deterministic Event.
 */
struct DtrmEvent{
	EventType eventType;

	/**
	 * Transition or place ID;
	 */
	int id;

	/**
	 * Occurrence time.
	 */
	double time;

	/**
	 * marking of region before this event happening.
	 */
	Marking* preRegionMarking;
	DtrmEvent* nextDtrmEvent;

	/*
	 * Set of next stochastic regions which could be entered after this deterministic event.
	 */
	std::vector<Region*>* nextRegions;

	/**
	 * marking of region before this event happening.
	 */
	Marking* postRegionMarking;

	DtrmEvent(EventType type) : eventType(type){nextRegions = new std::vector<Region*>();};
};

//bool lineGreaterSlopeFirst(Line & l1, Line & l2);
#endif
