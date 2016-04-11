/*
* TimedDiagram.cpp
*
*  Created on: 28 okt. 2011
*      Author: GhasemiehH
*/

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>

#include "TimedDiagram.h"

TimedDiagram* TimedDiagram::instance = NULL;

TimedDiagram::TimedDiagram() {
	currentTime = 0;
}

void TimedDiagram::clear(){
	dtrmEventList.clear();
	regionList.clear();
	currentTime = 0;
	instance == NULL;
}

void TimedDiagram::setModel(Model* model){
	clear();
	this->model = model;
}

TimedDiagram* TimedDiagram::getInstance() {
	if (instance == NULL)
		instance = new TimedDiagram();

	return instance;
}

void TimedDiagram::generateDiagram(Marking* initialMarking) {

	/*
	* TODO:
	* Generating the deterministic (under the \f$ t = s \f$ line) part of diagram.
	*/

	SegmentizeDtrmRegion(initialMarking);
//	for (DtrmEvent* e = dtrmEventList[0]; e != 0; e = e->nextDtrmEvent){
//		std::cout << "------------" << e->time << "--------------"<< std::endl;
//		printFluidLevels(model, e->postRegionMarking);
//	}


	/**
	* Generating the stochastic (top of \f$ t = s \f$ line) part of diagram.
	*/
	//TODO: consider the marking after the last dtrEvent.
	double prePoint = 0;
	gTrEnabledTime = 0;


//	std::cout << "Dtrm events: " << std::endl;
//	for (int j = 0; j < dtrmEventList.size() - 1; j++){
//		char* id;
//		if (dtrmEventList[j]->eventType == TRANSITION)
//			id = model->transitions[dtrmEventList[j]->id].id;
//		else if (dtrmEventList[j]->eventType == PLACE_LOWER_BOUNDRY || dtrmEventList[j]->eventType == PLACE_UPPER_BOUNDRY)
//			id = model->places[dtrmEventList[j]->id].id;
//		else if (dtrmEventList[j]->eventType == GAURD_ARC)
//			id = model->arcs[dtrmEventList[j]->id].id;
//		else
//			id = "max or first";
//
//		char* eventTypeMap [] = {"TRANSITION", "PLACE_LOWER_BOUNDRY", "PLACE_UPPER_BOUNDRY", "GAURD_ARC", "MAX_TIME_REACHED", "FIRST_NULL_EVENT"};
//
//		std::cout << "time: " << dtrmEventList[j]->time << "-- " <<  eventTypeMap[dtrmEventList[j]->eventType] << "-- name: " << id << std::endl;
//
//		printFluidRates(model, dtrmEventList[j]->postRegionMarking);
//		std::cout << std::endl;
//		std::cout << "---------------------------------------------------" << std::endl;
//
//	}
	// ignore the the last event.
	for (int j = 0; j < dtrmEventList.size() - 1; j++){

		if (!isGTransitionEnabled(model, dtrmEventList[j]->postRegionMarking)){
			prePoint = dtrmEventList[j]->nextDtrmEvent->time;
			gTrEnabledTime = dtrmEventList[j]->nextDtrmEvent->time;
			continue;
		}


		//in case that two consequent event are happening in the same time approximately.
		if (j < dtrmEventList.size() - 1 && IS_ZERO(dtrmEventList[j]->time - dtrmEventList[j+1]->time))
			continue;

		//moving to the frame with origin at (gTrEnabledTime, gTrEnabledTime)
		Point startPoint(prePoint /*- gTrEnabledTime*/, prePoint /*- gTrEnabledTime*/);	//******//
		Point endPoint(dtrmEventList[j]->nextDtrmEvent->time /*- gTrEnabledTime*/, dtrmEventList[j]->nextDtrmEvent->time /*- gTrEnabledTime*/);
		Marking* gtFiredMarking = copyMarking(model, dtrmEventList[j]->postRegionMarking);

		checkEnabled(model, gtFiredMarking);
		setActFluidRate(model, gtFiredMarking, 0);
		/**
		 * The reason for (- startPoint.Y) is that we should shift the time to the frame with origin at (startPoint.Y, startPoint.Y).
		 * (This is how the algorithm is designed, there are other ways too!!)
		 */
		gtFiredMarking = advanceMarking(model, gtFiredMarking, 1, -startPoint.Y, gtFiredMarking->fluidPlaceDeriv,gtFiredMarking->enabling);

		//TODO [IMPORTANT]: check whether the g-transition is enabled.
		fireGeneralTransition(model, gtFiredMarking);

		Segment* tsLine = new Segment(startPoint, endPoint);
		StochasticEvent* initEvent = new StochasticEvent(tsLine, TRANSITION);
		initEvent->id = gTransitionId(model);
		initEvent->preRegionMarking = dtrmEventList[j]->postRegionMarking;
		initEvent->preDtrmEvent = dtrmEventList[j];

//		std::cout << gtFiredMarking->fluid0[0]<< "**"<< gtFiredMarking->enabling[0] << "**" << gtFiredMarking->enabling[1]<< "**" << gtFiredMarking->fluidPlaceDeriv[0] <<  std::endl;


		segmentizeStochasticRegion(gtFiredMarking, initEvent, prePoint);

		prePoint = dtrmEventList[j]->nextDtrmEvent->time;
	}

	//std::ofstream regionfile("regions.out");

	//for (int j = 0; j < dtrmEventList.size(); j++){
	//	regionfile << "dtrm region(" << j << "): " << std::endl;
	//	regionfile << "time: " << dtrmEventList[j]->time << "-type: " << dtrmEventList[j]->eventType << "-id: " << dtrmEventList[j]->id << std::endl;
	//}

	//regionfile << "-------------------------------------------------" << std::endl;

	//for (int i = 0; i < regionList.size(); i++){
	//	regionfile << "stoch region(" << i << "): " << std::endl;
	//	regionList[i]->print(regionfile);


	//	regionfile << "-------------------------------------------------" << std::endl;
	//}


}

void TimedDiagram::SegmentizeDtrmRegion(Marking* initialMarking){

	//TODO: (IMPORTANT) This function does not take care of multiple events at the same time.
	// Do it later!!!!

	Marking *marking = copyMarking(model, initialMarking);

	double crntTime = 0;
	int *enabled;
	double firstT;

	int firstEventID = -1;
	EventType firstEventType = FIRST_NULL_EVENT;

	double *clock0, *fluid0, *drift;
	DtrmEvent* crntEvent = 0;
	DtrmEvent* preEvent = 0;

//	checkEnabled(model, marking);
//	setActFluidRate(model, marking, crntTime);

//	crntEvent = new DtrmEvent(FIRST_NULL_EVENT);
//	crntEvent->time = 0;
//	crntEvent->id = -1;

	Marking* preMarking = copyMarking(model, marking);

	while (crntTime < model->MaxTime){



		checkEnabled(model, marking);
		setActFluidRate(model, marking, crntTime);
		
		crntEvent= new DtrmEvent(firstEventType);
		crntEvent->time = crntTime;
		crntEvent->id = firstEventID;
		crntEvent->preRegionMarking = preMarking;
		crntEvent->postRegionMarking = marking;
		crntEvent->eventType = firstEventType;
		dtrmEventList.push_back(crntEvent);
		if (preEvent != NULL)
			preEvent->nextDtrmEvent = crntEvent;
		preEvent = crntEvent;


		preMarking = copyMarking(model, marking);

//		char* id;
//		if (crntEvent->eventType == TRANSITION)
//			id = model->transitions[crntEvent->id].id;
//		else if (crntEvent->eventType == PLACE_LOWER_BOUNDRY || crntEvent->eventType == PLACE_UPPER_BOUNDRY)
//			id = model->places[crntEvent->id].id;
//		else if (crntEvent->eventType == GAURD_ARC)
//			id = model->arcs[crntEvent->id].id;
//		else
//			id = "max or first";
//
//		char* eventTypeMap [] = {"TRANSITION", "PLACE_LOWER_BOUNDRY", "PLACE_UPPER_BOUNDRY", "GAURD_ARC", "MAX_TIME_REACHED", "FIRST_NULL_EVENT"};
//
//		std::cout << "***********************************************************************"<< std::endl;
//		std::cout << "time: " << crntEvent->time << "-- " <<  eventTypeMap[crntEvent->eventType] << "-- name: " << id << std::endl;
//		for (int t = 0; t < model->N_transitions; t++){
//			if (model->transitions[t].type == TT_FLUID || model->transitions[t].type == TT_FLUID_DYNA){
//				std::cout << model->transitions[t].id << "("<< marking->enabling[t] << "): " <<
//						marking->actFluidRate[model->transitions[t].idInMarking] << std::endl;
//			}
//		}



		enabled = marking->enabling;
		clock0 = marking->clock0;
		fluid0 = marking->fluid0;
		drift = marking->fluidPlaceDeriv;

		firstT = INF;

		int hasImmediateEnabled;
		std::vector<int> imTransVec;

		// count the number (and put in cache) the number of enabled transitions
		hasImmediateEnabled = 0;
		for (int i = 0; i < model->N_transitions; i++) {
//			std::cout << model->transitions[i].id <<  "(" << enabled[i] << ")" <<std::endl;

			if (enabled[i] != 1)
				continue;

			if ((model->transitions[i].type == TT_IMMEDIATE)) {
				firstT = 0.0;
				firstEventID = i;
				hasImmediateEnabled = 1;
				std::cout << "im Transition " << i << " is enabled (crnt time: " << currentTime <<  " ) -" << std::endl;
				printFluidLevels(model, marking);
				printFluidRates(model, marking);
				firstEventType = TRANSITION;
				imTransVec.push_back(firstEventID);
//				break;
			} else if (hasImmediateEnabled == 0 && (model->transitions[i].type == TT_DETERMINISTIC) &&
				(model->transitions[i].time - clock0[model->transitions[i].idInMarking] < firstT)){
					firstT = model->transitions[i].time - clock0[model->transitions[i].idInMarking];
					firstEventID = i;
					firstEventType = TRANSITION;
			}
		}

		if (hasImmediateEnabled) {
			for (int im = 0; im < imTransVec.size(); im++)
				fireTransition(model, marking, imTransVec[im]);
			continue;
		}

		for (int i = 0; i < model->N_places; i++) {
			if (model->places[i].type == PT_FLUID) {
				
				//if drift is close to zero nothing changes.
				if (IS_ZERO(drift[model->places[i].idInMarking]))
					continue;


				double f = fluid0[model->places[i].idInMarking];
				double d = drift[model->places[i].idInMarking];
				double fluidtime;

				if (d > ZERO_PREC) 
					fluidtime = (model->places[i].f_bound - f) / d;
				else if (d < -ZERO_PREC) 
					fluidtime = - (f) / d;

				if (fluidtime < firstT){
					firstT = fluidtime;
					firstEventID = i;
					if (d > ZERO_PREC) 
						firstEventType = PLACE_UPPER_BOUNDRY;
					else if (d < -ZERO_PREC)
						firstEventType = PLACE_LOWER_BOUNDRY;
				}
			}
		}


		for (int i = 0; i < model->N_arcs; i++) {
			// If it was a gaurd Arc.
			if ( (model->arcs[i].type == AT_TEST || model->arcs[i].type == AT_INHIBITOR) && model->places[model->arcs[i].placeId].type == PT_FLUID) {

				double f = fluid0[model->places[model->arcs[i].placeId].idInMarking];
				double d = drift[model->places[model->arcs[i].placeId].idInMarking];
				//if drift is close to zero nothing changes.

				if (IS_ZERO(d) )
					continue;
				double fluidtime = (model->arcs[i].weight - f) / d;

				if (fluidtime < ZERO_PREC) continue;

				if (fluidtime < firstT){
					firstT = fluidtime;
					firstEventID = i;
					firstEventType = GAURD_ARC;
				}
			}
		}


		//if there is no event before max time reached, we are done.
		if (crntTime + firstT >= model->MaxTime){
			crntEvent = new DtrmEvent(MAX_TIME_REACHED);
			crntEvent->time = model->MaxTime;
			crntEvent->id = -1;
			crntEvent->preRegionMarking = marking;
			crntEvent->postRegionMarking = NULL;
			crntEvent->nextDtrmEvent = 0;
			dtrmEventList.push_back(crntEvent);

			if (preEvent != NULL)
				preEvent->nextDtrmEvent = crntEvent;

			break;
		}


		marking = advanceMarking(model, marking, 0, firstT, drift, enabled);
		crntTime += firstT;

		if (firstEventType == TRANSITION)
			fireTransition(model, marking, firstEventID);

//		crntEvent= new DtrmEvent(firstEventType);
//		crntEvent->time = crntTime;
//		crntEvent->id = firstEventID;
//		crntEvent->preRegionMarking = preMarking;
//		crntEvent->postRegionMarking = marking;
//		crntEvent->eventType = firstEventType;
//		dtrmEventList.push_back(crntEvent);
//		if (preEvent != NULL)
//			preEvent->nextDtrmEvent = crntEvent;
//		preEvent = crntEvent;
	}
}
void TimedDiagram::segmentizeStochasticRegion(Marking* marking, StochasticEvent* eventSeg, double timeBias) {


	Segment* eventLine = eventSeg->timeSegment;
	if ((eventLine->p2.Y > model->MaxTime || eventLine->p2.X > model->MaxTime) && (eventLine->p1.X > model->MaxTime || eventLine->p1.Y > model->MaxTime)) {
		std::cout << ("\n Maximum time reached") << std::endl;
		return;
	}

	double start = eventLine->p1.X ;
	double end = eventLine->p2.X;

	int cntFirst;
	int *enabled;
	int *enabledTransitionCache = NULL; //TODO: it should be static? after each event a transition could be enabled or disabled.



//	if (eventSeg->preRegion != 0)
//		eventSeg->preRegion->print(std::cout);

	if (enabledTransitionCache == NULL)
		enabledTransitionCache = new int[model->N_transitions];


	checkEnabled(model, marking, eventSeg->timeSegment->p1.X, eventSeg->timeSegment->p2.X);
	setActFluidRate(model, marking, start /*- timeBias*/);

	std::cout << regionList.size() << std::endl;
//	printFluidLevels(model, marking);

	double* clock0, *clock1, *fluid0, *fluid1, *drift;

	enabled = marking->enabling;
	clock0 = marking->clock0;
	clock1 = marking->clock1;
	fluid0 = marking->fluid0;
	fluid1 = marking->fluid1;
	drift = marking->fluidPlaceDeriv;
	cntFirst = 0;

	double firstT = 1.0e100;

	int hasImmediateEnabled;

	std::vector<StochasticEvent*>* potentialEvents = new std::vector<StochasticEvent*>();

	cntFirst = 0;
	firstT = 1.0e100;

	// count the number (and put in cache) the number of enabled transitions
	hasImmediateEnabled = 0;
	for (int i = 0; i < model->N_transitions; i++) {
		if (enabled[i] != 1)
			continue;

		if ((model->transitions[i].type == TT_IMMEDIATE)) {
			if (firstT > 0.0) {
				firstT = 0.0;
				enabledTransitionCache[0] = i;
				cntFirst = 1;
				hasImmediateEnabled = 1;
			} else {
				enabledTransitionCache[cntFirst] = i;
				cntFirst++;
			}
		} else if ((model->transitions[i].type == TT_DETERMINISTIC) && (hasImmediateEnabled == 0)) {
				enabledTransitionCache[cntFirst] = i;
				cntFirst++;
		}
	}




	if (hasImmediateEnabled) {
		//TODO: Take care of immediate transitions!!
//		std::cout << "An immediate transition is enabled." << std::endl;
		Marking *newMarking;
		newMarking = copyMarking(model, marking);

		for (int i = 0; i < cntFirst; i++)
			fireTransition(model, newMarking, enabledTransitionCache[i]);

		StochasticEvent* imEvent = new StochasticEvent(eventSeg);
		imEvent->eventType = TRANSITION;
		imEvent->id = enabledTransitionCache[cntFirst-1];
		imEvent->postRegionMarking = newMarking;
		imEvent->preRegionMarking = marking;


		segmentizeStochasticRegion(newMarking, imEvent, timeBias);

	} else {
//		std::cout << "Other stuff!!!!." << std::endl;
		for (int i = 0; i < cntFirst; i++) {
			int trId = enabledTransitionCache[i];
			/**
			* A deterministic transition equation is determined by following:
			* \f$ \Delta t = T_D - (as + b), \Delta t = t - (\alpha s + \beta + [timeBias]) \f$
			* \f$ (as + b) \f$ --> is the current clock value of transition.
			* \f$ (\alpha s + \beta) \f$ --> is the event line equation passed as the argument.
			* \f$ t \f$ --> is the current time (vertical axis of the diagram).
			* Final equation would be :
			*  	\f$ t = (\alpha - a)s + (T_D + \beta - b + [timeBias]).
			* A timeBias should be subtracted from \f$ t \f$ in case this function is not called (in generate diagram) with the start from origin.
			*/
			double a = clock1[model->transitions[trId].idInMarking];
			double b = clock0[model->transitions[trId].idInMarking];

			//TODO:(IMPORTANT) fix the end of each new segment with coinciding the eventLine.
			Segment * detTransSeg = new Segment(eventLine->a - a,
				model->transitions[trId].time + eventLine->b - b /*+ timeBias*/,
				start, end);
			StochasticEvent* detTransEvent = new StochasticEvent(detTransSeg, TRANSITION);
			detTransEvent->id = enabledTransitionCache[i];
			detTransEvent->preRegionMarking = marking;
//			detTransSeg->print();
			potentialEvents->push_back(detTransEvent);

		}
		for (int i = 0; i < model->N_places; i++) {
			if (model->places[i].type == PT_FLUID) {

				//if drift is close to zero nothing changes.
				if (IS_ZERO(drift[model->places[i].idInMarking]))
					continue;


				/**
				* A fluid place equation is determined by following:
				* \f$ B - (as +b) = d\Delta t , \Delta t = t - (\alpha s + \beta) \f$
				* \f$ (as + b) \f$ --> is the current fluid level in the place.
				* \f$ (\alpha s + \beta) \f$ --> is the event line equation passed as the argument.
				* \f$ t \f$ --> is the current time (vertical axis of the diagram).
				* Final equation would be :
				* For d > 0 : \f$ t = (\alpha - \frac{a}{d})s + (\frac{-b + B}{d} + \beta).
				* For d < 0	: \f$ t = (\alpha + \frac{a}{|d|})s + (\frac{b}{|d|} + \beta). 
				*
				* Note: If the drift is negative the B is the lower bound and vice versa. 
				* and all d terms should be absolute.
				* ***********************************************************************************************************************
				* IMPORTANT NOTE:																										*
				* the s value used above is in the frame which its origin is at the point (gTrEnabledTime, gTrEnabledTime).				*
				* so in order to obtain the equation in the main frame we shoud replace all s and t should be minused by gTrEnabledTime.*
				* so the new equations would be:																						*
				* For d > 0 : \f$ t = (\alpha - \frac{a}{d})s + (\frac{-b + B}{d} + \beta) +											*
				*					  gTrEnabledTime(1 - \alpha + \frac{a}{d}).															*
				* For d < 0	: \f$ t = (\alpha + \frac{a}{|d|})s + (\frac{b}{|d|} + \beta) +												*
				*                     gTrEnabledTime(1 - \alpha - \frac{a}{|d|}).														*
				*************************************************************************************************************************
				*/

				double a = fluid1[model->places[i].idInMarking];
				double b = fluid0[model->places[i].idInMarking];
				double d = drift[model->places[i].idInMarking];

				StochasticEvent* placeEvent = new StochasticEvent();
				Segment* placeSeg;
				if (d > ZERO_PREC) {

					placeSeg = new Segment(eventLine->a - a / d,
						(model->places[i].f_bound - b) / d + eventLine->b,
						//+ gTrEnabledTime*(1- eventLine->a + a / d),
						/*+ timeBias*/
						start, end);
					placeEvent->eventType = PLACE_UPPER_BOUNDRY;
				} else if (d < -ZERO_PREC) {

//					std::cout << a << ", " << std::abs(d) + eventLine->b << ", " << start << ", " <<  end << std::endl;
					placeSeg = new Segment(eventLine->a + a / std::abs(d),
						b / std::abs(d) + eventLine->b,
						//+ gTrEnabledTime*(1- eventLine->a - a /abs(d))/*+ timeBias*/,
						start, end);
					placeEvent->eventType = PLACE_LOWER_BOUNDRY;
				}


				//regarding some 
				if (eventSeg->eventType ==  placeEvent->eventType && eventSeg->id == i){
					drift[model->places[i].idInMarking] = 0;
					continue;
				}

				placeEvent->timeSegment = placeSeg;
				
				placeEvent->id = i;
				placeEvent->preRegionMarking= marking;
				potentialEvents->push_back(placeEvent);
				cntFirst++;
			}
		}


		for (int i = 0; i < model->N_arcs; i++) {
			// If it was a gaurd Arc.
			if ( (model->arcs[i].type == AT_TEST || model->arcs[i].type == AT_INHIBITOR) && (model->places[model->arcs[i].placeId].type == PT_FLUID) ) {

				double a = fluid1[model->places[model->arcs[i].placeId].idInMarking];
				double b = fluid0[model->places[model->arcs[i].placeId].idInMarking];
				double d = drift[model->places[model->arcs[i].placeId].idInMarking];

				//if drift is close to zero nothing changes.

				if (IS_ZERO(d) )
					continue;
				/**
				* A fluid place equation is determined by following:
				* \f$ W - (as +b) = d\Delta t , \Delta t = t - (\alpha s + \beta) \f$
				* \f$ (as + b) \f$ --> is the current fluid level in the place.
				* \f$ (\alpha s + \beta) \f$ --> is the event line equation passed as the argument.
				* \f$ t \f$ --> is the current time (vertical axis of the diagram).
				* \f$ W \f$ --> is the gaurd associated with the arc.
				* Final equation would be :
				* \f$ t = (\alpha - \frac{a}{d})s + (\frac{-b + W}{d} + \beta).
				*
				**/


				StochasticEvent* arcEvent = new StochasticEvent();
				Segment* arcSeg;

				arcSeg = new Segment(eventLine->a - a / d, (model->arcs[i].weight - b) / d + eventLine->b, start, end);
				arcEvent->eventType = GAURD_ARC;


				//regarding some
				if (eventSeg->eventType ==  arcEvent->eventType && eventSeg->id == i){
					continue;
				}

				arcEvent->timeSegment = arcSeg;

				arcEvent->id = i;
				arcEvent->preRegionMarking= marking;
				potentialEvents->push_back(arcEvent);
				cntFirst++;

			}
		}


		if (cntFirst > 0) { // if there is at least an event that can happen

			std::vector<StochasticEvent*> *nextEvents = new std::vector<StochasticEvent*>();
//			computeNextEvents(potentialEvents, eventLine, nextEvents);
//
//			potentialEvents->clear();
//			delete potentialEvents;
//
//			//Region* region = new Region(nextEvents, eventLine);
//			//region->marking = marking;
//			//region->timeBias = timeBias;
//
//			//regionList.push_back(region);
//
//			createAddRegions(nextEvents, eventSeg, marking);

//			char* eventTypeMap [] = {"TRANSITION", "PLACE_LOWER_BOUNDRY", "PLACE_UPPER_BOUNDRY", "GAURD_ARC", "MAX_TIME_REACHED", "FIRST_NULL_EVENT"};
//			std::cout << "---------------------------" << std::endl;
//			printFluidRates(model, marking);
//			printFluidLevels(model, marking);
//			std::cout << "---------------------------" << std::endl;
//			std::cout << "Potential events: " << std::endl;
//			for(int ii = 0; ii < potentialEvents->size(); ii++){
//				std::cout << potentialEvents->at(ii)->id << ": " << eventTypeMap[potentialEvents->at(ii)->eventType] << std::endl;
//				potentialEvents->at(ii)->timeSegment->print();
//			}
//			std::cout << "---------------------------" << std::endl;

			computeNextEventCreatePoly(potentialEvents, eventSeg, marking, nextEvents);

//			std::cout << "---------------------------" << std::endl;
//			std::cout << "next events: " << std::endl;
//			for(int ii = 0; ii < nextEvents->size(); ii++)
//					nextEvents->at(ii)->timeSegment->print();
//			std::cout << "---------------------------" << std::endl;

			potentialEvents->clear();
			delete potentialEvents;


			for (unsigned int i = 0; i < nextEvents->size(); i++) {

				//if (IS_ZERO(nextEvents->at(i)->timeSegment->p1.X - nextEvents->at(i)->timeSegment->p2.X))
				//	continue;

				Marking *newMarking;
				//TODO: Here we should check for multiple possible events.

				double t1 = nextEvents->at(i)->timeSegment->a - eventLine->a;
				double t0 = nextEvents->at(i)->timeSegment->b - eventLine->b /*- timeBias*/;


				//in case that no time advancment is needed!
				//if (IS_ZERO(t1) && IS_ZERO(t0))
				//	continue;


				newMarking = advanceMarking(model, marking, t1, t0, drift, enabled);

				if (nextEvents->at(i)->eventType == TRANSITION)
					fireTransition(model, newMarking, nextEvents->at(i)->id);

				nextEvents->at(i)->postRegionMarking = newMarking;

				segmentizeStochasticRegion(newMarking, nextEvents->at(i), timeBias);
			}
		} else {

			Segment* maxSeg = new Segment(0, model->MaxTime, eventLine->p1.X, eventLine->p2.X);
			StochasticEvent* maxTimeReached = new StochasticEvent(maxSeg, MAX_TIME_REACHED);
			maxTimeReached->preRegionMarking = marking;
			maxTimeReached->postRegionMarking = marking;
			maxTimeReached->id = -1;
			std::vector<StochasticEvent*> *maxTimeEvents = new std::vector<StochasticEvent*>();
			maxTimeEvents->push_back(maxTimeReached);
			Region* region = new Region(maxTimeEvents, eventLine);
			region->marking = marking;
			region->timeBias = timeBias;

			if (eventSeg->preRegion != 0)
				eventSeg->preRegion->successors->push_back(region);
			else //if we have entered from the deterministic are to this new region.
				eventSeg->preDtrmEvent->nextRegions->push_back(region);

			regionList.push_back(region);
		}

	}

}

void TimedDiagram::computeNextEvents(std::vector<StochasticEvent*> * potentialEvents, Segment *uSegment, std::vector<StochasticEvent*> *nextEvents) {

	std::sort(potentialEvents->begin(), potentialEvents->end(),
		StochasticEvent::greaterSlopeFirst);


	std::vector<Segment*> workingSegVec;

	Point uPoint1 = uSegment->p1;
	Point uPoint2 = uSegment->p2;

	for (unsigned int i = 0; i < potentialEvents->size(); i++) {

		if (potentialEvents->at(i)->timeSegment->intersect(*uSegment, uPoint2) && uPoint2 != uSegment->p2 && uPoint2 != uSegment->p1){
			workingSegVec.push_back(new Segment(uPoint1, uPoint2));
			uPoint1 = uPoint2;
			uPoint2 = uSegment->p2;
		}
	}

	workingSegVec.push_back(new Segment(uPoint1, uPoint2));

	for (unsigned int i = 0; i < workingSegVec.size(); i++)
		minLines(potentialEvents, workingSegVec[i], nextEvents);
}

void TimedDiagram::minLines(std::vector<StochasticEvent*> * potentialEvents, Segment *uSegment, std::vector<StochasticEvent*> * nextEvents){

//	std::cout << "-------------------in min lines." << std::endl;
//	std::cout << "uSegment: ";
//	uSegment->print();
//	std::cout << "-------------------" << std::endl;

	double start = uSegment->p1.X;
	double end = uSegment->p2.X;

//	if (IS_ZERO(start - end)) return;


	double minStart = INF;
	int minIndex;
	for (unsigned int i = 0; i < potentialEvents->size(); i++) {
//		std::cout << "potential[" << i << "]";
//		potentialEvents->at(i)->timeSegment->print();

		if (potentialEvents->at(i)->timeSegment->getY(start) < uSegment->getY(start)) continue;
		if (potentialEvents->at(i)->timeSegment->getY(start) < minStart) {
			minStart = potentialEvents->at(i)->timeSegment->getY(start);
			minIndex = i;
		}
	}

	int crntIndex = minIndex;

	Point p1, p2;
	int nextIndex = crntIndex;

	p1.X = start;
	p1.Y = potentialEvents->at(crntIndex)->timeSegment->getY(start);

	while (crntIndex < potentialEvents->size()) {

//		std::cout << "----crntIndex: " << crntIndex << std::endl;

		p2 = potentialEvents->at(crntIndex)->timeSegment->p2;  //this is the first point that crntIndex line is intersecting with.

		for (unsigned int i = crntIndex + 1; i < potentialEvents->size(); i++) {
			Point p;
			if (potentialEvents->at(crntIndex)->timeSegment->intersect(*potentialEvents->at(i)->timeSegment, p))
				if (p.Y <= p2.Y) {
					p2 = p;
					nextIndex = i;
				}
		}
//		std::cout << "----nextIndex: " << nextIndex << std::endl;

		//if the intersection point after the end of current underlying segment we are done.
		if (p2.X >= end){
			p2.X = end;
			p2.Y = potentialEvents->at(crntIndex)->timeSegment->getY(end);
		}


		Segment* newSegment = new Segment(potentialEvents->at(crntIndex)->timeSegment->a, 
			potentialEvents->at(crntIndex)->timeSegment->b, p1.X, p2.X);



		StochasticEvent* sEvent = new StochasticEvent(potentialEvents->at(crntIndex));
		sEvent->timeSegment = newSegment;

		nextEvents->push_back(sEvent);

//		std::cout << "minIndex:"  << minIndex << std::endl;
//		std::cout << "crntIndex:" << crntIndex << std::endl;
//		std::cout << "potentialEvents->size():" << potentialEvents->size()<< std::endl;
//		uSegment->print();

		//if there is no intersection or if we have reached the end of the current underlying segment we are done.
		if (crntIndex == nextIndex || p2.X == end)
			break; 

		crntIndex = nextIndex;
		p1 = p2;
	}
}

void TimedDiagram::createAddRegions(std::vector<StochasticEvent*> * eventList, StochasticEvent * preEvent, Marking* marking){
	int index1 = 0; 
	int index2 = 0;
	Segment* lowerBoundray = preEvent->timeSegment;
	Point p1 = lowerBoundray->p1;
	Point p2 = lowerBoundray->p2;
	for (int i = 0; i < eventList->size(); i++){
		index2 = i;
		if (eventList->at(i)->timeSegment->p2.Y == lowerBoundray->getY(eventList->at(i)->timeSegment->p2.X) ||
			eventList->at(i)->timeSegment->p2.X == lowerBoundray->p2.X)
		{
			//in case we are here according to the second condition.
			p2.X = eventList->at(i)->timeSegment->p2.X;
			p2.Y = lowerBoundray->getY(eventList->at(i)->timeSegment->p2.X);

			std::vector<StochasticEvent*>* events = new std::vector<StochasticEvent*>();
			for (int j = index1; j <= index2; j++)
				events->push_back(eventList->at(j));

			Segment* lb = new Segment(p1, p2);
			Region* region = new Region(events, lb);
			region->marking = marking;

			regionList.push_back(region);

			//if we have entered from an stochastic region to this new region.
			if (preEvent->preRegion != 0)
				preEvent->preRegion->successors->push_back(region);
			else //if we have entered from the deterministic area to this new region.
				preEvent->preDtrmEvent->nextRegions->push_back(region);

			for (int j = index1; j <= index2; j++)
				eventList->at(j)->preRegion = region;

			p1 = p2;
			index1 = index2 + 1;
		}
	}
}

void TimedDiagram::computeNextEventCreatePoly(std::vector<StochasticEvent*> * potentialEvents, StochasticEvent * preEvent, Marking* marking, std::vector<StochasticEvent*> *nextEvents){

	std::sort(potentialEvents->begin(), potentialEvents->end(),
		StochasticEvent::greaterSlopeFirst);

	Segment* uSegment = preEvent->timeSegment;

	std::vector<Segment*> workingSegVec;

	Point uPoint1 = uSegment->p1;


	//Removing those lines under the preEvent segments.
	std::vector<int> toRemove(0);
	for (unsigned int i = 0; i < potentialEvents->size(); i++) {
//		if (uSegment->a >= potentialEvents->at(i)->timeSegment->a && uSegment->p1.Y >= potentialEvents->at(i)->timeSegment->p1.Y )
		if (uSegment->p1.Y >= potentialEvents->at(i)->timeSegment->p1.Y && uSegment->p2.Y >= potentialEvents->at(i)->timeSegment->p2.Y)
			toRemove.push_back(i);
	}
	for (unsigned int i = 0; i < toRemove.size(); i++){
		potentialEvents->erase(potentialEvents->begin() + toRemove[toRemove.size() - 1 - i]);
	}

//	std::cout << "-----------------------" << std::endl;
//	for(int ii = 0; ii < potentialEvents->size(); ii++)
//		potentialEvents->at(ii)->timeSegment->print();


	Point p;
	for (unsigned int i = 0; i < potentialEvents->size(); i++) {
//		potentialEvents->at(i)->timeSegment->print();
		if (potentialEvents->at(i)->timeSegment->intersect(*uSegment, p) && p != uSegment->p2 && p != uSegment->p1){
			if (!IS_ZERO(uPoint1.X - p.X)){
				workingSegVec.push_back(new Segment(uPoint1, p));
				uPoint1 = p;
			}
		}
	}

	if (!IS_ZERO(uPoint1.X - uSegment->p2.X))
		workingSegVec.push_back(new Segment(uPoint1, uSegment->p2));


	for (unsigned int i = 0; i < workingSegVec.size(); i++){
		std::vector<StochasticEvent*> *crntNextEvents = new std::vector<StochasticEvent*>();
		minLines(potentialEvents, workingSegVec[i], crntNextEvents);

		Region* region = new Region(crntNextEvents, workingSegVec[i]);
		region->marking = marking;

		regionList.push_back(region);

		//if we have entered from an stochastic region to this new region.
		if (preEvent->preRegion != 0)
			preEvent->preRegion->successors->push_back(region);
		else //if we have entered from the deterministic area to this new region.
			preEvent->preDtrmEvent->nextRegions->push_back(region);

		for (unsigned int j = 0; j < crntNextEvents->size(); j++){
			crntNextEvents->at(j)->preRegion = region;
			nextEvents->push_back(crntNextEvents->at(j));
		}

	}

//	std::cout << "new regions: " << std::endl;
//	for(int ii = 0; ii < workingSegVec.size(); ii++)
//			regionList[regionList.size() - 1 - ii]->print(std::cout);

}

double TimedDiagram::calProbAtTime(double time, double (*sPdfInt)(double), bool (*isPropHolds)(Model*, Marking*, double t0, double t1, double& , double& )){
	//TODO: by sorting region in more genuine way this function could be more effecient.

	double s1, s2;
	double prob = 0;
	Point p1, p2;
	
	if (time > gTrEnabledTime) {
		//stochastic part : only if at the given time g-transition could have been enabled.
		/**
		* IMPORTANT note: we have computed everything for the stochastic part in the frame with origin at (gTrEnabledTime, gTrEnabledTime)
		* So we should move to that frame.
		*/
		double sFrameTime = time - gTrEnabledTime; 

		Segment timeSeg(0, time , 0, model->MaxTime);
		for (int i = 0; i < regionList.size(); i++){
			//s1 and s2 are the validity interval for which the probability holds, returned by the function iPropHolds.
			//std::cout << "i=" << i << "-";
			if (regionList[i]->intersect(timeSeg, p1, p2)){
				s1 = p1.X; s2 = p2.X;
				double t0 = time - regionList[i]->lowerBoundry->b /*- regionList[i]->timeBias*/;
				double t1 = - regionList[i]->lowerBoundry->a;
				if (isPropHolds(model, regionList[i]->marking, t0, t1, s1, s2)){
					prob += fabs(sPdfInt(s1 - gTrEnabledTime) - sPdfInt(s2 - gTrEnabledTime));
				}
			}
		}
	}

	//determinestic part after g-transition firing
	int cc;
	for (cc = 0; cc < dtrmEventList.size() && time > dtrmEventList[cc]->time ; cc++);
	s1 = time > gTrEnabledTime ? time : 0;
	s2 = model->MaxTime;
	double t;
	if (cc == 0)
		t = time;
	else 
		t = time - dtrmEventList[cc - 1]->time;
	if (isPropHolds(model, dtrmEventList[cc]->preRegionMarking, t , 0, s1, s2))
		prob += fabs(sPdfInt(time > gTrEnabledTime ? s1 - gTrEnabledTime : s1) - sPdfInt(time > gTrEnabledTime ? s2 - gTrEnabledTime : s2));


	return prob;
}

//void TimedDiagram::saveDiagram(std::string filename){
//
//	const int row = model->MaxTime * scale;
//	const int col = model->MaxTime * scale;
//
//	debugImage = cv::Mat(row, col, CV_8UC3, cv::Scalar::all(255));
//
//
//
//	for (int i = 0; i < regionList.size(); i++){
//		std::cout << "-----------------------------------------" << std::endl;
////		regionList[i]->print(std::cout, model);
//		drawSegmet(debugImage, *regionList[i]->leftBoundry, scale);
//		drawSegmet(debugImage, *regionList[i]->rightBoundry, scale);
//		drawSegmet(debugImage, *regionList[i]->lowerBoundry, scale);
//		for (int j = 0; j < regionList[i]->eventSegments->size(); j++)
//			drawSegmet(debugImage, *regionList[i]->eventSegments->at(j)->timeSegment, scale);
//
//		cv::Point p(scale*(regionList[i]->lowerBoundry->p1.X + regionList[i]->lowerBoundry->p2.X)/2,
//			scale*(regionList[i]->lowerBoundry->p1.Y + regionList[i]->lowerBoundry->p2.Y)/2);
////		char rn[10];
////		sprintf(rn, "%d", i);
////		cv::putText(debugImage, std::string(rn), p, cv::FONT_HERSHEY_SIMPLEX, 1, cv::Scalar(0, 0, 255));
//	}
//
//	for (int i = 0; i < dtrmEventList.size(); i++){
//		cv::line(debugImage, cv::Point((dtrmEventList[i]->time/* - gTrEnabledTime*/) * scale, (dtrmEventList[i]->time /* - gTrEnabledTime*/)* scale),
//			cv::Point(col, (dtrmEventList[i]->time /* - gTrEnabledTime*/)* scale), cv::Scalar(0, 0, 0), 1);
//	}
//
//	cv::line(debugImage, cv::Point(0, 0), cv::Point(row, col), cv::Scalar(128, 128, 0), 2);
//
//	cv::Mat flipped;
//	cv::flip(debugImage, flipped, 0);
//	cv::imwrite(filename + ".jpg", flipped);
//
//	//cv::Mat im2show;
//	//cv::resize(diagram, im2show, cv::Size(1000, 1000));
//	//cv::imshow("regions", im2show);
//	//cv::waitKey();
//}

//void TimedDiagram::drawSegmet(cv::Mat & image, Segment& seg, const int scale){
////	double x, y;
////
////
////	y = std::min(seg.p1.Y * scale, (double)image.rows);
////	x = std::min(seg.p1.X * scale, (double)image.cols);
////	cv::Point p1(x, y);
////
////	y = std::min(seg.p2.Y * scale, (double)image.rows);
////	x =  std::min(seg.p2.X * scale, (double)image.cols);
////	cv::Point p2(x, y);
//
//	Point zero (0, 0);
//	Point mid (0, model->MaxTime);
//	Point end (model->MaxTime, model->MaxTime);
//
//	Segment upperBound(mid, end);
//	Segment leftBound(zero, mid);
//
//	Point intPoint;
//	cv::Point2d p1;
//	cv::Point2d p2;
//	if (seg.intersect(upperBound, intPoint)){
//		p1.x = intPoint.X; p1.y = intPoint.Y;
//		if (seg.p1.Y < seg.p2.Y){
//			p2.x = seg.p1.X; p2.y = seg.p1.Y;
//		}else{
//			p2.x = seg.p2.X; p2.y = seg.p2.Y;
//		}
//	}else if (seg.intersect(leftBound, intPoint)){
//		p1.x = intPoint.X; p1.y = intPoint.Y;
//		if (seg.p1.X > seg.p2.X){
//			p2.x = seg.p1.X; p2.y = seg.p1.Y;
//		}else{
//			p2.x = seg.p2.X; p2.y = seg.p2.Y;
//		}
//	}
//	else{
//		p1.x = seg.p1.X; p1.y = seg.p1.Y;
//		p2.x = seg.p2.X; p2.y = seg.p2.Y;
//	}
//
//	p1.x = p1.x* scale; p1.y = p1.y* scale;
//	p2.x = p2.x* scale; p2.y = p2.y* scale;
//
////
////
////	if (seg.p1.Y >= model->MaxTime){
////		y = (double)image.rows;
////		x = seg.getX(model->MaxTime) * scale;
////	}else{
////		y = seg.p1.Y * scale;
////		x = seg.p1.X * scale;
////	}
////	cv::Point p1(x, y);
////
////	if (seg.p2.Y >= model->MaxTime){
////		y = (double)image.rows;
////		x = seg.getX(model->MaxTime) * scale;
////	}else{
////		y = seg.p2.Y * scale;
////		x = seg.p2.X * scale;
////	}
////	cv::Point p2(x, y);
//
//
//	if (p1.x > p1.y || p2.x > p2.y)
//		cv::line(image, p1, p2, cv::Scalar(0, 0, 255), 2);
//	else
//		cv::line(image, p1, p2, cv::Scalar(0, 0, 0), 1);
//}
