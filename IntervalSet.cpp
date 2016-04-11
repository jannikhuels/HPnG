/*
 * IntervalSet.cpp
 *
 *  Created on: Mar 5, 2013
 *      Author: hamed
 */

#include "IntervalSet.h"


bool loweStartFirst(const Interval& i1,const Interval& i2) {
	return (i1.start < i2.start);
}

IntervalSet::IntervalSet() {
	// TODO Auto-generated constructor stub
}

IntervalSet::~IntervalSet() {
	// TODO Auto-generated destructor stub
}

IntervalSet* IntervalSet::intersect(const IntervalSet* intSet) {
	if (intSet == NULL) return NULL;
	IntervalSet *intersectionSet = new IntervalSet();
	for (unsigned int i = 0; i < intervals.size(); i++) {
		for (unsigned int j = 0; j < intSet->intervals.size(); j++) {
			Interval res;
			if (intervals.at(i).intersect(intSet->intervals.at(j), res))
				intersectionSet->intervals.push_back(res);
		}
	}

	return intersectionSet;

}


IntervalSet* IntervalSet::unionWith(const IntervalSet* uSet) {
	if (uSet == NULL) return this;
	IntervalSet *unionSet = new IntervalSet();
	for (unsigned int i = 0; i < intervals.size(); i++)
		unionSet->intervals.push_back(intervals.at(i));
	for (unsigned int i = 0; i < uSet->intervals.size(); i++)
		unionSet->intervals.push_back(uSet->intervals.at(i));
	unionSet->removeRedundency();
	return unionSet;
}

Interval::Interval(double _start, double _end) {
	start = std::min(_start, _end);
	end = std::max(_start, _end);
}

bool Interval::intersect(const Interval& in, Interval& result) {
	if (start >= in.end || end <= in.start)
		return false;

	if (start > in.start) {
		if (end < in.end) {
			result.start = start;
			result.end = end;
		} else {
			result.start = start;
			result.end = in.end;
		}
	} else {
		if (end > in.end) {
			result.start = in.start;
			result.end = in.end;
		} else {
			result.start = in.start;
			result.end = end;
		}
	}

	return true;
}

void IntervalSet::sort() {
	std::sort(intervals.begin(), intervals.end(), loweStartFirst);
}

void IntervalSet::removeRedundency() {
	sort();
	std::vector<Interval> newSet;
	for (unsigned int i = 0; i < intervals.size(); i++){
		Interval newInt(intervals.at(i).start, intervals.at(i).end);
		int j = i+1;
		while(j < intervals.size() && intervals.at(j).start < newInt.end){
			newInt.end = std::max(newInt.end, intervals.at(j).end);
			j++;
		}
		newSet.push_back(newInt);
		i = j-1;
	}

	intervals.clear();
	for (unsigned int i = 0; i < newSet.size(); i++){
		Interval toSave = newSet[i];
		while(toSave.end == newSet[i+1].start){
			toSave.end = newSet[i+1].end;
			i++;
		}
		intervals.push_back(toSave);
	}
}

IntervalSet* IntervalSet::minus(const IntervalSet* intSet) {
	if (intSet == NULL) return this;
	IntervalSet* result = new IntervalSet();
	IntervalSet* tmp = new IntervalSet();

	for (unsigned int i = 0; i < intervals.size(); i++)
		result->intervals.push_back(intervals.at(i));


	for (unsigned int i = 0; i < intSet->intervals.size(); i++){
		tmp->clear();
		for (unsigned int j = 0; j < result->intervals.size(); j++){
			IntervalSet* tmpInt = result->intervals.at(j).minus(intSet->intervals.at(i));
			if (tmpInt != NULL){
				for (unsigned int k = 0; k < tmpInt->intervals.size(); k++)
					tmp->intervals.push_back(tmpInt->intervals.at(k));
			}

		}

		result->clear();
		for (unsigned int k = 0; k < tmp->intervals.size(); k++)
			result->intervals.push_back(tmp->intervals.at(k));
	}

	return result;

}

IntervalSet* IntervalSet::intersect(const Interval inter) {
	IntervalSet* iSet = new IntervalSet();
	iSet->intervals.push_back(inter);
	return this->intersect(iSet);
}

IntervalSet* IntervalSet::unionWith(const Interval inter) {
	IntervalSet* iSet = new IntervalSet();
	iSet->intervals.push_back(inter);
	return this->unionWith(iSet);

}

IntervalSet* IntervalSet::minus(const Interval inter) {
	IntervalSet* iSet = new IntervalSet();
	iSet->intervals.push_back(inter);
	return this->minus(iSet);

}

void IntervalSet::print(std::ostream& out) {
	out << "{";
	for (unsigned int i = 0; i < intervals.size(); i++)
		out << "(" << intervals.at(i).start << "," << intervals.at(i).end << ") ";
	out << "}" << std::endl;
}

IntervalSet* Interval::minus(const Interval& tosub) {
	IntervalSet* ret = new IntervalSet();
	if (tosub.end <= this->start || tosub.start >= this->end)
		ret->intervals.push_back(*this);
	else if(tosub.end > this->start && tosub.end < this->end && tosub.start <= this->start)
		ret->intervals.push_back(Interval(tosub.end, this->end));
	else if(tosub.end >= this->end && tosub.start <= this->start)
		ret = NULL;
	else if(tosub.end < this->end && tosub.start > this->start){
		ret->intervals.push_back(Interval(this->start, tosub.start));
		ret->intervals.push_back(Interval(tosub.end, this->end));
	}
	else if(tosub.end >= this->end && tosub.start > this->start && tosub.start < this->end)
		ret->intervals.push_back(Interval(this->start, tosub.start));

	return ret;
}




