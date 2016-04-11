/*
 * IntervalSet.h
 *
 *  Created on: Mar 5, 2013
 *      Author: hamed
 */

#ifndef INTERVALSET_H_
#define INTERVALSET_H_

#include <vector>
#include <algorithm>
#include <iostream>

class IntervalSet;

class Interval{
public:
	double start;
	double end;

	Interval(double _start, double _end);
	Interval(): start(-1), end(-1){};

	Interval & operator =(const Interval &in){start = in.start; end = in.end; return *this;};
	bool intersect(const Interval & in, Interval & result);
	IntervalSet* minus(const Interval &in);
};


class IntervalSet {
public:
	std::vector<Interval> intervals;

	IntervalSet();
	virtual ~IntervalSet();

	IntervalSet* intersect(const IntervalSet* intSet);
	IntervalSet* unionWith(const IntervalSet* intSet);
	IntervalSet* minus(const IntervalSet* intSet);
	IntervalSet* intersect(const Interval intSet);
	IntervalSet* unionWith(const Interval intSet);
	IntervalSet* minus(const Interval intSet);

	/*
	 * Sorts based on the start point of intervals.
	 */
	void sort();

	/*
	 * merges the internal intersecting intervals with each others.
	 */
	void removeRedundency();

	void clear() {
		if (intervals.size() > 0)
			intervals.erase(intervals.begin(), intervals.end());
	};
	bool isEmpty() {return intervals.empty();};

	void print(std::ostream& out);
};

#endif /* INTERVALSET_H_ */
