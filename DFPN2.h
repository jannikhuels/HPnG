/*
 * DFPN.h
 *
 *  Created on: 10 nov. 2011
 *      Author: GhasemiehH
 */

#ifndef DFPN2_H_
#define DFPN2_H_

#define PT_DISCRETE   0
#define PT_FLUID      1

//#define IS_ZERO

typedef struct {
	int type; // Place type
	char *id; // Name of the place
	int d_mark; // Discrete marking (for discrete places)
	double f_level; // Fluid level (for fluid places)
	double f_bound; // Upper bound (for fluid places)

	//// computed field
	int idInMarking;

	int inputListSize; // number of input arcs arriving to this place
	int *inputList; // list  of input arcs arriving to this place
	int inhibListSize; // number of inhibitor arcs starting from this place
	int *inhibList; // list  of inhibitor arcs starting from this place
	int testListSize; // number of test arcs starting from this place
	int *testList; // list  of test arcs starting from this place
	int outputListSize; // number of output arcs exiting this place
	int *outputList; // list  of output arcs exiting to this place

} Place;

#define TT_DETERMINISTIC   0
#define TT_IMMEDIATE       1
#define TT_FLUID           2
#define TT_GENERAL         3
#define TT_FLUID_DYNA      4

typedef struct {
	///////// Model field
	int type; // Transition type
	char *id; // Name of the transition
	double time; // Firing time (deterministic)
	double weight; // Weight (for both deterministic and immediate)
	int priority; // Priority (for both deterministic and immediate)
	double flowRate; // Flow rate (for fluid transitions)
	char *distr; // distribution (for general)

	//////////For dynamic transitions only
	int dependencySize;
	int* dependencyList;  //list of transitions that this transition is dependent on. (only linear dependency is considered now.)
	double* coef;	//influence weight of each transition.
	//so for dynamic transitions: flowRate = \sume_i (coef[i]*dependencyList[i].actualFlowRate)

	////////// computed field
	int inputListSize; // number of input arcs arriving to this transition
	int *inputList; // list  of input arcs arriving to this transition
	int inhibListSize; // number of inhibitor arcs arriving to this transition
	int *inhibList; // list  of inhibitor arcs arriving to this transition
	int testListSize; // number of test arcs arriving to this transition
	int *testList; // list  of test arcs arriving to this transition
	int outputListSize; // number of output arcs exiting this transition
	int *outputList; // list  of output arcs exiting to this transition

	int idInMarking;
} Transition;

#define AT_DISCRETE_INPUT 0
#define AT_DISCRETE_OUTPUT 1
#define AT_FLUID_INPUT 2
#define AT_FLUID_OUTPUT 3
#define AT_INHIBITOR 4
#define AT_TEST 5

typedef struct {
	int type; // Transition type
	char *id; // Name of the transition
	int fromId; // id of the input (place or transition)
	int toId; // id of the output (place or transition)
	double weight; // arc weight
	double share; // fluid flow share
	int priority; // fluid flow priority

	////computed field
	int transId; // transition ID
	int placeId; // place ID
} Arc;

typedef struct {
	int *tokens; // Token for each discrete place
	double *fluid0; // Fluid level in fluid places
	double *fluid1; // expressed as: fluid1 * s + fluid0
	int *enabling; // Tells if a transition is enabled  - size equal to M->N_Transitions
	double *actFluidRate; // Actual fluid rate for fluid transitions
	double *clock0; // Remaining time to fire for deterministic transitions
	double *clock1; // expressed as: clock1 * s + clock0
	int *generalHasFired; // telss if a general transition has fired
	double *generalDisabled; // time that the general transition has been disbled
	int N_generalFired; // Number of general transitions that have fired

	double *fluidPlaceDeriv; // derivative of the fluid level in the place
} Marking;

struct StateTimeAlt_tag;

typedef struct State_tag {
	Marking *M; // Marking in the state
	double t0; // t0 and t1 determines the time T at which the system entered this state:
	double t1; // if s is the general transition sample, then
	// the state was entered at T = t1 * s + t0

	double leftInt; // Left bound of the validity region
	double rightInt; // Right bound of the validity region

	int stateID; // id of the state (for printing purposes);
	struct State_tag *next; // pointer to the "next" state (for printing purposes)

	struct StateTimeAlt_tag *nextStateList; // Next states

} State;

typedef struct StateProbAlt_tag {
	double p; // probability of choosing this alternative
	State *S; // next state that will be reached

	struct StateProbAlt_tag *next; // next element in the list
} StateProbAlt;

typedef struct StateTimeAlt_tag {
	StateProbAlt *sa; // State probability alternatives

	double leftInt; // the general transition duration should be greater than leftInt;
	double rightInt; // the general transition duration should be smaller than rightInt;

	struct StateTimeAlt_tag *next; // next element in the list
} StateTimeAlt;

typedef struct {
	int N_places; // number of places
	Place *places; // list of places

	int N_transitions; // number of transitions
	Transition *transitions; // list of transitions

	int N_arcs; // number of arcs
	Arc *arcs; // list of arcs

	double MaxTime;

	////computed field

	int N_discretePlaces;
	int N_fluidPlaces;
	int N_determTransitions;
	int N_fluidTransitions;
	int N_DynamicTransitions;
	int N_generalTransitions;

	State *initialState; // initial marking
} Model;

/*************** for minimum line algorithm *******************/

typedef struct {
	int lSize; // number of regions
	double *split; // array of the split points (size lSize + 1)
	int *nSeg; // array that contains the number of curve that have their minimum in that region (size lSize)
	int **min; // 2D array that contains the minimum of each region (1st size lSize, 2nd size nSeg[i])
} minList;

#define PRECISION 1.0e-12
#define BIGNUMBER 1.0e50

#define ZERO_PREC  1.0e-5

bool IS_ZERO(double x);

Model *ReadModel(const char *FileName);
void InitializeModel(Model *M);

Marking *allocMarking(Model *M);
void freeMarking(Model *M, Marking *K);
Marking *copyMarking(Model *M, Marking *Src);
void printState(Model *M, State *S);
Marking *createInitialMarking(Model *M);
void checkEnabled(Model *M, Marking *K, double s_min = 0, double s_max = 0);
void ShareFlow(Model *M, Marking *K, double Flux, int *_arcs, int N_arcs);
void setActFluidRate(Model *M, Marking *K, double s0);
StateTimeAlt *makeTimeAlt(double left, double right, StateTimeAlt *next);
StateProbAlt *makeProbAlt(State *S, double P);
Marking *advanceMarking(Model *M, Marking *m, double T1, double T0,
		double *drift, int *enabled);
void fireTransition(Model *M, Marking *NewM, int Tr);
bool isGTransitionEnabled(Model* model, Marking* marking);
void fireGeneralTransition(Model *M, Marking *NewM);

int gTransitionId(Model* M);

void printFluidRates(Model* M, Marking* K);
void printFluidLevels(Model* M, Marking* K);
int findPlace(Model *M, char *id);
int findTransition(Model *M, char *id);
int findArc(Model *M, char*id);


#endif /* DFPN_H_ */
