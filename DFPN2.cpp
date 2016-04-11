/*
 * DFPN.cpp
 *
 *  Created on: 10 nov. 2011
 *      Author: GhasemiehH
 */

#include "DFPN2.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "rapidxml/rapidxml.hpp"

using namespace rapidxml;
using namespace std;

char *transTypeName[] = { "Det. ", "Imm. ", "Fluid", "Gen. ", "Dyna." };
char *arcsTypeName[] = { "D.In ", "D.Out", "F.In ", "F.Out", "Inhib", "Test " };

/*********** Read Model (support function) *********************/
#define READ_BUF_DIM 16384

bool IS_ZERO(double x){
	return x < ZERO_PREC && x > -ZERO_PREC;
}

double ReadData(FILE *fp) {
	static char Buff[READ_BUF_DIM];
	double Val;

	fgets(Buff, READ_BUF_DIM, fp);
	while (!feof(fp) && (Buff[0] == '#'))
		fgets(Buff, READ_BUF_DIM, fp);
	sscanf(Buff, "%lf", &Val);

	return Val;
}

char *ReadString(FILE *fp) {
	static char Buff[READ_BUF_DIM];
	static char Val[READ_BUF_DIM];

	fgets(Buff, READ_BUF_DIM, fp);
	while (!feof(fp) && (Buff[0] == '#'))
		fgets(Buff, READ_BUF_DIM, fp);
	sscanf(Buff, "%s", Val);

	return strdup(Val);
}

char *ReadLine(FILE *fp) {
	static char Buff[READ_BUF_DIM];

	fgets(Buff, READ_BUF_DIM, fp);
	while (!feof(fp) && (Buff[0] == '#'))
		fgets(Buff, READ_BUF_DIM, fp);

	return strdup(Buff);
}

/************ Support functions ***************/

int findPlace(Model *M, char *id) {
	int i;

	for (i = 0; i < M->N_places; i++) {
		if (strcmp(M->places[i].id, id) == 0) {
			return i;
		}
	}
	//printf("Error: place %s not found\n\n", id);
	//exit(0);
	return -1;
}

int findTransition(Model *M, char *id) {
	int i;

	for (i = 0; i < M->N_transitions; i++) {
		if (strcmp(M->transitions[i].id, id) == 0) {
			return i;
		}
	}
	//printf("Error: transition %s not found\n\n", id);
	//exit(0);
	return -1;
}

int findArc(Model *M, char *id) {
    int i;
    
    for (i = 0; i < M->N_arcs; i++) {
        if (strcmp(M->arcs[i].id, id) == 0) {
            return i;
        }
    }
    //printf("Error: transition %s not found\n\n", id);
    //exit(0);
    return -1;
}

/*********** Read Model *********************/

#define ROOT_NODE_NAME "HPnG"
#define PLACES_NODE_NAME "places"
#define TRANSITIONS_NODE_NAME "transitions"
#define ARCS_NODE_NAME "arcs"

#define DISCRETE_PLACE_NAME "discretePlace"
#define CONTINUOUS_PLACE_NAME "continuousPlace"

#define FLUID_TRANSITION_NAME "fluidTransition"
#define DETERMINISTIC_TRANISTION_NAME "deterministicTransition"
#define IMMEDIATE_TRANSITION_NAME "immediateTransition"
#define DYNAMIC_TRANISITION_NAME "dynamicContinuousTransition"
#define PTT_TRANSTION_NAME "pttTransition"
#define GENERAL_TRANSITION_NAME "generalTransition"

#define DISCRETE_ARC_NAME "discreteArc"
#define CONTINUOUS_ARC_NAME "continuousArc"
#define GUARD_ARC_NAME "guardArc"

#define TT_PTT 5

#define MAX_ID_LEN 2048

// Function to set a place values. The type of the place must be set before.
void init_place(Model *M, xml_node<> *place_node, int index)
{
    Place *place = &M->places[index];
    // Set default values.
    place->d_mark = 0;
    place->f_level = 0;
    place->f_bound = 0;

    place->id = strdup(place_node->first_attribute("id")->value());


    if (((string)place_node->name()).compare(DISCRETE_PLACE_NAME) == 0) {
        place->type = PT_DISCRETE;
        place->d_mark = atoi(place_node->first_attribute("marking")->value());
        place->idInMarking = M->N_discretePlaces;
        M->N_discretePlaces++;
    } else {
        place->type = PT_FLUID;
        place->f_bound = (float) atof(place_node->first_attribute("capacity")->value());
        place->f_level = (float) atof(place_node->first_attribute("level")->value());
        place->idInMarking = M->N_fluidPlaces;
        M->N_fluidPlaces++;
    }
}

void init_transition(Model *M, xml_node<> *transition_node, int index)
{
    Transition *transition = &M->transitions[index];
    // Set default values.
    transition->time = 0;
    transition->weight = 0;
    transition->priority = 0;
    transition->flowRate = 0.0;
    transition->distr = strdup("na");

    transition->id = strdup(transition_node->first_attribute("id")->value());
    transition->priority = atoi(transition_node->first_attribute("priority")->value());
    transition->weight = (float) atof(transition_node->first_attribute("weight")->value());
    string transition_type = (string)transition_node->name();
    if (transition_type.compare(FLUID_TRANSITION_NAME) == 0) {
        transition->idInMarking = M->N_fluidTransitions;
        M->N_fluidTransitions += 1;
        transition->type = TT_FLUID;
        transition->flowRate = (float) atof(transition_node->first_attribute("rate")->value());
    } else if (transition_type.compare(DETERMINISTIC_TRANISTION_NAME) == 0) {
        transition->idInMarking = M->N_determTransitions;
        M->N_determTransitions += 1;
        transition->type = TT_DETERMINISTIC;
        transition->time = (float) atof(transition_node->first_attribute("discTime")->value());
    } else if (transition_type.compare(GENERAL_TRANSITION_NAME) == 0) {
        transition->idInMarking = M->N_generalTransitions;
        M->N_generalTransitions += 1;
        transition->type = TT_GENERAL;
        free(transition->distr);
        transition->distr = strdup(transition_node->first_attribute("cdf")->value());
        transition->time = (float) atof(transition_node->first_attribute("discTime")->value());
    } else if (transition_type.compare(DYNAMIC_TRANISITION_NAME) == 0) {
        transition->idInMarking = M->N_DynamicTransitions;
        M->N_DynamicTransitions += 1;
        transition->type = TT_FLUID_DYNA;
        // Read the additional data for dynamic transitions.
        int i = 0;
        int depNum = 0;
        for (xml_node<> *pid = transition_node->first_node(); pid; pid = pid->next_sibling()) {
        	depNum++;
        }
        if (depNum > 0)
        {
        	transition->dependencyList = new int[depNum];
        	transition->coef = new double[depNum];
        }
        for (xml_node<> *pid = transition_node->first_node(); pid; pid = pid->next_sibling()) {
        	transition->dependencyList[i] = findTransition(M, pid->value());
        	transition->coef[i] = (float) atof(pid->first_attribute("coef")->value());
            i++;
        }
        transition->dependencySize = i;

        transition->flowRate = 0.0;
        double flow=0;
        for (int ii = 0; ii < M->transitions[index].dependencySize; ii++) {
        	flow = flow + (float)transition->coef[ii] * (float)M->transitions[transition->dependencyList[ii]].flowRate;
        }

        transition->flowRate = (flow > 0) ? flow : 0;

    } else if (transition_type.compare(IMMEDIATE_TRANSITION_NAME) == 0) {
        transition->type = TT_IMMEDIATE;
        transition->weight = (double)atof(transition_node->first_attribute("weight")->value());
        transition->priority = (double) atoi(transition_node->first_attribute("priority")->value());
    } else {
        transition->type = TT_PTT;
    }
}

void init_arc(Model *M, xml_node<> *arc_node, int index)
{
    Arc *arc = &M->arcs[index];
    char* fromId = arc_node->first_attribute("fromNode")->value();
    char* toId = arc_node->first_attribute("toNode")->value();
    char isInput = 1;
    int placeId = findPlace(M, fromId);
    int transitionId = findTransition(M, toId);
    if (placeId == -1) {
        placeId = findPlace(M, toId);
        transitionId = findTransition(M, fromId);
        isInput = 0;
    }
    if (((string)arc_node->name()).compare(DISCRETE_ARC_NAME) == 0) {
        if (isInput) {
            arc->type = AT_DISCRETE_INPUT;
        } else {
            arc->type = AT_DISCRETE_OUTPUT;
        }
    } else if(((string)arc_node->name()).compare(CONTINUOUS_ARC_NAME) == 0) {
        if (isInput) {
            arc->type = AT_FLUID_INPUT;
        } else {
            arc->type = AT_FLUID_OUTPUT;
        }
    } else {
        arc->type = AT_TEST;
        int inhibitor = atoi(arc_node->first_attribute("isInhibitor")->value());
        if(inhibitor == 1)
        	arc->type = AT_INHIBITOR;
    }

    arc->id = arc_node->first_attribute("id")->value();
    arc->priority = atoi(arc_node->first_attribute("priority")->value());
    arc->weight = (double) atof(arc_node->first_attribute("weight")->value());
    arc->share = (double) atof(arc_node->first_attribute("share")->value());
    if (isInput || arc->type == AT_TEST || arc->type == AT_INHIBITOR) {
        arc->fromId = placeId;
        arc->toId = transitionId;
    } else {
        arc->fromId = transitionId;
        arc->toId = placeId;
    }
    arc->transId = transitionId;
    arc->placeId = placeId;

}

Model *ReadModel(const char *FileName) {
    Model *M;
    int i;

    xml_document<> doc;
    xml_node<> * root_node;
    // Read the xml file into a vector
    cout << "Reading model " << FileName << endl;
    ifstream file(FileName);
    if (!file) {
        cout << "\n\n Error: cannot find model file " << FileName << endl << endl;
        //exit(0);
    }

    M = (Model *) malloc(sizeof(Model));

    M->N_discretePlaces = 0;
    M->N_fluidPlaces = 0;
    M->N_determTransitions = 0;
    M->N_fluidTransitions = 0;
    M->N_generalTransitions = 0;
    M->N_DynamicTransitions = 0;

    vector<char> buffer((istreambuf_iterator<char>(file)), istreambuf_iterator<char>());
    buffer.push_back('\0');
    // Parse the buffer using the xml file parsing library into doc
    doc.parse<0>(&buffer[0]);
    // Find our root node
    root_node = doc.first_node(ROOT_NODE_NAME);

    // Read all places
    xml_node<> *places_node;
    places_node = root_node->first_node(PLACES_NODE_NAME);
    int number_of_places = 0;
    i = 0;
    M->places = NULL;
    for (xml_node<> *place = places_node->first_node(); place; place = place->next_sibling()) {
        number_of_places++;
        M->places = (Place *)realloc(M->places, sizeof(Place)*number_of_places);
        init_place(M, place, i);

        printf("Place %d: %s %s %d %g %g\n", i,
               (M->places[i].type == 0 ? "disc." : "fluid"),
               M->places[i].id, M->places[i].d_mark, M->places[i].f_level,
               M->places[i].f_bound);
        i++;
    }
    M->N_places = number_of_places;
    printf("Number of places:                                        : %d\n\n",
           M->N_places);

    //Read all transitions
    xml_node<> *transitions_node;
    transitions_node = root_node->first_node(TRANSITIONS_NODE_NAME);
    int number_of_transitions = 0;
    i = 0;
    M->transitions = NULL;
    for (xml_node<> *transition = transitions_node->first_node(); transition; transition = transition->next_sibling()) {
        number_of_transitions++;
        M->N_transitions = number_of_transitions;
        M->transitions = (Transition *)realloc(M->transitions, sizeof(Transition)*number_of_transitions);
        init_transition(M, transition, i);

        printf("Transition %d: %s %s %lg %lg %d %lg %s\n", i,
               transTypeName[M->transitions[i].type], M->transitions[i].id,
               M->transitions[i].time, M->transitions[i].weight,
               M->transitions[i].priority, M->transitions[i].flowRate,
               M->transitions[i].distr);
        if (M->transitions[i].type == TT_FLUID_DYNA) {
            for (int ii = 0; ii < M->transitions[i].dependencySize; ii++) {
                printf("   Id: %d  Coef: %g\n", M->transitions[i].dependencyList[ii], M->transitions[i].coef[ii]);
            }
        } else if (M->transitions[i].type == TT_PTT) {
            // Currently not implemented.
        }

        i++;
    }
    printf("Number of transitions:                                   : %d\n\n",
           M->N_transitions);
    for (int i = 0; i < M->N_transitions; i++) {
        if (M->transitions[i].type == TT_FLUID_DYNA)
            M->transitions[i].idInMarking += M->N_fluidTransitions;
    }

    //Read all arcs
    xml_node<> *arcs_node;
    arcs_node = root_node->first_node(ARCS_NODE_NAME);
    i = 0;
    M->arcs = NULL;
    M->N_arcs = 0;
    for (xml_node<> *arc = arcs_node->first_node(); arc; arc = arc->next_sibling()) {
        M->N_arcs++;
        M->arcs = (Arc *)realloc(M->arcs, sizeof(Arc) * M->N_arcs);
        init_arc(M, arc, i);
        printf("Arc %d: %s %s %d %d %lg %lg %d\n", i,
               arcsTypeName[M->arcs[i].type], M->arcs[i].id,
               M->arcs[i].fromId, M->arcs[i].toId,
               M->arcs[i].weight, M->arcs[i].share, M->arcs[i].priority);
        i++;
    }
    printf("Number of arcs:                                          : %d\n\n",
           M->N_arcs);
    return M;
}


void InitializeModel(Model *M) {
	int i, j, k;
	int tmp;
	for (i = 0; i < M->N_transitions; i++) {
		M->transitions[i].inputListSize = 0;
		M->transitions[i].inhibListSize = 0;
		M->transitions[i].testListSize = 0;
		M->transitions[i].outputListSize = 0;
	}
	for (i = 0; i < M->N_places; i++) {
		M->places[i].inputListSize = 0;
		M->places[i].inhibListSize = 0;
		M->places[i].testListSize = 0;
		M->places[i].outputListSize = 0;
	}

	for (i = 0; i < M->N_arcs; i++) {
		// Depending on the arc type, it increases the counter (of the corresponding transition)
		switch (M->arcs[i].type) {
		case AT_DISCRETE_OUTPUT:
		case AT_FLUID_OUTPUT:
			// check which transition it is
			M->transitions[M->arcs[i].fromId].outputListSize += 1;
			M->places[M->arcs[i].toId].inputListSize += 1;
			break;
		case AT_DISCRETE_INPUT:
		case AT_FLUID_INPUT:
			M->transitions[M->arcs[i].toId].inputListSize += 1;
			M->places[M->arcs[i].fromId].outputListSize += 1;
			break;
		case AT_INHIBITOR:
			M->transitions[M->arcs[i].toId].inhibListSize += 1;
			M->places[M->arcs[i].fromId].inhibListSize += 1;
			break;
		case AT_TEST:
			M->transitions[M->arcs[i].toId].testListSize += 1;
			M->places[M->arcs[i].fromId].testListSize += 1;
			break;
		default:
			printf("error");
			break;

		}
	}

	for (i = 0; i < M->N_transitions; i++) {
		// Now that we know the sizes, we can allocate the memory for the lists
		M->transitions[i].inputList = (int *) calloc(
				M->transitions[i].inputListSize, sizeof(int));
		M->transitions[i].inhibList = (int *) calloc(
				M->transitions[i].inhibListSize, sizeof(int));
		M->transitions[i].testList = (int *) calloc(
				M->transitions[i].testListSize, sizeof(int));
		M->transitions[i].outputList = (int *) calloc(
				M->transitions[i].outputListSize, sizeof(int));

		//TODO: horrible variable change!!!!
		M->transitions[i].inputListSize = 0; // now this variable becomes the index of the next element of the list
		M->transitions[i].outputListSize = 0;
		M->transitions[i].inhibListSize = 0;
		M->transitions[i].testListSize = 0;
	}

	for (i = 0; i < M->N_places; i++) {
		M->places[i].inputList = (int *) calloc(M->places[i].inputListSize,
				sizeof(int));
		M->places[i].inhibList = (int *) calloc(M->places[i].inhibListSize,
				sizeof(int));
		M->places[i].testList = (int *) calloc(M->places[i].testListSize,
				sizeof(int));
		M->places[i].outputList = (int *) calloc(M->places[i].outputListSize,
				sizeof(int));
		M->places[i].inputListSize = 0;
		M->places[i].inhibListSize = 0;
		M->places[i].testListSize = 0;
		M->places[i].outputListSize = 0;
	}

	for (i = 0; i < M->N_arcs; i++) {
		// Depending on the arc type, it puts the elements in the array

		switch (M->arcs[i].type) {
		case AT_DISCRETE_OUTPUT:
		case AT_FLUID_OUTPUT:
			// check which transition it is
			M->transitions[M->arcs[i].fromId].outputList[M->transitions[M->arcs[i].fromId].outputListSize] =
					i;
			M->transitions[M->arcs[i].fromId].outputListSize += 1;
			M->places[M->arcs[i].toId].inputList[M->places[M->arcs[i].toId].inputListSize] =
					i;
			M->places[M->arcs[i].toId].inputListSize += 1;
			break;
		case AT_DISCRETE_INPUT:
		case AT_FLUID_INPUT:
			M->transitions[M->arcs[i].toId].inputList[M->transitions[M->arcs[i].toId].inputListSize] =
					i;
			M->transitions[M->arcs[i].toId].inputListSize += 1;
			M->places[M->arcs[i].fromId].outputList[M->places[M->arcs[i].fromId].outputListSize] =
					i;
			M->places[M->arcs[i].fromId].outputListSize += 1;
			break;
		case AT_INHIBITOR:
			M->transitions[M->arcs[i].toId].inhibList[M->transitions[M->arcs[i].toId].inhibListSize] =
					i;
			M->transitions[M->arcs[i].toId].inhibListSize += 1;
			M->places[M->arcs[i].fromId].inhibList[M->places[M->arcs[i].fromId].inhibListSize] =
					i;
			M->places[M->arcs[i].fromId].inhibListSize += 1;
			break;
		case AT_TEST:
			M->transitions[M->arcs[i].toId].testList[M->transitions[M->arcs[i].toId].testListSize] =
					i;
			M->transitions[M->arcs[i].toId].testListSize += 1;
			M->places[M->arcs[i].fromId].testList[M->places[M->arcs[i].fromId].testListSize] =
					i;
			M->places[M->arcs[i].fromId].testListSize += 1;
			break;
		default:
			printf("error");
			break;
		}
	}

	///// Sorts arcs by priority and virtual saturation time
	for (i = 0; i < M->N_places; i++) {
		if (M->places[i].type == PT_FLUID) {
			// Only for fluid places

			// input bag:
			for (j = 0; j < M->places[i].inputListSize - 1; j++) {
				for (k = j + 1; k < M->places[i].inputListSize; k++) {
					if ((M->arcs[M->places[i].inputList[k]].priority
							< M->arcs[M->places[i].inputList[j]].priority)
							|| ((M->arcs[M->places[i].inputList[k]].priority
									== M->arcs[M->places[i].inputList[j]].priority)
									&& (M->arcs[M->places[i].inputList[k]].share
											> M->arcs[M->places[i].inputList[j]].share))) {
						//printf("Is");
						tmp = M->places[i].inputList[k];
						M->places[i].inputList[k] = M->places[i].inputList[j];
						M->places[i].inputList[j] = tmp;
					}
				}
			}

			// output bag:
			for (j = 0; j < M->places[i].outputListSize - 1; j++) {
				for (k = j + 1; k < M->places[i].outputListSize; k++) {
					if ((M->arcs[M->places[i].outputList[k]].priority
							< M->arcs[M->places[i].outputList[j]].priority)
							|| ((M->arcs[M->places[i].outputList[k]].priority
									== M->arcs[M->places[i].outputList[j]].priority)
									&& (M->arcs[M->places[i].outputList[k]].share
											> M->arcs[M->places[i].outputList[j]].share))) {
						//printf("Os [%d]%d [%d]%d,", j, M->places[i].outputList[j], k, M->places[i].outputList[k]);
						tmp = M->places[i].outputList[k];
						M->places[i].outputList[k] = M->places[i].outputList[j];
						M->places[i].outputList[j] = tmp;
					}
				}
			}

			printf("\n Place %s: \n input  -> ", M->places[i].id);
			for (j = 0; j < M->places[i].inputListSize; j++) {
				printf("%d=P:%d %s (%lg), ", M->places[i].inputList[j],
						M->arcs[M->places[i].inputList[j]].priority,
						M->transitions[M->arcs[M->places[i].inputList[j]].fromId].id,
						M->arcs[M->places[i].inputList[j]].share);
			}
			printf("\n output -> ");
			for (j = 0; j < M->places[i].outputListSize; j++) {
				printf("%d=P:%d %s (%lg), ", M->places[i].outputList[j],
						M->arcs[M->places[i].outputList[j]].priority,
						M->transitions[M->arcs[M->places[i].outputList[j]].toId].id,
						M->arcs[M->places[i].outputList[j]].share);
			}
			printf("\n\n");
		}
	}
	/*
	 // read out the transition bags for every transition
	 for(i = 0; i < M->N_transitions; i++) {
	 printf("Transition %d: inputListSize %d, arcs: ",
	 i,
	 M->transitions[i].inputListSize);

	 for(j =0; j < M->transitions[i].inputListSize; j++){

	 printf(" %d ",
	 M->transitions[i].inputList[j]);
	 }
	 printf("\n");
	 }
	 printf("\n");
	 for(i = 0; i < M->N_transitions; i++) {
	 printf("Transition %d: outputListSize %d, arcs: ",
	 i,
	 M->transitions[i].outputListSize);

	 for(j =0; j < M->transitions[i].outputListSize; j++){

	 printf(" %d ",
	 M->transitions[i].outputList[j]);
	 }
	 printf("\n");
	 }
	 printf("\n");
	 for(i = 0; i < M->N_transitions; i++) {
	 printf("Transition %d: inhibListSize %d, arcs: ",
	 i,
	 M->transitions[i].inhibListSize);

	 for(j =0; j < M->transitions[i].inhibListSize; j++){

	 printf(" %d ",
	 M->transitions[i].inhibList[j]);
	 }
	 printf("\n");
	 }
	 printf("\n");
	 for(i = 0; i < M->N_transitions; i++) {
	 printf("Transition %d: testListSize %d, arcs: ",
	 i,
	 M->transitions[i].testListSize);

	 for(j =0; j < M->transitions[i].testListSize; j++){

	 printf(" %d ",
	 M->transitions[i].testList[j]);
	 }
	 printf("\n");
	 }
	 */

}

Marking *allocMarking(Model *M) {
	Marking *out;
	out = (Marking *) malloc(sizeof(Marking));
	out->tokens = (int *) calloc(M->N_discretePlaces, sizeof(int));
	out->fluid0 = (double *) calloc(M->N_fluidPlaces, sizeof(double));
	out->fluid1 = (double *) calloc(M->N_fluidPlaces, sizeof(double));
	out->enabling = (int *) calloc(M->N_transitions, sizeof(int));
	out->actFluidRate = (double *) calloc(
			M->N_fluidTransitions + M->N_DynamicTransitions, sizeof(double));
	out->clock0 = (double *) calloc(M->N_determTransitions, sizeof(double));
	out->clock1 = (double *) calloc(M->N_determTransitions, sizeof(double));
	out->generalHasFired = (int *) calloc(M->N_generalTransitions, sizeof(int));
	out->generalDisabled = (double *) calloc(M->N_generalTransitions,
			sizeof(double));

	out->fluidPlaceDeriv = (double *) calloc(M->N_fluidPlaces, sizeof(double));
	out->N_generalFired = 0;

	return out;
}

void freeMarking(Model *M, Marking *K) {
	free(K->tokens);
	free(K->fluid0);
	free(K->fluid1);
	free(K->enabling);
	free(K->actFluidRate);
	free(K->clock0);
	free(K->clock1);
	free(K->generalHasFired);
	free(K->generalDisabled);

	free(K->fluidPlaceDeriv);

	free(K);
}

Marking *copyMarking(Model *M, Marking *Src) {
	Marking *out;
	int i;

	out = allocMarking(M);

	for (i = 0; i < M->N_discretePlaces; i++)
		out->tokens[i] = Src->tokens[i];
	for (i = 0; i < M->N_fluidPlaces; i++)
		out->fluid0[i] = Src->fluid0[i];
	for (i = 0; i < M->N_fluidPlaces; i++)
		out->fluid1[i] = Src->fluid1[i];
	for (i = 0; i < M->N_fluidPlaces; i++)
		out->fluidPlaceDeriv[i] = Src->fluidPlaceDeriv[i];

	for (i = 0; i < M->N_transitions; i++)
		out->enabling[i] = Src->enabling[i];
	for (i = 0; i < M->N_determTransitions; i++)
		out->clock0[i] = Src->clock0[i];
	for (i = 0; i < M->N_determTransitions; i++)
		out->clock1[i] = Src->clock1[i];
	//for(i = 0; i < M->N_fluidTransitions; i++) out->actFluidRate[i] = Src->actFluidRate[i];
	for (i = 0; i < M->N_transitions; i++) {
		if (M->transitions[i].type == TT_FLUID
				|| M->transitions[i].type == TT_FLUID_DYNA) {
			out->actFluidRate[M->transitions[i].idInMarking] =
					Src->actFluidRate[M->transitions[i].idInMarking];
			//M->transitions[i].flowRate;
		}
	}
	for (i = 0; i < M->N_generalTransitions; i++)
		out->generalHasFired[i] = Src->generalHasFired[i];
	for (i = 0; i < M->N_generalTransitions; i++)
		out->generalDisabled[i] = Src->generalDisabled[i];

	//for(i = 0; i < M->N_fluidPlaces; i++) out->fluidPlaceDeriv[i] = Src->fluidPlaceDeriv[i];
	out->N_generalFired = Src->N_generalFired;

	return out;
}

Marking *createInitialMarking(Model *M) {
	int i;
	Marking *Mrk;

	Mrk = allocMarking(M); // allocate the initial marking
	// for all discrete places
	for (i = 0; i < M->N_places; i++) {
		if (M->places[i].type == 0) //discrete place
				{
			Mrk->tokens[M->places[i].idInMarking] = M->places[i].d_mark;
		} else //fluid place
		{
			Mrk->fluid0[M->places[i].idInMarking] = M->places[i].f_level;
			Mrk->fluid1[M->places[i].idInMarking] = 0.0;
		}
	}

	// set the enabling for all transitions to zero
	// for fluid transitions set actFluidRate to fluid rate
	for (i = 0; i < M->N_transitions; i++) {
		Mrk->enabling[i] = 0;
		if (M->transitions[i].type == TT_FLUID
				|| M->transitions[i].type == TT_FLUID_DYNA) {
			Mrk->actFluidRate[M->transitions[i].idInMarking] =
					M->transitions[i].flowRate;
		}
		if (M->transitions[i].type == TT_DETERMINISTIC) {
			Mrk->clock0[M->transitions[i].idInMarking] = 0.0;
			Mrk->clock1[M->transitions[i].idInMarking] = 0.0;
		}
		if (M->transitions[i].type == TT_GENERAL) {
			Mrk->generalHasFired[M->transitions[i].idInMarking] = 0;
		}
	}

	return Mrk;
}

void checkEnabled(Model *M, Marking *K, double s_min, double s_max) {
	int i, j;
	int enable;

	// fills the K->enabling field for all the transitions
	K->N_generalFired = 0.0;

	for (i = 0; i < M->N_transitions; i++) {

		enable = 1;


//		printf("%s\n", M->transitions[i].id);
		// deterministic transition------Hamed: Anything except Fluid !?
		if (M->transitions[i].type != TT_FLUID && M->transitions[i].type != TT_FLUID_DYNA) {
			// check all input places


			for (j = 0; j < M->transitions[i].inputListSize; j++) {
				// disable if the input place has at least as many tokens in the current marking as the weight of the corresponding arc
				// Hamed: If at least exists an arc with weight greater than input place marking.


//				printf("place: %s, arc weight = %f\n",
//						M->places[M->arcs[M->transitions[i].inputList[j]].placeId].id, M->arcs[M->transitions[i].inputList[j]].weight);
				if (K->tokens[M->places[M->arcs[M->transitions[i].inputList[j]].placeId].idInMarking]
						< M->arcs[M->transitions[i].inputList[j]].weight) {
					enable = 0;
					//printf(
					//		"disable transition %d due to inputarc from place %d \n",
					//		i, M->arcs[M->transitions[i].inputList[j]].fromId);
					break;

				}
			}
		}



		// check all test arcs
		if (enable == 1) {
			for (j = 0; j < M->transitions[i].testListSize; j++) {
				// disable if one test condition is not valid
				// if the place is discrete
				if (M->places[M->arcs[M->transitions[i].testList[j]].placeId].type == PT_DISCRETE) {
					if (K->tokens[M->places[M->arcs[M->transitions[i].testList[j]].placeId].idInMarking]
							< M->arcs[M->transitions[i].testList[j]].weight) {
						enable = 0;
						//printf(
						//		"disable transition %d due to testarc from place %d \n",
						//		i,
						//		M->arcs[M->transitions[i].testList[j]].fromId);
						break;
					}
				} else {
					float fContent =
							K->fluid0[M->places[M->arcs[M->transitions[i].testList[j]].placeId].idInMarking]
									+ (s_max + s_min) / 2
											* K->fluid1[M->places[M->arcs[M->transitions[i].testList[j]].placeId].idInMarking];

//					printf("fcontnet: %f\n", fContent);
//					printf("weight %f\n ", M->arcs[M->transitions[i].testList[j]].weight);
//					printf("fluidPlaceDeriv:%f\n", K->fluidPlaceDeriv[M->places[M->arcs[M->transitions[i].testList[j]].placeId].idInMarking]);
//
//
//					if (!IS_ZERO(M->arcs[M->transitions[i].testList[j]].weight- fContent))
//						printf("1\n");
//					if (fContent < M->arcs[M->transitions[i].testList[j]].weight)
//						printf("2\n");
//					if (IS_ZERO(M->arcs[M->transitions[i].testList[j]].weight- fContent))
//						printf("3\n");
//					if (K->fluidPlaceDeriv[M->places[M->arcs[M->transitions[i].testList[j]].placeId].idInMarking] < 0)
//						printf("4\n");

					//the second condition considers if we have entered a new region because of a guard arc. In this case the place content is
					// equal to the weight of the arc and if the (previous) derivative is negative the place should be disabled,
					//because we had a guard arc event and the place is about to violate the condition.
					//(note that in this case the the content of the place is independent of parameter 's').
					// For the FIRST condition we have to check that the if fContent is not equal to the arc weight (for percition issues)!
					if ((!IS_ZERO(M->arcs[M->transitions[i].testList[j]].weight- fContent)
							&& fContent < M->arcs[M->transitions[i].testList[j]].weight)
							|| (IS_ZERO(M->arcs[M->transitions[i].testList[j]].weight- fContent)
									&& K->fluidPlaceDeriv[M->places[M->arcs[M->transitions[i].testList[j]].placeId].idInMarking]
											< 0)) {
						enable = 0;
						break;
					}
				}
			}
		}

		// check all inhib arcs

		if (enable == 1) {
			for (j = 0; j < M->transitions[i].inhibListSize; j++) {
				// disable if one inhibitor condition is not valid
				// if the place is discrete
				if (M->places[M->arcs[M->transitions[i].inhibList[j]].placeId].type
						== PT_DISCRETE) {
					if (K->tokens[M->places[M->arcs[M->transitions[i].inhibList[j]].placeId].idInMarking]
							>= M->arcs[M->transitions[i].inhibList[j]].weight) {
						enable = 0;
						//printf(
						//		"disable transition %d due to inhibarc from place %d \n",
						//		i,
						//		M->arcs[M->transitions[i].inhibList[j]].fromId);
						break;
					}
				} else {
					float fContent =
							K->fluid0[M->places[M->arcs[M->transitions[i].inhibList[j]].placeId].idInMarking]
									+ (s_max + s_min) / 2
											* K->fluid1[M->places[M->arcs[M->transitions[i].inhibList[j]].placeId].idInMarking];
					if ((!IS_ZERO(M->arcs[M->transitions[i].inhibList[j]].weight- fContent)
							&& fContent > M->arcs[M->transitions[i].inhibList[j]].weight)
							|| (IS_ZERO( fContent - M->arcs[M->transitions[i].inhibList[j]].weight)
									&& K->fluidPlaceDeriv[M->places[M->arcs[M->transitions[i].inhibList[j]].placeId].idInMarking] >= 0)) {
						enable = 0;
						break;
					}
				}
			}
		}
		K->enabling[i] = enable;

		if (M->transitions[i].type == TT_GENERAL) {
			if (K->generalHasFired[M->transitions[i].idInMarking]) {
				K->N_generalFired++;
			}
		}

		// printf("transition %d enabled %d \n",i,enable);
	}

}

void ShareFlow(Model *M, Marking *K, double Flux, int *_arcs, int N_arcs) {
	int i, j;
	double Dsum;
	int currentPri;
	double fluxReq, sharedFlux;

	i = 0;
	while (i < N_arcs) {
		//printf("Entering share at %d (Flux = %lg)\n", i, Flux);
		while ((K->enabling[M->arcs[_arcs[i]].transId] == 0) && (i < N_arcs)) {
			i++; // skip all the transitions that are not enabled
		}
		//printf("Skipped disabled until %d\n", i);

		// first arc at this priority
		currentPri = M->arcs[_arcs[i]].priority;
//		printf("%d\n", currentPri);

		fluxReq = M->transitions[M->arcs[_arcs[i]].transId].flowRate;
		Dsum = M->arcs[_arcs[i]].share
				* M->transitions[M->arcs[_arcs[i]].transId].flowRate;
		j = i + 1;
		while ((j < N_arcs) && (M->arcs[_arcs[j]].priority == currentPri)) { // look for other enabled transitions at the same priority level
			//TODO: (very IMPORTANT) the "i" here should be "j".
			//commented by HAMED.
			//if (K->enabling[M->arcs[_arcs[i]].transId] == 1) { // consider only the enabled transitions
			if (K->enabling[M->arcs[_arcs[j]].transId] == 1) { // consider only the enabled transitions
				fluxReq += M->transitions[M->arcs[_arcs[j]].transId].flowRate;
				Dsum += M->arcs[_arcs[j]].share
						* M->transitions[M->arcs[_arcs[j]].transId].flowRate;
			}
			j++;
		}
		//printf("Found transitions with same priority until %d\n", j);

		if (Flux > fluxReq) { // Here we have enough for this priority level
			//printf("Enough flux to satisfy priority level\n");
			Flux -= fluxReq; // We give every transition what it requires
			i = j; // j is the first element with the next lower priority
		} else {
			if (j == i + 1) { // there is only one transition in this priority level
			//printf("Only one transition gets the whole share %g to %d\n", Flux, M->transitions[M->arcs[_arcs[i]].transId].idInMarking);
			// We give what we have
				K->actFluidRate[M->transitions[M->arcs[_arcs[i]].transId].idInMarking] =
						Flux;
				i = j;
				break;
			} else { // There are more than one transition in this priority level
				//printf("Starting priority sharing algorithm\n");
				j = i;
				while ((j < N_arcs)
						&& (M->arcs[_arcs[j]].priority == currentPri)) { // look for other enabled transitions at the same priority level
					//TODO: (very IMPORTANT) the "i" here should be "j".
					//commented by HAMED
					//if (K->enabling[M->arcs[_arcs[i]].transId] == 1) { // consider only the enabled transitions
					if (K->enabling[M->arcs[_arcs[j]].transId] == 1) { // consider only the enabled transitions
						sharedFlux =
								Flux * M->arcs[_arcs[j]].share
										* M->transitions[M->arcs[_arcs[j]].transId].flowRate
										/ Dsum;
						if (sharedFlux
								> M->transitions[M->arcs[_arcs[j]].transId].flowRate) {
							Flux -=
									M->transitions[M->arcs[_arcs[j]].transId].flowRate;
							Dsum -=
									M->transitions[M->arcs[_arcs[j]].transId].flowRate
											* M->arcs[_arcs[j]].share;
						} else {
							K->actFluidRate[M->transitions[M->arcs[_arcs[j]].transId].idInMarking] =
									sharedFlux;
						}
					}
					j++;
				}
				i = j;
				break;
			}
		}
	}
	while (i < N_arcs) { // Put to 0 all the remaining transitions
		//TODO: I think we should do this only to enabled transitions. However, I don't know if there is a significant difference.
		K->actFluidRate[M->transitions[M->arcs[_arcs[i]].transId].idInMarking] =
				0.0;
		i++;
	}
}

void updateDynamicRates(Model *M, Marking *K, int *_arcs, int N_arcs){
	//I know this piece of code is really nasty, by there is no escape by this stupid design, that I have in hand.
	//we iterate over all transitions that had a rate change. if they were NOT dynamic,we change the rate of all dynamic transitions which were dependent on them.
		for (int i = 0; i < N_arcs; i++) {
			int stId = M->arcs[_arcs[i]].transId;
			if (M->transitions[stId].type == TT_FLUID) {
				for (int d = 0; d < M->N_transitions; d++){
					if (M->transitions[d].type == TT_FLUID_DYNA) {
						bool updateRate = false;
						for (int jj = 0; jj < M->transitions[d].dependencySize; jj++){
							if (M->transitions[d].dependencyList[jj] == stId){
								updateRate = true;
								break;
							}
						}
						if (updateRate){
							M->transitions[d].flowRate = 0;
							for (int jj = 0; jj < M->transitions[d].dependencySize; jj++)
								if (K->enabling[M->transitions[d].dependencyList[jj]])
									M->transitions[d].flowRate += (float)M->transitions[d].coef[jj]* (float) M->transitions[M->transitions[d].dependencyList[jj]].flowRate;

							M->transitions[d].flowRate = M->transitions[d].flowRate >= 0 ? M->transitions[d].flowRate : 0;
							K->actFluidRate[M->transitions[d].idInMarking] = M->transitions[d].flowRate;
						}
					}
				}
			}
		}
}

#define MAX_F_CONF_SOL_ITER  5000
/**
 * Hamed: In this function I see a significant difference between implementation and what is described in the paper!!!
 **/
void setActFluidRate(Model *M, Marking *K, double s0) {

	static int callNum = 0;
//	printf("setActFluidRate # call: %d \n", callNum++);
//	printFluidLevels(M, K);
//	printFluidRates(M, K);

	//	static double *lastIterVal = NULL;
	int i, j, trId;
	double inFlux, outFlux;
	// for all places
	int noChange = 0;
	int nIter = 0;

	/*
	 ****************************************************************************************************
	 * I am pessimistic about this algorithm. investigate weather this stuck in an infinite loop or not.*
	 ****************************************************************************************************
	 */

	//Assuming nominal rates for all  the fluid transitions.
	//TODO:  very important! is this correct?
	for (int d = 0; d < M->N_transitions; d++) {
		if (M->transitions[d].type == TT_FLUID) {
			if (K->enabling[d] > 0 )
				K->actFluidRate[M->transitions[d].idInMarking] = M->transitions[d].flowRate;
			else
				K->actFluidRate[M->transitions[d].idInMarking] = 0;
		}
	}

	for (int d = 0; d < M->N_transitions; d++) {
		if (M->transitions[d].type == TT_FLUID_DYNA) {
			M->transitions[d].flowRate = 0;
			for (int jj = 0; jj < M->transitions[d].dependencySize; jj++){
				if (K->enabling[M->transitions[d].dependencyList[jj]]){
					M->transitions[d].flowRate +=
						(float)M->transitions[d].coef[jj]
								* (float)M->transitions[M->transitions[d].dependencyList[jj]].flowRate;
				}
			}
			M->transitions[d].flowRate  = (M->transitions[d].flowRate > 0) ?M->transitions[d].flowRate:0 ;

			if (K->enabling[d] > 0 ){
				K->actFluidRate[M->transitions[d].idInMarking] = M->transitions[d].flowRate;
			}
			else
				K->actFluidRate[M->transitions[d].idInMarking] = 0;
		}
	}


//	printf("-----------------in setActualFlowRate------------------\n");
//
//	printf("***At the start: \n");
//	printFluidRates(M, K);


	while (noChange == 0) {
		noChange = 1;
		nIter++;

		for (i = 0; i < M->N_places; i++) {
			// check all arcs whether they are fluid input / ouput arcs to that place
			if (M->places[i].type == PT_FLUID) {
				inFlux = outFlux = 0.0;


				for (j = 0; j < M->places[i].inputListSize; j++) {
					trId = M->arcs[M->places[i].inputList[j]].fromId;
					if (K->enabling[trId] > 0) {
						inFlux +=
								K->actFluidRate[M->transitions[trId].idInMarking];
					}
				}

				for (j = 0; j < M->places[i].outputListSize; j++) {
					trId = M->arcs[M->places[i].outputList[j]].toId;
					if (K->enabling[trId] > 0) {
						outFlux +=
								K->actFluidRate[M->transitions[trId].idInMarking];
					}
				}

				K->fluidPlaceDeriv[M->places[i].idInMarking] = inFlux - outFlux;

//				printf("%s(%f + %fs)-->%f: in Flux: %f, outFlux: %f, drift: %f\n", M->places[i].id,
//						K->fluid0[M->places[i].idInMarking], K->fluid1[M->places[i].idInMarking], M->places[i].f_bound,
//						inFlux, outFlux, inFlux - outFlux);

				if ((inFlux - outFlux > ZERO_PREC)
						/*
						 * Note that this equation is independent of any value of s \in [s_min, s_max].
						 * because, this equation will be zero only if we have entered this new region,
						 * because of this fluid place reaching its upper/lower boundary.
						 */
						&& (IS_ZERO(K->fluid0[M->places[i].idInMarking] + s0* K->fluid1[M->places[i].idInMarking]
										- M->places[i].f_bound))
						&& (K->fluid1[M->places[i].idInMarking] >= -ZERO_PREC)) {

					if (nIter> 2000){
						printf("sth is wrong here!!!!\n");
						printFluidLevels(M, K);
						printFluidRates(M, K);
					}
					// Here, since inFlux > outFlux, we know that some transition will not receive its full amount of fluid
					//					printf("\n\n Problem at the upper bound at place %s \n\n", M->places[i].id);
					ShareFlow(M, K, outFlux, M->places[i].inputList, M->places[i].inputListSize);
//					printFluidRates(M, K);

					updateDynamicRates(M, K, M->places[i].inputList, M->places[i].inputListSize);
//					printFluidRates(M, K);

					K->fluidPlaceDeriv[M->places[i].idInMarking] = 0;
					noChange = 0;
				}

				if ((inFlux - outFlux < -ZERO_PREC)
						&& (IS_ZERO(K->fluid0[M->places[i].idInMarking] + s0* K->fluid1[M->places[i].idInMarking]))
						&& (K->fluid1[M->places[i].idInMarking] <= ZERO_PREC)) {
					//					printf("\n\n Problem at the lower bound at place %s \n\n", M->places[i].id);
					ShareFlow(M, K, inFlux, M->places[i].outputList, M->places[i].outputListSize);

					updateDynamicRates(M, K, M->places[i].outputList, M->places[i].outputListSize);
					K->fluidPlaceDeriv[M->places[i].idInMarking] = 0;
					noChange = 0;
				}


//
//				printf("----------------------------\n");
//				printFluidRates(M, K);
//

			}
		}

		if (nIter > MAX_F_CONF_SOL_ITER) {
			printf("\n ERROR: Loop in conflict resolution\n\n");
			//exit(0);
		}
	}

//	printf("***At the END: \n");
//	for (int t = 0; t < M->N_transitions; t++){
//		if (M->transitions[t].type == TT_FLUID || M->transitions[t].type == TT_FLUID_DYNA){
//				printf("%s(%d) --> Nominal rate = %f\n",  M->transitions[t].id, K->enabling[t], K->actFluidRate[M->transitions[t].idInMarking]);
//		}
//	}


}

StateTimeAlt *makeTimeAlt(double left, double right, StateTimeAlt *next) {
	StateTimeAlt *A;
	A = (StateTimeAlt *) malloc(sizeof(StateTimeAlt));
	A->leftInt = left;
	A->rightInt = right;
	A->next = next;
	return A;
}

StateProbAlt *makeProbAlt(State *S, double P) {
	StateProbAlt *A;

	A = (StateProbAlt *) malloc(sizeof(StateProbAlt));
	A->p = P;
	A->next = NULL;
	A->S = S;
	return A;
}

#define POS_INFINITY 1.0e100

Marking *advanceMarking(Model *M, Marking *m, double T1, double T0,
		double *drift, int *enabled) {
	Marking *NewM;
	int i;

	NewM = copyMarking(M, m);

	// Advance fluid
	for (i = 0; i < M->N_places; i++) {
		if (M->places[i].type == PT_FLUID) {
			//printf("Updating fluid level f1:%g, f0:%g, T1:%g, T0:%g\n",
			//		NewM->fluid1[M->places[i].idInMarking],
			//		NewM->fluid0[M->places[i].idInMarking], T1, T0);

				NewM->fluid0[M->places[i].idInMarking] +=
					drift[M->places[i].idInMarking] * T0;
				NewM->fluid1[M->places[i].idInMarking] +=
					drift[M->places[i].idInMarking] * T1;
			//printf("%d => %g s + %g\n", i,
			//		NewM->fluid1[M->places[i].idInMarking],
			//		NewM->fluid0[M->places[i].idInMarking]);
		}
	}
	// Advance clocks
	for (i = 0; i < M->N_transitions; i++) {
		if ((M->transitions[i].type == TT_DETERMINISTIC) && (enabled[i] == 1)) {
//			printf("Updating clock with fluid Dep.\n");
			NewM->clock0[M->transitions[i].idInMarking] += T0;
			NewM->clock1[M->transitions[i].idInMarking] += T1;
			//printf("%d => %g s + %g\n", i,
			//		NewM->clock1[M->transitions[i].idInMarking],
			//		NewM->clock0[M->transitions[i].idInMarking]);
		}
		if ((M->transitions[i].type == TT_GENERAL) && (enabled[i] == 0)
				&& (NewM->generalHasFired[M->transitions[i].idInMarking] == 0)) {
			NewM->generalDisabled[M->transitions[i].idInMarking] += T0;
			//printf(
			//		"AAAAAARRRRRRGGGGHHHHHHHHHH! This should not happen now!\n\n\n\n\n\n");
		}
	}

	return NewM;
}

void fireTransition(Model *M, Marking *NewM, int Tr) {
	int j;
	for (j = 0; j < M->transitions[Tr].inputListSize; j++) {
		NewM->tokens[M->places[M->arcs[M->transitions[Tr].inputList[j]].fromId].idInMarking] -=
				M->arcs[M->transitions[Tr].inputList[j]].weight;
		//printf("remove one token from place %s \n",  M->places[M->arcs[Tr].inputList[j]].fromId].id);
	} // add tokens to each output place according to outputarc.weight  !!!!! is this correct

	for (j = 0; j < M->transitions[Tr].outputListSize; j++) {
		NewM->tokens[M->places[M->arcs[M->transitions[Tr].outputList[j]].toId].idInMarking] +=
				M->arcs[M->transitions[Tr].outputList[j]].weight;
		//printf("add one token to place %s \n",  M->places[M->arcs[M->transitions[Tr].outputList[j]].toId].id);
	}
	if (M->transitions[Tr].type == TT_DETERMINISTIC) {
		NewM->clock0[M->transitions[Tr].idInMarking] = 0.0;
		NewM->clock1[M->transitions[Tr].idInMarking] = 0.0;
	}
}

bool isGTransitionEnabled(Model* model, Marking* marking) {
	checkEnabled(model, marking);

	int gtId;
	for (gtId = 0; gtId < model->N_transitions; gtId++) {
		if (model->transitions[gtId].type == TT_GENERAL)
			break;
	}

	return marking->enabling[gtId];

}

void fireGeneralTransition(Model *M, Marking *NewM) {
	int gtId;
	for (gtId = 0; gtId < M->N_transitions; gtId++) {
		if (M->transitions[gtId].type == TT_GENERAL)
			break;
	}

	fireTransition(M, NewM, gtId);
}

int gTransitionId(Model* M) {
	int gtId;
	for (gtId = 0; gtId < M->N_transitions; gtId++) {
		if (M->transitions[gtId].type == TT_GENERAL)
			return gtId;
	}
    return -1;
}


void printFluidRates(Model* M, Marking* K){
	for (int t = 0; t < M->N_transitions; t++){
		if (M->transitions[t].type == TT_FLUID || M->transitions[t].type == TT_FLUID_DYNA){
			printf("*%s (%d): %f \t* \n", M->transitions[t].id, K->enabling[t], K->actFluidRate[M->transitions[t].idInMarking]);
		}
	}
}

void printFluidLevels(Model* M, Marking* K){
	for (int p = 0; p < M->N_places; p++){
		if (M->places[p].type == PT_FLUID){
			printf ("*%s (%f + %fs) - drift: %f \t* \n",
					M->places[p].id, K->fluid0[M->places[p].idInMarking],
					K->fluid1[M->places[p].idInMarking],
					K->fluidPlaceDeriv[M->places[p].idInMarking]);
		}
	}

}

