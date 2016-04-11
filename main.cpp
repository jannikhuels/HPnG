/*
 * main.cpp
 *
 *  Created on: 10 nov. 2011
 *      Author: GhasemiehH
 */

#include <iostream>
#include <sstream>
#include <cstring>
#include <fstream>
#include <sys/time.h>
#include <unistd.h>
#include <stdlib.h>

#include <cmath>
//#include <direct.h>


#include "DFPN2.h"
#include "Line.h"
#include "TimedDiagram.h"
#include "IntervalSet.h"

#include "ModelChecker.h"

using namespace std;

double expPrarameter = 2;

double expo(double s){
    return (1 - exp(-s/expPrarameter));
}

int K = 4;
double teta = .5;

int fact(int n){
    if (n == 0)
        return 1;
    else
        return n*fact(n-1);
}

double gamma(double s){
    double sum = 0;
    for (int i = 0; i < K; i++){
        sum += (1/fact(i))*pow((s/teta), i)*exp(-s/teta);
    }
    return 1 - sum;
}


double sigma = 1;
double mu = 2;
double normal(double s){
    return .5*(1 + erf((s - mu)/(sigma*1.414213562)));
}

double foldedNormal(double s){
    return .5*(erf((s + mu)/(sigma*1.414213562)) + erf((s - mu)/(sigma*1.414213562)));
}

//the property: if Pf1 has more than 5 units of fluid.

int pIndex = 16;
double amount = 0;
bool propertyTest(Model* model, Marking* marking, double t0, double t1, double &s1, double &s2){
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
    if (p <= s2 && p >= s1){
        if (a < -ZERO_PREC) s1 = p;
        if (a > +ZERO_PREC) s2 = p;
        return true;
    } else if (p > s2){
        if (a * s2 + b <= amount) return true;
    } else if (p < s1){
        if (a * s1 + b <= amount) return true;
    }
    
    return false;
}

int mainIntegrity(int argc, char* argv[]){
    Model *model;
    
    model = ReadModel(argv[1]);
    
    double time = 0; //model->transitions[tIndex].time;
    mu =  9;
    sigma = 1;
    
    model->MaxTime = 24;
    
    long mtime;
    long seconds = 0;
    long useconds = 0;
    timeval gt1, gt2;
    
    char* EventNames[] = { "TRANSITION", "PLACE_LOWER_BOUNDRY", "PLACE_UPPER_BOUNDRY", "GAURD_ARC", "MAX_TIME_REACHED", "FIRST_NULL_EVENT"};
    
    
    stringstream ss;
    ss << "./" << argv[0] << "-" << argv[1] << ".result";
    stringstream ss3;
    ss3 << "./" << argv[0] << "-" << argv[1] << ".minimum";
    ofstream oFile(ss.str().c_str());
    ofstream mFile(ss3.str().c_str());
    
    Interval bound(0, 14);
    
    int batteryId = findPlace(model, "battery");
    int fullId = findPlace(model, "full");
    int stopId = findPlace(model, "stop");
    int enabledId = findPlace(model, "enabled");
    int gArcId = findArc(model, "0.17");
    int loadingId = findPlace(model, "loading");
    
    double bCapacity = model->places[batteryId].f_bound;
    double min_prob = 1.0;
    
    AtomDisFormula *trueFML = new AtomDisFormula(enabledId, 1);
    AtomDisFormula *stopFML = new AtomDisFormula(stopId, 1);
    AtomDisFormula *fullFML = new AtomDisFormula(fullId, 1);
    AtomDisFormula *loadingFML = new AtomDisFormula(loadingId, 1);
    //AtomContNegFormula *fullFML = new AtomContNegFormula(batteryId, bCapacity * 0.99);
    
    
    //Here we cant parametrize things!!
    
    //int capacities[5] = {670, 1100};
    //for (int i = 1; i < 2; i++) {
    //    bCapacity = capacities[i];
    //for (bCapacity = 300; bCapacity <= 3000; bCapacity += 300){
        
        
        //    model->places[0].f_bound = bCapacity;
        //    model->places[0].f_level = bCapacity;
        //    model->arcs[4].weight = bCapacity - .01;
        //    model->arcs[5].weight = bCapacity;
        
        //		model->arcs[6].weight = bCapacity*.3;
        //		model->arcs[7].weight = bCapacity*.3;
        //
        //		model->arcs[6].weight = bCapacity*.5;
        //		model->arcs[7].weight = bCapacity*.3;
        
        //    model->arcs[6].weight = 0;
        //    model->arcs[7].weight = 0;
        
        //for (double portion = 0.1; portion <= 0.9; portion += 0.2){
            
            //model->arcs[6].weight = round(bCapacity*portion/100);
            //model->arcs[7].weight = round(bCapacity*portion/100);
            //model->arcs[6].weight = bCapacity;
            //model->arcs[7].weight = bCapacity;
            //model->arcs[6].weight = bCapacity < 550?bCapacity:550;
            //model->arcs[7].weight = bCapacity < 550?bCapacity:550;
            
            
            //		model->transitions[model->N_transitions-1].flowRate = ((.5-.3)/8)*bCapacity;
            //model->arcs[9].weight = bCapacity*0.2;
            //model->arcs[8].weight = bCapacity*0.5;
            
            //for (double failTime = 0; failTime < model->MaxTime; failTime += 1){
                            //model->transitions[6].time = failTime;
                //int transId = findTransition(model, "grid_failed");
                //model->transitions[transId].time = failTime;
                //time = model->transitions[6].time;
                //time = failTime;
    //for (double load = 1; load >= 0.5; load -= .1) {
        for (sigma=0.25; sigma<=1.5; sigma+=0.25) {
        
                InitializeModel(model);
            //model->arcs[gArcId].weight = bCapacity * load;
            
                
                Marking* initialMarking = createInitialMarking(model);
                checkEnabled(model, initialMarking, 0 ,0);
                
                TimedDiagram::getInstance()->clear();
                TimedDiagram::getInstance()->setModel(model);
                
                gettimeofday(&gt1, NULL);
                TimedDiagram::getInstance()->generateDiagram(initialMarking);
                gettimeofday(&gt2, NULL);
                seconds  += gt2.tv_sec  - gt1.tv_sec;
                useconds += gt2.tv_usec - gt1.tv_usec;
                
                //Model checking started
                ModelChecker modelChecker(model, TimedDiagram::getInstance());
                //modelChecker.ttc = time;
    
                modelChecker.scale = TimedDiagram::getInstance()->scale;
        
                //AtomContNegFormula *fullFML = new AtomContNegFormula(batteryId, bCapacity);
    IntervalSet *full;
    IntervalSet *stop;
    //stop = modelChecker.until(trueFML, stopFML, bound, time);
    //full = modelChecker.until(trueFML, fullFML, bound, time);
    //modelChecker.ttc = 16;
    //modelChecker.iSetAtomDis(stop, stopFML);
    //modelChecker.iSetAtomCont(full, fullFML);
    
    
                IntervalSet *res;
            res = modelChecker.until(loadingFML, fullFML, bound, time);
   // modelChecker.iSetAnd(res, stop, full);
    

                //modelChecker.iSetAnd(res, notFull, stop);
                //modelChecker.iSetNeg(notFull, res);
                double prob = modelChecker.calcProb(res, normal, 0);
    /*pIndex = batteryId;
    amount = 26999;
            modelChecker.ttc = 17;
    prob = TimedDiagram::getInstance()->calProbAtTime(17, normal, propertyTest);*/
    
                oFile << prob << "	";
                /*if (prob < min_prob && failTime <= 48) {
                    min_prob=prob;
                }*/
                
                //res->clear();
            //}
            //mFile << min_prob << "	";
            //min_prob = 1.0;
            // oInt << endl;
            // oFile << endl;
            
            //}
        }
        oFile << endl;
    //}
//}

    
    
    return 0;
}


/******************************************This where you have to start **********************************************/
int mainBattery(int argc, char* argv[]){
    Model *model;
    
    model = ReadModel(argv[1]);
    
    unsigned int bIndex = findPlace(model, "battery");
    unsigned int gOnIndex = findPlace(model, "grid_on");
    
    
    double time = 0; //model->transitions[tIndex].time;
    mu =  .5;
    
    model->MaxTime = 2*24 + 10;
    
    long mtime;
    long seconds = 0;
    long useconds = 0;
    timeval gt1, gt2;
    
    char* EventNames[] = { "TRANSITION", "PLACE_LOWER_BOUNDRY", "PLACE_UPPER_BOUNDRY", "GAURD_ARC", "MAX_TIME_REACHED", "FIRST_NULL_EVENT"};
    
    
    AtomContPosFormula* phi1 = new AtomContPosFormula(bIndex, 0.01 );
    AtomDisFormula* psi = new AtomDisFormula(gOnIndex, 1);

    
    stringstream ss;
    ss << "./" << argv[0] << "-" << argv[1] << ".result";
    stringstream ss2;
    ss2 << "./" << argv[0] << "-" << argv[1] << ".time";
    stringstream ss3;
    ss3 << "./" << argv[0] << "-" << argv[1] << ".minimum";
    ofstream oFile(ss.str().c_str());
    ofstream oInt("./intervals.out");
    ofstream tFile(ss2.str().c_str());
    ofstream mFile(ss3.str().c_str());
    
    Interval bound(0, model->MaxTime);
    IntervalSet* res;
    
    int maxRgNum = 0;
    int avgRgNum = 0;
    int ptNum = 0;
    
    //Here we cant parametrize things!!
    double bCapacity = 3000;
    double min_prob = 1.0;
    //int capacities[5] = {670, 1100};
    //for (int i = 1; i < 2; i++) {
    //    bCapacity = capacities[i];
    int t3a_id = findArc(model, "0.18");
    int t4a_id = findArc(model, "0.19");
    int t1a_id = findArc(model, "0.16");
    int t2a_id = findArc(model, "0.17");
    int b_id = findPlace(model, "battery");
        //for (bCapacity = 2500; bCapacity <= 3000; bCapacity += 300){
            
            //model->places[b_id].f_bound = bCapacity;
            //model->places[b_id].f_level = bCapacity;
            //model->arcs[t1a_id].weight = bCapacity;
            //model->arcs[t2a_id].weight = bCapacity;
            
            
        
    //    model->places[0].f_bound = bCapacity;
    //    model->places[0].f_level = bCapacity;
    //    model->arcs[4].weight = bCapacity - .01;
    //    model->arcs[5].weight = bCapacity;
        
        //		model->arcs[6].weight = bCapacity*.3;
        //		model->arcs[7].weight = bCapacity*.3;
        //
        //		model->arcs[6].weight = bCapacity*.5;
        //		model->arcs[7].weight = bCapacity*.3;
        
    //    model->arcs[6].weight = 0;
    //    model->arcs[7].weight = 0;
        
        //for (double portion = 0.1; portion <= 0.9; portion += 0.2){
            
            //model->arcs[t3a_id].weight = bCapacity * portion;
            //model->arcs[t4a_id].weight = bCapacity * portion;
            //model->arcs[6].weight = round(bCapacity*portion/100);
            //model->arcs[7].weight = round(bCapacity*portion/100);
            //model->arcs[6].weight = bCapacity;
            //model->arcs[7].weight = bCapacity;
            //model->arcs[6].weight = bCapacity < 550?bCapacity:550;
            //model->arcs[7].weight = bCapacity < 550?bCapacity:550;
            
            
            //		model->transitions[model->N_transitions-1].flowRate = ((.5-.3)/8)*bCapacity;
            //model->arcs[9].weight = bCapacity*0.2;
            //model->arcs[8].weight = bCapacity*0.5;
            
            for (double failTime = 0; failTime < 56; failTime += 1){
                //model->transitions[6].time = failTime;
                //int transId = findTransition(model, "grid_failed");
                //model->transitions[transId].time = failTime;
                //time = model->transitions[6].time;
                time = failTime;
                
                InitializeModel(model);
                
                Marking* initialMarking = createInitialMarking(model);
                checkEnabled(model, initialMarking, 0 ,0);
                
                TimedDiagram::getInstance()->clear();
                TimedDiagram::getInstance()->setModel(model);
                
                gettimeofday(&gt1, NULL);
                TimedDiagram::getInstance()->generateDiagram(initialMarking);
                gettimeofday(&gt2, NULL);
                seconds  += gt2.tv_sec  - gt1.tv_sec;
                useconds += gt2.tv_usec - gt1.tv_usec;
                /******************************* YOU MAY NOT NEED THIS *********************************************/
                int regNum = TimedDiagram::getInstance()->regionList.size();
                maxRgNum = (maxRgNum < regNum) ? regNum : maxRgNum;
                avgRgNum += regNum;
                ptNum++;
                double rgTime = ((gt2.tv_sec  - gt1.tv_sec) * 1000 + (gt2.tv_usec - gt1.tv_usec)/1000.0) + 0.5;
                
                tFile << ptNum << ": " << bCapacity << " - " << failTime <<  "\t" <<  regNum << "\t" << rgTime << "ms\t" ;
                /******************************* YOU MAY NOT NEED THIS *********************************************/
                
                //For debugging purposes.
                //			TimedDiagram::getInstance()->scale = 5;
                //			if (argc == 3){
                //				std::cout << "Writing the debug region diagram...." << std::endl;
                //				std::stringstream ss;
                //				ss << argv[1] << "_rd";
                //				TimedDiagram::getInstance()->saveDiagram(ss.str());
                //				cv::Mat flipped;
                //				cv::flip(TimedDiagram::getInstance()->debugImage, flipped, 0);
                //				cv::imshow("test", flipped);
                //				cv::waitKey(0);
                //			}
                
                
                //Model checking started
                ModelChecker modelChecker(model, TimedDiagram::getInstance());
                
                amount = 1500;
                pIndex = bIndex;
                double prob = TimedDiagram::getInstance()->calProbAtTime(failTime + .1, foldedNormal, propertyTest);
                
                /*modelChecker.ttc = failTime;
                modelChecker.scale = TimedDiagram::getInstance()->scale;
                AtomContNegFormula *test = new AtomContNegFormula(bIndex, 3000);
                int id = findPlace(model, "P1");
                AtomDisFormula *tt = new AtomDisFormula(id,0);*/
                
                //IntervalSet *result;
                //Interval b(0,24);
                //result = modelChecker.until(tt, test, b, failTime);
                //result = modelChecker.calcAtomContISetAtTime(failTime, bIndex, 3000);
                
                //			modelChecker.debugImage = TimedDiagram::getInstance()->debugImage;
                //modelChecker.scale = TimedDiagram::getInstance()->scale;
                
                //			gettimeofday(&gt1, NULL);
                //res = modelChecker.until(phi1, psi, bound, time);
                //			gettimeofday(&gt2, NULL);
                //			seconds  += gt2.tv_sec  - gt1.tv_sec;
                //			useconds += gt2.tv_usec - gt1.tv_usec;
                
                //			double mcTime = ((gt2.tv_sec  - gt1.tv_sec) * 1000 + (gt2.tv_usec - gt1.tv_usec)/1000.0) + 0.5;
                //			tFile << mcTime << "ms"<< std::endl;
                
                /******************************* YOU MAY NOT NEED THIS *********************************************/
                /*cout << "-------------------------------------" << endl;
                cout << "The satisfaction set of Until formula: " << endl;
                for (unsigned int i = 0; i < res->intervals.size(); i++){
                    res->intervals[i].end = res->intervals[i].end - time;
                    res->intervals[i].start = res->intervals[i].start - time;
                }
                
                res->print(cout);
                res->print(oInt);*/
                /******************************* YOU MAY NOT NEED THIS *********************************************/
                
                //double prob =modelChecker.calcProb(result, foldedNormal, 0);
                oFile << prob << "	";
                if (prob < min_prob && failTime <= 48) {
                    min_prob=prob;
                }
                
                
                //			if (argc == 3){
                //				cv::Mat flipped;
                //				cv::flip(modelChecker.debugImage, flipped, 0);
                //				cv::imshow("test", flipped);
                //				cv::waitKey(0);
                //			}
                
                
                //res->clear();
                
            }
            mFile << min_prob << "	";
            min_prob = 1.0;
           // oInt << endl;
           // oFile << endl;
            
        //}
        //}
        mFile << endl;
    //}
    //}
    
    ofstream timeFile("./ctime.txt", ios::app);
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cout << "Time : " << mtime<< "ms" << std::endl;
    timeFile << argv[0] << "\t" << argv[1] << "\t" << maxRgNum << "\t" << avgRgNum / ptNum << "\t" << mtime<< "ms" << std::endl;
    
    
    
    //	FILE* gnuplotPipe;
    //	if (gnuplotPipe = popen("gnuplot -persist","w"))
    //	{
    //		fprintf(gnuplotPipe,
    //				"set pm3d \n unset surface\n "
    //				"set xlabel \"Failure Time\" font \"Times-Roman, 16\"\n"
    ////				"set xlabel \"Buffer Capacity\" font \"Times-Roman, 20\"\n"
    //				"set zlabel \"Probability\" font \"Times-Roman, 16\"\n"
    //				"set ylabel \"Capacity\" font \"Times-Roman, 16\"\n"
    ////				"set xtics (\"0.5\" 0, \"5\" 9) \n "
    ////				"set xtics (\"0\" 0,\"5\" 5,\"10\" 10,\"15\" 15, \"25\" 20,\"30\" 25) font \"Times-Roman, 20\"\n "
    //
    ////				"set zrange [0:1] \n"
    ////				"set xrange [0:48] \n"
    ////				"set xtics (\"0\" 0, \"4\" 4, \"8\" 8, \"12\" 12, \"16\" 16, \"20\" 20, \"24\" 24) font \"Times-Roman, 18\"\n "
    ////				"set ytics font \"Times-Roman, 18\"\n "
    //				"set palette rgbformulae 33,13,10 \n"
    ////				"set terminal wxt size 1300,800 \n"
    //				);
    //
    //
    //
    //		fflush(gnuplotPipe);
    //		fprintf(gnuplotPipe,"splot \"%s\" matrix with lines\n", "./results.out");
    //		fflush(gnuplotPipe);
    //		fprintf(gnuplotPipe,"exit \n");
    //		fclose(gnuplotPipe);
    //	}
    //
    
    
    return 0;
}

int mainBattery1(int argc, char* argv[]){
    Model *model;
    
    model = ReadModel(argv[1]);
    
    unsigned int bIndex = findPlace(model, "battery");
    unsigned int gOnIndex = findPlace(model, "grid_on");
    int t3a_id = findArc(model, "0.18");
    int t4a_id = findArc(model, "0.19");
    int t1a_id = findArc(model, "0.16");
    int t2a_id = findArc(model, "0.17");
    int fail =findTransition(model, "grid_failed");
    
    
    double time = 0; //model->transitions[tIndex].time;
    mu =  .5;
    
    model->MaxTime = 2*24 + 10;
    
    long seconds = 0;
    long useconds = 0;
    timeval gt1, gt2;
    
    char* EventNames[] = { "TRANSITION", "PLACE_LOWER_BOUNDRY", "PLACE_UPPER_BOUNDRY", "GAURD_ARC", "MAX_TIME_REACHED", "FIRST_NULL_EVENT"};
    
    
    AtomContPosFormula* phi1 = new AtomContPosFormula(bIndex, 0.01 );
    AtomDisFormula* psi = new AtomDisFormula(gOnIndex, 1);
    
    stringstream ss;
    ss << "./" << argv[0] << "-" << argv[1] << ".result";
    ofstream oFile(ss.str().c_str());
    
    Interval bound(0, model->MaxTime);
    IntervalSet* res;
    
    int maxRgNum = 0;
    int avgRgNum = 0;
    int ptNum = 0;
    
    //Here we can parametrize things!!
    //double bCapacity = 350;
    //double bCapacity = model->places[bIndex].f_bound;
    //	for (double bCapacity = 500; bCapacity <  3001; bCapacity += 500){
    //model->places[bIndex].f_bound = bCapacity;
    //model->places[bIndex].f_level = bCapacity;
    //model->arcs[t1a_id].weight = bCapacity;
    //model->arcs[t2a_id].weight = bCapacity;
    
    //for (double portion = 0; portion <= 0.9; portion += .1){
      //  model->arcs[t3a_id].weight = bCapacity*portion;
        //model->arcs[t4a_id].weight = bCapacity*portion;
    //for (double error = 1.01; error >= 0.2; error -= 0.2) {
    double error = 0.8;
        model->transitions[findTransition(model, "p1_rate")].flowRate = model->transitions[findTransition(model, "p1_rate")].flowRate * error;
        model->transitions[findTransition(model, "p2_rate")].flowRate = model->transitions[findTransition(model, "p2_rate")].flowRate * error;
        model->transitions[findTransition(model, "p3_rate")].flowRate = model->transitions[findTransition(model, "p3_rate")].flowRate * error;
        model->transitions[findTransition(model, "p4_rate")].flowRate = model->transitions[findTransition(model, "p4_rate")].flowRate * error;
        model->transitions[findTransition(model, "p5_rate")].flowRate = model->transitions[findTransition(model, "p5_rate")].flowRate * error;
        model->transitions[findTransition(model, "p6_rate")].flowRate = model->transitions[findTransition(model, "p6_rate")].flowRate * error;
        model->transitions[findTransition(model, "p7_rate")].flowRate = model->transitions[findTransition(model, "p7_rate")].flowRate * error;
        model->transitions[findTransition(model, "p8_rate")].flowRate = model->transitions[findTransition(model, "p8_rate")].flowRate * error;
        model->transitions[findTransition(model, "p9_rate")].flowRate = model->transitions[findTransition(model, "p9_rate")].flowRate * error;
        
        
        for (double failTime = 0; failTime < model->MaxTime; failTime += 1){
            model->transitions[fail].time = failTime;
            time = failTime;
            
            InitializeModel(model);
            
            Marking* initialMarking = createInitialMarking(model);
            checkEnabled(model, initialMarking, 0 ,0);
            
            TimedDiagram::getInstance()->clear();
            TimedDiagram::getInstance()->setModel(model);
            
            TimedDiagram::getInstance()->generateDiagram(initialMarking);
            
            
            //Model checking started
            ModelChecker modelChecker(model, TimedDiagram::getInstance());
            
            modelChecker.scale = TimedDiagram::getInstance()->scale;
            
            res = modelChecker.until(phi1, psi, bound, time + 0.01);
            for (unsigned int i = 0; i < res->intervals.size(); i++){
                res->intervals[i].end = res->intervals[i].end - time;
                res->intervals[i].start = res->intervals[i].start - time;
            }
            
            res->print(cout);
            
            oFile << modelChecker.calcProb(res, foldedNormal, 0) << "	";
            
            
            
            res->clear();
        }
        
        oFile << endl;
        
    //}
    
    
    return 0;
}

int mainPriorities(int argc, char* argv[]){
    Model *model;
    
    model = ReadModel(argv[1]);
    
    unsigned int bIndex = findPlace(model, "battery");
    unsigned int gOnIndex = findPlace(model, "grid_on");
    int t3a_id = findArc(model, "0.18");
    int t4a_id = findArc(model, "0.19");
    int t1a_id = findArc(model, "0.16");
    int t2a_id = findArc(model, "0.17");
    int fail =findTransition(model, "grid_failed");
    
    
    double time = 0; //model->transitions[tIndex].time;
    mu =  .5;
    
    model->MaxTime = 2*24 + 10;
    
    long seconds = 0;
    long useconds = 0;
    timeval gt1, gt2;
    
    char* EventNames[] = { "TRANSITION", "PLACE_LOWER_BOUNDRY", "PLACE_UPPER_BOUNDRY", "GAURD_ARC", "MAX_TIME_REACHED", "FIRST_NULL_EVENT"};
    
    
    AtomContPosFormula* phi1 = new AtomContPosFormula(bIndex, 0.01 );
    AtomDisFormula* psi = new AtomDisFormula(gOnIndex, 1);
    
    stringstream ss;
    ss << "./" << argv[0] << "-" << argv[1] << ".result";
    ofstream oFile(ss.str().c_str());
    
    Interval bound(0, model->MaxTime);
    IntervalSet* res;
    
    int maxRgNum = 0;
    int avgRgNum = 0;
    int ptNum = 0;
    
    //Here we can parametrize things!!
    //double bCapacity = 350;
    //double bCapacity = model->places[bIndex].f_bound;
    //	for (double bCapacity = 500; bCapacity <  3001; bCapacity += 500){
    //model->places[bIndex].f_bound = bCapacity;
    //model->places[bIndex].f_level = bCapacity;
    //model->arcs[t1a_id].weight = bCapacity;
    //model->arcs[t2a_id].weight = bCapacity;
    
    //for (double portion = 0; portion <= 0.9; portion += .1){
    //  model->arcs[t3a_id].weight = bCapacity*portion;
    //model->arcs[t4a_id].weight = bCapacity*portion;
    //for (double error = 1.01; error >= 0.2; error -= 0.2) {
    /*double error = 0.8;
    model->transitions[findTransition(model, "p1_rate")].flowRate = model->transitions[findTransition(model, "p1_rate")].flowRate * error;
    model->transitions[findTransition(model, "p2_rate")].flowRate = model->transitions[findTransition(model, "p2_rate")].flowRate * error;
    model->transitions[findTransition(model, "p3_rate")].flowRate = model->transitions[findTransition(model, "p3_rate")].flowRate * error;
    model->transitions[findTransition(model, "p4_rate")].flowRate = model->transitions[findTransition(model, "p4_rate")].flowRate * error;
    model->transitions[findTransition(model, "p5_rate")].flowRate = model->transitions[findTransition(model, "p5_rate")].flowRate * error;
    model->transitions[findTransition(model, "p6_rate")].flowRate = model->transitions[findTransition(model, "p6_rate")].flowRate * error;
    model->transitions[findTransition(model, "p7_rate")].flowRate = model->transitions[findTransition(model, "p7_rate")].flowRate * error;
    model->transitions[findTransition(model, "p8_rate")].flowRate = model->transitions[findTransition(model, "p8_rate")].flowRate * error;
    model->transitions[findTransition(model, "p9_rate")].flowRate = model->transitions[findTransition(model, "p9_rate")].flowRate * error;*/
    
    
    for (double failTime = 0; failTime < model->MaxTime; failTime += 1){
        model->transitions[fail].time = failTime;
        time = failTime;
        
        InitializeModel(model);
        
        Marking* initialMarking = createInitialMarking(model);
        checkEnabled(model, initialMarking, 0 ,0);
        
        TimedDiagram::getInstance()->clear();
        TimedDiagram::getInstance()->setModel(model);
        
        TimedDiagram::getInstance()->generateDiagram(initialMarking);
        
        
        //Model checking started
        ModelChecker modelChecker(model, TimedDiagram::getInstance());
        
        modelChecker.scale = TimedDiagram::getInstance()->scale;
        
        res = modelChecker.until(phi1, psi, bound, time + 0.01);
        for (unsigned int i = 0; i < res->intervals.size(); i++){
            res->intervals[i].end = res->intervals[i].end - time;
            res->intervals[i].start = res->intervals[i].start - time;
        }
        
        res->print(cout);
        
        oFile << modelChecker.calcProb(res, foldedNormal, 0) << "	";
        
        
        
        res->clear();
    }
    
    oFile << endl;
    
    //}
    
    
    return 0;
}


int mainStocRainTime(int argc, char *argv[]) {
    Model *model;
    
    if (argc < 3) {
        //std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName STD(optional, should be present if you want the STD as a JPG file.)" << std::endl;
        std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName" << std::endl;
        return 0;
    }
    
    model = ReadModel(argv[1]);
    
    unsigned int pIndex, disPIndex, NormTIndex;
    
    for (int i = 0; i < model->N_places; i++) {
        //if (model->places[i].type != PT_FLUID) continue;
        if (model->places[i].type == PT_FLUID && strncmp(model->places[i].id, argv[2], strlen(argv[2])) == 0)  pIndex = i;
        if (model->places[i].type == PT_DISCRETE && strncmp(model->places[i].id, "normalAgain_P", strlen("normalAgain_P")) == 0)  disPIndex = i;
    }
    
    for (int i = 0; i < model->N_transitions; i++) {
        if (strncmp(model->transitions[i].id, "startNormal_T", strlen("startNormal_T")) == 0)  NormTIndex= i;
    }
    
    double time = 0; //model->transitions[tIndex].time;
    mu =  .5;
    
    model->MaxTime = 40.0;
    
    long mtime;
    long seconds = 0;
    long useconds = 0;
    timeval gt1, gt2;
    
    
    
    AtomContNegFormula* phi1 = new AtomContNegFormula(pIndex, 0.01 );
    AtomDisFormula* psi = new AtomDisFormula(disPIndex, 1);
    
    ofstream oFile("./results.out");
    ofstream oInt("./intervals.out");
    for (double T = 30; T < 31; T += 1){
        Interval bound(0, T);
        IntervalSet* res;
        for (mu = .1; mu <= 4; mu += .1){
            //			for (double rainTime = .5; rainTime <= 5 ; rainTime += .5){
            for (double bufferCapacity = 5; bufferCapacity <= 30 ; bufferCapacity += 1){
                model->places[4].f_bound = bufferCapacity;
                model->arcs[14].weight = bufferCapacity;
                model->arcs[15].weight = bufferCapacity;
                time = model->transitions[2].time;
                
                InitializeModel(model);
                
                Marking* initialMarking = createInitialMarking(model);
                checkEnabled(model, initialMarking, 0 ,0);
                
                TimedDiagram::getInstance()->clear();
                TimedDiagram::getInstance()->setModel(model);
                
                gettimeofday(&gt1, NULL);
                TimedDiagram::getInstance()->generateDiagram(initialMarking);
                gettimeofday(&gt2, NULL);
                seconds  += gt2.tv_sec  - gt1.tv_sec;
                useconds += gt2.tv_usec - gt1.tv_usec;
                
                
                //				TimedDiagram::getInstance()->scale = 20;
                //				if (argc == 4){
                //					std::cout << "Writing the debug region diagram...." << std::endl;
                //					std::stringstream ss;
                //					ss << argv[1] << "_rd";
                //					TimedDiagram::getInstance()->saveDiagram(ss.str());
                //					cv::Mat flipped;
                //					cv::flip(TimedDiagram::getInstance()->debugImage, flipped, 0);
                //					cv::imshow("test", flipped);
                //					cv::waitKey(0);
                //				}
                
                
                ModelChecker modelChecker(model, TimedDiagram::getInstance());
                //				modelChecker.debugImage = TimedDiagram::getInstance()->debugImage;
                modelChecker.scale = TimedDiagram::getInstance()->scale;
                
                gettimeofday(&gt1, NULL);
                res = modelChecker.until(phi1, psi, bound, 0);
                gettimeofday(&gt2, NULL);
                seconds  += gt2.tv_sec  - gt1.tv_sec;
                useconds += gt2.tv_usec - gt1.tv_usec;
                
                
                cout << "-------------------------------------" << endl;
                cout << "The satisfaction set of Until formula: " << endl;
                for (unsigned int i = 0; i < res->intervals.size(); i++){
                    res->intervals[i].end = res->intervals[i].end - time;
                    res->intervals[i].start = res->intervals[i].start - time;
                }
                
                res->print(cout);
                res->print(oInt);
                
                oFile << modelChecker.calcProb(res, foldedNormal, 0) << "	";
                
                
                //				if (argc == 4){
                //					cv::Mat flipped;
                //					cv::flip(modelChecker.debugImage, flipped, 0);
                //					cv::imshow("test", flipped);
                //					cv::waitKey(0);
                //				}
                
                
                res->clear();
            }
            oInt << endl;
            oFile << endl;
        }
    }
    
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cout << "Time : " << mtime<< "ms" << std::endl;
    
    FILE* gnuplotPipe;
    if (gnuplotPipe = popen("gnuplot -persist","w"))
    {
        fprintf(gnuplotPipe,"set pm3d \n unset surface\n "
                //				"set xlabel \"Start Rain Time\" \n"
                "set xlabel \"Buffer Capacity\" font \"Times-Roman, 20\"\n"
                "set zlabel \"Probability\" font \"Times-Roman, 20\"\n"
                "set ylabel \"mu \"\n"
                //				"set xtics (\"0.5\" 0, \"5\" 9) \n "
                "set xtics (\"5\" 0,\"10\" 5,\"15\" 10,\"20\" 15, \"25\" 20,\"30\" 25) font \"Times-Roman, 20\"\n "
                "set ytics (\"0.1\" 0, \"4\" 38) font \"Times-Roman, 20\"\n "
                "set ztics font \"Times-Roman, 20\"\n "
                "set zrange [0:1] \n"
                "set palette rgbformulae 33,13,10 \n"
                "set terminal wxt size 1300,800 \n");
        
        fflush(gnuplotPipe);
        fprintf(gnuplotPipe,"splot \"%s\" matrix with lines\n", "./results.out");
        fflush(gnuplotPipe);
        fprintf(gnuplotPipe,"exit \n");
        fclose(gnuplotPipe);
    }
    
    
    return 0;
}

int mainFailureBufffer(int argc, char *argv[]) {
    
    Model *model;
    
    if (argc < 3) {
        //std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName STD(optional, should be present if you want the STD as a JPG file.)" << std::endl;
        std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName" << std::endl;
        return 0;
    }
    
    model = ReadModel(argv[1]);
    
    unsigned int pIndex, disPIndex, NormTIndex;
    
    for (int i = 0; i < model->N_places; i++) {
        //if (model->places[i].type != PT_FLUID) continue;
        if (model->places[i].type == PT_FLUID && strncmp(model->places[i].id, argv[2], strlen(argv[2])) == 0)  pIndex = i;
        if (model->places[i].type == PT_DISCRETE && strncmp(model->places[i].id, "repaired_P", strlen("repaired_P")) == 0)  disPIndex = i;
    }
    
    for (int i = 0; i < model->N_transitions; i++) {
        if (strncmp(model->transitions[i].id, "startNormal_T", strlen("startNormal_T")) == 0)  NormTIndex= i;
    }
    
    double time = 0; //model->transitions[tIndex].time;
    mu =  .5;
    
    model->MaxTime = 40.0;
    
    long mtime, seconds, useconds;
    expPrarameter = 1;
    
    
    AtomContNegFormula* phi1 = new AtomContNegFormula(pIndex, 0.01 );
    AtomDisFormula* psi = new AtomDisFormula(disPIndex, 1);
    
    ofstream oFile("./results.out");
    ofstream oInt("./intervals.out");
    for (double T = 30; T < 31; T += 1){
        Interval bound(0, T);
        IntervalSet* res;
        model->transitions[0].flowRate = 12.2;
        
        for (double bufferCapacity = 5; bufferCapacity <= 30 ; bufferCapacity += 1){
            model->places[4].f_bound = bufferCapacity;
            model->arcs[14].weight = bufferCapacity;
            model->arcs[15].weight = bufferCapacity;
            
            for (double failureTime = .5; failureTime <= 5 ; failureTime += .1){
                model->transitions[5].time = failureTime;
                time = failureTime;
                InitializeModel(model);
                
                Marking* initialMarking = createInitialMarking(model);
                checkEnabled(model, initialMarking, 0 ,0);
                
                TimedDiagram::getInstance()->clear();
                TimedDiagram::getInstance()->setModel(model);
                TimedDiagram::getInstance()->generateDiagram(initialMarking);
                
                TimedDiagram::getInstance()->scale = 20;
                //				if (argc == 4){
                //					std::cout << "Writing the debug region diagram...." << std::endl;
                //					std::stringstream ss;
                //					ss << argv[1] << "_rd";
                //					TimedDiagram::getInstance()->saveDiagram(ss.str());
                //					cv::Mat flipped;
                //					cv::flip(TimedDiagram::getInstance()->debugImage, flipped, 0);
                //					cv::imshow("test", flipped);
                //					cv::waitKey(0);
                //				}
                
                
                ModelChecker modelChecker(model, TimedDiagram::getInstance());
                //				modelChecker.debugImage = TimedDiagram::getInstance()->debugImage;
                modelChecker.scale = TimedDiagram::getInstance()->scale;
                
                res = modelChecker.until(phi1, psi, bound, 0);
                
                cout << "-------------------------------------" << endl;
                cout << "The satisfaction set of Until formula: " << endl;
                for (unsigned int i = 0; i < res->intervals.size(); i++){
                    res->intervals[i].end = res->intervals[i].end - time;
                    res->intervals[i].start = res->intervals[i].start - time;
                }
                
                res->print(cout);
                res->print(oInt);
                
                oFile << modelChecker.calcProb(res, expo, 0) << "	";
                
                
                //				if (argc == 4){
                //					cv::Mat flipped;
                //					cv::flip(modelChecker.debugImage, flipped, 0);
                //					cv::imshow("test", flipped);
                //					cv::waitKey(0);
                //				}
                
                
                res->clear();
            }
            oInt << endl;
            oFile << endl;
        }
    }
    
    FILE* gnuplotPipe;
    if (gnuplotPipe = popen("gnuplot -persist","w"))
    {
        fprintf(gnuplotPipe,"set pm3d \n unset surface\n "
                "set xlabel \"Failure Time\" \n"
                "set ylabel \"Rain Rate \"\n"
                "set zlabel \"Probability\"\n"
                "set xtics (\"0.5\" 0, \"5\" 44) \n "
                "set ytics (\"6\" 0, \"13\" 33) \n "
                "set palette rgbformulae 33,13,10 \n"
                "set terminal wxt size 1300,800 \n");
        
        fflush(gnuplotPipe);
        fprintf(gnuplotPipe,"splot \"%s\" matrix with lines\n", "./results.out");
        fflush(gnuplotPipe);
        fprintf(gnuplotPipe,"exit \n");
        fclose(gnuplotPipe);
    }
    
    
    return 0;
}


int mainFailure(int argc, char *argv[]) {
    
    Model *model;
    
    if (argc < 3) {
        //std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName STD(optional, should be present if you want the STD as a JPG file.)" << std::endl;
        std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName" << std::endl;
        return 0;
    }
    
    model = ReadModel(argv[1]);
    
    unsigned int pIndex, disPIndex, NormTIndex;
    
    for (int i = 0; i < model->N_places; i++) {
        //if (model->places[i].type != PT_FLUID) continue;
        if (model->places[i].type == PT_FLUID && strncmp(model->places[i].id, argv[2], strlen(argv[2])) == 0)  pIndex = i;
        if (model->places[i].type == PT_DISCRETE && strncmp(model->places[i].id, "repaired_P", strlen("repaired_P")) == 0)  disPIndex = i;
    }
    
    for (int i = 0; i < model->N_transitions; i++) {
        if (strncmp(model->transitions[i].id, "startNormal_T", strlen("startNormal_T")) == 0)  NormTIndex= i;
    }
    
    double time = 0; //model->transitions[tIndex].time;
    mu =  .5;
    
    model->MaxTime = 40.0;
    
    long mtime;
    long seconds = 0;
    long useconds = 0;
    timeval gt1, gt2;
    
    expPrarameter = 2;
    
    
    AtomContNegFormula* phi1 = new AtomContNegFormula(pIndex, 0.01 );
    AtomDisFormula* psi = new AtomDisFormula(disPIndex, 1);
    
    ofstream oFile("./results.out");
    ofstream oInt("./intervals.out");
    for (double T = 30; T < 31; T += 1){
        Interval bound(0, T);
        IntervalSet* res;
        //		TODO: [very importante] debug this case: for (double rainRate = 3; rainRate <= 13; rainRate+= .2)--> noisy behavior.
        //		double rainRate = 10;
        //		double failureTime = 3.5;
        for (double rainRate = 6; rainRate <= 13; rainRate+= .5){
            model->transitions[0].flowRate = rainRate;
            for (double failureTime = .5; failureTime <= 5 ; failureTime += .2){
                model->transitions[5].time = failureTime;
                time = failureTime;
                bound.end = time + T;
                InitializeModel(model);
                
                Marking* initialMarking = createInitialMarking(model);
                checkEnabled(model, initialMarking, 0 ,0);
                
                TimedDiagram::getInstance()->clear();
                TimedDiagram::getInstance()->setModel(model);
                TimedDiagram::getInstance()->generateDiagram(initialMarking);
                
                //				TimedDiagram::getInstance()->scale = 20;
                //				if (argc == 4){
                //					std::cout << "Writing the debug region diagram...." << std::endl;
                //					std::stringstream ss;
                //					ss << argv[1] << "_rd";
                //					TimedDiagram::getInstance()->saveDiagram(ss.str());
                //					cv::Mat flipped;
                //					cv::flip(TimedDiagram::getInstance()->debugImage, flipped, 0);
                //					cv::imshow("test", flipped);
                //					cv::waitKey(0);
                //				}
                
                
                ModelChecker modelChecker(model, TimedDiagram::getInstance());
                //				modelChecker.debugImage = TimedDiagram::getInstance()->debugImage;
                modelChecker.scale = TimedDiagram::getInstance()->scale;
                
                res = modelChecker.until(phi1, psi, bound, 0);
                
                cout << "-------------------------------------" << endl;
                cout << "The satisfaction set of Until formula: " << endl;
                for (unsigned int i = 0; i < res->intervals.size(); i++){
                    res->intervals[i].end = res->intervals[i].end - time;
                    res->intervals[i].start = res->intervals[i].start - time;
                }
                
                res->print(cout);
                res->print(oInt);
                
                oFile << modelChecker.calcProb(res, expo, 0) << "	";
                
                
                //				if (argc == 4){
                //					cv::Mat flipped;
                //					cv::flip(modelChecker.debugImage, flipped, 0);
                //					cv::imshow("test", flipped);
                //					cv::waitKey(0);
                //				}
                
                
                res->clear();
            }
            oInt << endl;
            oFile << endl;
        }
    }
    
    FILE* gnuplotPipe;
    if (gnuplotPipe = popen("gnuplot -persist","w"))
    {
        fprintf(gnuplotPipe,"set pm3d \n unset surface\n "
                "set xlabel \"Failure Time\" font \"Times-Roman, 20\"\n"
                "set ylabel \"Intake Rate \" font \"Times-Roman, 20\"\n"
                "set zlabel \"Probability\" font \"Times-Roman, 20\"\n"
                "set xtics (\"0.5\" 0, \"5\" 22) font \"Times-Roman, 20\"\n "
                "set ytics (\"6\" 0, \"13\" 14) font \"Times-Roman, 20\"\n "
                "set ztics font \"Times-Roman, 20\"\n "
                "set zrange [0:1] \n"
                "set palette rgbformulae 33,13,10 \n"
                "set terminal wxt size 1300,800 \n");
        
        fflush(gnuplotPipe);
        fprintf(gnuplotPipe,"splot \"%s\" matrix with lines\n", "./results.out");
        fflush(gnuplotPipe);
        fprintf(gnuplotPipe,"exit \n");
        fclose(gnuplotPipe);
    }
    
    
    return 0;
}

int mainStreetRate(int argc, char *argv[]) {
    Model *model;
    
    if (argc < 3) {
        //std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName STD(optional, should be present if you want the STD as a JPG file.)" << std::endl;
        std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName" << std::endl;
        return 0;
    }
    
    model = ReadModel(argv[1]);
    
    unsigned int pIndex, disPIndex, NormTIndex;
    
    for (int i = 0; i < model->N_places; i++) {
        //if (model->places[i].type != PT_FLUID) continue;
        if (model->places[i].type == PT_FLUID && strncmp(model->places[i].id, argv[2], strlen(argv[2])) == 0)  pIndex = i;
        if (model->places[i].type == PT_DISCRETE && strncmp(model->places[i].id, "normalAgain_P", strlen("normalAgain_P")) == 0)  disPIndex = i;
    }
    
    for (int i = 0; i < model->N_transitions; i++) {
        if (strncmp(model->transitions[i].id, "startNormal_T", strlen("startNormal_T")) == 0)  NormTIndex= i;
    }
    
    double time = 0; //model->transitions[tIndex].time;
    mu =  .5;
    
    model->MaxTime = 40.0;
    
    long mtime;
    long seconds = 0;
    long useconds = 0;
    timeval gt1, gt2;
    
    
    AtomContNegFormula* phi1 = new AtomContNegFormula(pIndex, 0.01 );
    AtomDisFormula* psi = new AtomDisFormula(disPIndex, 1);
    
    ofstream oFile("./results.out");
    ofstream oInt("./intervals.out");
    for (double T = 30; T < 31; T += 1){
        Interval bound(0, T);
        IntervalSet* res;
        for (mu = .1; mu <= 4; mu += .1){
            //			for (double rainTime = .5; rainTime <= 5 ; rainTime += .5){
            for (double streetRate = 2; streetRate <= 5 ; streetRate += .2){
                model->transitions[11].flowRate = streetRate;
                model->transitions[12].flowRate = streetRate;
                model->transitions[13].flowRate = streetRate;
                
                time = model->transitions[2].time;
                
                InitializeModel(model);
                
                Marking* initialMarking = createInitialMarking(model);
                checkEnabled(model, initialMarking, 0 ,0);
                
                TimedDiagram::getInstance()->clear();
                TimedDiagram::getInstance()->setModel(model);
                
                gettimeofday(&gt1, NULL);
                TimedDiagram::getInstance()->generateDiagram(initialMarking);
                gettimeofday(&gt2, NULL);
                seconds  += gt2.tv_sec  - gt1.tv_sec;
                useconds += gt2.tv_usec - gt1.tv_usec;
                
                
                //				TimedDiagram::getInstance()->scale = 20;
                //				if (argc == 4){
                //					std::cout << "Writing the debug region diagram...." << std::endl;
                //					std::stringstream ss;
                //					ss << argv[1] << "_rd";
                //					TimedDiagram::getInstance()->saveDiagram(ss.str());
                //					cv::Mat flipped;
                //					cv::flip(TimedDiagram::getInstance()->debugImage, flipped, 0);
                //					cv::imshow("test", flipped);
                //					cv::waitKey(0);
                //				}
                
                
                ModelChecker modelChecker(model, TimedDiagram::getInstance());
                //				modelChecker.debugImage = TimedDiagram::getInstance()->debugImage;
                modelChecker.scale = TimedDiagram::getInstance()->scale;
                
                gettimeofday(&gt1, NULL);
                res = modelChecker.until(phi1, psi, bound, 0);
                gettimeofday(&gt2, NULL);
                seconds  += gt2.tv_sec  - gt1.tv_sec;
                useconds += gt2.tv_usec - gt1.tv_usec;
                
                
                cout << "-------------------------------------" << endl;
                cout << "The satisfaction set of Until formula: " << endl;
                for (unsigned int i = 0; i < res->intervals.size(); i++){
                    res->intervals[i].end = res->intervals[i].end - time;
                    res->intervals[i].start = res->intervals[i].start - time;
                }
                
                res->print(cout);
                res->print(oInt);
                
                oFile << modelChecker.calcProb(res, foldedNormal, 0) << "	";
                
                
                //				if (argc == 4){
                //					cv::Mat flipped;
                //					cv::flip(modelChecker.debugImage, flipped, 0);
                //					cv::imshow("test", flipped);
                //					cv::waitKey(0);
                //				}
                
                
                res->clear();
            }
            oInt << endl;
            oFile << endl;
        }
    }
    
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cout << "Time : " << mtime<< "ms" << std::endl;
    
    FILE* gnuplotPipe;
    if (gnuplotPipe = popen("gnuplot -persist","w"))
    {
        fprintf(gnuplotPipe,"set pm3d \n unset surface\n "
                //				"set xlabel \"Start Rain Time\" \n"
                "set xlabel \"Street Rate\" font \"Times-Roman, 20\"\n"
                "set zlabel \"Probability\" font \"Times-Roman, 20\"\n"
                "set ylabel \"mu \"\n"
                //				"set xtics (\"0.5\" 0, \"5\" 9) \n "
                "set xtics (\"2\" 0, \"5\" 14) font \"Times-Roman, 20\"\n "
                "set ytics (\"0.1\" 0, \"4\" 38) font \"Times-Roman, 20\"\n "
                "set ztics font \"Times-Roman, 20\"\n "
                "set zrange [0:1] \n"
                "set palette rgbformulae 33,13,10 \n"
                "set terminal wxt size 1300,800 \n");
        
        fflush(gnuplotPipe);
        fprintf(gnuplotPipe,"splot \"%s\" matrix with lines\n", "./results.out");
        fflush(gnuplotPipe);
        fprintf(gnuplotPipe,"exit \n");
        fclose(gnuplotPipe);
    }
    
    
    return 0;
}

int mainVersPomp(int argc, char *argv[]) {
    Model *model;
    
    if (argc < 3) {
        //std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName STD(optional, should be present if you want the STD as a JPG file.)" << std::endl;
        std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName" << std::endl;
        return 0;
    }
    
    model = ReadModel(argv[1]);
    
    unsigned int pIndex, disPIndex, NormTIndex;
    
    for (int i = 0; i < model->N_places; i++) {
        //if (model->places[i].type != PT_FLUID) continue;
        if (model->places[i].type == PT_FLUID && strncmp(model->places[i].id, argv[2], strlen(argv[2])) == 0)  pIndex = i;
        if (model->places[i].type == PT_DISCRETE && strncmp(model->places[i].id, "normalAgain_P", strlen("normalAgain_P")) == 0)  disPIndex = i;
    }
    
    for (int i = 0; i < model->N_transitions; i++) {
        if (strncmp(model->transitions[i].id, "startNormal_T", strlen("startNormal_T")) == 0)  NormTIndex= i;
    }
    
    double time = 0; //model->transitions[tIndex].time;
    mu =  .5;
    
    model->MaxTime = 40.0;
    
    long mtime;
    long seconds = 0;
    long useconds = 0;
    timeval gt1, gt2;
    
    
    AtomContNegFormula* phi1 = new AtomContNegFormula(pIndex, 0.01 );
    AtomDisFormula* psi = new AtomDisFormula(disPIndex, 1);
    
    ofstream oFile("./results.out");
    ofstream oInt("./intervals.out");
    for (double T = 30; T < 31; T += 1){
        Interval bound(0, T);
        IntervalSet* res;
        for (mu = .1; mu <= 4; mu += .1){
            //			for (double rainTime = .5; rainTime <= 5 ; rainTime += .5){
            for (double verSlipRate = .25; verSlipRate  <= 5 ; verSlipRate += .25){
                model->transitions[17].flowRate = verSlipRate;
                
                time = model->transitions[2].time;
                
                InitializeModel(model);
                
                Marking* initialMarking = createInitialMarking(model);
                checkEnabled(model, initialMarking, 0 ,0);
                
                TimedDiagram::getInstance()->clear();
                TimedDiagram::getInstance()->setModel(model);
                
                gettimeofday(&gt1, NULL);
                TimedDiagram::getInstance()->generateDiagram(initialMarking);
                gettimeofday(&gt2, NULL);
                seconds  += gt2.tv_sec  - gt1.tv_sec;
                useconds += gt2.tv_usec - gt1.tv_usec;
                
                //				TimedDiagram::getInstance()->scale = 20;
                //				if (argc == 4){
                //					std::cout << "Writing the debug region diagram...." << std::endl;
                //					std::stringstream ss;
                //					ss << argv[1] << "_rd";
                //					TimedDiagram::getInstance()->saveDiagram(ss.str());
                //					cv::Mat flipped;
                //					cv::flip(TimedDiagram::getInstance()->debugImage, flipped, 0);
                //					cv::imshow("test", flipped);
                //					cv::waitKey(0);
                //				}
                
                
                ModelChecker modelChecker(model, TimedDiagram::getInstance());
                //				modelChecker.debugImage = TimedDiagram::getInstance()->debugImage;
                modelChecker.scale = TimedDiagram::getInstance()->scale;
                
                
                gettimeofday(&gt1, NULL);
                res = modelChecker.until(phi1, psi, bound, 0);
                gettimeofday(&gt2, NULL);
                seconds  += gt2.tv_sec  - gt1.tv_sec;
                useconds += gt2.tv_usec - gt1.tv_usec;
                
                
                
                cout << "-------------------------------------" << endl;
                cout << "The satisfaction set of Until formula: " << endl;
                for (unsigned int i = 0; i < res->intervals.size(); i++){
                    res->intervals[i].end = res->intervals[i].end - time;
                    res->intervals[i].start = res->intervals[i].start - time;
                }
                
                res->print(cout);
                res->print(oInt);
                
                oFile << modelChecker.calcProb(res, foldedNormal, 0) << "	";
                
                
                //				if (argc == 4){
                //					cv::Mat flipped;
                //					cv::flip(modelChecker.debugImage, flipped, 0);
                //					cv::imshow("test", flipped);
                //					cv::waitKey(0);
                //				}
                
                
                res->clear();
            }
            oInt << endl;
            oFile << endl;
        }
    }
    
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cout << "Time : " << mtime<< "ms" << std::endl;
    
    FILE* gnuplotPipe;
    if (gnuplotPipe = popen("gnuplot -persist","w"))
    {
        fprintf(gnuplotPipe,"set pm3d \n unset surface\n "
                //				"set xlabel \"Start Rain Time\" \n"
                "set xlabel \"Fresh Sludge Rate\" font \"Times-Roman, 20\"\n"
                "set zlabel \"Probability\" font \"Times-Roman, 20\"\n"
                "set ylabel \"mu \"\n"
                //				"set xtics (\"0.5\" 0, \"5\" 9) \n "
                "set xtics (\"0.25\" 0, \"1\" 3,\"2\" 7,\"3\" 11,\"4\" 15,\"5\" 19) font \"Times-Roman, 20\"\n "
                "set ytics (\"0.1\" 0, \"4\" 38) font \"Times-Roman, 20\"\n "
                "set ztics font \"Times-Roman, 20\"\n "
                "set zrange [0:1] \n"
                "set palette rgbformulae 33,13,10 \n"
                "set terminal wxt size 1300,800 \n");
        
        fflush(gnuplotPipe);
        fprintf(gnuplotPipe,"splot \"%s\" matrix with lines\n", "./results.out");
        fflush(gnuplotPipe);
        fprintf(gnuplotPipe,"exit \n");
        fclose(gnuplotPipe);
    }
    
    
    return 0;
}


int mainSurPlusPomp(int argc, char *argv[]) {
    Model *model;
    
    if (argc < 3) {
        //std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName STD(optional, should be present if you want the STD as a JPG file.)" << std::endl;
        std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName" << std::endl;
        return 0;
    }
    
    model = ReadModel(argv[1]);
    
    unsigned int pIndex, disPIndex, NormTIndex;
    
    for (int i = 0; i < model->N_places; i++) {
        //if (model->places[i].type != PT_FLUID) continue;
        if (model->places[i].type == PT_FLUID && strncmp(model->places[i].id, argv[2], strlen(argv[2])) == 0)  pIndex = i;
        if (model->places[i].type == PT_DISCRETE && strncmp(model->places[i].id, "normalAgain_P", strlen("normalAgain_P")) == 0)  disPIndex = i;
    }
    
    for (int i = 0; i < model->N_transitions; i++) {
        if (strncmp(model->transitions[i].id, "startNormal_T", strlen("startNormal_T")) == 0)  NormTIndex= i;
    }
    
    double time = 0; //model->transitions[tIndex].time;
    mu =  .5;
    
    model->MaxTime = 40.0;
    
    
    
    AtomContNegFormula* phi1 = new AtomContNegFormula(pIndex, 0.01 );
    AtomDisFormula* psi = new AtomDisFormula(disPIndex, 1);
    
    ofstream oFile("./results.out");
    ofstream oInt("./intervals.out");
    
    long mtime;
    long seconds = 0;
    long useconds = 0;
    timeval gt1, gt2;
    
    for (double T = 30; T < 31; T += 1){
        Interval bound(0, T);
        IntervalSet* res;
        for (mu = .1; mu <= 4; mu += .1){
            //			for (double rainTime = .5; rainTime <= 5 ; rainTime += .5){
            for (double surPlusRate = .25; surPlusRate   <= 5 ; surPlusRate  += .25){
                model->transitions[18].flowRate = surPlusRate ;
                
                time = model->transitions[2].time;
                
                InitializeModel(model);
                
                Marking* initialMarking = createInitialMarking(model);
                checkEnabled(model, initialMarking, 0 ,0);
                
                TimedDiagram::getInstance()->clear();
                TimedDiagram::getInstance()->setModel(model);
                
                gettimeofday(&gt1, NULL);
                TimedDiagram::getInstance()->generateDiagram(initialMarking);
                gettimeofday(&gt2, NULL);
                seconds  += gt2.tv_sec  - gt1.tv_sec;
                useconds += gt2.tv_usec - gt1.tv_usec;
                
                
                //				TimedDiagram::getInstance()->scale = 20;
                //				if (argc == 4){
                //					std::cout << "Writing the debug region diagram...." << std::endl;
                //					std::stringstream ss;
                //					ss << argv[1] << "_rd";
                //					TimedDiagram::getInstance()->saveDiagram(ss.str());
                //					cv::Mat flipped;
                //					cv::flip(TimedDiagram::getInstance()->debugImage, flipped, 0);
                //					cv::imshow("test", flipped);
                //					cv::waitKey(0);
                //				}
                
                
                ModelChecker modelChecker(model, TimedDiagram::getInstance());
                //				modelChecker.debugImage = TimedDiagram::getInstance()->debugImage;
                modelChecker.scale = TimedDiagram::getInstance()->scale;
                
                gettimeofday(&gt1, NULL);
                res = modelChecker.until(phi1, psi, bound, 0);
                gettimeofday(&gt2, NULL);
                seconds  += gt2.tv_sec  - gt1.tv_sec;
                useconds += gt2.tv_usec - gt1.tv_usec;
                
                
                cout << "-------------------------------------" << endl;
                cout << "The satisfaction set of Until formula: " << endl;
                for (unsigned int i = 0; i < res->intervals.size(); i++){
                    res->intervals[i].end = res->intervals[i].end - time;
                    res->intervals[i].start = res->intervals[i].start - time;
                }
                
                res->print(cout);
                res->print(oInt);
                
                oFile << modelChecker.calcProb(res, foldedNormal, 0) << "	";
                
                //				if (argc == 4){
                //					cv::Mat flipped;
                //					cv::flip(modelChecker.debugImage, flipped, 0);
                //					cv::imshow("test", flipped);
                //					cv::waitKey(0);
                //				}
                //
                
                res->clear();
            }
            oInt << endl;
            oFile << endl;
        }
    }
    
    
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cout << "Time : " << mtime<< "ms" << std::endl;
    
    
    FILE* gnuplotPipe;
    if (gnuplotPipe = popen("gnuplot -persist","w"))
    {
        fprintf(gnuplotPipe,"set pm3d \n unset surface\n "
                //				"set xlabel \"Start Rain Time\" \n"
                "set xlabel \"Surplus Sludge Rate\" font \"Times-Roman, 20\"\n"
                "set zlabel \"Probability\" font \"Times-Roman, 20\"\n"
                "set ylabel \"mu \"\n"
                //				"set xtics (\"0.5\" 0, \"5\" 9) \n "
                "set xtics (\"0.25\" 0, \"5\" 19) font \"Times-Roman, 20\"\n "
                "set ytics (\"0.1\" 0, \"4\" 38) font \"Times-Roman, 20\"\n "
                "set ztics font \"Times-Roman, 20\"\n "
                "set zrange [0:1] \n"
                "set palette rgbformulae 33,13,10 \n"
                "set terminal wxt size 1300,800 \n");
        
        fflush(gnuplotPipe);
        fprintf(gnuplotPipe,"splot \"%s\" matrix with lines\n", "./results.out");
        fflush(gnuplotPipe);
        fprintf(gnuplotPipe,"exit \n");
        fclose(gnuplotPipe);
    }
    
    
    return 0;
}



int main(int argc, char *argv[]) {
    //	long mtime, seconds, useconds;
    //	timeval gt1, gt2;
    //	gettimeofday(&gt1, NULL);
    
    //	gettimeofday(&gt2, NULL);
    //	seconds  = gt2.tv_sec  - gt1.tv_sec;
    //	useconds = gt2.tv_usec - gt1.tv_usec;
    //	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    //	std::cout << "Time mainStocRainTime: " << mtime<< "ms" << std::endl;
    
    mainPriorities(argc, argv);
    //mainIntegrity(argc, argv);
    
    //	mainStocRainTime(argc, argv);
    //	mainStreetRate(argc, argv);
    //	mainVersPomp(argc, argv);
    //	mainSurPlusPomp(argc, argv);
    //	argv[1] = "./models/WT-zandFailure.m";
    //	mainFailure(argc, argv);
}


int main2(int argc, char *argv[]) {
    
    Model *model;
    
    if (argc < 3) {
        //std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName STD(optional, should be present if you want the STD as a JPG file.)" << std::endl;
        std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName" << std::endl;
        return 0;
    }
    
    model = ReadModel(argv[1]);
    
    unsigned int pIndex, disPIndex;
    
    for (int i = 0; i < model->N_places; i++) {
        //if (model->places[i].type != PT_FLUID) continue;
        if (model->places[i].type == PT_FLUID && strncmp(model->places[i].id, argv[2], strlen(argv[2])) == 0)  pIndex = i;;
        if (model->places[i].type == PT_DISCRETE && strncmp(model->places[i].id, "Input1On", strlen("Input1On")) == 0)  disPIndex = i;;
        
    }
    
    unsigned int tIndex = 0;
    for (int i = 0; i < model->N_transitions; i++) {
        if (model->transitions[i].type != TT_DETERMINISTIC) continue;
        if (strncmp(model->transitions[i].id, "failure", strlen("failure")) != 0) continue;
        
        tIndex = i;
    }
    
    double time = model->transitions[tIndex].time;
    cout << "time: " <<  time << endl;
    cout << "pIndex: " <<  pIndex << endl;
    cout << "disPIndex: " <<  disPIndex << endl;
    
    
    InitializeModel(model);
    //Input and output list of places and transitions are created and sorted wrt to their priority and share.
    
    //expPrarameter = atoi(model->transitions[gTransitionId(model)].distr);
    
    mu =  2;
    sigma = 1;
    expPrarameter = 2;
    
    model->MaxTime = 20.0;
    Marking* initialMarking = createInitialMarking(model);
    checkEnabled(model, initialMarking, 0 ,0);
    
    
    long mtime, seconds, useconds;
    
    TimedDiagram::getInstance()->setModel(model);
    
    
    timeval gt1, gt2;
    gettimeofday(&gt1, NULL);
    TimedDiagram::getInstance()->generateDiagram(initialMarking);
    gettimeofday(&gt2, NULL);
    
    seconds  = gt2.tv_sec  - gt1.tv_sec;
    useconds = gt2.tv_usec - gt1.tv_usec;
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    
    std::cout << "Number of regions: " <<TimedDiagram::getInstance()->getNumberOfRegions() << std::endl;
    std::cout << "Time to generate STD: " << mtime<< "ms" << std::endl;
    
    
    //	TimedDiagram::getInstance()->scale = 40;
    //	if (argc == 4){
    //		std::cout << "Writing the debug region diagram...." << std::endl;
    //		std::stringstream ss;
    //		ss << argv[1] << "_rd";
    //		TimedDiagram::getInstance()->saveDiagram(ss.str());
    //		cv::Mat flipped;
    //		cv::flip(TimedDiagram::getInstance()->debugImage, flipped, 0);
    //		cv::imshow("test", flipped);
    //		cv::waitKey(0);
    //	}
    
    
    
    
    ModelChecker modelChecker(model, TimedDiagram::getInstance());
    //cout << "------------------------------"<<TimedDiagram::getInstance()->dtrmEventList[0]->postRegionMarking->tokens[7] <<"--------------------------" << endl;
    
    
    AtomContNegFormula* phi1 = new AtomContNegFormula(pIndex, 1);
    AtomContPosFormula* phi2 = new AtomContPosFormula(pIndex, 10);
    
    AtomContPosFormula* tt = new AtomContPosFormula(pIndex, -1);
    
    //	AtomDisFormula* psi2_dis = new AtomDisFormula(disPIndex, 1);
    //	AtomContFormula* psi2_cont = new AtomContPosFormula(pIndex, 10);
    //	ADFormula* psi2 = new ADFormula(psi2_cont, psi2_dis);
    
    
    //	modelChecker.debugImage = TimedDiagram::getInstance()->debugImage;
    modelChecker.scale = TimedDiagram::getInstance()->scale;
    
    ofstream oFile("./results.out");
    ofstream oInt("./intervals.out");
    for (double T = time + 10; T <= time + 10; T += 1){
        Interval bound(0, T);
        IntervalSet* i1 = modelChecker.until(tt, phi1, bound, time);
        //		IntervalSet* i2 = modelChecker.until(tt, phi2, bound, time);
        
        //		i1->print(cout);
        //		i2->print(cout);
        
        IntervalSet* res = i1; //->unionWith(i2);
        
        cout << "The satisfaction set of Until formula: " << endl;
        for (unsigned int i = 0; i < res->intervals.size(); i++){
            res->intervals[i].end = res->intervals[i].end - time;
            res->intervals[i].start = res->intervals[i].start - time;
        }
        
        res->print(cout);
        res->print(oInt);
        
        oFile << modelChecker.calcProb(res, gamma, 0) << endl;
    }
    
    return 0;
}




int main3(int argc, char *argv[]) {
    
    Model *model;
    
    if (argc < 3) {
        //std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName STD(optional, should be present if you want the STD as a JPG file.)" << std::endl;
        std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName" << std::endl;
        return 0;
    }
    
    model = ReadModel(argv[1]);
    InitializeModel(model);
    //Input and output list of places and transitions are created and sorted wrt to their priority and share.
    
    expPrarameter = 3; //atoi(model->transitions[gTransitionId(model)].distr);
    
    model->MaxTime = 10.0;
    Marking* initialMarking = createInitialMarking(model);
    
    long mtime, seconds, useconds;
    
    TimedDiagram::getInstance()->setModel(model);
    
    timeval gt1, gt2;
    gettimeofday(&gt1, NULL);
    TimedDiagram::getInstance()->generateDiagram(initialMarking);
    gettimeofday(&gt2, NULL);
    
    seconds  = gt2.tv_sec  - gt1.tv_sec;
    useconds = gt2.tv_usec - gt1.tv_usec;
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    
    std::cout << "Number of regions: " <<TimedDiagram::getInstance()->getNumberOfRegions() << std::endl;
    std::cout << "Time to generate STD: " << mtime<< "ms" << std::endl;
    
    
    //	TimedDiagram::getInstance()->scale = 20;
    //	if (argc == 4){
    //		std::cout << "Writing the debug region diagram...." << std::endl;
    //		std::stringstream ss;
    //		ss << argv[1] << "_rd";
    //		TimedDiagram::getInstance()->saveDiagram(ss.str());
    //		cv::Mat flipped;
    //		cv::flip(TimedDiagram::getInstance()->debugImage, flipped, 0);
    //		cv::imshow("test", flipped);
    //		cv::waitKey(0);
    //	}
    
    std::cout << "starting measure computation..." << std::endl;
    
    //Model* model = TimedDiagram::getInstance()->model;
    for (int i = 0; i < model->N_places; i++) {
        if (model->places[i].type != PT_FLUID) continue;
        if (strncmp(model->places[i].id, argv[2], strlen(argv[2])) != 0) continue;
        
        pIndex = i;
    }
    
    
    std::ofstream oFile;
    oFile.open(argv[2], std::ios::out);
    
    mtime = 0; seconds = 0; useconds = 0;
    timeval t0, t1;
    
    for (double t = .02; t <= model->MaxTime + .01; t += .5){
        for (amount = .05; amount <= model->places[pIndex].f_bound; amount += .2){
            //std::cout << "t = " << t << "--amount = " << amount << std::endl;
            gettimeofday(&t0, NULL);
            double p = TimedDiagram::getInstance()->calProbAtTime(t, expo, propertyTest);
            gettimeofday(&t1, NULL);
            
            seconds  += t1.tv_sec  - t0.tv_sec;
            useconds += t1.tv_usec - t0.tv_usec;
            
            oFile << " "<< p;
            //std::cout << " "<< p;
        }
        oFile << std::endl;
        //std::cout << std::endl;
    }
    
    mtime = ((seconds * 1000 + useconds/1000.0) + 0.5);
    std::cout << "total time computing measures:  " << mtime << "ms" << std::endl;
    
    oFile.close();
    FILE* gnuplotPipe;
    if (gnuplotPipe = popen("gnuplot -persist","w"))
    {
        fprintf(gnuplotPipe,"set pm3d \n unset surface\n "
                "set xlabel \"Storage Content [m^3]\" \n"
                "set zlabel \"Probability\"\n"
                "set ylabel \"time [hours]\"\n"
                "set xtics (\"0\" 0, \"2\" 10, \"4\" 20,\"6\" 30, \"8\" 40) \n "
                "set ytics (\"0\" 0, \"20\" 40, \"40\" 80,\"60\" 120, \"80\" 160, \"100\" 200) \n "
                "set palette rgbformulae 33,13,10 \n");
        
        fflush(gnuplotPipe);
        fprintf(gnuplotPipe,"splot \"%s\" matrix with lines\n",argv[2]);
        fflush(gnuplotPipe);
        fprintf(gnuplotPipe,"exit \n");
        fclose(gnuplotPipe);
    }
    
    
    
    
    //getch();
    
    
    return 0;
}




//******************************************************************************************************************************************
//for (unsigned int i = 0; i <  TimedDiagram::getInstance()->regionList.size(); i++){
//	Region* reg =  TimedDiagram::getInstance()->regionList[i];
//	std::cout <<"-----------------------------"<< std::endl;
//	std::cout << "Region Number: " << i << ": ";
//	reg->print(std::cout);
//	std::cout << std::endl;
//
//	for (unsigned int j = 0; j < reg->successors->size(); j++)
//		reg->successors->at(j)->print(std::cout);
//}
//
//for (unsigned int i = 0; i <  TimedDiagram::getInstance()->dtrmEventList.size(); i++){
//	DtrmEvent* dtrmEvent =  TimedDiagram::getInstance()->dtrmEventList[i];
//	std::cout <<"-----------------------------"<< std::endl;
//	std::cout << "DtrmEvent time: " << i << ": " << dtrmEvent->time << std::endl;
//
//
//	for (unsigned int j = 0; j < dtrmEvent->nextRegions->size(); j++)
//		dtrmEvent->nextRegions->at(j)->print(std::cout);
//
//
//	std::cout << "Next DtrmEvent time: " << ((dtrmEvent->nextDtrmEvent!=0)? dtrmEvent->nextDtrmEvent->time : -1) << std::endl;
//
//}
//******************************************************************************************************************************************

//int main(int argc, char *argv[]){
//	IntervalSet* set1 = new IntervalSet();
//	set1 = set1->unionWith(Interval(3, 5));
//	set1->intervals.push_back(Interval(4, 5));
//	set1->intervals.push_back(Interval(7, 8));
//	set1->intervals.push_back(Interval(9, 10));
//	set1->intervals.push_back(Interval(12, 15));
//
//	IntervalSet* set2 = new IntervalSet();
//	set2->intervals.push_back(Interval(1, 2));
//	set2->intervals.push_back(Interval(4, 11));
//	set2->intervals.push_back(Interval(13, 16));
//
//	set1->print();
//	set2->print();
//	IntervalSet* intersection = set1->intersect(*set2);
//	intersection->print();
//
//	IntervalSet* unionres= set1->unionWith(*set2);
//	unionres->print();

//	set1->print();
//	set2->print();
//
//	std::cout << "--------------------------" << std::endl;
//	IntervalSet* minusres= set1->minus(*set2);
//	minusres->print();
//	minusres= set2->minus(*set1);
//	minusres->print();

//}


//int main(int argc, char *argv[]) {
//
//	Model *model;
//
//	if (argc < 3) {
//		//std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName STD(optional, should be present if you want the STD as a JPG file.)" << std::endl;
//		std::cout << "\n\n Error! Usage: DFPN ModelFile.m FluidPlaceName" << std::endl;
//		return 0;
//	}
//
//	model = ReadModel(argv[1]);
//	InitializeModel(model);
//	//Input and output list of places and transitions are created and sorted wrt to their priority and share.
//
//	expPrarameter = atoi(model->transitions[gTransitionId(model)].distr);
//
//	model->MaxTime = 100.0;
//	Marking* initialMarking = createInitialMarking(model);
//
//	long mtime, seconds, useconds;
//
//	TimedDiagram::getInstance()->setModel(model);
//
//	timeval gt1, gt2;
//	gettimeofday(&gt1, NULL);
//	TimedDiagram::getInstance()->generateDiagram(initialMarking);
//	gettimeofday(&gt2, NULL);
//
//	seconds  = gt2.tv_sec  - gt1.tv_sec;
//	useconds = gt2.tv_usec - gt1.tv_usec;
//	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
//
//	std::cout << "Number of regions: " <<TimedDiagram::getInstance()->getNumberOfRegions() << std::endl;
//	std::cout << "Time to generate STD: " << mtime<< "ms" << std::endl;
//
//
//	if (argc == 4){
//		std::cout << "Writing the debug region diagram...." << std::endl;
//		std::stringstream ss;
//		ss << argv[1] << "_rd";
//		TimedDiagram::getInstance()->saveDiagram(ss.str());
//	}
//
//	std::cout << "starting measure computation..." << std::endl;
//
//	//Model* model = TimedDiagram::getInstance()->model;
//	for (int i = 0; i < model->N_places; i++) {
//		if (model->places[i].type != PT_FLUID) continue;
//		if (strncmp(model->places[i].id, argv[2], strlen(argv[2])) != 0) continue;
//
//		pIndex = i;
//	}
//
//
//	std::ofstream oFile;
//	oFile.open(argv[2], std::ios::out);
//
//	mtime = 0; seconds = 0; useconds = 0;
//	timeval t0, t1;
//
//	for (double t = .02; t <= model->MaxTime + .01; t += .5){
//		for (amount = .05; amount <= model->places[pIndex].f_bound; amount += .2){
//			//std::cout << "t = " << t << "--amount = " << amount << std::endl;
//			gettimeofday(&t0, NULL);
//			double p = TimedDiagram::getInstance()->calProbAtTime(t, expo, propertyTest);
//			gettimeofday(&t1, NULL);
//
//			seconds  += t1.tv_sec  - t0.tv_sec;
//			useconds += t1.tv_usec - t0.tv_usec;
//
//			oFile << " "<< p;
//			//std::cout << " "<< p;
//		}
//		oFile << std::endl;
//		//std::cout << std::endl;
//	}
//
//	mtime = ((seconds * 1000 + useconds/1000.0) + 0.5);
//	std::cout << "total time computing measures:  " << mtime << "ms" << std::endl;
//
//	oFile.close();
//	FILE* gnuplotPipe;
//	if (gnuplotPipe = popen("gnuplot -persist","w"))
//	{
//		fprintf(gnuplotPipe,"set pm3d \n unset surface\n "
//				"set xlabel \"Storage Content [m^3]\" \n"
//				"set zlabel \"Probability\"\n"
//				"set ylabel \"time [hours]\"\n"
//				"set xtics (\"0\" 0, \"2\" 10, \"4\" 20,\"6\" 30, \"8\" 40) \n "
//				"set ytics (\"0\" 0, \"20\" 40, \"40\" 80,\"60\" 120, \"80\" 160, \"100\" 200) \n "
//				"set palette rgbformulae 33,13,10 \n");
//
//		fflush(gnuplotPipe);
//		fprintf(gnuplotPipe,"splot \"%s\" matrix with lines\n",argv[2]);
//		fflush(gnuplotPipe);
//		fprintf(gnuplotPipe,"exit \n");
//		fclose(gnuplotPipe);
//	}
//
//
//
//
//	//getch();
//
//
//	return 0;
//}

//int main(int argc, char *argv[]) {
//	Model *model;
//
//	if (argc != 2) {
//		std::cout << "\n\n Error! Usage: DFPN ModelFile.m " << std::endl;
//		return 0;
//	}
//
//	model = ReadModel(argv[1]);
//	InitializeModel(model);
//	//Input and output list of places and transitions are created and sorted wrt to their priority and sahre.
//
//	model->MaxTime = 100.0;
//	Marking* initialMarking = createInitialMarking(model);
//
//
//
//	TimedDiagram::getInstance()->setModel(model);
//
//	clock_t gt1, gt2;
//	gt1= clock();
//	TimedDiagram::getInstance()->generateDiagram(initialMarking);
//	gt2 = clock();
//
//	std::cout << "Number of regions: " <<TimedDiagram::getInstance()->getNumberOfRegions() << std::endl;
//	std::cout << "Cmp of regions: " <<TimedDiagram::getInstance()->getNumberOfRegions() << std::endl;
//
//
//	std::cout << "Writting the debug region diagram...." << std::endl;
//	std::stringstream ss;
//	ss << argv[1] << "_rd";
//	//TimedDiagram::getInstance()->saveDiagram(ss.str());
//
//	clock_t time = 0;
//	clock_t t0, t1;
//
//	std::cout << "starting measure computation..." << std::endl;
//
//	//Model* model = TimedDiagram::getInstance()->model;
//	for (int i = 0; i < model->N_places; i++) {
//		if (model->places[i].type != PT_FLUID) continue;
//		if (strncmp(model->places[i].id, "stor", strlen("stor")) != 0) continue;
//
//		pIndex = i;
//		std::stringstream name;
//		name << argv[1] << "_" << model->places[i].id << ".out";
//		std::ofstream oFile;
//		oFile.open(name.str().c_str(), std::ios::out);
//		//pIndex = 2;
//		for (double t = .02; t <= model->MaxTime + .01; t += .5){
//			for (amount = .05; amount <= model->places[i].f_bound; amount += .2){
//				//std::cout << "t = " << t << "--amount = " << amount << std::endl;
//				t0 = clock();
//				double p = TimedDiagram::getInstance()->calProbAtTime(t, spdf, propertyTest);
//				t1 = clock();
//
//				time += (t1 - t0);
//
//				oFile << " "<< p;
//				//std::cout << " "<< p;
//			}
//			oFile << std::endl;
//			//std::cout << std::endl;
//		}
//
//		oFile.close();
//		FILE* gnuplotPipe;
//		if (gnuplotPipe = popen("gnuplot -persist","w"))
//		{
//			fprintf(gnuplotPipe,"set pm3d \n unset surface\n "
//					"set xlabel \"Storage Content [m^3]\" \n"
//					"set zlabel \"Probability\"\n"
//					"set ylabel \"time [hours]\"\n"
//					"set xtics (\"0\" 0, \"2\" 10, \"4\" 20,\"6\" 30, \"8\" 40) \n "
//					"set ytics (\"0\" 0, \"20\" 40, \"40\" 80,\"60\" 120, \"80\" 160, \"100\" 200) \n "
//					"set palette rgbformulae 33,13,10 \n");
//
//			fflush(gnuplotPipe);
//			fprintf(gnuplotPipe,"splot \"%s\" matrix with lines\n",name.str().c_str());
//			fflush(gnuplotPipe);
//			fprintf(gnuplotPipe,"exit \n");
//			fclose(gnuplotPipe);
//		}
//	}
//
//
//
//	std::cout << "total time computing measures:  " << 100*time << "ms" << std::endl;
//
//	//getch();
//
//
//	return 0;
//}



//int main(int argc, char *argv[]) {
//	Model *model;
//
//	std::string dir = "D:\\Profiles\\GhasemiehH\\Scientific\\UTwente\\Research\\DFPN\\DFPN\\Test\\rate\\";
//	std::string fpreffix = "dsn_";
//	std::string fSuffix = "rates.m";
//	std::stringstream file;
//	std::stringstream ss;
//	std::stringstream name;
//
//	file << dir << "result.txt";
//	std::ofstream outFile (file.str());
//
//	for (int l = 2; l < 12; l++){
//
//		outFile << l << "\t";
//
//		file.str("");
//		if ( l==11)
//			file << dir << fpreffix << 2 << fSuffix;
//		else
//			file << dir << fpreffix << l << fSuffix;
//		model = ReadModel(file.str().c_str());
//		InitializeModel(model);
//		//Input and output list of places and transitions are created and sorted wrt to their priority and sahre.
//
//		model->MaxTime = 100.0;
//		Marking* initialMarking = createInitialMarking(model);
//
//
//
//		TimedDiagram::getInstance()->setModel(model);
//
//		clock_t gt1, gt2;
//		gt1= clock();
//		TimedDiagram::getInstance()->generateDiagram(initialMarking);
//		gt2 = clock();
//
//		outFile << TimedDiagram::getInstance()->getNumberOfRegions() << "\t" << gt2 - gt1 << "\t";
//
//
//		std::cout << "Writting the debug region diagram...." << std::endl;
//		ss.str("");
//		ss << file.str() << "_rd";
//		//TimedDiagram::getInstance()->saveDiagram(ss.str());
//
//		clock_t time = 0;
//		clock_t t0, t1;
//
//		std::cout << "starting measure computation..." << std::endl;
//
//		//Model* model = TimedDiagram::getInstance()->model;
//		for (int i = 0; i < model->N_places; i++) {
//			if (model->places[i].type != PT_FLUID) continue;
//			if (strncmp(model->places[i].id, "stor", strlen("stor")) != 0) continue;
//
//			pIndex = i;
//			name.str("");
//			name << file.str() << "_" << model->places[i].id << ".out";
//			std::ofstream oFile = std::ofstream(name.str());
//			//pIndex = 2;
//			for (double t = 0.01; t < model->MaxTime + .001 ; t += .5){
//				for (amount = 0; amount <= model->places[i].f_bound; amount += .2){
//					//std::cout << "t = " << t << "--amount = " << amount << std::endl;
//					t0 = clock();
//					double p = TimedDiagram::getInstance()->calProbAtTime(t, spdf, propertyTest);
//					t1 = clock();
//
//					time += (t1 - t0);
//
//					oFile << " "<< p;
//					//std::cout << " "<< p;
//				}
//				oFile << std::endl;
//				//std::cout << std::endl;
//			}
//		}
//
//		outFile << time << std::endl;
//	}
//
//	outFile.close();
//	getch();
//	return 0;
//}
//

