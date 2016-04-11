/*
 * Formula.h
 *
 *  Created on: Mar 1, 2013
 *      Author: hamed
 */

#ifndef FORMULA_H_
#define FORMULA_H_

enum Formula_type {DISCRITE, CONT_POS, COMBI, CONT_NEG};

class Formula {
    
public:
    virtual Formula_type getType() = 0;
    Formula(){};
    virtual ~Formula(){};
};


class AtomDisFormula : public Formula {
    
private:
    unsigned int placeIndex;
    int n;
    
public:
    AtomDisFormula(unsigned int _placeIndex, int _n): placeIndex(_placeIndex), n(_n){};
    Formula_type getType() {return DISCRITE;};
    
    int getN() const {
        return n;
    }
    
    unsigned int getPlaceIndex() const {
        return placeIndex;
    }
};

class AtomContFormula : public Formula{
private:
    unsigned int placeIndex;
    double c;
    
public:
    /**
     * @param placeIndex index of place \f$ x \f$.
     * @param c camparison value.
     */
    AtomContFormula(unsigned int _placeIndex, double _c): placeIndex(_placeIndex), c(_c){};
    
    
    double getC() const {
        return c;
    }
    
    unsigned int getPlaceIndex() const {
        return placeIndex;
    }
    
};


/*
 * \brief Contains atomic formula of type \f$ x < c \f$.
 */
class AtomContNegFormula : public AtomContFormula{
public:
    AtomContNegFormula(unsigned int _placeIndex, double _c): AtomContFormula(_placeIndex, _c){}
    Formula_type getType() {return CONT_NEG;}
};


/*
 * \brief Contains atomic formula of type \f$x \geq c \f$.
 */
class AtomContPosFormula : public AtomContFormula{
public:
    AtomContPosFormula(unsigned int _placeIndex, double _c): AtomContFormula(_placeIndex, _c){}
    Formula_type getType() {return CONT_POS;}
};



/*
 * Temporary class! because of deadline!
 */
class ADFormula: public Formula{
public:
    ADFormula(	AtomContFormula* cont, AtomDisFormula* dis): contFormula(cont), disFormula(dis){};
    
    Formula_type getType() {return COMBI;};
    AtomContFormula* contFormula;
    AtomDisFormula* disFormula;
    
};

#endif /* FORMULA_H_ */