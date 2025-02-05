
#ifndef MYMODELCURVE1_H
#define MYMODELCURVE1_H


#include "../../gmlib-master/modules/parametrics/gmpcurve.h"

using namespace GMlib;


// Model curve 1
// Lissajous curve
template <typename T>
class MyModelCurve1 : public PCurve<T,3> {
    GM_SCENEOBJECT(MyModelCurve1)

public:
    MyModelCurve1(T a, T b, T k_x, T k_y);
    MyModelCurve1( const MyModelCurve1<T>& copy );
    virtual ~MyModelCurve1();

    bool                isClosed() const override;

protected:
    // Virtual functions from PCurve, which have to be implemented locally
    void                eval(T t, int d, bool l) const override;
    T                   getStartP() const override;
    T                   getEndP()   const override;

    // Protected data for the curve
    T                   _a;
    T                   _b;
    T                   _k_x;
    T                   _k_y;

}; // END class MyModelCurve1



// Constructor
template <typename T>
inline
    MyModelCurve1<T>::MyModelCurve1(T a, T b, T k_x, T k_y) : PCurve<T,3>(20, 0, 0) {
    _a = a;
    _b = b;
    _k_x = k_x;
    _k_y = k_y;
}



// Copy constructor
template <typename T>
inline
    MyModelCurve1<T>::MyModelCurve1( const MyModelCurve1<T>& copy ) : PCurve<T,3>(copy) {
    _a = copy._a;
    _b = copy._b;
    _k_x = copy._k_x;
    _k_y = copy._k_y;
}



// The destructor
template <typename T>
MyModelCurve1<T>::~MyModelCurve1() {}



template <typename T>
bool MyModelCurve1<T>::isClosed() const {
    return true;
}



/*!
   *  Evaluation of the curve at a given parameter value
   *  To compute position and d derivatives at parameter value t on the curve.
   *  7 derivatives are implemented
   *
   *  \param  t[in]  The parameter value to evaluate at
   *  \param  d[in]  The number of derivatives to compute
   *  \param  l[in]  (dummy) because left and right are always equal
   */
template <typename T>
void MyModelCurve1<T>::eval( T t, int d, bool /*l*/ ) const {
    this->_p.setDim( d + 1 );

    this->_p[0][0] = _a * cos(_k_x * t);            // x
    this->_p[0][1] = _b * sin(_k_y * t);            // y
    this->_p[0][2] = 0.5 * (_a+_b) / 2 * sin(t);    // z
}



template <typename T>
T MyModelCurve1<T>::getStartP() const {
    return T(0);
}



template <typename T>
T MyModelCurve1<T>::getEndP() const {
    return T(M_2PI);
}



#endif // MYMODELCURVE1_H
