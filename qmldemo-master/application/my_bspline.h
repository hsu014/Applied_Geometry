
#ifndef MYBSPLINE_H
#define MYBSPLINE_H


#include "../../gmlib-master/modules/parametrics/gmpcurve.h"

using namespace GMlib;


// B-Spline
// Hardcoded degree 2
template <typename T>
class MyBSpline : public PCurve<T,3> {
    GM_SCENEOBJECT(MyBSpline)

public:
    MyBSpline(const DVector<Vector<T,3>>& c);
    MyBSpline(const DVector<Vector<T,3>>& p, int n);
    MyBSpline( const MyBSpline<T>& copy );
    virtual ~MyBSpline();

    bool                isClosed() const override;

protected:
    // Virtual functions from PCurve, which have to be implemented locally
    void                eval( T t, int d, bool /*l*/ ) const override;
    T                   getStartP() const override;
    T                   getEndP()   const override;

    // Protected data for the curve

private:
    std::vector<T>          _t;                     // knot vector
    DVector<Vector<T,3>>    _c;                     // control points

    T           W(int d, int i, T t) const;         // Map domain to 0-1
    Vector<T,3> B(int i, T t) const;                // Basis functions for given t
    int         find_i(T t) const;                  // Find index i for given t
    void        make_knot_vector(int n, int d=2);   // Create the knot vector _t

}; // END class MyBSpline



// Constructor 1
template <typename T>
inline
    MyBSpline<T>::MyBSpline(const DVector<Vector<T,3>>& c) : PCurve<T,3>(20, 0, 0), _c(c) {
    std::cout << "Control points _c: " << _c << std::endl;

    make_knot_vector(c.getDim());
}



// Constructor 2
template <typename T>
inline
    MyBSpline<T>::MyBSpline(const DVector<Vector<T,3>>& p, int n) : PCurve<T,3>(20, 0, 0) {
    // make the c-vector, then same as above
}



// Copy constructor
template <typename T>
inline
    MyBSpline<T>::MyBSpline( const MyBSpline<T>& copy ) : PCurve<T,3>(copy) {
    _t = copy._t;
    _c = copy._c;
}



// Destructor
template <typename T>
MyBSpline<T>::~MyBSpline() {}



template <typename T>
bool MyBSpline<T>::isClosed() const {
    return false;
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
void MyBSpline<T>::eval( T t, int d, bool /*l*/ ) const {
    this->_p.setDim( d + 1 );

    int i = find_i(t);
    Vector<T,3> b = B(i, t);

    this->_p[0] = b[0]*_c[i-2] + b[1]*_c[i-1] + b[2]*_c[i];
}



template <typename T>
T MyBSpline<T>::getStartP() const {
    return _t[2];
}



template <typename T>
T MyBSpline<T>::getEndP() const {
    return _t[_c.getDim()];
}



template <typename T>
T MyBSpline<T>::W(int d, int i, T t) const {
    return (t - _t[i])/(_t[i+d] - _t[i]);
}



template <typename T>
Vector<T,3> MyBSpline<T>::B(int i, T t) const {
    Vector<T,3> b;

    b[0] = (1 - W(1, i, t))*(1 - W(2, i-1, t)); // b_i-2;

    b[1] = (1 - W(1, i, t))*(W(2, i-1, t)) +    // b_i-1;
           (W(1, i, t))*(1 - W(2, i, t));

    b[2] = (W(1, i, t))*(W(2, i, t));           // b_i;

    return b;
}



template <typename T>
int MyBSpline<T>::find_i(T t) const {
    int i = 2; // = d
    int n = _c.getDim();
    for (; i<n; i++) {
        if (_t[i+1] > t) {
            //std::cout << "Index: " << i << std::endl;
            return i;
        }
    }
    //std::cout << "Index: " << n-1 << std::endl;
    return n-1;
}



template <typename T>
void MyBSpline<T>::make_knot_vector(int n, int d) {
    for (int i = 0; i<=d; i++) _t.push_back(0);

    for (int i = d+1; i<n; i++) _t.push_back(i-d);

    for (int i = 0; i<=d; i++) _t.push_back(n-d);

    std::cout << "Knot vector _t: " << _t << std::endl;

}



#endif // MYBSPLINE_H
