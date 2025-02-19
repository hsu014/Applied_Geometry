
#ifndef MYBSPLINE_H
#define MYBSPLINE_H


#include "../gmlib-master/modules/parametrics/gmpcurve.h"

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

    bool                    isClosed() const override;

protected:
    // Virtual functions from PCurve, which have to be implemented locally
    void                    eval( T t, int d, bool /*l*/ ) const override;
    T                       getStartP() const override;
    T                       getEndP()   const override;

private:
    std::vector<T>          _t;                                         // knot vector
    DVector<Vector<T,3>>    _c;                                         // control points
    int                     _n;                                         // Number of control points

    T                       W(int d, int i, T t) const;                 // Map domain to 0-1
    Vector<T,3>             B(int i, T t) const;                        // Basis functions for given t
    int                     find_i(T t) const;                          // Find index i for given t
    void                    make_knot_vector(int n, int d, T s, T e);   // Create the knot vector _t

}; // END class MyBSpline



// Constructor 1
template <typename T>
inline
    MyBSpline<T>::MyBSpline(const DVector<Vector<T,3>>& c) : PCurve<T,3>(20, 0, 0), _c(c) {

    _n = c.getDim();

    int d = 2;                                      // Degree of b-spline hardcoded as 2
    make_knot_vector(_n, d, 0, _n-d);
}



// Constructor 2
template <typename T>
inline
    MyBSpline<T>::MyBSpline(const DVector<Vector<T,3>>& p, int n) : PCurve<T,3>(20, 0, 0), _n(n) {

    int m = p.getDim();                             // Number of points
    std::vector<T> x(m);                            // Total distance traveled in each point

    x[0] = 0;
    for (int i = 1; i < m; i++) {
        x[i] = x[i-1] + (p[i] - p[i-1]).getLength();
    }

    make_knot_vector(_n, 2, x[0], x.back());

    DMatrix<T> A(m, _n, T(0));
    for (int j=0; j<m; j++) {
        int i = find_i(x[j]);

        auto b = B(i, x[j]);
        A[j][i-2] = b[0];
        A[j][i-1] = b[1];
        A[j][i] = b[2];
    }

    DMatrix<T> A_t(A);
    A_t.transpose();

    DMatrix<T> B = A_t*A;

    DVector<Vector<T,3>> y = A_t*p;

    DMatrix<T> B_inv(B);
    B_inv.invert();

    DVector<Vector<T,3>> c = B_inv * y;
    std::cout << "c: " << c << std::endl;

    _c = c;
}



// Copy constructor
template <typename T>
inline
    MyBSpline<T>::MyBSpline(const MyBSpline<T>& copy ) : PCurve<T,3>(copy) {
    _t = copy._t;
    _c = copy._c;
    _n = copy._n;
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

    // std::cout << "t: " << t << ". b: " << b << std::endl;

    // const int d = 2;
    // Vector<T,(d+1)> b2;

    // b2[0] = 1;
    // for (int j=1; j<=d; j++) {
    //     b2[j] = W(j, i, t) * b2[j-1];

    //     for (int k=j-1; k>0; k--) {
    //         b2[k] = W(j, i-j+k, t) * b2[k-1] +
    //                 (1 - W(j, i-j+k+1, t)) * b2[k];
    //     }
    //     b2[0] = (1 - W(j, i-j, t)) * b2[0];
    // }

    // std::cout << "t: " << t << ". b: " << b2 << std::endl;


    return b;
}



template <typename T>
int MyBSpline<T>::find_i(T t) const {
    int i = 2;

    for (; i<_n; i++) {
        if (_t[i+1] > t) {
            return i;
        }
    }

    return _n-1;
}



template <typename T>
void MyBSpline<T>::make_knot_vector(int n, int d, T s, T e) {
    for (int i = 0; i<=d; i++) _t.push_back(s);

    for (int i = d+1; i<n; i++) _t.push_back((e-s)/(n-d)*(i-d));

    for (int i = 0; i<=d; i++) _t.push_back(e);

    // std::cout << "Knot vector _t: " << _t << std::endl;
}



#endif // MYBSPLINE_H
