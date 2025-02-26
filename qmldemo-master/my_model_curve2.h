
#ifndef MYMODELCURVE2_H
#define MYMODELCURVE2_H


#include "../gmlib-master/modules/parametrics/gmpcurve.h"

using namespace GMlib;


// Model curve 2
// Hypotrochoid
template <typename T>
class MyModelCurve2 : public PCurve<T,3> {
    GM_SCENEOBJECT(MyModelCurve2)

public:
    MyModelCurve2(int R, int r, T d);
    MyModelCurve2( const MyModelCurve2<T>& copy );
    virtual ~MyModelCurve2();

    bool                isClosed() const override;

protected:
    // Virtual functions from PCurve, which have to be implemented locally
    void                eval(T t, int d, bool l) const override;
    T                   getStartP() const override;
    T                   getEndP()   const override;

private:
    int                 lcm(int a, int b) const;
    T                   _R;
    T                   _r;
    T                   _d;

}; // END class MyModelCurve2



// Constructor
template <typename T>
inline
    MyModelCurve2<T>::MyModelCurve2(int R, int r, T d) : PCurve<T,3>(20, 0, 0) {

    _R = R;
    _r = r;
    _d = d;
}



// Copy constructor
template <typename T>
inline
    MyModelCurve2<T>::MyModelCurve2( const MyModelCurve2<T>& copy ) : PCurve<T,3>(copy) {

    _R = copy._R;
    _r = copy._r;
    _d = copy._d;
}



// Destructor
template <typename T>
MyModelCurve2<T>::~MyModelCurve2() {}



template <typename T>
bool MyModelCurve2<T>::isClosed() const {
    return true;
}



template <typename T>
void MyModelCurve2<T>::eval( T t, int d, bool /*l*/ ) const {

    this->_p.setDim( d + 1 );

    this->_p[0][0] = (_R-_r) * cos(t) + _d * cos(((_R-_r)/_r) * t);     // x
    this->_p[0][1] = (_R-_r) * sin(t) - _d * sin(((_R-_r)/_r) * t);     // y
    this->_p[0][2] = T(0);                                              // z
}



template <typename T>
T MyModelCurve2<T>::getStartP() const {
    return T(0);
}



template <typename T>
T MyModelCurve2<T>::getEndP() const {

    int num = lcm(int(_r), int(_R));

    return T(M_2PI * T(num) / _R);
}



template <typename T>
int MyModelCurve2<T>::lcm(int a, int b) const {

    /*
     * Least common multiple function from GeeksforGeeks:
     * https://www.geeksforgeeks.org/program-to-find-lcm-of-two-numbers/
    */
    int greater = std::max(a, b);
    int smallest = std::min(a, b);
    for (int i = greater; ; i += greater) {
        if (i % smallest  == 0)
            return i;
    }
}



#endif // MYMODELCURVE2_H
