#ifndef MY_BLENDING_SPLINE_SURFACE_H
#define MY_BLENDING_SPLINE_SURFACE_H


#include "../gmlib-master/modules/parametrics/gmpsurf.h"
#include "simplesubsurf.h"

using namespace GMlib;


// Blending Spline Surface
template <typename T>
class MyBlendingSplineSurface : public PSurf<T,3> {
    GM_SCENEOBJECT(MyBlendingSplineSurface)

public:
    MyBlendingSplineSurface(PSurf<T, 3>* su, int nu, int nv);
    MyBlendingSplineSurface(const MyBlendingSplineSurface<T>& copy);
    virtual ~MyBlendingSplineSurface();

    void                    sample(int m1, int m2, int d1, int d2);
    bool                    isClosedU() const override;
    bool                    isClosedV() const override;

protected:
    void                    eval(T u, T v, int d1, int d2, bool lu = true, bool lv = true ) const override;
    T                       getStartPU() const override;
    T                       getEndPU()   const override;
    T                       getStartPV() const override;
    T                       getEndPV()   const override;
    void                    localSimulate(double dt) override;

private:
    int                     _m1;
    int                     _m2;

    PSurf<T, 3>*            _su;
    DMatrix<PSurf<T,3>*>    _c;
    int                     _nu;
    int                     _nv;
    std::vector<T>          _u;         // Knot vector u
    std::vector<T>          _v;         // Knot vector v

    T                       W(const std::vector<T>& tt, int d, int i, T t) const;
    Vector<T,3>             B(const std::vector<T>& tt, int i, T t) const;
    T                       Bl(T x) const;
    T                       dW(const std::vector<T>& tt, int d, int i) const;       // Derivative of W
    Vector<T,3>             dB(const std::vector<T>& tt, int i, T t) const;         // Derivative of B
    T                       dBl(T x) const;                                         // Derivative of Bl

    int                     find_i(const std::vector<T>& tt, int d, int n, T t) const;                  // Find index i for given tt
    void                    make_knot_vector(std::vector<T>& tt, int n, int d, T s, T e, bool closed);  // Populate knot vector tt
    void                    make_local_surfaces(int nu, int nv, bool closedU, bool closedV);            // Create sub surfaces


};



// Constructor 1
template <typename T>
inline
    MyBlendingSplineSurface<T>::MyBlendingSplineSurface(PSurf<T, 3>* su, int nu, int nv) {

    _su = su;
    _nu = nu;
    _nv = nv;

    make_knot_vector(_u, _nu, 1, _su->getParStartU(), _su->getParEndU(), isClosedU());  // u
    make_knot_vector(_v, _nv, 1, _su->getParStartV(), _su->getParEndV(), isClosedV());  // v

    make_local_surfaces(_nu, _nv, isClosedU(), isClosedV());


}



// Copy constructor
template <typename T>
inline
    MyBlendingSplineSurface<T>::MyBlendingSplineSurface(const MyBlendingSplineSurface<T>& copy) {

    _su = copy._su;
    _c = copy._c;
    _nu = copy._nu;
    _nv = copy._nv;
    _u = copy._u;
    _v = copy._v;

}



// Destructor
template <typename T>
MyBlendingSplineSurface<T>::~MyBlendingSplineSurface() {}



template <typename T>
void MyBlendingSplineSurface<T>::sample(int m1, int m2, int d1, int d2) {
    _m1 = m1;
    _m2 = m2;

    PSurf<T, 3>::sample(m1, m2, d1, d2);
}



template <typename T>
bool MyBlendingSplineSurface<T>::isClosedU() const {
    return _su->isClosedU();
}


template <typename T>
bool MyBlendingSplineSurface<T>::isClosedV() const {
    return _su->isClosedV();
}



template <typename T>
void MyBlendingSplineSurface<T>::eval( T u, T v, int d1, int d2, bool /*lu*/, bool /*lv*/ ) const {

    this->_p.setDim( d1+1, d2+1 );

    int i = find_i(_u, _nu, 1, u);
    int j = find_i(_v, _nv, 1, v);

    DMatrix<Vector<T,3>> a00 = _c(i-1)(j-1)->evaluateParent(u, v, d1, d2);  // S i-1, j-1
    DMatrix<Vector<T,3>> a01 = _c(i-1)(j)->evaluateParent(u, v, d1, d2);    // S i-1, j
    DMatrix<Vector<T,3>> a10 = _c(i)(j-1)->evaluateParent(u, v, d1, d2);    // S i,   j-1
    DMatrix<Vector<T,3>> a11 = _c(i)(j)->evaluateParent(u, v, d1, d2);      // S i,   j

    Vector<T,3> bu = B(_u, i, u);
    Vector<T,3> bv = B(_v, j, v);
    Vector<T,3> dbu = dB(_u, i, u);
    Vector<T,3> dbv = dB(_v, j, v);

    auto c0 = a00[0][0] + bu[1] * (a10[0][0] - a00[0][0]);
    auto c1 = a01[0][0] + bu[1] * (a11[0][0] - a01[0][0]);
    this->_p[0][0] = c0 + bv[1] *  (c1 - c0);

    // Du: [1][0]
    // Dv: [0][1]
    auto du0 = a00[1][0] + dbu[1] * (a10[0][0] - a00[0][0]) + bu[1] * (a10[1][0] - a00[0][0]);
    auto dv0 = a00[0][1] +  bu[1] * (a10[0][1] - a00[0][1]);

    auto du1 = a01[1][0] + dbu[1] * (a11[0][0] - a01[0][0]) + bu[1] * (a11[1][0] - a01[0][0]);
    auto dv1 = a01[0][1] +  bu[1] * (a11[0][1] - a01[0][1]);

    this->_p[1][0] = du0 +  bv[1] * (du1 - du0);
    this->_p[0][1] = dv0 + dbv[1] * (c1 - c0) + bv[1] * (dv1 - dv0);
}



template <typename T>
T MyBlendingSplineSurface<T>::getStartPU() const {
    return _su->getParStartU();
}



template <typename T>
T MyBlendingSplineSurface<T>::getEndPU() const {
    return _su->getParEndU();
}



template <typename T>
T MyBlendingSplineSurface<T>::getStartPV() const {
    return _su->getParStartV();
}



template <typename T>
T MyBlendingSplineSurface<T>::getEndPV() const {
    return _su->getParEndV();
}



template <typename T>
void MyBlendingSplineSurface<T>::localSimulate(double dt) {
    this->sample(_m1, _m2, 1, 1);
    this->setEditDone(true);

    SceneObject::localSimulate(dt);
}



template <typename T>
T MyBlendingSplineSurface<T>::W(const std::vector<T>& tt, int d, int i, T t) const {
    return (t - tt[i])/(tt[i+d] - tt[i]);
}



template <typename T>
Vector<T,3> MyBlendingSplineSurface<T>::B(const std::vector<T>& tt, int i, T t) const {

    Vector<T,3> b;

    b[0] = 1 - Bl(W(tt, 1, i, t));
    b[1] = Bl(W(tt, 1, i, t));

    return b;
}



template <typename T>
T MyBlendingSplineSurface<T>::Bl(T x) const {
    return 0.5 - 0.5 * cos(M_PI * x);                   // same as pow(sin(M_Pi * x / 2), 2);
}



template <typename T>
T MyBlendingSplineSurface<T>::dW(const std::vector<T>& tt, int d, int i) const {
    return 1/(tt[i+d] - tt[i]);
}



template <typename T>
Vector<T,3> MyBlendingSplineSurface<T>::dB(const std::vector<T>& tt, int i, T t) const {

    Vector<T,3> b;

    b[0] = -dBl(W(tt, 1, i, t)) * dW(tt, 1, i);
    b[1] = -b[0];

    return b;
}



template <typename T>
T MyBlendingSplineSurface<T>::dBl(T x) const {
    return 0.5 * sin(M_PI * x) * M_PI;
}



template <typename T>
int MyBlendingSplineSurface<T>::find_i(const std::vector<T>& tt, int n ,int d, T t) const {

    int i = d;
    for (; i<n; i++) {
        if (tt[i+1] > t) {
            return i;
        }
    }

    return n-1;
}



template <typename T>
void MyBlendingSplineSurface<T>::make_knot_vector(std::vector<T>& tt, int n, int d, T s, T e, bool closed) {

    if (closed) {
        for (int i = -d; i<n+d; i++) tt.push_back(s + (e-s)/(n-d)*(i));
    }
    else{
        for (int i = 0; i<=d; i++) tt.push_back(s);

        for (int i = d+1; i<n; i++) tt.push_back(s + (e-s)/(n-d)*(i-d));

        for (int i = 0; i<=d; i++) tt.push_back(e);
    }
}



template <typename T>
void MyBlendingSplineSurface<T>::make_local_surfaces(int nu, int nv, bool closedU, bool closedV) {

    _c = DMatrix<PSurf<T,3>*>(nu, nv);

    if (closedU) nu -= 1;
    if (closedV) nv -= 1;

    for (int i = 0; i<nu; i++) {
        for (int j = 0; j<nv; j++){
            auto sub_surf = new PSimpleSubSurf<T>(_su,
                                                  _u[i], _u[i+2], _u[i+1],
                                                  _v[j], _v[j+2], _v[j+1]);
            _c[i][j] = sub_surf;

            // Visualize subcurves:
            sub_surf->toggleDefaultVisualizer();
            sub_surf->sample(10, 10, 1, 1);
            sub_surf->setColor(GMlib::Color(200, 50, 0));
            sub_surf->setCollapsed(true);

            this->insert(sub_surf);
        }
    }
    if (closedU){
        _c[nu] = _c[0];
    }
    if (closedV){
        for (int i = 0; i<nu; i++) {
            _c[i][nv] = _c[i][0];
        }
    }
    if (closedU and closedV){
        _c[nu][nv] = _c[0][0];
    }
}



#endif // MY_BLENDING_SPLINE_SURFACE_H
