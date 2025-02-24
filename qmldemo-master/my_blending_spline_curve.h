
#ifndef MYBLENDINGSPLINECURVE_H
#define MYBLENDINGSPLINECURVE_H


#include "../gmlib-master/modules/parametrics/gmpcurve.h"
#include "../gmlib-master/modules/parametrics/curves/gmpsubcurve.h"


using namespace GMlib;


// Blending Spline Curve
// Hardcoded degree 1
template <typename T>
class MyBlendingSplineCurve : public PCurve<T,3> {
    GM_SCENEOBJECT(MyBlendingSplineCurve)

public:
    MyBlendingSplineCurve(PCurve<T,3>* cu, int n);
    MyBlendingSplineCurve(const MyBlendingSplineCurve<T>& copy);
    virtual ~MyBlendingSplineCurve();

    void                    sample(int m, int d = 0) override;
    bool                    isClosed() const override;
    void                    toggleAnimate();

protected:
    void                    eval( T t, int d, bool /*l*/ ) const override;
    T                       getStartP() const override;
    T                       getEndP()   const override;
    void                    localSimulate(double dt) override;
    void                    animate(double dt);

private:
    std::vector<T>          _t;                                         // knot vector
    DVector<PCurve<T,3>*>   _c;                                         // subcurves
    int                     _n;                                         // Number of control points
    int                     _d;                                         // Degree
    PCurve<T,3>*            _cu;                                        // Model curve
    int                     _sample;

    T                       W(int d, int i, T t) const;                 // Map domain to 0-1
    Vector<T,3>             B(int i, T t) const;                        // Basis functions for given t
    T                       Bl(T x) const;                              // Blending(B)-function
    int                     find_i(T t) const;                          // Find index i for given t
    void                    make_knot_vector(int n, int d, T s, T e, bool closed);
    void                    make_local_curves(int n, bool closed);

    // For animation
    bool                    _animate;
    float                   _rotation_val;
    int                     _rotation_dir;
    float                   _color_val;
    int                     _color_dir;
    float                   _trans_val;
    int                     _trans_dir;
    float                   _scale;
    float                   _scale_val;
    bool                    _scale_increase;

}; // END class MyBlendingSplineCurve



// Constructor 1
template <typename T>
inline
    MyBlendingSplineCurve<T>::MyBlendingSplineCurve(PCurve<T,3>* cu, int n) : PCurve<T,3>(20, 0, 0),
    _cu(cu), _n(n) {

    _d = 1;                       // Degree of blending spline curve hardcoded as 1

    _animate = false;
    _rotation_val = 0.0;
    _rotation_dir = 1;
    _color_val = 0.0;
    _color_dir = 1;
    _trans_val = 0.0;
    _trans_dir = 1;
    _scale = 1.0;
    _scale_val = 1.002;
    _scale_increase = true;

    auto s = _cu->getParStart();
    auto e = _cu->getParEnd();

    make_knot_vector(_n, _d, s, e, isClosed());
    make_local_curves(_n, isClosed());
}



// Copy constructor
template <typename T>
inline
    MyBlendingSplineCurve<T>::MyBlendingSplineCurve(const MyBlendingSplineCurve<T>& copy) : PCurve<T,3>(copy) {

    _cu = copy._cu;
    _n = copy._n;
    _t = copy._t;
}



// Destructor
template <typename T>
MyBlendingSplineCurve<T>::~MyBlendingSplineCurve() {}



template <typename T>
void MyBlendingSplineCurve<T>::sample(int m, int d) {

    _sample = m;

    PCurve<T,3>::sample(m, d);
}



template <typename T>
bool MyBlendingSplineCurve<T>::isClosed() const {
    return _cu->isClosed();
}



template <typename T>
void MyBlendingSplineCurve<T>::toggleAnimate() {
    _animate = not _animate;
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
void MyBlendingSplineCurve<T>::eval( T t, int d, bool /*l*/ ) const {

    this->_p.setDim( d + 1 );

    int i = find_i(t);
    Vector<T,3> b = B(i, t);

    auto c0 = _c[i-1]->evaluateParent(t, 0);
    auto c1 = _c[i]->evaluateParent(t, 0);

    this->_p[0] = b[0]*c0[0] + b[1]*c1[0];
}



template <typename T>
T MyBlendingSplineCurve<T>::getStartP() const {
    return _cu->getParStart();
}



template <typename T>
T MyBlendingSplineCurve<T>::getEndP() const {
    return _cu->getParEnd();
}



template <typename T>
void MyBlendingSplineCurve<T>::localSimulate(double dt) {

    if (_animate) animate(dt);


    this->sample(_sample, 0);
    this->setEditDone(true);

    SceneObject::localSimulate(dt);
}



template <typename T>
void MyBlendingSplineCurve<T>::animate(double dt) {

    // Rotate segments
    float lim_rot = 0.7;
    if (_rotation_val > lim_rot)    _rotation_dir = -1;
    if (_rotation_val < -1*lim_rot) _rotation_dir = 1;

    _rotation_val += _rotation_dir * dt;

    // Translate curve
    float lim_trans = 2.6;
    float trans_speed = 0.4;
    if (_trans_val > lim_trans)     _trans_dir = -1;
    if (_trans_val < -1*lim_trans)  _trans_dir = 1;

    _trans_val += _trans_dir * dt;

    // Color
    float col_scale = 100.0;
    if (_color_val >= 255.0) {
        _color_dir = -1;
        _color_val = 255.0;
    }
    if (_color_val <= 150.0) {
        _color_dir = 1;
        _color_val = 150.0;
    }

    _color_val += _color_dir * dt * col_scale;

    // Scale curve
    float scale_min = 0.8;
    float scale_max = 1.3;
    if (_scale > scale_max) _scale_increase = false;
    if (_scale < scale_min) _scale_increase = true;

    if (_scale_increase)    _scale = _scale*_scale_val;
    else                    _scale = _scale*(1.0/_scale_val);


    for (int i = 0; i<_c.getDim()-1; i++) {
        // Position values used as s vector that points away from center of parrent.
        auto pos = _c[i]->getPos();

        if (i%2 == 0) {
            _c[i]->rotate(dt * _rotation_dir, Vector<float,3>(pos[0], pos[1], pos[2])*-1);
        }
        else {
            _c[i]->rotate(dt * _rotation_dir, Vector<float,3>(pos[0], pos[1], pos[2]));
        }
    }

    this->translate(Vector<float,3>(0, 0, 1)*dt*_trans_dir*trans_speed);
    this->setColor(GMlib::Color(int(_color_val), 150, 0));
    if (_scale_increase) this->scale(_scale_val);
    else this->scale(1 / _scale_val);
}



template <typename T>
T MyBlendingSplineCurve<T>::W(int d, int i, T t) const {
    return (t - _t[i])/(_t[i+d] - _t[i]);
}



template <typename T>
Vector<T,3> MyBlendingSplineCurve<T>::B(int i, T t) const {

    Vector<T,3> b;

    b[0] = 1 - Bl(W(1, i, t));
    b[1] = Bl(W(1, i, t));

    return b;
}



template <typename T>
T MyBlendingSplineCurve<T>::Bl(T x) const {
    return 0.5 - 0.5 * cos(M_PI * x);                   // same as pow(sin(M_Pi * x / 2), 2);
}



template <typename T>
int MyBlendingSplineCurve<T>::find_i(T t) const {

    int i = _d;
    for (; i<_n; i++) {
        if (_t[i+1] > t) {
            return i;
        }
    }

    return _n-1;
}



template <typename T>
void MyBlendingSplineCurve<T>::make_knot_vector(int n, int d, T s, T e, bool closed) {

    if (closed) {
        for (int i = -d; i<n+d; i++) _t.push_back(s + (e-s)/(n-d)*(i));
    }
    else {
        for (int i = 0; i<=d; i++) _t.push_back(s);

        for (int i = d+1; i<n; i++) _t.push_back(s + (e-s)/(n-d)*(i-d));

        for (int i = 0; i<=d; i++) _t.push_back(e);
    }
}



template <typename T>
void MyBlendingSplineCurve<T>::make_local_curves(int n, bool closed) {

    if (closed) n -= 1;

    for (int i = 0; i<n; i++){
        _c.append(new PSubCurve<T>(_cu, _t[i], _t[i+2], _t[i+1]));

        // Visualize subcurves:
        _c.back()->toggleDefaultVisualizer();
        _c.back()->sample(10,0);
        _c.back()->setColor(GMlib::Color(200, 50, 0, 150));
        _c.back()->setCollapsed(true);

        this->insert(_c.back());
    }

    if (closed) _c.append(_c[0]);
}



#endif // MYBLENDINGSPLINECURVE_H
