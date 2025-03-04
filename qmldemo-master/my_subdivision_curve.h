
#ifndef MYSUBDIVISIONCURVE_H
#define MYSUBDIVISIONCURVE_H


#include "../gmlib-master/modules/parametrics/gmpcurve.h"
#include "../gmlib-master/modules/parametrics/curves/gmpcircle.h"

using namespace GMlib;

/*
 * Closed Subdivision Curve
 * Lane Riesenfeld
*/
template <typename T>
class MySubdivisionCurve : public PCurve<T,3> {
    GM_SCENEOBJECT(MySubdivisionCurve)

public:
    MySubdivisionCurve(const DVector<Vector<T,3>>& P, bool isClosed=true);
    MySubdivisionCurve(const MySubdivisionCurve<T>& copy );
    virtual ~MySubdivisionCurve();

    void                    sample(int m, int d = 0) override;
    bool                    isClosed() const override;

protected:
    // Virtual functions from PCurve, which have to be implemented locally
    void                    eval( T t, int d, bool /*l*/ ) const override;
    T                       getStartP() const override;
    T                       getEndP()   const override;

private:
    void                    LaneRiesenfeldClosed(std::vector<DVector<Vector<T,3>>>& ph, int k, int d) const;
    int                     doublePart(std::vector<DVector<Vector<T,3>>>& ph, int n) const;
    void                    smoothPartClosed(std::vector<DVector<Vector<T,3>>>& ph, int n, int d) const;

    DVector<Vector<T,3>>    _P;                                     // points
    bool                    _isClosed;

}; // END class MySubdivisionCurve



// Constructor
template <typename T>
inline
    MySubdivisionCurve<T>::MySubdivisionCurve(const DVector<Vector<T,3>>& P, bool isClosed)
    : PCurve<T,3>(20, 0, 0), _P(P), _isClosed(isClosed) {}



// Copy constructor
template <typename T>
inline
    MySubdivisionCurve<T>::MySubdivisionCurve(const MySubdivisionCurve<T>& copy ) : PCurve<T,3>(copy) {

    _P = copy._P;
    _isClosed = copy._isClosed;
}



// Destructor
template <typename T>
MySubdivisionCurve<T>::~MySubdivisionCurve() {}



template <typename T>
bool MySubdivisionCurve<T>::isClosed() const {
    return _isClosed;
}



template <typename T>
void MySubdivisionCurve<T>::eval( T t, int d, bool /*l*/ ) const {
    this->_p.setDim( d + 1 );
}



template <typename T>
T MySubdivisionCurve<T>::getStartP() const {
    return 0;
}



template <typename T>
T MySubdivisionCurve<T>::getEndP() const {
    return 0;
}



template <typename T>
void MySubdivisionCurve<T>::LaneRiesenfeldClosed(std::vector<DVector<Vector<T,3>>>& ph, int k, int d) const {

    // k - number of doublings of points
    // d - degree of smoothing

    int n = _P.getDim();
    int m = pow(2, k) * n + 1;

    ph.resize(m);

    for (int i=0; i<n; i++) {
        ph[i][0] = _P[i];
    }

    ph[n++][0] = _P[0];

    for (int i=0; i<k; i++) {
        n = doublePart(ph, n);
        smoothPartClosed(ph, n, d);
    }


}



template <typename T>
int MySubdivisionCurve<T>::doublePart(std::vector<DVector<Vector<T,3>>>& ph, int n) const {

    for (int i=n-1; i>0; i--) {
        ph[2*i][0] = ph[i][0];
        ph[2*i-1][0] = 0.5 * (ph[i][0] + ph[i-1][0]);
    }
    return 2 * n - 1;
}



template <typename T>
void MySubdivisionCurve<T>::smoothPartClosed(std::vector<DVector<Vector<T,3>>>& ph, int n, int d) const {

    for (int j=1; j < d; j++) {
        for (int i=0; i < n-1; i++) {
            ph[i][0] = 0.5 * (ph[i][0] + ph[i+1][0]);
        }
        ph[n-1][0] = ph[0][0];
    }
}



template <typename T>
void MySubdivisionCurve<T>::sample( int m, int d ) {

    _visu.resize(1);

    LaneRiesenfeldClosed(_visu[0].sample_val, m, d);

    computeSurroundingSphere(_visu[0].sample_val, _visu[0].sur_sphere);

    // Visualize points
    int b = pow(2, m) * _P.getDim() + 1;
    for (int i=0; i<b; i++) {
        auto point = new PCircle<T>(T(0));

        point->toggleDefaultVisualizer();
        point->sample(1,0);
        point->setCollapsed(true);
        point->translate(_visu[0].sample_val[i][0]);

        this->insert(point);
    }

    this->setEditDone();
}



#endif // MYSUBDIVISIONCURVE_H
