
#include <iostream>

#include "scenario.h"
#include "testtorus.h"
#include "my_model_curve1.h"
#include "my_model_curve2.h"
#include "my_bspline.h"
#include "my_subdivision_curve.h"
#include "my_blending_spline_curve.h"
#include "my_blending_spline_surface.h"
#include "../gmlib-master/modules/parametrics/surfaces/gmpplane.h"
#include "../gmlib-master/modules/parametrics/visualizers/gmpsurfnormalsvisualizer.h"


// hidmanager
#include "hidmanager/defaulthidmanager.h"

// gmlib
#include <scene/light/gmpointlight.h>
#include <scene/sceneobjects/gmpathtrack.h>
#include <scene/sceneobjects/gmpathtrackarrows.h>

// qt
#include <QQuickItem>


template <typename T>
inline
    std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
    out << v.size() << std::endl;
    for(uint i=0; i<v.size(); i++) out << " " << v[i];
    out << std::endl;
    return out;
}


void Scenario::initializeScenario() {
    // Insert a light
    GMlib::Point<GLfloat,3> init_light_pos( 2.0, 4.0, 10 );
    GMlib::PointLight *light = new GMlib::PointLight(  GMlib::GMcolor::white(), GMlib::GMcolor::white(),
                                                     GMlib::GMcolor::white(), init_light_pos );
    light->setAttenuation(0.8f, 0.002f, 0.0008f);
    this->scene()->insertLight( light, false );

    // Insert Sun
    this->scene()->insertSun();

    // Default camera parameters
    int init_viewport_size = 600;
    GMlib::Point<float,3>  init_cam_pos( 0.0f, 0.0f, 0.0f );
    GMlib::Vector<float,3> init_cam_dir( 0.0f, 1.0f, 0.0f );
    GMlib::Vector<float,3> init_cam_up(  1.0f, 0.0f, 0.0f );

    // Projection cam
    auto proj_rcpair = createRCPair("Projection");
    proj_rcpair.camera->set(init_cam_pos,init_cam_dir,init_cam_up);
    proj_rcpair.camera->setCuttingPlanes( 1.0f, 8000.0f );
    proj_rcpair.camera->rotateGlobal( GMlib::Angle(-45), GMlib::Vector<float,3>( 1.0f, 0.0f, 0.0f ) );
    proj_rcpair.camera->translateGlobal( GMlib::Vector<float,3>( 0.0f, -20.0f, 20.0f ) );
    scene()->insertCamera( proj_rcpair.camera.get() );
    proj_rcpair.renderer->reshape( GMlib::Vector<int,2>(init_viewport_size, init_viewport_size) );


    /***************************************************************************
   *                                                                         *
   * Standar example, including path track and path track arrows             *
   *                                                                         *
   ***************************************************************************/

    GMlib::Material mm(GMlib::GMmaterial::polishedBronze());
    mm.set(45.0);

    GMlib::Material mm2(GMlib::GMmaterial::jade());
    mm2.set(45.0);

    bool drawTorus = false;
    bool drawModelCurve1 = false;
    bool drawModelCurve2 = false;
    bool drawBSpline1 = false;
    bool drawBSpline2 = false;
    bool drawSubdivisionCurve = false;
    bool drawBlendingSplineCurve = false;
    bool drawBlendingSplineSurface = true;


    /* Test torus */
    if (drawTorus) {
        auto ptom = new TestTorus(1.0f, 0.4f, 0.6f);
        ptom->toggleDefaultVisualizer();
        ptom->sample(60,60,1,1);
        this->scene()->insert(ptom);
        auto ptrack = new GMlib::PathTrack();
        ptrack->setLineWidth(2);
        ptom->insert(ptrack);
        auto ptrack2 = new GMlib::PathTrackArrows();
        ptrack2->setArrowLength(2);
        ptom->insert(ptrack2);
    }

    /* Test model curve */
    if (drawModelCurve1) {
        auto tcurve = new MyModelCurve1<float>(3, 3, 1, 2);
        tcurve->toggleDefaultVisualizer();
        tcurve->sample(200,0);
        tcurve->setLineWidth(4);
        tcurve->translate(Vector<float,3>(0,0,0));
        this->scene()->insert(tcurve);

        auto tcurve2 = new MyModelCurve1<float>(3, 3, 5, 3);
        tcurve2->toggleDefaultVisualizer();
        tcurve2->sample(200,0);
        tcurve2->setLineWidth(4);
        tcurve2->translate(Vector<float,3>(-7,0,0));
        this->scene()->insert(tcurve2);
    }

    if (drawModelCurve2) {
        auto tcurve = new MyModelCurve2<float>(8, 3, 2);
        tcurve->toggleDefaultVisualizer();
        tcurve->sample(200,0);
        tcurve->setLineWidth(4);
        tcurve->translate(Vector<float,3>(0,-5,0));
        this->scene()->insert(tcurve);

        auto tcurve2 = new MyModelCurve2<float>(5, 2, 1);
        tcurve2->toggleDefaultVisualizer();
        tcurve2->sample(200,0);
        tcurve2->setLineWidth(4);
        tcurve2->translate(Vector<float,3>(0, 5, 0));
        this->scene()->insert(tcurve2);
    }

    /* Test b-spline, Constructor 1 */
    if (drawBSpline1) {
        GMlib::DVector<GMlib::Vector<float,3>> points(10);
        points[0] = GMlib::Vector<float,3>(0, 0, 0);
        points[1] = GMlib::Vector<float,3>(0, 2, 0);
        points[2] = GMlib::Vector<float,3>(2, 2, 0);
        points[3] = GMlib::Vector<float,3>(4, 1.5, 0);
        points[4] = GMlib::Vector<float,3>(4, 0, 0);
        points[5] = GMlib::Vector<float,3>(3, 0, 1);
        points[6] = GMlib::Vector<float,3>(3, 1, 0);
        points[7] = GMlib::Vector<float,3>(1, 1, 0);
        points[8] = GMlib::Vector<float,3>(1, 0.5, 0);
        points[9] = GMlib::Vector<float,3>(2, 0, 0);

        auto test_bspline1 = new MyBSpline<float>(points);

        test_bspline1->toggleDefaultVisualizer();
        test_bspline1->sample(100,0);
        test_bspline1->setLineWidth(4);
        test_bspline1->translate(Vector<float,3>(1,1,1));
        this->scene()->insert(test_bspline1);
    }

    /* Test b-spline, Constructor 2 */
    if (drawBSpline2) {
        int p = 100;
        float a = 4;
        float b = 6;

        GMlib::DVector<GMlib::Vector<float,3>> points(p);
        for (int i = 0; i < p; i++) {
            float t = i/float(p) * M_2PI;
            points[i] = GMlib::Vector<float,3>(a*sin(t), b*cos(t), 0);
        }

        auto test_bspline2 = new MyBSpline<float>(points, 6);

        test_bspline2->toggleDefaultVisualizer();
        test_bspline2->sample(100,0);
        test_bspline2->setLineWidth(4);
        this->scene()->insert(test_bspline2);
    }

    /* Test Subdivision curve */
    if (drawSubdivisionCurve) {
        GMlib::DVector<GMlib::Vector<float,3>> points(4);
        points[0] = GMlib::Vector<float,3>(0, 0, 0);
        points[1] = GMlib::Vector<float,3>(0, 5, 0);
        points[2] = GMlib::Vector<float,3>(7, 5, 0);
        points[3] = GMlib::Vector<float,3>(7, 0, 0);

        auto test_sub = new MySubdivisionCurve<float>(points);
        test_sub->toggleDefaultVisualizer();
        test_sub->sample(4,2);
        this->scene()->insert(test_sub);
    }

    /* Test blending spline curve */
    if (drawBlendingSplineCurve) {
        // auto model_curve = new MyModelCurve1<float>(5, 5, 3, 2);
        auto model_curve = new MyModelCurve2<float>(8, 3, 1.8);
        auto blending_spline = new MyBlendingSplineCurve<float>(model_curve, 9);

        blending_spline->toggleDefaultVisualizer();
        blending_spline->sample(200,0);
        blending_spline->toggleAnimate();

        blending_spline->translate(Vector<float,3>(-3,0,0));
        blending_spline->rotate(0.6, Vector<float,3>(0, 1, 0));
        blending_spline->scale(0.7);
        blending_spline->setLineWidth(5);
        blending_spline->setColor(GMlib::Color(255, 0, 0));

        this->scene()->insert(blending_spline);
    }

    /* Test blending spline surface */
    if (drawBlendingSplineSurface) {
        auto n_vis = new PSurfNormalsVisualizer<float, 3>();

        auto plane = new PPlane<float>(
            Point<float,3>(-4,-4,0),
            Vector<float,3>(8,0,0),
            Vector<float,3>(0,8,0));

        // auto plane2 = new PPlane<float>(
        //     Point<float,3>(-4,-4,0),
        //     Vector<float,3>(8,0,0),
        //     Vector<float,3>(0,8,0));
        // plane2->toggleDefaultVisualizer();
        // plane2->sample(20,20,1,1);
        // plane2->translate(Vector<float,3>(0,0,-3));
        // plane2->setColor(GMlib::Color(255, 255, 0));
        // plane2->setMaterial(mm2);
        // this->scene()->insert(plane2);

        auto blending_spline = new MyBlendingSplineSurface<float>(
            plane, 3, 3);
        blending_spline->toggleDefaultVisualizer();
        blending_spline->sample(10,10,1,1);
        blending_spline->setMaterial(mm2);
        blending_spline->insertVisualizer(n_vis);
        this->scene()->insert(blending_spline);
    }
}




void Scenario::cleanupScenario() {

}




void Scenario::callDefferedGL() {

    GMlib::Array< const GMlib::SceneObject*> e_obj;
    this->scene()->getEditedObjects(e_obj);

    for(int i=0; i < e_obj.getSize(); i++)
        if(e_obj(i)->isVisible()) e_obj[i]->replot();
}

