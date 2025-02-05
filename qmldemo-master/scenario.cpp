
#include <iostream>

#include "scenario.h"
#include "testtorus.h"
#include "application/my_model_curve1.h"
#include "application/my_bspline.h".h"


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

  /* Test torus */
  // auto ptom = new TestTorus(1.0f, 0.4f, 0.6f);
  // ptom->toggleDefaultVisualizer();
  // ptom->sample(60,60,1,1);
  // this->scene()->insert(ptom);
  // auto ptrack = new GMlib::PathTrack();
  // ptrack->setLineWidth(2);
  // ptom->insert(ptrack);
  // auto ptrack2 = new GMlib::PathTrackArrows();
  // ptrack2->setArrowLength(2);
  // ptom->insert(ptrack2);

  /* Test model curve */
  // auto tcurve = new MyModelCurve1<float>(5, 5, 5, 4);
  // tcurve->toggleDefaultVisualizer();
  // tcurve->sample(200,0);
  // tcurve->setLineWidth(4);
  // this->scene()->insert(tcurve);

  /* Test b-spline */
  GMlib::DVector<GMlib::Vector<float,3>> points(10);
  points[0] = GMlib::Vector<float,3>(0, 0, 0);
  points[1] = GMlib::Vector<float,3>(0, 2, 0);
  points[2] = GMlib::Vector<float,3>(2, 2, 0);
  points[3] = GMlib::Vector<float,3>(4, 1.5, 0);
  points[4] = GMlib::Vector<float,3>(4, 0, 0);
  points[5] = GMlib::Vector<float,3>(3, 0, 0);
  points[6] = GMlib::Vector<float,3>(3, 1, 0);
  points[7] = GMlib::Vector<float,3>(1, 1, 0);
  points[8] = GMlib::Vector<float,3>(1, 0.5, 0);
  points[9] = GMlib::Vector<float,3>(2, 0, 0);



  auto test_bspline = new MyBSpline<float>(points);
  test_bspline->toggleDefaultVisualizer();
  test_bspline->sample(100,0);
  test_bspline->setLineWidth(4);
  this->scene()->insert(test_bspline);

}




void Scenario::cleanupScenario() {

}




void Scenario::callDefferedGL() {

  GMlib::Array< const GMlib::SceneObject*> e_obj;
  this->scene()->getEditedObjects(e_obj);

  for(int i=0; i < e_obj.getSize(); i++)
    if(e_obj(i)->isVisible()) e_obj[i]->replot();
}

