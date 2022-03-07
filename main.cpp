#include <dart/dart.hpp>
#include <dart/gui/osg/osg.hpp>
#include <dart/utils/utils.hpp>
#include <dart/utils/urdf/urdf.hpp>
#include <dart/gui/gui.hpp>
#include "HRP4WorldNode.hpp"
#include "HRP4EventHandler.hpp"

int main(int argc, char* argv[])
{

  // Create a world
  dart::simulation::WorldPtr world(new dart::simulation::World);

  // Load ground and HRP4 robot and add them to the world
  dart::utils::DartLoader urdfLoader;
  auto ground = urdfLoader.parseSkeleton(realpath("../ground.urdf", NULL));
  auto hrp4 = urdfLoader.parseSkeleton(realpath("../urdf/hrp4.urdf", NULL));
  Eigen::Isometry3d tf = Eigen::Isometry3d::Identity();
  tf.translation() += Eigen::Vector3d(3.0, 0.0, 0.0);
  
  ground->getJoint(0)->setTransformFromParentBodyNode(tf);  
  world->addSkeleton(ground);
  world->addSkeleton(hrp4);

  // set joint actuator type and force limits
  double forceLimit = 100;

  for (int i = 0; i < hrp4->getNumJoints(); i++) {
	  size_t  dim   = hrp4->getJoint(i)->getNumDofs();
	  if(dim==6) {
		  hrp4->getJoint(i)->setActuatorType(dart::dynamics::Joint::PASSIVE);
	  }
	  if(dim==1) {
		  hrp4->getJoint(i)->setActuatorType(dart::dynamics::Joint::ACCELERATION);
		  hrp4->getJoint(i)->setForceUpperLimit(0,  forceLimit);
		  hrp4->getJoint(i)->setForceLowerLimit(0, -forceLimit);
		  hrp4->getJoint(i)->setPositionLimitEnforced(true);
	  }
  }

  // Set gravity of the world
  world->setGravity(Eigen::Vector3d(0.0, 0.0, -9.81));
  world->setTimeStep(1.0/100.0);

  // Wrap a WorldNode around it
  osg::ref_ptr<HRP4WorldNode> node = new HRP4WorldNode(world, hrp4);
  node->setNumStepsPerCycle(1);

  // Create a Viewer and set it up with the WorldNode
  dart::gui::osg::ImGuiViewer viewer;
  viewer.addWorldNode(node);
  viewer.switchHeadlights(false);

  // Set recording
  //viewer.record("../data","frame"); 
  // use this to generate the video from the frames
  //ffmpeg -framerate 100 -i frame%06d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p sim_video.mp4  
  
  // Pass in the custom event handler
  viewer.addEventHandler(new HRP4EventHandler(node));

  // Set the dimensions for the window
  viewer.setUpViewInWindow(0, 0, 1280, 960);

  // Set the window name
  viewer.realize();
  osgViewer::Viewer::Windows windows;
  viewer.getWindows(windows);
  windows.front()->setWindowName("HRP4 with IS-MPC");

  // Adjust the viewpoint of the Viewer
  // read .txt files to store trajectory if needed - first store the file size
  std::ifstream inFile; inFile.open("../view.txt");
  if (!inFile) {
       std::cout << "Unable to open File" << std::endl;
       exit(1);
  }
  int file_size = 0; double buffer;
  while (inFile >> buffer) {
      file_size = file_size + 1;
  }
  inFile.close();
  Eigen::VectorXd viewpoint = Eigen::VectorXd::Zero(10); 
  inFile.open("../view.txt"); 
  if (!inFile) {
       std::cout << "Unable to open File" << std::endl;
       exit(1);
  }
  for (int i = 0; i < 10; i++) {
       inFile >> viewpoint(i);
  }
  inFile.close(); 

  viewer.getCameraManipulator()->setHomePosition(
        ::osg::Vec3d( viewpoint(0),  viewpoint(1), viewpoint(2))*viewpoint(3),
        ::osg::Vec3d( viewpoint(4),  viewpoint(5), viewpoint(6)),
        ::osg::Vec3d( viewpoint(7),  viewpoint(8), viewpoint(9)));  

  viewer.setCameraManipulator(viewer.getCameraManipulator());

  // Run the simulation!
  viewer.run();
}
