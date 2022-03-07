#include "HRP4WorldNode.hpp"

//==============================================================================
HRP4WorldNode::HRP4WorldNode(
    const dart::simulation::WorldPtr world,
    const dart::dynamics::SkeletonPtr hrp4)
  : dart::gui::osg::WorldNode(world),
    mExternalForce(Eigen::Vector3d::Zero()),
    mForceDuration(0.0)
{
  assert(world);
  assert(hrp4);

  mWorld = world;
  mRobot = hrp4;

  mController.reset(new Controller(hrp4, world));
  mController->setInitialConfiguration();


    // Create red ball
    /*Eigen::Isometry3d tf =mRobot->getBodyNode("base_link")->getWorldTransform();
    mTarget = std::make_shared<dart::dynamics::SimpleFrame>(dart::dynamics::Frame::World(), "target", tf);
    dart::dynamics::ShapePtr ball(new dart::dynamics::SphereShape(0.025));
    mTarget->setShape(ball);
    mTarget->getVisualAspect(true)->setColor(Eigen::Vector3d(0.9,0,0));
    mWorld->addSimpleFrame(mTarget);*/
}

//==============================================================================
void HRP4WorldNode::customPreStep()
{
  mController->update();
}

//==============================================================================
void HRP4WorldNode::reset()
{
  //mExternalForce.setZero();
  mController.reset(new Controller(mRobot, mWorld));
  mController->setInitialConfiguration();

}

//==============================================================================
void HRP4WorldNode::pushForwardAtlas(double force, int frames)
{/*
  mExternalForce.x() = force;
  mForceDuration = frames;*/
}

//==============================================================================
void HRP4WorldNode::pushBackwardAtlas(double force, int frames)
{/*
  mExternalForce.x() = -force;
  mForceDuration = frames;*/
}

//==============================================================================
void HRP4WorldNode::pushLeftAtlas(double force, int frames)
{/*
  mExternalForce.z() = force;
  mForceDuration = frames;*/
}

//==============================================================================
void HRP4WorldNode::pushRightAtlas(double force, int frames)
{/*
  mExternalForce.z() = -force;
  mForceDuration = frames;*/
}

//==============================================================================
void HRP4WorldNode::switchToNormalStrideWalking()
{/*
  mController->changeStateMachine("walking", mWorld->getTime());*/
}

//==============================================================================
void HRP4WorldNode::switchToShortStrideWalking()
{/*
  mController->changeStateMachine("running", mWorld->getTime());*/
}

//==============================================================================
void HRP4WorldNode::switchToNoControl()
{/*
  mController->changeStateMachine("standing", mWorld->getTime());*/
}

dart::dynamics::SkeletonPtr HRP4WorldNode::getRobot() {
	return mRobot;
}

dart::simulation::WorldPtr HRP4WorldNode::getWorld() {
	return mWorld;
}

std::shared_ptr<Controller> HRP4WorldNode::getController() {
	return mController;
}

