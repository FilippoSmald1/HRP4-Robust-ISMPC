#ifndef DART_EXAMPLE_OSG_OSGATLASSIMBICON_ATLASSIMBICONEVENTHANDLER_HPP_
#define DART_EXAMPLE_OSG_OSGATLASSIMBICON_ATLASSIMBICONEVENTHANDLER_HPP_

#include <dart/dart.hpp>
#include <dart/utils/utils.hpp>
#include <dart/gui/osg/osg.hpp>

#include "HRP4WorldNode.hpp"

class HRP4EventHandler : public osgGA::GUIEventHandler
{
public:

  HRP4EventHandler(HRP4WorldNode* node);

  bool handle(const osgGA::GUIEventAdapter& ea,
              osgGA::GUIActionAdapter&) override;

protected:

  HRP4WorldNode* mNode;

};

#endif // DART_EXAMPLE_OSG_OSGATLASSIMBICON_ATLASSIMBICONEVENTHANDLER_HPP_
