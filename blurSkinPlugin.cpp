#include <maya/MFnPlugin.h>

#include "blurSkinCmd.h"
#include "blurSkinEdit.h"
#include "pointsDisplay.h"

MStatus initializePlugin(MObject obj) {
    MStatus status;
    MFnPlugin plugin(obj, "blur studios", "1.0", "Any");

    status = plugin.registerCommand("blurSkinCmd", blurSkinCmd::creator, blurSkinCmd::newSyntax);
    CHECK_MSTATUS_AND_RETURN_IT(status);

    status = plugin.registerNode("blurSkinDisplay", blurSkinDisplay::id, blurSkinDisplay::creator,
                                 blurSkinDisplay::initialize);

    if (!status) {
        status.perror("registerNode");
        return (status);
    }

    status = plugin.registerNode("pointsDisplay", pointsDisplay::id, &pointsDisplay::creator,
                                 &pointsDisplay::initialize, MPxNode::kLocatorNode,
                                 &pointsDisplay::drawDbClassification);
    // sUseLegacyDraw ? NULL : &pointsDisplay::drawDbClassification);
    if (!status) {
        status.perror("registerNode");
        return status;
    }

    status = MHWRender::MDrawRegistry::registerDrawOverrideCreator(
        pointsDisplay::drawDbClassification, pointsDisplay::drawRegistrantId,
        PointsDisplayDrawOverride::Creator);
    if (!status) {
        status.perror("registerDrawOverrideCreator");
        return status;
    }

    return (status);
}

MStatus uninitializePlugin(MObject obj) {
    MStatus status;
    MFnPlugin plugin(obj);

    status = plugin.deregisterCommand("blurSkinCmd");
    CHECK_MSTATUS_AND_RETURN_IT(status);

    status = plugin.deregisterNode(blurSkinDisplay::id);
    if (!status) {
        status.perror("deregisterNode");
        return (status);
    }

    status = MHWRender::MDrawRegistry::deregisterDrawOverrideCreator(
        pointsDisplay::drawDbClassification, pointsDisplay::drawRegistrantId);
    if (!status) {
        status.perror("deregisterDrawOverrideCreator");
        return status;
    }

    status = plugin.deregisterNode(pointsDisplay::id);
    if (!status) {
        status.perror("deregisterNode");
        return status;
    }

    return status;
}