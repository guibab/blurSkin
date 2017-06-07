#include <maya/MFnPlugin.h>

#include "blurSkinCmd.h"

MStatus initializePlugin(MObject obj) {
    MStatus result;
    MFnPlugin plugin(obj, "blur studios", "1.0", "Any");

    result = plugin.registerCommand("blurSkinCmd", blurSkinCmd::creator, blurSkinCmd::newSyntax);
    CHECK_MSTATUS_AND_RETURN_IT(result);

    return result;
}

MStatus uninitializePlugin(MObject obj) {
    MStatus status;
    MFnPlugin plugin(obj);

    status = plugin.deregisterCommand("blurSkinCmd");
    CHECK_MSTATUS_AND_RETURN_IT(status);
    return status;
}