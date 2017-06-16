#include "functions.h"

unsigned int getMIntArrayIndex(MIntArray& myArray, int searching) {
    unsigned int toReturn = -1;
    for (unsigned int element = 0; element < myArray.length(); ++element) {
        if (myArray[element] == searching) {
            toReturn = element;
            break;
        }
    }
    return toReturn;
}

void CVsAround(int storedU, int storedV, int numCVsInU, int numCVsInV, bool UIsPeriodic,
               bool VIsPeriodic, MIntArray& vertices) {
    int resCV;
    // plus U
    int UNext = storedU + 1;
    if (UNext < numCVsInU) {
        resCV = numCVsInV * UNext + storedV;
        if (getMIntArrayIndex(vertices, resCV) == -1) vertices.append(resCV);
    } else if (UIsPeriodic) {
        UNext -= numCVsInU;
        resCV = numCVsInV * UNext + storedV;
        if (getMIntArrayIndex(vertices, resCV) == -1) vertices.append(resCV);
    }
    // minus U
    int UPrev = storedU - 1;
    if (UPrev >= 0) {
        resCV = numCVsInV * UPrev + storedV;
        if (getMIntArrayIndex(vertices, resCV) == -1) vertices.append(resCV);
    } else if (UIsPeriodic) {
        UPrev += numCVsInU;
        resCV = numCVsInV * UPrev + storedV;
        if (getMIntArrayIndex(vertices, resCV) == -1) vertices.append(resCV);
    }
    // plus V
    int VNext = storedV + 1;
    if (VNext < numCVsInV) {
        resCV = numCVsInV * storedU + VNext;
        if (getMIntArrayIndex(vertices, resCV) == -1) vertices.append(resCV);
    } else if (VIsPeriodic) {
        VNext -= numCVsInV;
        resCV = numCVsInV * storedU + VNext;
        if (getMIntArrayIndex(vertices, resCV) == -1) vertices.append(resCV);
    }
    // minus V
    int VPrev = storedV - 1;
    if (VPrev >= 0) {
        resCV = numCVsInV * storedU + VPrev;
        if (getMIntArrayIndex(vertices, resCV) == -1) vertices.append(resCV);
    } else if (VIsPeriodic) {
        VPrev += numCVsInV;
        resCV = numCVsInV * storedU + VPrev;
        if (getMIntArrayIndex(vertices, resCV) == -1) vertices.append(resCV);
    }
    // vertInd = numCVsInV * indexU + indexV;
}

// from the mesh retrieves the skinCluster
MStatus findSkinCluster(MDagPath MeshPath, MObject& theSkinCluster, int indSkinCluster,
                        bool verbose) {
    if (verbose) MGlobal::displayInfo(MString(" ---- findSkinCluster ----"));
    MStatus stat;

    MFnDagNode dagNode(MeshPath);  // path to the visible mesh
    // MFnMesh meshFn(MeshPath, &stat);     // this is the visible mesh
    MObject inObj;
    MObject dataObj1;

    MObjectArray listSkinClusters;
    // the deformed mesh comes into the visible mesh
    // through its "inmesh" plug
    MPlug inMeshPlug;
    if (MeshPath.apiType() == MFn::kMesh)
        inMeshPlug = dagNode.findPlug("inMesh", &stat);
    else if (MeshPath.apiType() == MFn::kNurbsSurface)
        inMeshPlug = dagNode.findPlug("create", &stat);

    if (stat == MS::kSuccess && inMeshPlug.isConnected()) {
        // walk the tree of stuff upstream from this plug
        MItDependencyGraph dgIt(inMeshPlug, MFn::kInvalid, MItDependencyGraph::kUpstream,
                                MItDependencyGraph::kDepthFirst, MItDependencyGraph::kPlugLevel,
                                &stat);
        if (MS::kSuccess == stat) {
            dgIt.disablePruningOnFilter();
            int count = 0;

            for (; !dgIt.isDone(); dgIt.next()) {
                MObject thisNode = dgIt.thisNode();
                // go until we find a skinCluster
                if (thisNode.apiType() == MFn::kSkinClusterFilter) {
                    listSkinClusters.append(thisNode);
                    // return MS::kSuccess;
                }
            }
        }
        int listSkinClustersLength = listSkinClusters.length();
        if (verbose)
            MGlobal::displayInfo(MString("    nb skinClusters is ") + listSkinClustersLength);
        if (listSkinClustersLength > indSkinCluster) {
            theSkinCluster = listSkinClusters[indSkinCluster];

            MFnDependencyNode nodeFn(theSkinCluster);
            if (verbose)
                MGlobal::displayInfo(MString("    returned skinCluster: ") + nodeFn.name());

            return MS::kSuccess;
        }
        // std::cout << "skinCluster: " << returnedSkinCluster.name().asChar() << "\n";
    }
    return MS::kFailure;
}

MStatus findMesh(MObject& skinCluster, MDagPath& theMeshPath, bool verbose) {
    if (verbose) MGlobal::displayInfo(MString(" ---- findMesh ----"));
    MFnSkinCluster theSkinCluster(skinCluster);
    MObjectArray objectsDeformed;
    theSkinCluster.getOutputGeometry(objectsDeformed);
    int objectsDeformedCount = objectsDeformed.length();
    bool doContinue = false;
    if (objectsDeformedCount == 0) {
        int j = 0;
        // for (int j = 0; j < objectsDeformedCount; j++) {
        theMeshPath.getAPathTo(objectsDeformed[j]);
        if (verbose) {
            MFnDependencyNode deformedNameMesh(objectsDeformed[j]);
            MString deformedNameMeshSTR = deformedNameMesh.name();
            if (verbose) MGlobal::displayInfo("     -> DEFORMING : " + deformedNameMeshSTR + "\n");
        }
        //}
        return MS::kSuccess;
    }
    return MS::kFailure;
}
