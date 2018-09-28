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

MStatus findOrigMesh(MObject& skinCluster, MObject& origMesh, bool verbose) {
    if (verbose) MGlobal::displayInfo(MString(" ---- find Orig Mesh ----"));
    MFnSkinCluster theSkinCluster(skinCluster);
    MObjectArray objectsDeformed;
    theSkinCluster.getInputGeometry(objectsDeformed);
    origMesh = objectsDeformed[0];
    if (verbose) {
        MFnDependencyNode deformedNameMesh(origMesh);
        MGlobal::displayInfo("     -> DEFORMING : " + deformedNameMesh.name() + "\n");
    }
    return MS::kSuccess;
}

MStatus getListColorsJoints(MObject& skinCluster, MColorArray& jointsColors) {
    MStatus stat;
    /*
    MDagPathArray  listOfJoints;
    MFnSkinCluster theSkinCluster(skinCluster);
    theSkinCluster.influenceObjects(listOfJoints, &stat);
    int nbJoints = listOfJoints.length();

    // now get the colors per joint ----------------------------------------
    MStringArray allJointsNames;
    float color[3];
    MColorArray jointsColors(nbJoints);
    for (int i = 0; i < nbJoints; i++) {
    MFnDagNode jnt(listOfJoints[i]);
    //jnt = MFnDagNode (listOfJoints[iterator]);
    MPlug lockInfluenceWeightsPlug = jnt.findPlug("wireColorRGB");
    //MPlug lockInfluenceWeightsPlug(jnt, MFnDagNode::wireColorRGB);

    color[0] = lockInfluenceWeightsPlug.child(0).asFloat();
    color[1] = lockInfluenceWeightsPlug.child(1).asFloat();
    color[2] = lockInfluenceWeightsPlug.child(2).asFloat();

    jointsColors[i] = MColor(color);
    //MString jointName = jnt.name();
    }
    */

    MFnDependencyNode skinClusterDep(skinCluster);
    MPlug influenceColor_plug = skinClusterDep.findPlug("influenceColor");
    int nbElements = influenceColor_plug.numElements();
    // MGlobal::displayInfo(influenceColor_plug.name() + " - " +nbJoints);
    jointsColors.clear();
    MIntArray plugIndices;
    influenceColor_plug.getExistingArrayAttributeIndices(plugIndices);
    int nbJoints = plugIndices[plugIndices.length() - 1] + 1;
    jointsColors.setLength(nbJoints);
    float black[4] = {0, 0, 0, 1};
    for (int i = 0; i < nbElements; ++i) {
        // weightList[i]
        MPlug colorPlug = influenceColor_plug.elementByPhysicalIndex(i);
        int logicalInd = colorPlug.logicalIndex();
        if (colorPlug.isConnected()) {
            MPlugArray connections;
            colorPlug.connectedTo(connections, true, false);
            if (connections.length() > 0) {
                MPlug theConn = connections[0];
                float element[4] = {theConn.child(0).asFloat(), theConn.child(1).asFloat(),
                                    theConn.child(2).asFloat(), 1};
                jointsColors.set(element, logicalInd);
            } else
                jointsColors.set(black, logicalInd);
        } else {
            // MGlobal::displayInfo(colorPlug.name());
            float element[4] = {colorPlug.child(0).asFloat(), colorPlug.child(1).asFloat(),
                                colorPlug.child(2).asFloat(), 1};
            // MGlobal::displayInfo(colorPlug.name()+ " "+ element[0] + " " + element[1] + " " +
            // element[2]);
            jointsColors.set(element, logicalInd);
        }
    }
    return stat;
}

MStatus getListLockJoints(MObject& skinCluster, MIntArray& jointsLocks) {
    MStatus stat;

    MFnDependencyNode skinClusterDep(skinCluster);
    MPlug influenceColor_plug = skinClusterDep.findPlug("lockWeights");

    int nbJoints = influenceColor_plug.numElements();
    jointsLocks.clear();
    jointsLocks.setLength(nbJoints);

    for (int i = 0; i < nbJoints; ++i) {
        // weightList[i]
        MPlug lockPlug = influenceColor_plug.elementByPhysicalIndex(i);
        if (lockPlug.isConnected()) {
            MPlugArray connections;
            lockPlug.connectedTo(connections, true, false);
            if (connections.length() > 0) {
                MPlug theConn = connections[0];
                jointsLocks.set(theConn.asInt(), i);
            } else
                jointsLocks.set(0, i);
        } else {
            jointsLocks.set(lockPlug.asInt(), i);
        }
        // MGlobal::displayInfo(lockPlug.name() + " " + jointsLocks [i]);
    }
    return stat;
}

MStatus getListLockVertices(MObject& skinCluster, MIntArray& vertsLocks) {
    MStatus stat;

    MFnSkinCluster theSkinCluster(skinCluster);
    MObjectArray objectsDeformed;
    theSkinCluster.getOutputGeometry(objectsDeformed);
    MFnDependencyNode deformedNameMesh(objectsDeformed[0]);
    MPlug lockedVerticesPlug = deformedNameMesh.findPlug("lockedVertices", &stat);
    if (MS::kSuccess != stat) {
        MGlobal::displayError(MString("cant find lockerdVertices plug"));
        return stat;
    }

    MFnDependencyNode skinClusterDep(skinCluster);
    MPlug weight_list_plug = skinClusterDep.findPlug("weightList");

    int nbVertices = weight_list_plug.numElements();

    // vertsLocks.clear();
    MObject Data;
    stat = lockedVerticesPlug.getValue(Data);  // to get the attribute

    MFnIntArrayData intData(Data);
    MIntArray vertsLocksIndices = intData.array(&stat);
    vertsLocks.clear();
    vertsLocks = MIntArray(nbVertices, 0);
    for (int i = 0; i < vertsLocksIndices.length(); ++i) vertsLocks[vertsLocksIndices[i]] = 1;
    // MGlobal::displayInfo(MString(" getListLockVertices | ") + currentColorSet.name () + MString("
    // ") + vertsLocks.length());
    return stat;
}

MStatus getSymetryAttributes(MObject& skinCluster, MIntArray& symetryList) {
    MStatus stat;

    MFnSkinCluster theSkinCluster(skinCluster);
    MObjectArray objectsDeformed;
    theSkinCluster.getOutputGeometry(objectsDeformed);

    MFnDagNode deformedNameMesh(objectsDeformed[0]);
    MObject prt = deformedNameMesh.parent(0);

    MFnDependencyNode prtDep(prt);
    MPlug symVerticesPlug = prtDep.findPlug("symmetricVertices", &stat);
    if (MS::kSuccess != stat) {
        MGlobal::displayError(MString("cant find symmetricVertices plug"));
        return stat;
    }

    MObject Data;
    stat = symVerticesPlug.getValue(Data);  // to get the attribute

    MFnIntArrayData intData(Data);
    symetryList = intData.array(&stat);
}

MStatus getMirrorVertices(MIntArray mirrorVertices, MIntArray& theEditVerts,
                          MIntArray& theMirrorVerts, MIntArray& editAndMirrorVerts) {
    MStatus status;

    MIntArray vertExists(mirrorVertices.length(), 0);
    editAndMirrorVerts.copy(theEditVerts);

    theMirrorVerts.setLength(theEditVerts.length());  // to change
    for (int i = 0; i < theEditVerts.length(); ++i) vertExists[theEditVerts[i]] = 1;
    for (int i = 0; i < theEditVerts.length(); ++i) {
        int theVert = theEditVerts[i];
        int theMirroredVert = mirrorVertices[theVert];

        theMirrorVerts[i] = theMirroredVert;     // to change
        if (vertExists[theMirroredVert] == 0) {  // not in first array
            // theMirrorVerts.append(theMirroredVert);
            editAndMirrorVerts.append(theMirroredVert);
        }
    }
    // MGlobal::displayError(MString("theEditVerts ") + theEditVerts.length() + MString("
    // theMirrorVerts ") + theMirrorVerts.length() + MString(" editAndMirrorVerts ") +
    // editAndMirrorVerts.length() );

    return status;
}

MStatus editLocks(MObject& skinCluster, MIntArray& inputVertsToLock, bool addToLock,
                  MIntArray& vertsLocks) {
    MStatus stat;

    MFnSkinCluster theSkinCluster(skinCluster);
    MObjectArray objectsDeformed;
    theSkinCluster.getOutputGeometry(objectsDeformed);
    MFnDependencyNode deformedNameMesh(objectsDeformed[0]);
    MPlug lockedVerticesPlug = deformedNameMesh.findPlug("lockedVertices", &stat);
    if (MS::kSuccess != stat) {
        MGlobal::displayError(MString("cant find lockerdVertices plug"));
        return stat;
    }

    // now expand the array -----------------------
    int val = 0;
    if (addToLock) val = 1;
    for (unsigned int i = 0; i < inputVertsToLock.length(); ++i) {
        int vtx = inputVertsToLock[i];
        vertsLocks[vtx] = val;
    }

    MIntArray theArrayValues;
    for (unsigned int vtx = 0; vtx < vertsLocks.length(); ++vtx) {
        if (vertsLocks[vtx] == 1) theArrayValues.append(vtx);
    }
    // now set the value ---------------------------
    MFnIntArrayData tmpIntArray;
    stat = lockedVerticesPlug.setValue(tmpIntArray.create(theArrayValues));  // to set the attribute
    return stat;
}

MStatus getListColors(MObject& skinCluster, int nbVertices, MColorArray& currColors,
                      bool useMPlug) {
    MStatus stat;
    MFnDagNode skinClusterDag(skinCluster);
    MFnDependencyNode skinClusterDep(skinCluster);
    MColorArray jointsColors;
    getListColorsJoints(skinCluster, jointsColors);

    if (!useMPlug) {
        MFnSkinCluster theSkinCluster(skinCluster);
        int nbJoints = jointsColors.length();

        // now get the weights ----------------------------------------
        MObject allVerticesObj;

        MFnSingleIndexedComponent allVertices;
        MDoubleArray fullOrigWeights;

        allVertices.setCompleteData(nbVertices);
        allVerticesObj = allVertices.create(MFn::kMeshVertComponent);
        unsigned int infCount;

        MDagPath path;
        theSkinCluster.getPathAtIndex(0, path);
        theSkinCluster.getWeights(path, allVerticesObj, fullOrigWeights, infCount);

        // now get the colors per vertices ----------------------------------------
        int currentWeightsLength = fullOrigWeights.length();
        int indexInfluence, i, indexWeight;
        double theWeight, maxWeight;

        currColors.clear();
        currColors.setLength(nbVertices);

        for (i = 0; i < nbVertices; ++i) {
            MColor theColor;
            maxWeight = 0.;
            for (indexInfluence = 0; indexInfluence < nbJoints; indexInfluence++) {
                indexWeight = i * nbJoints + indexInfluence;
                theWeight = fullOrigWeights[indexWeight];
                if (theWeight > 0.) {
                    theColor += jointsColors[indexInfluence] * theWeight;
                }
                currColors[i] = theColor;
            }
        }
    } else {
        MPlug weight_list_plug = skinClusterDep.findPlug("weightList");
        // MGlobal::displayInfo(weight_list_plug.name());
        int nbElements = weight_list_plug.numElements();
        currColors.clear();
        currColors.setLength(nbVertices);

        for (int i = 0; i < nbVertices; ++i) {
            // weightList[i]
            // if (i > 50) break;
            MPlug ith_weights_plug = weight_list_plug.elementByPhysicalIndex(i);
            int vertexIndex = ith_weights_plug.logicalIndex();
            // MGlobal::displayInfo(ith_weights_plug.name());

            // weightList[i].weight
            MPlug plug_weights = ith_weights_plug.child(0);  // access first compound child
            int nb_weights = plug_weights.numElements();
            // MGlobal::displayInfo(plug_weights.name() + nb_weights);
            MColor theColor;
            for (int j = 0; j < nb_weights; j++) {  // for each joint
                MPlug weight_plug = plug_weights.elementByPhysicalIndex(j);
                // weightList[i].weight[j]
                int indexInfluence = weight_plug.logicalIndex();
                double theWeight = weight_plug.asDouble();
                // MGlobal::displayInfo(weight_plug.name() + " " + indexInfluence + " " +
                // theWeight);
                if (theWeight > 0.01) {
                    theColor += jointsColors[indexInfluence] * theWeight;
                }
            }
            currColors[vertexIndex] = theColor;
        }
    }

    return MS::kSuccess;
}

MStatus editArray(int command, int influence, int nbJoints, MIntArray& lockJoints,
                  MDoubleArray& fullWeightArray, MIntArray& vertices, MDoubleArray& verticesWeight,
                  MDoubleArray& theWeights, bool normalize) {
    MStatus stat;
    // 0 Add - 1 Remove - 2 AddPercent - 3 Absolute - 4 Smooth - 5 Sharpen - 6 LockVertices - 7
    // UnLockVertices
    //
    if (command == 5) {  // sharpen  -----------------------
        for (int i = 0; i < vertices.length(); ++i) {
            int theVert = vertices[i];
            double theVal = verticesWeight[i] + 1.0;
            double substract = theVal / nbJoints;
            double sum = 0.0;
            for (int j = 0; j < nbJoints; ++j) {
                // check the zero val ----------
                double jntVal = (fullWeightArray[theVert * nbJoints + j] * theVal) - substract;
                jntVal = std::max(0.0, std::min(jntVal, 1.0));  // clamp
                theWeights[i * nbJoints + j] = jntVal;
                if (normalize) sum += jntVal;
            }
            // now normalize
            if ((normalize) && (sum != 1.0))
                for (int j = 0; j < nbJoints; ++j) theWeights[i * nbJoints + j] /= sum;
        }
    } else {
        // do the command --------------------------
        for (int i = 0; i < vertices.length(); ++i) {
            int theVert = vertices[i];
            double theVal = verticesWeight[i];

            double currentW = fullWeightArray[theVert * nbJoints + influence];
            if (((command == 1) || (command == 3)) &&
                (currentW > 0.99999)) {  // value is 1 we cant do anything
                for (int j = 0; j < nbJoints; ++j)
                    theWeights[i * nbJoints + j] = fullWeightArray[theVert * nbJoints + j];
                continue;
            }
            double newW = currentW;
            if (command == 0)
                newW += theVal;  // ADD
            else if (command == 1)
                newW -= theVal;  // Remove
            else if (command == 2)
                newW += theVal * newW;  // AddPercent
            else if (command == 3)
                newW = theVal;  // Absolute

            newW = std::max(0.0, std::min(newW, 1.0));  // clamp
            double newRest = 1.0 - newW;
            double oldRest = 1.0 - currentW;
            double div = oldRest / newRest;

            double sum = 0.0;
            for (int j = 0; j < nbJoints; ++j) {
                // check the zero val ----------
                if (j != influence) {
                    if (newW == 1.0)
                        theWeights[i * nbJoints + j] = 0.0;
                    else
                        theWeights[i * nbJoints + j] =
                            fullWeightArray[theVert * nbJoints + j] / div;
                } else
                    theWeights[i * nbJoints + j] = newW;

                if (normalize)
                    theWeights[i * nbJoints + j] =
                        std::max(0.0, std::min(theWeights[i * nbJoints + j], 1.0));  // clamp
                sum += theWeights[i * nbJoints + j];
            }
            if (sum < .01)  // zero problem revert weights ----------------------
                for (int j = 0; j < nbJoints; ++j)
                    theWeights[i * nbJoints + j] = fullWeightArray[theVert * nbJoints + j];
            else if (normalize && (sum != 1.0))  // normalize ---------------
                for (int j = 0; j < nbJoints; ++j) theWeights[i * nbJoints + j] /= sum;
        }
    }
    return stat;
}

MStatus setAverageWeight(MIntArray& verticesAround, int currentVertex, int indexCurrVert,
                         int nbJoints, MIntArray& lockJoints, MDoubleArray& fullWeightArray,
                         MDoubleArray& theWeights) {
    MStatus stat;
    int sizeVertices = verticesAround.length();
    int i, j, posi;
    // MGlobal::displayInfo(MString(" paint smooth vtx [")+ currentVertex+ MString("] index - ") +
    // indexCurrVert + MString(" aroundCount ") + sizeVertices);

    MDoubleArray sumWeigths(nbJoints, 0.0);
    // compute sum weights
    for (i = 0; i < sizeVertices; i++) {
        for (j = 0; j < nbJoints; j++) {
            posi = verticesAround[i] * nbJoints + j;
            sumWeigths[j] += fullWeightArray[posi];
        }
    }
    double total = 0.0;
    double totalBaseVtx = 0.0;
    for (j = 0; j < nbJoints; j++) {
        int posi = currentVertex * nbJoints + j;
        sumWeigths[j] /= sizeVertices;

        total += sumWeigths[j];
        totalBaseVtx += fullWeightArray[posi];
    }
    if (total > 0. && totalBaseVtx > 0.) {
        double mult = totalBaseVtx / total;
        for (j = 0; j < nbJoints; j++) {
            int posiToSet = indexCurrVert * nbJoints + j;
            sumWeigths[j] *= mult;  // normalement divide par 1
            theWeights[posiToSet] = sumWeigths[j];
        }
    }
    return MS::kSuccess;
}
