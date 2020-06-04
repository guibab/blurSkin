
#include "blurSkinCmd.h"

#include "functions.h"

const char* blurSkinCmd::kQueryFlagShort = "-q";
const char* blurSkinCmd::kQueryFlagLong = "-query";

const char* blurSkinCmd::kSkinClusterNameFlagShort = "-skn";
const char* blurSkinCmd::kSkinClusterNameFlagLong = "-skinCluster";

const char* blurSkinCmd::kMeshNameFlagShort = "-mn";
const char* blurSkinCmd::kMeshNameFlagLong = "-meshName";

const char* blurSkinCmd::kIndexSkinClusterFlagShort = "-si";
const char* blurSkinCmd::kIndexSkinClusterFlagLong = "-skinClusterIndex";

const char* blurSkinCmd::kPercentMovementFlagShort = "-pc";
const char* blurSkinCmd::kPercentMovementFlagLong = "-percentMvt";

const char* blurSkinCmd::kListJointsFlagShort = "-lj";
const char* blurSkinCmd::kListJointsFlagLong = "-listJoints";

const char* blurSkinCmd::kListJointsValuesFlagShort = "-ljv";
const char* blurSkinCmd::kListJointsValuesFlagLong = "-listJointsValues";

const char* blurSkinCmd::kVerboseFlagShort = "-vrb";
const char* blurSkinCmd::kVerboseFlagLong = "-verbose";

const char* blurSkinCmd::kListVerticesIndicesFlagShort = "-li";
const char* blurSkinCmd::kListVerticesIndicesFlagLong = "-listVerticesIndices";

const char* blurSkinCmd::kListCVsIndicesFlagShort = "-lcv";
const char* blurSkinCmd::kListCVsIndicesFlagLong = "-listCVsIndices";

const char* blurSkinCmd::kListVerticesWeightFlagShort = "-lw";
const char* blurSkinCmd::kListVerticesWeightFlagLong = "-listVerticesWeight";

const char* blurSkinCmd::kRepeatFlagShort = "-rp";
const char* blurSkinCmd::kRepeatFlagLong = "-repeat";

const char* blurSkinCmd::kDepthFlagShort = "-d";
const char* blurSkinCmd::kDepthFlagLong = "-depth";

const char* blurSkinCmd::kRespectLocksFlagShort = "-rl";
const char* blurSkinCmd::kRespectLocksFlagLong = "-respectLocks";

const char* blurSkinCmd::kThresholdFlagShort = "-th";
const char* blurSkinCmd::kThresholdFlagLong = "-threshold";

const char* blurSkinCmd::kZeroInfluencesFlagShort = "-zi";
const char* blurSkinCmd::kZeroInfluencesFlagLong = "-zeroInfluences";

const char* blurSkinCmd::kCommandFlagShort = "-c";
const char* blurSkinCmd::kCommandFlagLong = "-command";

const char* blurSkinCmd::kHelpFlagShort = "-h";
const char* blurSkinCmd::kHelpFlagLong = "-help";

/**
Displays command instructions.
*/
void DisplayHelp() {
    MString help;
    help += "Flags:\n";
    help += "-skinCluster         -skn   String     Name of the skinCluster\n";
    help +=
        "-meshName            -mn    String	 Name of the mesh if skincluster is not passed\n";
    help +=
        "                                        If -skn and -mn are not passed uses selection\n";
    help += "-skinClusterIndex    -si    Int        Index of SkinCluster if -mn passed  \n";
    help += "                                        Default 0\n";
    help +=
        "-percentMvt          -pc    Float      Between 0. and 1. percent of value set default "
        "1.\n";
    help += "-verbose             -vrb   N/A        Verbose print\n";
    help += "-listCVsIndices      -lcv   Strings    List cvs indices [(u1, v1), (u2, v2), ...]\n";
    help += "-listVerticesIndices -li    Strings    List vertices indices\n";
    help += "-listVerticesWeight  -lw    Doubles    List vertices weights\n";
    help += "-listJoints          -lj    String     List joints names\n";
    help += "-listJointsValues    -ljw   Doubles    List joints weights\n";
    help += "-repeat              -rp    Int        Repeat the calculation             default 1\n";
    help += "-depth               -d     Int        Depth for the smooth               default 1\n";
    help +=
        "-respectLocks        -rl    N/A        Respect locks                      default True \n";
    help += "-zeroInfluences      -zi    N/A        get zero columns \n";
    help += "-command             -c     N/A        the command action correct inputs are :\n";
    help +=
        "                                        smooth - add - absolute - percentage - average - "
        "colors - prune\n";
    help += "-threshold			  -th   Double      threshold for the prune weights\n";
    help += "-help                -h     N/A        Display this text.\n";
    MGlobal::displayInfo(help);
}

MSyntax blurSkinCmd::newSyntax() {
    MSyntax syntax;
    syntax.addFlag(kQueryFlagShort, kQueryFlagLong);
    syntax.addFlag(kSkinClusterNameFlagShort, kSkinClusterNameFlagLong, MSyntax::kString);
    syntax.addFlag(kMeshNameFlagShort, kMeshNameFlagLong, MSyntax::kString);
    syntax.addFlag(kPercentMovementFlagShort, kPercentMovementFlagLong, MSyntax::kDouble);
    syntax.addFlag(kIndexSkinClusterFlagShort, kIndexSkinClusterFlagLong, MSyntax::kLong);
    syntax.addFlag(kVerboseFlagShort, kVerboseFlagLong, MSyntax::kBoolean);

    syntax.addFlag(kListCVsIndicesFlagShort, kListCVsIndicesFlagLong, MSyntax::kLong,
                   MSyntax::kLong);
    syntax.makeFlagMultiUse(kListCVsIndicesFlagShort);

    syntax.addFlag(kListVerticesIndicesFlagShort, kListVerticesIndicesFlagLong, MSyntax::kLong);
    syntax.makeFlagMultiUse(kListVerticesIndicesFlagShort);

    syntax.addFlag(kListVerticesWeightFlagShort, kListVerticesWeightFlagLong, MSyntax::kDouble);
    syntax.makeFlagMultiUse(kListVerticesWeightFlagShort);

    syntax.addFlag(kListJointsFlagShort, kListJointsFlagLong, MSyntax::kString);
    syntax.makeFlagMultiUse(kListJointsFlagShort);

    syntax.addFlag(kListJointsValuesFlagShort, kListJointsValuesFlagLong, MSyntax::kDouble);
    syntax.makeFlagMultiUse(kListJointsValuesFlagShort);

    syntax.addFlag(kRepeatFlagShort, kRepeatFlagLong, MSyntax::kLong);
    syntax.addFlag(kDepthFlagShort, kDepthFlagLong, MSyntax::kLong);
    syntax.addFlag(kRespectLocksFlagShort, kRespectLocksFlagLong);
    syntax.addFlag(kThresholdFlagShort, kThresholdFlagLong, MSyntax::kDouble);
    syntax.addFlag(kZeroInfluencesFlagShort, kZeroInfluencesFlagLong);

    syntax.addFlag(kCommandFlagShort, kCommandFlagLong, MSyntax::kString);

    syntax.addFlag(kHelpFlagShort, kHelpFlagLong);
    // syntax.setObjectType(MSyntax::kSelectionList, 0, 255);
    // syntax.useSelectionAsDefault(true);
    return syntax;
}

blurSkinCmd::blurSkinCmd()
    : command_(kCommandSmooth),
      repeat_(1),
      depth_(1),
      percentMvt_(1.0),
      indSkinCluster_(0),
      threshold_(0.0001),
      respectLocks_(true),
      verbose(false) {}

MStatus blurSkinCmd::getListLockJoints() {
    if (verbose) MGlobal::displayInfo(MString(" ---- get list lock joints ----"));
    MStatus stat;
    MDagPathArray listOfJoints;
    MFnSkinCluster theSkinCluster(skinCluster_);
    theSkinCluster.influenceObjects(listOfJoints, &stat);
    nbJoints = listOfJoints.length();
    int nbJntsInput = listJoints_.length();

    if (verbose) MGlobal::displayInfo(MString("    nbJoints : ") + nbJoints);

    // preSet the operation per joint at zero
    for (int i = 0; i < nbJoints; i++) perJointAddingValues_.append(0.0);
    MStringArray allJointsNames;

    for (int i = 0; i < nbJoints; i++) {
        if (nbJntsInput == 0) {
            jointsInputIndices_.append(i);
        };
        MFnDagNode jnt(listOfJoints[i]);
        // jnt = MFnDagNode (listOfJoints[iterator]);
        MPlug lockInfluenceWeightsPlug = jnt.findPlug("lockInfluenceWeights");
        bool isLockInfluenceWeights = lockInfluenceWeightsPlug.asBool();

        if (isLockInfluenceWeights)
            lockJoints.append(1);
        else
            lockJoints.append(0);

        // get the index from the name for our list of joints
        MString jointName = jnt.name();
        int jntIndex = listJoints_.indexOf(jointName);
        if ((jntIndex != -1) &&
            listJointsValues_.length() > i) {  // if the joint is in the list of joints name
            perJointAddingValues_.set(listJointsValues_[jntIndex], i);
        }
        if (verbose) {
            if (isLockInfluenceWeights)
                MGlobal::displayInfo("    " + jnt.name() + " is Locked ");
            else
                MGlobal::displayInfo("    " + jnt.name() + " is not Locked ");
        }
        allJointsNames.append(jointName);
    }

    // set the list of joints -----------------------------------------------
    if (verbose)
        MGlobal::displayInfo(MString(" set the list of joints   ... nbJntsInput : ") + nbJntsInput +
                             MString(" "));
    for (int j = 0; j < nbJntsInput; ++j) {
        MString jntInput = listJoints_[j];
        int indexInRealJoints = allJointsNames.indexOf(jntInput);
        jointsInputIndices_.append(indexInRealJoints);
    }
    if (verbose) MGlobal::displayInfo(" end set the list of joints   ...");

    return MS::kSuccess;
}

MStatus blurSkinCmd::printWeigth(int vertex, int u, int v) {
    if (verbose) MGlobal::displayInfo(MString(" ---- printWeigth ----"));
    MFnSkinCluster theSkinCluster(skinCluster_);

    // 3 get the weights
    MDoubleArray weightsVertex;
    unsigned int infCount;
    MString toDisplay;
    MObject tmpComponent;
    if (isMeshSurface_) {
        MFnSingleIndexedComponent theVertex;
        tmpComponent = theVertex.create(MFn::kMeshVertComponent);
        theVertex.create(MFn::kMeshVertComponent);
        theVertex.addElement(vertex);
        toDisplay = MString("weigth of vtx (") + vertex + MString(") : ");
    } else if (isNurbsSurface_) {
        MFnDoubleIndexedComponent doubleFn;
        // MFnSingleIndexedComponent theVertex;
        tmpComponent = doubleFn.create(MFn::kSurfaceCVComponent);

        // MItGeometry gIter(meshPath_);
        // gIter.;
        doubleFn.addElement(u, v);
        toDisplay = MString("weigth of CV (") + u + MString(",") + v + MString(") : ");
        ;
    }
    theSkinCluster.getWeights(meshPath_, tmpComponent, weightsVertex, infCount);
    for (int j = 0; j < nbJoints; j++) {
        toDisplay += MString("[") + weightsVertex[j] + MString("]");
    }
    toDisplay += MString("\n stored value is ");
    for (int j = 0; j < nbJoints; j++) {
        int posiToSet = vertex * nbJoints + j;
        toDisplay += MString("[") + newWeights[posiToSet] + MString("]");
    }

    MGlobal::displayInfo(toDisplay);
    return MS::kSuccess;
}

MStatus blurSkinCmd::getAverageWeight(MIntArray vertices, int currentVertex) {
    if (verbose) MGlobal::displayInfo(MString(" ---- getAverageWeight ----"));
    if (verbose) MGlobal::displayInfo(MString("nbJoints ") + nbJoints);
    MStatus stat;
    int sizeVertices = vertices.length();
    int i, j, posi;

    MDoubleArray sumWeigths(nbJoints, 0.0);
    // compute sum weights
    for (i = 0; i < sizeVertices; ++i) {
        for (j = 0; j < nbJoints; ++j) {
            posi = vertices[i] * nbJoints + j;
            sumWeigths[j] += currentWeights[posi];
        }
    }

    double totalBaseVtxUnlock = 0.0, totalBaseVtxLock = 0.0;
    ;
    double totalVtxUnlock = 0.0, totalVtxLock = 0.0;

    for (j = 0; j < nbJoints; j++) {
        sumWeigths[j] /= sizeVertices;

        // get if jnt is locked
        bool isLockJnt = lockJoints[j] == 1;
        // get currentWeight of currentVtx
        int posiToSet = currentVertex * nbJoints + j;
        double currentW = currentWeights[posiToSet];
        double targetW = sumWeigths[j];

        // sum it all
        if (!isLockJnt) {
            totalBaseVtxUnlock += currentW;
            totalVtxUnlock += targetW;
        } else {
            totalBaseVtxLock += currentW;
            totalVtxLock += targetW;
        }
    }
    // setting part ---------------
    double normalizedValueAvailable = 1.0 - totalBaseVtxLock;

    if (normalizedValueAvailable > 0.0 && totalVtxUnlock > 0.0) {  // we have room to set weights
        double mult = normalizedValueAvailable / totalVtxUnlock;
        for (j = 0; j < nbJoints; j++) {
            bool isLockJnt = lockJoints[j] == 1;
            int posiToSet = currentVertex * nbJoints + j;
            double currentW = currentWeights[posiToSet];
            double targetW = sumWeigths[j];

            if (isLockJnt) {
                newWeights.set(currentW, posiToSet);
            } else {
                targetW *= mult;  // normalement divide par 1, sauf cas lock joints
                newWeights.set(targetW, posiToSet);
            }
        }
    }
    return MS::kSuccess;
}

void blurSkinCmd::verboseSetWeights(int currentVertex) {
    if (verbose) {
        MString toDisplay("new weigth of vtx (");
        toDisplay += currentVertex;
        toDisplay += MString(") : ");
        double total = 0.0;
        for (int j = 0; j < nbJoints; j++) {
            int posiToSet = currentVertex * nbJoints + j;
            toDisplay += MString("[") + newWeights[posiToSet] + MString("]");
            total += newWeights[posiToSet];
        }
        toDisplay += MString(" = ") + total;
        MGlobal::displayInfo(toDisplay);
    }
}

MStatus blurSkinCmd::addWeights(int currentVertex) {
    if (verbose) {
        if (command_ == kCommandAdd)
            MGlobal::displayInfo(MString(" command is Add"));
        else if (command_ == kCommandAbsolute)
            MGlobal::displayInfo(MString(" command is Absolute"));
        else if (command_ == kCommandPercentage)
            MGlobal::displayInfo(MString(" command is Percentage"));
    }

    double totalWithLocks = 0.0;
    double totalOfAdding = 0.0;
    double otherJointsTotal = 0.0;
    int posi, j;

    int nbAbsolute = 0;
    // get the totals of the weights ------------------------------
    // ------------------------------------------------------------
    for (j = 0; j < nbJoints; j++) {
        posi = currentVertex * nbJoints + j;
        double value = currentWeights[posi];
        if (!lockJoints[j]) {
            totalWithLocks += value;
            if (perJointAddingValues_[j] != 0) {
                if (command_ == kCommandAdd) {
                    totalOfAdding += value + perJointAddingValues_[j];
                } else if (command_ == kCommandAbsolute) {
                    totalOfAdding += perJointAddingValues_[j];
                    nbAbsolute++;
                } else if (command_ == kCommandPercentage) {
                    totalOfAdding += value * (1 + perJointAddingValues_[j]);
                }
            } else {
                otherJointsTotal += value;
            }
        }
    }
    // it we're trying to set an absolute on only one joint  and nowhere to set the rest!!
    if (command_ == kCommandAbsolute) {
        if (otherJointsTotal == 0. && nbAbsolute == 1)  // we can't set
            return MS::kSuccess;
        /*
        if (totalOfAdding+ otherJointsTotal != totalWithLocks)
                totalOfAdding = totalWithLocks - otherJointsTotal;
        */
    }

    // now do the setting of the values ---------------------------
    // ------------------------------------------------------------
    double newVal;
    if ((totalOfAdding >= totalWithLocks) ||
        (command_ == kCommandAbsolute && otherJointsTotal == 0.)) {
        // normalize, rest of joints get zero
        double multFactor = totalWithLocks / totalOfAdding;
        for (j = 0; j < nbJoints; j++) {
            posi = currentVertex * nbJoints + j;
            double value = currentWeights[posi];
            if (!lockJoints[j]) {
                if (perJointAddingValues_[j] != 0) {
                    if (command_ == kCommandAdd) {
                        newVal = (value + perJointAddingValues_[j]) * multFactor;
                    } else if (command_ == kCommandAbsolute) {
                        newVal = perJointAddingValues_[j] * multFactor;
                    } else if (command_ == kCommandPercentage) {
                        newVal = (value * (1 + perJointAddingValues_[j])) * multFactor;
                    }
                } else {
                    newVal = 0;
                }
                newWeights.set(newVal, posi);
            }
        }
    } else {
        // rest is the value to set to the other joints
        double rest = totalWithLocks - totalOfAdding;
        double multFactor = rest / otherJointsTotal;
        for (j = 0; j < nbJoints; j++) {
            posi = currentVertex * nbJoints + j;
            double value = currentWeights[posi];
            if (!lockJoints[j]) {
                if (perJointAddingValues_[j] != 0) {
                    if (command_ == kCommandAdd) {
                        newVal = value + perJointAddingValues_[j];
                    } else if (command_ == kCommandAbsolute) {
                        newVal = perJointAddingValues_[j];
                    } else if (command_ == kCommandPercentage) {
                        newVal = value * (1 + perJointAddingValues_[j]);
                    }
                } else {
                    newVal = value * multFactor;
                }
                newWeights.set(newVal, posi);
            }
        }
    }
    // print of result computation
    verboseSetWeights(currentVertex);

    return MS::kSuccess;
}

MStatus blurSkinCmd::getSoftSelection(bool getSoft) {
    if (verbose) MGlobal::displayInfo(MString("---- getSoftSelection ----"));
    MStatus stat;

    MSelectionList richSelList;
    if (getSoft) {
        MRichSelection richSel;
        MGlobal::getRichSelection(richSel);
        richSel.getSelection(richSelList);
    } else {
        MGlobal::getActiveSelectionList(richSelList);
    }

    for (MItSelectionList iter(richSelList, MFn::kInvalid); !iter.isDone(); iter.next()) {
        // MDagPath meshPath_;
        iter.getDagPath(meshPath_, component);

        if (verbose) {
            MGlobal::displayInfo("     softSelection on mesh   " + meshPath_.fullPathName());
            if (component.apiType() == MFn::kMeshVertComponent) {
                MGlobal::displayInfo("     vertices are selected");
            } else if (component.apiType() == MFn::kSurfaceCVComponent) {
                MGlobal::displayInfo("     CVs are selected");
            }
        }
        if (component.isNull() || !(component.apiType() == MFn::kMeshVertComponent ||
                                    component.apiType() == MFn::kSurfaceCVComponent)) {
            // MFn::kSurfaceCVComponent
            // do on all vertices
            if (verbose) MGlobal::displayInfo(MString("     do on all vertices"));
            meshPath_.extendToShape();
            getTypeOfSurface();
            useAllVertices();
        } else {
            getTypeOfSurface();
            MFnComponent componentFn(component);
            int count = componentFn.elementCount();
            MStatus stat;
            if (verbose) MGlobal::displayInfo(MString("     check selection and weights"));

            if (component.apiType() == MFn::kMeshVertComponent) {
                MFnSingleIndexedComponent singleFn(component, &stat);
                if (MS::kSuccess == stat) {
                    for (int i = 0; i < count; i++) {
                        MWeight weight = componentFn.weight(i);
                        int vertInd = singleFn.element(i);
                        float influence = weight.influence();
                        if (verbose)
                            MGlobal::displayInfo(
                                MString("      Vertex[") + vertInd +
                                MString("] has influence weight ") +
                                influence);  //+ MString(" and seam weight ") + weight.seam());
                        indicesVertices_.append(vertInd);
                        weightVertices_.append(influence);
                    }
                }
            } else if (component.apiType() == MFn::kSurfaceCVComponent) {
                MFnDoubleIndexedComponent doubleFn(component, &stat);
                if (MS::kSuccess == stat) {
                    // MFnNurbsSurface nurbsFn(meshPath_, &stat);
                    int sizeInV = numCVsInV_;  // nurbsFn.numCVsInV();
                    int vertInd;
                    for (int i = 0; i < count; i++) {
                        MWeight weight = componentFn.weight(i);
                        float influence = weight.influence();
                        int u, v;
                        // from help
                        // index = numCVsInV() * uIndex + vIndex

                        doubleFn.getElement(i, u, v);
                        vertInd = sizeInV * u + v;
                        indicesVertices_.append(vertInd);
                        weightVertices_.append(influence);
                        indicesU_.append(u);
                        indicesV_.append(v);
                        if (verbose)
                            MGlobal::displayInfo(
                                MString("      Component vertInd -") + vertInd + MString("- [") +
                                u + MString(",") + v + MString("] has influence weight ") +
                                weight.influence() + MString(" and seam weight ") + weight.seam());
                    }
                }
            }
        }
    }
    return MS::kSuccess;
}

MStatus blurSkinCmd::useAllVertices() {
    if (verbose) MGlobal::displayInfo(MString("---- useAllVertices ----"));

    MStatus stat;
    indicesVertices_.clear();
    weightVertices_.clear();

    if (meshPath_.apiType() == MFn::kMesh) {  // if is mesh
        MFnMesh meshFn(meshPath_, &stat);     // this is the visible mesh
        MIntArray ObjVertices;
        int nbVertices = meshFn.numVertices(&stat);
        for (int i = 0; i < nbVertices; i++) {
            ObjVertices.append(i);
            indicesVertices_.append(i);
            weightVertices_.append(1.0);
        }
        MFnSingleIndexedComponent allVertices;
        component = allVertices.create(MFn::kMeshVertComponent);
        allVertices.addElements(ObjVertices);
    } else if (meshPath_.apiType() == MFn::kNurbsSurface) {  // if is nurbs
        // cmds.select("nurbsPlane1.cv[*][*]")
        if (verbose)
            MGlobal::displayInfo(MString("useAllVertices  :    numCVsInU ") + numCVsInU_ +
                                 MString(" numCVsInV ") + numCVsInV_);

        int indexU, indexV;
        int vertInd;
        indicesU_.clear();
        indicesV_.clear();
        for (indexU = 0; indexU < numCVsInU_; ++indexU) {
            for (indexV = 0; indexV < numCVsInV_; ++indexV) {
                vertInd = numCVsInV_ * indexU + indexV;

                indicesVertices_.append(vertInd);
                indicesU_.append(indexU);
                indicesV_.append(indexV);
                weightVertices_.append(1.0);
            }
        }
    }
    return MS::kSuccess;
}

MStatus blurSkinCmd::setColors() {
    MStatus stat;

    MFnMesh meshFn(meshPath_, &stat);  // this is the visible mesh
    int nbVertices = meshFn.numVertices();

    MColorArray theColors;
    getListColors(skinCluster_, nbVertices, theColors, true);

    MIntArray vertexIndices;
    vertexIndices.setLength(nbVertices);
    for (unsigned int i = 0; i < nbVertices; i++) {
        vertexIndices[i] = i;
    }
    MObject origMeshObj;
    findOrigMesh(skinCluster_, origMeshObj, false);

    MFnMesh meshFnOrig(origMeshObj, &stat);  // this is the orig mesh
    meshFnOrig.setVertexColors(theColors, vertexIndices);
    meshFnOrig.setDisplayColors(true);
    meshFn.setDisplayColors(true);
    return stat;
}

MIntArray blurSkinCmd::getZeroInfluences() {
    MStatus status = MS::kSuccess;
    // This plug is an array (one element for each vertex in your mesh
    MFnDependencyNode skinClusterDep(skinCluster_);

    MPlug weight_list_plug = skinClusterDep.findPlug("weightList");
    MIntArray jointUsed(nbJoints, 0);
    int nbAt1 = 0;
    for (int i = 0; i < indicesVertices_.length(); ++i) {
        if (nbAt1 == nbJoints) break;  // allFound

        int vertexIndex = indicesVertices_[i];
        // weightList[i]
        MPlug ith_weights_plug = weight_list_plug.elementByLogicalIndex(vertexIndex);

        // weightList[i].weight
        MPlug plug_weights = ith_weights_plug.child(0);  // access first compound child
        int nb_weights = plug_weights.numElements();

        for (int j = 0; j < nb_weights; j++) {  // for each joint
            MPlug weight_plug = plug_weights.elementByPhysicalIndex(j);
            // weightList[i].weight[j]
            int indexInfluence = weight_plug.logicalIndex();
            if (jointUsed[indexInfluence] == 0) {  // check the value if zero or not
                double theWeight = weight_plug.asDouble();
                if (theWeight != 0) {
                    jointUsed[indexInfluence] = 1;
                    nbAt1++;
                }
            }
        }
    }
    if (verbose) MGlobal::displayInfo(" get zero Columns ");

    return jointUsed;
}

MStatus blurSkinCmd::getAllWeights() {
    if (verbose) MGlobal::displayInfo(MString(" ---- getAllWeights ----"));
    MFnSkinCluster theSkinCluster(skinCluster_);

    MStatus stat;
    MObject allVerticesObj;
    if (meshPath_.apiType() == MFn::kMesh) {  // if is mesh

        MFnMesh meshFn(meshPath_, &stat);  // this is the visible mesh
        MIntArray ObjVertices;
        if (verbose) MGlobal::displayInfo(MString("    mesh is") + meshPath_.fullPathName());

        int nbVertices = meshFn.numVertices(&stat);
        for (int i = 0; i < nbVertices; i++) ObjVertices.append(i);

        MFnSingleIndexedComponent allVertices;
        allVertices.addElements(ObjVertices);
        allVerticesObj = allVertices.create(MFn::kMeshVertComponent);
        unsigned int infCount;
        stat = theSkinCluster.getWeights(meshPath_, allVerticesObj, fullOrigWeights, infCount);
        if (stat == MS::kFailure)
            MGlobal::displayError(MString(" can not get the skin weights for mesh ") +
                                  meshPath_.fullPathName());
    } else if (meshPath_.apiType() == MFn::kNurbsSurface) {  // if is nurbs
        /*
        MFnDoubleIndexedComponent allCVs;
        MFnNurbsSurface MfnSurface(meshPath_);
        int sizeInV = MfnSurface.numCVsInV();
        int sizeInU = MfnSurface.numCVsInU();

        for (int indexU = 0; indexU < sizeInU; sizeInU++) {
                for (int indexV = 0; indexV < sizeInV; sizeInV++) {
                        allCVs.addElement(indexU, indexV);
                }
        }
        allVerticesObj = allCVs.create(MFn::kSurfaceCVComponent);
        allCVs.setCompleteData( sizeInU,  sizeInV);

        */
        MItGeometry gIter(meshPath_);
        MDoubleArray wts;

        for (; !gIter.isDone(); gIter.next()) {
            MObject comp = gIter.component(&stat);
            // Get the weights for this vertex (one per influence object)
            //
            unsigned int infCount;
            stat = theSkinCluster.getWeights(meshPath_, comp, wts, infCount);
            for (unsigned int i = 0; i < infCount; ++i) fullOrigWeights.append(wts[i]);
        }
        weigthsForUndo.copy(fullOrigWeights);
    }
    if (verbose)
        MGlobal::displayInfo(MString("    full weights nb weights ") + fullOrigWeights.length());

    currentWeights.copy(fullOrigWeights);
    newWeights.copy(fullOrigWeights);
    return MS::kSuccess;
}

void blurSkinCmd::getTypeOfSurface() {
    if (verbose) MGlobal::displayInfo(MString("---- getTypeOfSurface ----"));

    MFn::Type fType = meshPath_.apiType();
    isNurbsSurface_ = fType == MFn::kNurbsSurface;
    isMeshSurface_ = fType == MFn::kMesh;
    isNurbsCurve_ = fType == MFn::kNurbsCurve;
    isBezierCurve_ = fType == MFn::kBezierCurve;
    MString info = MString("     type of shape  is  : ");
    if (isNurbsSurface_)
        info += MString("Nurbs Surface");
    else if (isMeshSurface_)
        info += MString("Mesh");
    else if (isNurbsCurve_)
        info += MString("Nurbs Curve");
    else if (isBezierCurve_)
        info += MString("Bezier Curve");
    if (verbose) MGlobal::displayInfo(info);

    if (isNurbsSurface_) {
        MFnNurbsSurface MfnSurface(meshPath_);
        numCVsInV_ = MfnSurface.numCVsInV();
        numCVsInU_ = MfnSurface.numCVsInU();
        MFnNurbsSurface::Form formInU = MfnSurface.formInU();
        UIsPeriodic_ = formInU == MFnNurbsSurface::kPeriodic;
        MFnNurbsSurface::Form formInV = MfnSurface.formInV();
        VIsPeriodic_ = formInV == MFnNurbsSurface::kPeriodic;
        UDeg_ = MfnSurface.degreeU();
        VDeg_ = MfnSurface.degreeV();
        // int vertInd;
        if (VIsPeriodic_) numCVsInV_ -= VDeg_;
        if (UIsPeriodic_) numCVsInU_ -= UDeg_;
    }
}

MStatus blurSkinCmd::GatherCommandArguments(const MArgList& args) {
    if (verbose) MGlobal::displayInfo(MString(" ---- GatherCommandArguments ----"));
    MStatus status;
    MArgDatabase argData(syntax(), args);
    // argData.getObjects(selectionList_);
    // display help  --------------------------------------------------------------------------
    if (argData.isFlagSet(kHelpFlagShort)) {
        command_ = kCommandHelp;
        return MS::kSuccess;
    }

    // get the type of the command -------------------------------------------------------
    if (argData.isFlagSet(kCommandFlagShort)) {
        MString commandStringName = argData.flagArgumentString(kCommandFlagShort, 0, &status);
        if (commandStringName == "smooth")
            command_ = kCommandSmooth;
        else if (commandStringName == "add")
            command_ = kCommandAdd;
        else if (commandStringName == "absolute")
            command_ = kCommandAbsolute;
        else if (commandStringName == "percentage")
            command_ = kCommandPercentage;
        else if (commandStringName == "average")
            command_ = kCommandAverage;
        else if (commandStringName == "colors")
            command_ = kCommandSetColors;
        else if (commandStringName == "prune")
            command_ = kCommandPruneWeights;
    }
    // overwrites commands
    // --------------------------------------------------------------------------
    if (argData.isFlagSet(kQueryFlagShort)) {
        command_ = kCommandQuery;
    }
    // get basic arguments ----------------------------------------------------------------------
    if (argData.isFlagSet(kVerboseFlagShort))
        verbose = argData.flagArgumentBool(kVerboseFlagShort, 0, &status);

    if (argData.isFlagSet(kRepeatFlagShort))
        repeat_ = argData.flagArgumentInt(kRepeatFlagShort, 0, &status);

    if (argData.isFlagSet(kPercentMovementFlagShort))
        percentMvt_ = argData.flagArgumentDouble(kPercentMovementFlagShort, 0, &status);

    if (argData.isFlagSet(kThresholdFlagShort))
        threshold_ = argData.flagArgumentDouble(kThresholdFlagShort, 0, &status);

    if (argData.isFlagSet(kDepthFlagShort))
        depth_ = argData.flagArgumentInt(kDepthFlagShort, 0, &status);

    if (argData.isFlagSet(kRespectLocksFlagShort))
        respectLocks_ = argData.flagArgumentBool(kRespectLocksFlagShort, 0, &status);

    if (argData.isFlagSet(kZeroInfluencesFlagShort)) {
        getZeroInfluences_ = argData.flagArgumentBool(kZeroInfluencesFlagShort, 0, &status);
        command_ = kCommandGetZeroInfluences;
    }

    // get list input joints  --------------------------------------------------------------------
    if (argData.isFlagSet(kListJointsFlagShort)) {
        int nbUse = argData.numberOfFlagUses(kListJointsFlagShort);
        MString toDisplay("List Joints : ");
        for (int i = 0; i < nbUse; i++) {
            MArgList flagArgs;
            argData.getFlagArgumentList(kListJointsFlagShort, i, flagArgs);
            MString jointName;
            flagArgs.get(0, jointName);
            listJoints_.append(jointName);
            toDisplay += jointName;
            toDisplay += MString(" - ");
        }
        if (verbose) MGlobal::displayInfo(toDisplay);
    }
    if (argData.isFlagSet(kListJointsValuesFlagShort)) {
        int nbUse = argData.numberOfFlagUses(kListJointsValuesFlagShort);
        MString toDisplay("Joints values : ");
        for (int i = 0; i < nbUse; i++) {
            MArgList flagArgs;
            argData.getFlagArgumentList(kListJointsValuesFlagShort, i, flagArgs);
            double value;
            flagArgs.get(0, value);
            listJointsValues_.append(value);
            toDisplay += value;
            toDisplay += MString(" - ");
        }
        if (verbose) MGlobal::displayInfo(toDisplay);
    }

    // get list input vertices --------------------------------------------------------------------
    bool foundListVerticesIndices = false;
    int nbVerts = 0;
    if (argData.isFlagSet(kListCVsIndicesFlagShort)) {
        int nbUse = argData.numberOfFlagUses(kListCVsIndicesFlagShort);
        MString toDisplay("List CVs : ");
        for (int i = 0; i < nbUse; i++) {
            MArgList flagArgs;
            argData.getFlagArgumentList(kListCVsIndicesFlagShort, i, flagArgs);
            int uIndex, vIndex;
            flagArgs.get(0, uIndex);
            flagArgs.get(1, vIndex);
            indicesU_.append(uIndex);
            indicesV_.append(vIndex);
            toDisplay += uIndex;
            toDisplay += MString("_");
            toDisplay += vIndex;
            toDisplay += MString(" - ");
        }
        if (verbose) MGlobal::displayInfo(toDisplay);
        foundListVerticesIndices = true;
        nbVerts = indicesU_.length();
    } else if (argData.isFlagSet(kListVerticesIndicesFlagShort)) {
        foundListVerticesIndices = true;
        int nbUse = argData.numberOfFlagUses(kListVerticesIndicesFlagShort);
        MString toDisplay("List Vertices : ");
        for (int i = 0; i < nbUse; i++) {
            MArgList flagArgs;
            argData.getFlagArgumentList(kListVerticesIndicesFlagShort, i, flagArgs);
            int vtxIndex;
            flagArgs.get(0, vtxIndex);
            indicesVertices_.append(vtxIndex);
            toDisplay += vtxIndex;
            toDisplay += MString(" - ");
        }
        if (verbose) MGlobal::displayInfo(toDisplay);
        nbVerts = indicesVertices_.length();
    }
    if (argData.isFlagSet(kListVerticesWeightFlagShort)) {
        int nbUse = argData.numberOfFlagUses(kListVerticesWeightFlagShort);
        MString toDisplay("Weight Vertices : ");
        for (int i = 0; i < nbUse; i++) {
            MArgList flagArgs;
            argData.getFlagArgumentList(kListVerticesWeightFlagShort, i, flagArgs);
            double vtxWeight;
            flagArgs.get(0, vtxWeight);
            weightVertices_.append(float(vtxWeight));
            toDisplay += vtxWeight;
            toDisplay += MString(" - ");
        }
        if (verbose) MGlobal::displayInfo(toDisplay);
    } else if (foundListVerticesIndices) {
        for (int i = 0; i < nbVerts; i++) weightVertices_.append(1.0);
    }

    // get index skinCluster (default 0) -------------------------------------------------------
    if (argData.isFlagSet(kIndexSkinClusterFlagShort)) {
        indSkinCluster_ = argData.flagArgumentInt(kIndexSkinClusterFlagShort, 0, &status);
    }

    // get input skinCluster name -------------------------------------------------------
    bool foundSkinCluster = false;
    useSelection = true;

    if (argData.isFlagSet(kSkinClusterNameFlagShort)) {
        MString skinClusterName = argData.flagArgumentString(kSkinClusterNameFlagShort, 0, &status);
        MSelectionList selList;
        MGlobal::getSelectionListByName(skinClusterName, selList);
        selList.getDependNode(0, skinCluster_);

        MFnDependencyNode nodeFn(skinCluster_);
        if (verbose) MGlobal::displayInfo(MString("    input skin name: ") + nodeFn.name());

        selList.clear();
        foundSkinCluster = true;
        // now get the mesh .... seems to CRASH
        useSelection = false;

    }  // OR get input mesh name
    if (argData.isFlagSet(kMeshNameFlagShort)) {
        MString meshName = argData.flagArgumentString(kMeshNameFlagShort, 0, &status);
        MSelectionList selList;
        MGlobal::getSelectionListByName(meshName, selList);
        selList.getDagPath(0, meshPath_);
        selList.clear();
        //
        if (!foundSkinCluster)
            status = findSkinCluster(meshPath_, skinCluster_, indSkinCluster_, verbose);
        useSelection = false;

    } else if (foundSkinCluster) {
        status = findMesh(skinCluster_, meshPath_, verbose);
    }
    // now get input of vertices indices (if selection the return already done previously)
    //-------------------------------------------------------

    return MS::kSuccess;
}

MStatus blurSkinCmd::doIt(const MArgList& args) {
    MStatus stat;
    MString info;
    info += "\n";

    stat = GatherCommandArguments(args);
    CHECK_MSTATUS_AND_RETURN_IT(stat);
    if (verbose) {
        info += MString("\n    percentMvt : ") + percentMvt_;
        info += MString("\n    repeat : ") + repeat_;
        info += MString("\n    depth : ") + depth_;
        info += MString("\n    verbose :") + verbose;
        info += MString("\n    respectLocks_ :") + respectLocks_;
        info += MString("\n    useSelection :") + useSelection;
        info += MString("\n    indSkinCluster :") + indSkinCluster_ + MString("\n");

        MGlobal::displayInfo(info);
    }
    if (command_ == kCommandHelp) {
        DisplayHelp();
        return MS::kSuccess;
    }

    // create the component
    if (useSelection) {
        getSoftSelection();  // dont get soft selection
        stat = findSkinCluster(meshPath_, skinCluster_, indSkinCluster_, verbose);
        if (stat == MS::kFailure) {
            MGlobal::displayError("cant find skin Cluster");
            return MS::kFailure;
        }
    } else {
        getTypeOfSurface();
        if (((isMeshSurface_) && (indicesVertices_.length() == 0)) ||
            ((isNurbsSurface_) && (indicesU_.length() == 0))) {
            useAllVertices();
        } else {
            if (isMeshSurface_) {
                // build array of vertices with weight
                MFnSingleIndexedComponent theVertices;
                component = theVertices.create(MFn::kMeshVertComponent);
                theVertices.addElements(indicesVertices_);

                MFnComponent componentFn(component);
                for (unsigned int i = 0; i < indicesVertices_.length(); i++) {
                    MWeight wght = componentFn.weight(i);
                    wght.setInfluence(weightVertices_[i]);
                    componentFn.setWeight(i, wght);
                }
            } else if (isNurbsSurface_) {
                // MFnNurbsSurface nurbsFn(meshPath_, &stat);
                CHECK_MSTATUS_AND_RETURN_IT(stat);
                int sizeInV = numCVsInV_;  // nurbsFn.numCVsInV();
                int vertInd;
                int u, v;
                float influence;

                for (int i = 0; i < indicesU_.length(); ++i) {
                    influence = weightVertices_[i];
                    u = indicesU_[i];
                    v = indicesV_[i];
                    vertInd = sizeInV * u + v;
                    indicesVertices_.append(vertInd);
                    weightVertices_.append(influence);
                }
            } else {
                return MS::kSuccess;
            }
        }
    }

    executeAction();

    return MS::kSuccess;
}

MStatus blurSkinCmd::executeAction() {
    if (verbose) MGlobal::displayInfo(MString(" ---- executeAction ----"));
    MStatus stat;
    if (verbose) {
        if (component.isNull())
            MGlobal::displayInfo(MString(" component is null "));
        else
            MGlobal::displayInfo(MString(" component is NOT null "));
        if (component.apiType() == MFn::kMeshVertComponent)
            MGlobal::displayInfo(MString(" component is vertices"));
        else
            MGlobal::displayInfo(MString(" component is NOT vertices"));
    }

    if (command_ == kCommandGetZeroInfluences) {
        MIntArray usedInfluences = getZeroInfluences();
        for (int i = 0; i < usedInfluences.length(); ++i) appendToResult(usedInfluences[i]);
        return MS::kSuccess;
    }

    // 1 first get list locked joints
    getListLockJoints();
    // 2 get all weights of vertices
    getAllWeights();
    // 3 get list of locked vertices
    getListLockVertices(skinCluster_, lockVertices_);
    if (command_ == kCommandSetColors) {
        if (verbose) MGlobal::displayInfo(MString(" ---- set Colors ----"));
        setColors();
    } else if (command_ == kCommandQuery) {
        if (verbose) MGlobal::displayInfo(MString(" ---- QUERY RETURN ----"));
        int nbVertices = indicesVertices_.length();
        int index, j, posiToSet;
        int nbJntsInput = jointsInputIndices_.length();
        if (verbose) {
            MString toDisplay = MString("joints indices ");
            for (int k = 0; k < nbJntsInput; ++k) {
                j = jointsInputIndices_[k];
                toDisplay += MString("  ") + j;
            }
            MGlobal::displayInfo(toDisplay);
        }
        //------- now do the setting ----------
        for (int i = 0; i < nbVertices; ++i) {
            index = indicesVertices_[i];
            // for (int j = 0; j<nbJoints; ++j) {
            for (int k = 0; k < nbJntsInput; ++k) {
                j = jointsInputIndices_[k];
                if (j == -1) {
                    appendToResult(-1.);
                } else {
                    posiToSet = index * nbJoints + j;
                    appendToResult(newWeights[posiToSet]);
                    // if (perJointAddingValues_[j] != 0) {
                }
            }
        }
        // appendToResult(1.0);
        return MS::kSuccess;
    } else if (command_ == kCommandPruneWeights) {
        MGlobal::displayInfo(MString("Prune Weights , threshold is  ") + threshold_ +
                             MString("currentWeights  length ") + currentWeights.length());
        doPruneWeight(currentWeights, nbJoints, threshold_);
        // return MS::kSuccess;
    } else if (command_ == kCommandAverage) {
        MDoubleArray averageWeights;
        int index;
        int nbVertices = indicesVertices_.length();
        for (int j = 0; j < nbJoints; j++) {
            double sumJointWeights = 0;
            for (int i = 0; i < nbVertices; ++i) {
                index = indicesVertices_[i];
                int posiToSet = index * nbJoints + j;
                sumJointWeights += newWeights[posiToSet];
            }
            averageWeights.append(sumJointWeights / nbVertices);
        }
        //------- now do the setting ----------
        for (int i = 0; i < nbVertices; ++i) {
            for (int j = 0; j < nbJoints; j++) {
                index = indicesVertices_[i];
                int posiToSet = index * nbJoints + j;
                newWeights.set(averageWeights[j], posiToSet);
            }
        }
        currentWeights.copy(newWeights);
    } else if (isNurbsSurface_) {
        int index, storedU, storedV;
        MIntArray vertices;

        for (int r = 0; r < repeat_; r++) {
            if (verbose) MGlobal::displayInfo(MString("repeat nb :") + r);
            for (int i = 0; i < indicesVertices_.length(); ++i) {
                index = indicesVertices_[i];
                if (lockVertices_[index] != 1) {  // if not locked

                    storedU = indicesU_[i];
                    storedV = indicesV_[i];
                    if (verbose)
                        MGlobal::displayInfo(MString("stored are          U :") + storedU +
                                             MString(" V :") + storedV);

                    if (verbose) stat = printWeigth(index, storedU, storedV);
                    if (command_ != kCommandSmooth) {
                        addWeights(index);
                    } else {
                        vertices.clear();
                        /*
                        for (int d = 0; d < depth_; d++) { // <= to add one more
                        CVsAround(storedU, storedV, numCVsInU, numCVsInV, UIsPeriodic, VIsPeriodic,
                        vertices);
                        }
                        */
                        CVsAround(storedU, storedV, numCVsInU_, numCVsInV_, UIsPeriodic_,
                                  VIsPeriodic_, vertices);
                        stat = getAverageWeight(vertices, index);
                    }
                }
            }
            currentWeights.copy(newWeights);
        }
    } else if (isMeshSurface_) {  // component.apiType() == MFn::kMeshVertComponent) {
                                  // 3 cycle through all vertices
        MItMeshVertex itVertex(meshPath_, component, &stat);
        if (stat == MS::kFailure) {
            MGlobal::displayError(
                MString(" MItMeshVertex itVertex(meshPath_, component, &stat); "));
            return MS::kFailure;
        }
        MIntArray vertices, repeatVertices, tmpVertices;
        MItMeshVertex itTempVertex(meshPath_);
        int prevIndex;
        // repeat the function
        for (int r = 0; r < repeat_; r++) {
            while (!itVertex.isDone()) {
                int currentVertex = itVertex.index();
                if (lockVertices_[currentVertex] != 1) {  // if not locked
                    if (verbose) MGlobal::displayInfo(MString(" vtx :") + currentVertex);
                    if (verbose) stat = printWeigth(currentVertex);

                    if (command_ == kCommandSmooth) {
                        // here depth of get vertices connected
                        itVertex.getConnectedVertices(vertices);
                        // check the depth
                        if (depth_ > 1) {
                            // init the array
                            std::unordered_set<int> setOfVerts;
                            for (unsigned int itVtx = 0; itVtx < vertices.length(); itVtx++) {
                                setOfVerts.insert(vertices[itVtx]);
                            }
                            // for the repeats
                            std::unordered_set<int> setOfVertsTmp;
                            for (int d = 1; d < depth_; d++) {  // <= to add one more
                                setOfVertsTmp.clear();
                                for (int vtx : setOfVerts) {
                                    itTempVertex.setIndex(vtx, prevIndex);
                                    itTempVertex.getConnectedVertices(repeatVertices);
                                    for (unsigned int itVtx = 0; itVtx < repeatVertices.length();
                                         itVtx++)
                                        setOfVertsTmp.insert(repeatVertices[itVtx]);
                                }
                                for (int vtx : setOfVertsTmp) setOfVerts.insert(vtx);
                            }
                            // now set the MIntArray
                            vertices.clear();
                            for (int vtx : setOfVerts) vertices.append(vtx);
                        }
                        stat = getAverageWeight(vertices, currentVertex);
                        if (stat == MS::kFailure) {
                            MGlobal::displayError(
                                MString("something is failing, select and try again"));
                            return MS::kFailure;
                        }
                    } else
                        addWeights(currentVertex);
                }
                itVertex.next();
            }
            currentWeights.copy(newWeights);
            itVertex.reset();
        }
    }
    // FINAL PART SETTING THE WEIGHTS --------------------------------------------
    // we make it short first
    weightsForSetting.clear();
    double oppositePercentMvt = 1.0 - percentMvt_;
    // for all the vertices
    int i = 0, currentVertex, start, end, index;
    double theNewWeight;

    float influence;
    int currentWeightsLength = currentWeights.length();
    for (int i = 0; i < indicesVertices_.length(); ++i) {
        currentVertex = indicesVertices_[i];
        influence = weightVertices_[i];
        start = currentVertex * nbJoints;
        end = start + nbJoints;
        for (index = start; index < end; index++) {
            // set percentage
            if (index >= currentWeightsLength) {
                if (verbose)
                    MGlobal::displayInfo(MString(" BREAK currentVertex :") + currentVertex +
                                         MString(" index  :") + index);
                break;
            }
            theNewWeight =
                currentWeights[index] * percentMvt_ + fullOrigWeights[index] * oppositePercentMvt;
            // set the influence value
            theNewWeight = theNewWeight * influence + fullOrigWeights[index] * (1. - influence);

            if (isMeshSurface_)
                weightsForSetting.append(theNewWeight);
            else if (isNurbsSurface_)
                newWeights.set(theNewWeight, index);
        }
    }

    // now set the values
    redoIt();
    return MS::kFailure;
}

MStatus blurSkinCmd::redoIt() {
    MStatus stat;

    MIntArray influenceIndices;
    for (int i = 0; i < nbJoints; i++) influenceIndices.append(i);
    MFnSkinCluster theSkinCluster(skinCluster_);
    if (isMeshSurface_)
        theSkinCluster.setWeights(meshPath_, component, influenceIndices, weightsForSetting, false,
                                  &weigthsForUndo);
    else {
        MDoubleArray wts;
        // MFnSingleIndexedComponent theVertex;
        int index, storedU, storedV;
        for (int i = 0; i < indicesVertices_.length(); ++i) {
            MFnDoubleIndexedComponent doubleFn;
            MObject tmpComponent = doubleFn.create(MFn::kSurfaceCVComponent);
            MDoubleArray tmpWeightsUndo;

            index = indicesVertices_[i];
            storedU = indicesU_[i];
            storedV = indicesV_[i];
            doubleFn.addElement(storedU, storedV);
            //------- build array -----------
            wts.clear();
            for (int j = 0; j < nbJoints; j++) {
                int posiToSet = index * nbJoints + j;
                wts.append(newWeights[posiToSet]);
            }
            stat = theSkinCluster.setWeights(meshPath_, tmpComponent, influenceIndices, wts, false,
                                             &tmpWeightsUndo);
        }
        /*
        MItGeometry gIter(meshPath_);
        MDoubleArray wts;
        MFnNurbsSurface MfnSurface(meshPath_);
        int index = 0;// indexU, indexV;
        int sizeInV = MfnSurface.numCVsInV();
        int resArrayInd;
        for (; !gIter.isDone(); gIter.next()) {
                MObject comp = gIter.component(&stat);
                resArrayInd = getMIntArrayIndex(indicesVertices_, index);
                if (resArrayInd != -1) {
                        MDoubleArray tmpWeightsUndo;
                        //------- build array -----------
                        wts.clear();
                        for (int j = 0; j<nbJoints; j++) {
                                int posiToSet = index*nbJoints + j;
                                wts.append(newWeights[posiToSet]);
                        }
                        stat = theSkinCluster.setWeights(meshPath_, comp, influenceIndices, wts,
        false, &tmpWeightsUndo);
                }
                index++;
        }
        */
    }
    return MS::kSuccess;
}

MStatus blurSkinCmd::undoIt() {
    MStatus stat;

    MIntArray influenceIndices;
    for (int i = 0; i < nbJoints; i++) influenceIndices.append(i);
    MFnSkinCluster theSkinCluster(skinCluster_);
    if (isMeshSurface_)
        theSkinCluster.setWeights(meshPath_, component, influenceIndices, weigthsForUndo, false,
                                  &weightsForSetting);
    else if (isNurbsSurface_) {
        MDoubleArray wts;
        // MFnSingleIndexedComponent theVertex;
        int index, storedU, storedV;
        for (int i = 0; i < indicesVertices_.length(); ++i) {
            MFnDoubleIndexedComponent doubleFn;
            MObject tmpComponent = doubleFn.create(MFn::kSurfaceCVComponent);
            MDoubleArray tmpWeightsUndo;

            index = indicesVertices_[i];
            storedU = indicesU_[i];
            storedV = indicesV_[i];
            doubleFn.addElement(storedU, storedV);
            //------- build array -----------
            wts.clear();
            for (int j = 0; j < nbJoints; j++) {
                int posiToSet = index * nbJoints + j;
                wts.append(fullOrigWeights[posiToSet]);
            }
            stat = theSkinCluster.setWeights(meshPath_, tmpComponent, influenceIndices, wts, false,
                                             &tmpWeightsUndo);
        }
        /*
        MItGeometry gIter(meshPath_);
        MDoubleArray wts;
        MFnNurbsSurface MfnSurface(meshPath_);
        int index = 0;// , indexU, indexV;
        int sizeInV = MfnSurface.numCVsInV();
        int resArrayInd;
        for (; !gIter.isDone(); gIter.next()) {
                MObject comp = gIter.component(&stat);
                resArrayInd = getMIntArrayIndex(indicesVertices_, index);
                if (resArrayInd != -1) {
                        MDoubleArray tmpWeightsUndo;
                        //------- build array -----------
                        wts.clear();
                        for (int j = 0; j<nbJoints; j++) {
                                int posiToSet = index*nbJoints + j;
                                wts.append(fullOrigWeights[posiToSet]);
                        }
                        stat = theSkinCluster.setWeights(meshPath_, comp, influenceIndices, wts,
        false, &tmpWeightsUndo);
                }
                index++;
        }
        */
    }
    return MS::kSuccess;
}

void* blurSkinCmd::creator() { return new blurSkinCmd(); }

blurSkinCmd::~blurSkinCmd() {
    // Note that we do nothing with fComponent which is owned by Maya.
    /*
            meshOrigin.~MObject ();
            symetriqueFaces.~MIntArray ();
            symetriqueVertices.~MIntArray ();
            symetriqueEdges.~MIntArray ();
            rightLeftVertices.~MIntArray ();
            rightLeftEdges.~MIntArray ();
            rightLeftFaces.~MIntArray ();
            unSymFaceArray.~MIntArray ();
            unSymEdgeArray.~MIntArray ();
    */
}

bool blurSkinCmd::isUndoable() const { return true; }
