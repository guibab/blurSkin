
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

const char* blurSkinCmd::kListVerticesWeightFlagShort = "-lw";
const char* blurSkinCmd::kListVerticesWeightFlagLong = "-listVerticesWeight";

const char* blurSkinCmd::kRepeatFlagShort = "-rp";
const char* blurSkinCmd::kRepeatFlagLong = "-repeat";

const char* blurSkinCmd::kDepthFlagShort = "-d";
const char* blurSkinCmd::kDepthFlagLong = "-depth";

const char* blurSkinCmd::kRespectLocksFlagShort = "-rl";
const char* blurSkinCmd::kRespectLocksFlagLong = "-respectLocks";

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
    help += "-skinCluster (-skn):         String     Name of the skinCluster\n";
    help +=
        "-meshName (-mn):             String	 Name of the mesh if skincluster is not passed\n";
    help +=
        "                                        If -skn and -mn are not passed uses selection\n";
    help += "-skinClusterIndex (-si):     Int        Index of SkinCluster if -mn passed  \n";
    help += "                                        Default 0\n";
    help +=
        "-percentMvt (-pc):           Float      Between 0. and 1. percent of value set default "
        "1.\n";
    help += "-verbose (-v):               N/A        Verbose print\n";
    help += "-listVerticesIndices (-li)   Strings    List vertices indices\n";
    help += "-listVerticesWeight (-lw)    Doubles    List vertices weights\n";
    help += "-listJoints         (-lj)    String     List joints names\n";
    help += "-listJointsValues   (-ljw)   Doubles    List joints weights\n";
    help +=
        "-repeat (-rp)                Int        Repeat the calculation             default 1\n";
    help +=
        "-depth (-d):                 Int        Depth for the smooth               default 1\n";
    help +=
        "-respectLocks (-rl):         N/A        Respect locks                      default True "
        "\n";
    help += "-command (-c):               N/A        the command action correct inputs are :\n";
    help += "                                        smooth - add - absolute - percentage\n";
    help += "-help (-h)                   N/A        Display this text.\n";
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
      respectLocks_(true),
      verbose(false) {}

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
    }

    // get basic arguments ----------------------------------------------------------------------
    if (argData.isFlagSet(kVerboseFlagShort))
        verbose = argData.flagArgumentBool(kVerboseFlagShort, 0, &status);

    if (argData.isFlagSet(kRepeatFlagShort))
        repeat_ = argData.flagArgumentInt(kRepeatFlagShort, 0, &status);

    if (argData.isFlagSet(kPercentMovementFlagShort))
        percentMvt_ = argData.flagArgumentDouble(kPercentMovementFlagShort, 0, &status);

    if (argData.isFlagSet(kDepthFlagShort))
        depth_ = argData.flagArgumentInt(kDepthFlagShort, 0, &status);

    if (argData.isFlagSet(kRespectLocksFlagShort))
        respectLocks_ = argData.flagArgumentBool(kRespectLocksFlagShort, 0, &status);

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
    if (argData.isFlagSet(kListVerticesIndicesFlagShort)) {
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
        for (unsigned int i = 0; i < indicesVertices_.length(); i++) weightVertices_.append(1.0);
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
MStatus blurSkinCmd::getListLockJoints() {
    if (verbose) MGlobal::displayInfo(MString(" ---- get list lock joints ----"));
    MStatus stat;
    MDagPathArray listOfJoints;
    MFnSkinCluster theSkinCluster(skinCluster_);
    theSkinCluster.influenceObjects(listOfJoints, &stat);
    nbJoints = listOfJoints.length();
    if (verbose) MGlobal::displayInfo(MString("    nbJoints : ") + nbJoints);

    // preSet the operation per joint at zero
    for (int i = 0; i < nbJoints; i++) perJointAddingValues_.append(0.0);
    for (int i = 0; i < nbJoints; i++) {
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
        if (jntIndex != -1) {  // if the joint is in the list of joints name
            perJointAddingValues_.set(listJointsValues_[jntIndex], i);
        }

        if (verbose) {
            if (isLockInfluenceWeights)
                MGlobal::displayInfo("    " + jnt.name() + " is Locked ");
            else
                MGlobal::displayInfo("    " + jnt.name() + " is not Locked ");
        }
    }

    return MS::kSuccess;
}

MStatus blurSkinCmd::printWeigth(int vertex) {
    if (verbose) MGlobal::displayInfo(MString(" ---- printWeigth ----"));
    MFnSkinCluster theSkinCluster(skinCluster_);

    // 3 get the weights
    MDoubleArray weightsVertex;
    unsigned int infCount;

    MFnSingleIndexedComponent theVertex;
    MObject tmpComponent = theVertex.create(MFn::kMeshVertComponent);
    theVertex.addElement(vertex);

    theSkinCluster.getWeights(meshPath_, tmpComponent, weightsVertex, infCount);
    int j;
    MString toDisplay("weigth of vtx (");
    toDisplay += vertex;
    toDisplay += MString(") : ");
    for (j = 0; j < nbJoints; j++) {
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

    MDoubleArray sumWeigths;
    for (j = 0; j < nbJoints; j++) sumWeigths.append(0.0);

    // std::cout << " allWeigths NEW\t";
    for (i = 0; i < sizeVertices; i++) {
        // std::cout << "\nindexVTX: " <<vertices[i]<<" | ";
        for (j = 0; j < nbJoints; j++) {
            posi = vertices[i] * nbJoints + j;
            // std::cout << currentWeights[posi] <<" - ";
            sumWeigths[j] += currentWeights[posi];
        }
    }
    double total = 0.0;
    double totalBaseVtx = 0.0;
    int assign = 0;
    for (j = 0; j < nbJoints; j++) {
        int posi = currentVertex * nbJoints + j;
        sumWeigths[j] /= sizeVertices;

        if (!lockJoints[j]) {
            total += sumWeigths[j];
            totalBaseVtx += currentWeights[posi];
            assign += 1;
        }
    }
    if (total > 0. && totalBaseVtx > 0.) {
        double mult = totalBaseVtx / total;

        // std::cout << "\n total :\t" <<total ;
        // std::cout << "\n vtx :" << currentVertex << "\n";
        for (j = 0; j < nbJoints; j++) {
            int posiToSet = currentVertex * nbJoints + j;
            if (!lockJoints[j]) {
                sumWeigths[j] *= mult;  // normalement divide par 1, sauf cas lock joints
                                        // std::cout << " " << sumWeigths[j];
                newWeights.set(sumWeigths[j], posiToSet);
                // newWeights[posiToSet]=sumWeigths[j];
            }
        }
    }
    // print of result computation
    verboseSetWeights(currentVertex);
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
                } else if (command_ == kCommandPercentage) {
                    totalOfAdding += value * (1 + perJointAddingValues_[j]);
                }
            } else {
                otherJointsTotal += value;
            }
        }
    }
    // now do the setting of the values ---------------------------
    // ------------------------------------------------------------
    double newVal;
    if (totalOfAdding >= totalWithLocks) {
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

MStatus blurSkinCmd::getSoftSelection() {
    if (verbose) MGlobal::displayInfo(MString("---- getSoftSelection ----"));
    MStatus stat;

    MRichSelection richSel;
    MGlobal::getRichSelection(richSel);
    MSelectionList richSelList;
    richSel.getSelection(richSelList);

    for (MItSelectionList iter(richSelList, MFn::kInvalid); !iter.isDone(); iter.next()) {
        // MDagPath meshPath_;
        iter.getDagPath(meshPath_, component);
        if (verbose) MGlobal::displayInfo(" softSelection on mesh   " + meshPath_.fullPathName());

        if (component.isNull() || (component.apiType() != MFn::kMeshVertComponent)) {
            // do on all vertices
            meshPath_.extendToShape();

            MFnMesh meshFn(meshPath_, &stat);  // this is the visible mesh
            MIntArray ObjVertices;
            int nbVertices = meshFn.numVertices(&stat);
            for (int i = 0; i < nbVertices; i++) ObjVertices.append(i);
            MFnSingleIndexedComponent allVertices;
            component = allVertices.create(MFn::kMeshVertComponent);
            allVertices.addElements(ObjVertices);
        }
    }
    return MS::kSuccess;
}

MStatus blurSkinCmd::executeAction() {
    if (verbose) MGlobal::displayInfo(MString(" ---- executeAction ----"));
    MStatus stat;
    // 1 first get list locked joints
    getListLockJoints();
    // 2 get all weights of vertices
    getAllWeights();
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
    if (!component.isNull() && (component.apiType() == MFn::kMeshVertComponent)) {
        // 3 cycle through all vertices
        MItMeshVertex itVertex(meshPath_, component, &stat);
        MIntArray vertices, repeatVertices, tmpVertices;
        MItMeshVertex itTempVertex(meshPath_);
        int prevIndex;
        // repeat the function
        for (int r = 0; r < repeat_; r++) {
            while (!itVertex.isDone()) {
                int currentVertex = itVertex.index();
                if (verbose) MGlobal::displayInfo(MString(" vtx :") + currentVertex);
                if (verbose) stat = printWeigth(currentVertex);

                if (command_ == kCommandSmooth) {
                    // here depth of get vertices connected
                    itVertex.getConnectedVertices(vertices);
                    // check the depth
                    if (depth_ > 1) {
                        // init the array
                        std::unordered_set<int> setOfVerts;
                        for (unsigned int itVtx = 0; itVtx < vertices.length(); itVtx++)
                            setOfVerts.insert(vertices[itVtx]);
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
                } else
                    addWeights(currentVertex);

                itVertex.next();
            }
            currentWeights.copy(newWeights);
            itVertex.reset();
        }
        // FINAL PART SETTING THE WEIGHTS --------------------------------------------
        // we make it short first
        newWeights.clear();
        double oppositePercentMvt = 1.0 - percentMvt_;
        // for all the vertices
        int i = 0, currentVertex, start, end, index;
        double theNewWeight;

        // prepare the array
        MFnComponent componentFn(component);
        int count = componentFn.elementCount();
        MStatus stat;
        MFnSingleIndexedComponent singleFn(component, &stat);
        if (MS::kSuccess == stat) {
            for (int i = 0; i < count; i++) {
                currentVertex = singleFn.element(i);
                float influence = componentFn.weight(i).influence();

                start = currentVertex * nbJoints;
                end = start + nbJoints;
                for (index = start; index < end; index++) {
                    // set percentage
                    theNewWeight = currentWeights[index] * percentMvt_ +
                                   fullOrigWeights[index] * oppositePercentMvt;
                    // set the influence value
                    theNewWeight =
                        theNewWeight * influence + fullOrigWeights[index] * (1. - influence);
                    newWeights.append(theNewWeight);
                }
            }
        }
        // now set the values
        redoIt();
    }
    return MS::kFailure;
}

MStatus blurSkinCmd::getAllWeights() {
    if (verbose) MGlobal::displayInfo(MString(" ---- getAllWeights ----"));
    MFnSkinCluster theSkinCluster(skinCluster_);

    MStatus stat;
    MFnMesh meshFn(meshPath_, &stat);  // this is the visible mesh
    MIntArray ObjVertices;
    if (verbose) MGlobal::displayInfo(MString("    mesh is") + meshPath_.fullPathName());

    int nbVertices = meshFn.numVertices(&stat);
    for (int i = 0; i < nbVertices; i++) ObjVertices.append(i);

    MFnSingleIndexedComponent allVertices;
    allVertices.addElements(ObjVertices);
    MObject allVerticesObj = allVertices.create(MFn::kMeshVertComponent);

    unsigned int infCount;
    theSkinCluster.getWeights(meshPath_, allVerticesObj, fullOrigWeights, infCount);
    if (verbose)
        MGlobal::displayInfo(MString("    full weights nb weights ") + fullOrigWeights.length());

    currentWeights.copy(fullOrigWeights);
    newWeights.copy(fullOrigWeights);
    return MS::kSuccess;
}

MStatus blurSkinCmd::doIt(const MArgList& args) {
    MStatus status;
    MString info;
    info += "\n";

    status = GatherCommandArguments(args);
    CHECK_MSTATUS_AND_RETURN_IT(status);
    if (verbose) {
        info += MString("\n    percentMvt : ") + percentMvt_;
        info += MString("\n    repeat : ") + repeat_;
        info += MString("\n    depth : ") + depth_;
        info += MString("\n    verbose :") + verbose;
        info += MString("\n    respectLocks_ :") + respectLocks_;
        info += MString("\n    indSkinCluster :") + indSkinCluster_ + MString("\n");
        MGlobal::displayInfo(info);
    }
    if (command_ == kCommandHelp) {
        DisplayHelp();
        return MS::kSuccess;
    }

    // create the component
    if (useSelection) {
        getSoftSelection();
        status = findSkinCluster(meshPath_, skinCluster_, indSkinCluster_, verbose);
    } else {
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
    }
    executeAction();

    return MS::kSuccess;
}
MStatus blurSkinCmd::redoIt() {
    MStatus status;

    MIntArray influenceIndices;
    for (int i = 0; i < nbJoints; i++) influenceIndices.append(i);
    MFnSkinCluster theSkinCluster(skinCluster_);
    theSkinCluster.setWeights(meshPath_, component, influenceIndices, newWeights, false,
                              &weigthsForUndo);
    return MS::kSuccess;
}

MStatus blurSkinCmd::undoIt() {
    MStatus status;

    MIntArray influenceIndices;
    for (int i = 0; i < nbJoints; i++) influenceIndices.append(i);
    MFnSkinCluster theSkinCluster(skinCluster_);
    theSkinCluster.setWeights(meshPath_, component, influenceIndices, weigthsForUndo, false,
                              &newWeights);
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
