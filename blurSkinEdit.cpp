#include "blurSkinEdit.h"

#include "functions.h"

MTypeId blurSkinDisplay::id(0x001226F9);

MObject blurSkinDisplay::_inMesh;
MObject blurSkinDisplay::_outMesh;
MObject blurSkinDisplay::_paintableAttr;
MObject blurSkinDisplay::_clearArray;
MObject blurSkinDisplay::_callUndo;
MObject blurSkinDisplay::_postSetting;
MObject blurSkinDisplay::_commandAttr;
MObject blurSkinDisplay::_smoothRepeat;
MObject blurSkinDisplay::_influenceAttr;
MObject blurSkinDisplay::_s_per_joint_weights;
MObject blurSkinDisplay::_s_skin_weights;

blurSkinDisplay::blurSkinDisplay() {}
blurSkinDisplay::~blurSkinDisplay() {}

MStatus blurSkinDisplay::getAttributes(MDataBlock& dataBlock) {
    MStatus status;

    if (verbose) MGlobal::displayInfo(" reloadCommand  ");
    MDataHandle commandData = dataBlock.inputValue(_commandAttr);
    MDataHandle influenceData = dataBlock.inputValue(_influenceAttr);
    MDataHandle smoothRepeatData = dataBlock.inputValue(_smoothRepeat);
    MDataHandle postSettingData = dataBlock.inputValue(_postSetting);

    this->influenceIndex = influenceData.asInt();
    this->commandIndex = commandData.asInt();
    this->smoothRepeat = smoothRepeatData.asInt();
    this->postSetting = postSettingData.asBool();

    if (verbose) MGlobal::displayInfo(MString(" commandeIndex  ") + this->commandIndex + " - ");
    if (verbose) MGlobal::displayInfo(MString(" influenceIndex ") + this->influenceIndex + " - ");
    if (verbose) MGlobal::displayInfo(MString(" smoothRepeat   ") + this->smoothRepeat + " - ");

    this->reloadCommand = false;
    this->applyPaint = false;

    return status;
}

MStatus blurSkinDisplay::compute(const MPlug& plug, MDataBlock& dataBlock) {
    MStatus status;

    if (plug.attribute() == blurSkinDisplay::_outMesh) {
        if (verbose) MGlobal::displayInfo(" _outMesh CALL ");  // beginning opening of node

        MDataHandle inMeshData = dataBlock.inputValue(blurSkinDisplay::_inMesh);
        MDataHandle outMeshData = dataBlock.outputValue(blurSkinDisplay::_outMesh);

        if (this->reloadCommand) getAttributes(dataBlock);

        if (skinCluster_ == MObject::kNullObj) {
            outMeshData.copy(inMeshData);  // copy the mesh

            getConnectedSkinCluster();                              // get the skinCluster
            getListColorsJoints(skinCluster_, this->jointsColors);  // get the joints colors
            status = fillArrayValues(true);   // get the skin data and all the colors
            set_skinning_weights(dataBlock);  // set the skin data and all the colors
            getAttributes(dataBlock);         // for too fast a start
        }

        if (skinCluster_ != MObject::kNullObj) {
            // 1 get the colors
            MObject outMesh = outMeshData.asMesh();
            MFnMesh meshFn(outMesh);
            int nbVertices = meshFn.numVertices();
            if (this->init) {  // init of node /////////////////
                this->init = false;
                this->paintedValues = MDoubleArray(nbVertices, 0);
                if (verbose) MGlobal::displayInfo(" set COLORS ");  // beginning opening of node
                MIntArray VertexCountPerPolygon, fullVvertexList;
                meshFn.getVertices(VertexCountPerPolygon, fullVvertexList);
                // now setting colors
                meshFn.createColorSetDataMesh(MString("paintColors"));
                // meshFn.createColorSetWithName("paintColorsSet");
                meshFn.setCurrentColorSetName(MString("paintColorsSet"));
                meshFn.setColors(this->currColors);
                meshFn.assignColors(fullVvertexList);
                // locks joints // all unlock
                this->lockJoints = MIntArray(nbVertices, 0);

                // get connected vertices --------------------
                getConnectedVertices(outMesh, nbVertices);
            } else if (this->applyPaint) {
                if (verbose) MGlobal::displayInfo("  -- > applyPaint  ");
                this->applyPaint = false;

                /////////////////////////////////////////////////////////////////
                // consider painted attribute --------
                /////////////////////////////////////////////////////////////////
                MColor multColor(1, 1, 1);
                float intensity = 1.0;
                // 0 Add - 1 Remove - 2 AddPercent - 3 Absolute - 4 Smooth - 5 Sharpen - 6 Colors

                if ((commandIndex < 4) && (influenceIndex > -1))
                    multColor = this->jointsColors[influenceIndex];
                else
                    intensity = .7;

                MColorArray theEditColors;
                MIntArray theEditVerts;
                MDoubleArray verticesWeight;

                if (this->clearTheArray) {
                    this->clearTheArray = false;
                    MDataHandle clearArrayData = dataBlock.inputValue(_clearArray);
                    bool clearArrayVal = clearArrayData.asBool();
                    if (verbose) {
                        MString strVal = "False";
                        if (clearArrayVal) strVal = "True";
                        MGlobal::displayInfo("  -- > clearArrayVal   " + strVal);
                    }
                    if (clearArrayVal) {
                        for (unsigned int i = 0; i < this->paintedValues.length(); i++) {
                            if (this->paintedValues[i] != 0) {
                                theEditColors.append(this->currColors[i]);
                                theEditVerts.append(i);
                                verticesWeight.append(this->paintedValues[i]);
                                this->paintedValues[i] = 0.0;
                            }
                        }
                        // post brushing apply values
                        // ----------------------------------------------------
                        if (this->postSetting)
                            applyCommand(dataBlock, theEditVerts, verticesWeight);
                        // refresh the colors with real values
                        // -------------------------------------------
                        refreshColors(theEditVerts, theEditColors);

                        MPlug clearArrayPlug(thisMObject(), _clearArray);
                        clearArrayPlug.setBool(false);
                    }
                } else if (this->callUndo) {
                    this->callUndo = false;
                    MDataHandle callUndoData = dataBlock.inputValue(_callUndo);
                    bool callUndoVal = callUndoData.asBool();
                    if (verbose) {
                        MString strVal = "False";
                        if (callUndoVal) strVal = "True";
                        MGlobal::displayInfo("  -- > CALL UNDO" + strVal);
                    }
                    if (callUndoVal) {                             // do the undo
                        if (this->undoVertsIndices_.size() > 0) {  // if stack is more than zero
                            std::vector<int> undoVerts = this->undoVertsIndices_.back();
                            std::vector<double> undoWeight = this->undoVertsValues_.back();

                            MDoubleArray undoWeight_MArr;  // (undoWeight.size());
                            MIntArray undoVerts_MArr;      // (undoWeight.size());
                            for (int vtx : undoVerts) undoVerts_MArr.append(vtx);
                            for (double val : undoWeight) undoWeight_MArr.append(val);

                            applyCommand(dataBlock, undoVerts_MArr, undoWeight_MArr,
                                         false);  // not storing undos - should prepare redos

                            this->undoVertsIndices_.pop_back();
                            this->undoVertsValues_.pop_back();
                            refreshColors(undoVerts_MArr, theEditColors);
                        } else {
                            MGlobal::displayInfo("  NO MORE UNDOS ");
                        }
                    }
                    // now set Attr false
                    MPlug callUndoPlug(thisMObject(), _callUndo);
                    callUndoPlug.setBool(false);
                } else {
                    // read paint values ---------------------------
                    MFnDoubleArrayData arrayData;
                    MObject dataObj = dataBlock.inputValue(_paintableAttr).data();
                    arrayData.setObject(dataObj);

                    unsigned int length = arrayData.length();
                    for (unsigned int i = 0; i < length; i++) {
                        double val = arrayData[i];
                        if (val > 0) {
                            if (val !=
                                this->paintedValues[i]) {  // only if other zone painted ----------
                                val = std::max(0.0, std::min(val, 1.0));  // clamp
                                theEditVerts.append(i);
                                verticesWeight.append(val);
                                // MGlobal::displayInfo(MString(" paint value ") + i + MString(" -
                                // ") + val);
                                this->paintedValues[i] = val;  // store to not repaint

                                val *= intensity;
                                theEditColors.append(val * multColor +
                                                     (1.0 - val) * this->currColors[i]);
                            }
                        }
                    }
                    // during brushing apply values
                    // ---------------------------------------------------
                    if (!this->postSetting) applyCommand(dataBlock, theEditVerts, verticesWeight);
                }
                // meshFn.setVertexColors(theEditColors, theEditVerts);
                meshFn.setSomeColors(theEditVerts, theEditColors);
                dataBlock.setClean(plug);
            }
        }
    }
    /*
    else if (plug.attribute() == blurSkinDisplay::_s_skin_weights) {
            //MGlobal::displayInfo(" _s_skin_weights CALL ");
    }*/
    dataBlock.setClean(plug);

    return status;
}

MStatus blurSkinDisplay::applyCommand(MDataBlock& dataBlock, MIntArray& theEditVerts,
                                      MDoubleArray& verticesWeight, bool storeUndo) {
    MStatus status;
    if (verbose) MGlobal::displayInfo(" applyCommand ");

    if ((commandIndex == 0) || (commandIndex == 4)) {
        // MDoubleArray previousWeights(this->nbJoints*theEditVerts.length(), 0.0);
        std::vector<double> previousWeights;
        std::vector<int> undoVerts;
        undoVerts.resize(theEditVerts.length());
        previousWeights.resize(this->nbJoints * theEditVerts.length());

        MDoubleArray theWeights(this->nbJoints * theEditVerts.length(), 0.0);
        int repeatLimit = 1;
        if (commandIndex == 4) repeatLimit = this->smoothRepeat;
        for (int repeat = 0; repeat < repeatLimit; ++repeat) {
            if (commandIndex == 4) {
                for (int i = 0; i < theEditVerts.length(); ++i) {
                    int theVert = theEditVerts[i];
                    double theVal = verticesWeight[i];

                    MIntArray vertsAround = this->connectedVertices[theVert];
                    status = setAverageWeight(vertsAround, theVert, i, this->nbJoints,
                                              this->lockJoints, this->skinWeightList, theWeights);
                }
            } else {
                status = editArray(commandIndex, influenceIndex, this->nbJoints,
                                   this->skinWeightList, theEditVerts, verticesWeight, theWeights);
            }
            // now set the weights
            for (int i = 0; i < theEditVerts.length(); ++i) {
                int theVert = theEditVerts[i];
                for (int j = 0; j < this->nbJoints; ++j) {
                    if (repeat == 0 && storeUndo)
                        previousWeights[i * this->nbJoints + j] =
                            this->skinWeightList[theVert * this->nbJoints + j];
                    this->skinWeightList[theVert * this->nbJoints + j] =
                        verticesWeight[i] * theWeights[i * this->nbJoints + j] +
                        (1.0 - verticesWeight[i]) *
                            this->skinWeightList[theVert * this->nbJoints + j];
                }
            }
        }
        if (storeUndo) {
            // MIntArray undoVerts;
            // undoVerts.copy(theEditVerts);
            // now store the undo ----------------
            for (int i = 0; i < theEditVerts.length(); ++i) undoVerts[i] = theEditVerts[i];
            this->undoVertsIndices_.push_back(undoVerts);
            this->undoVertsValues_.push_back(previousWeights);
        }
        replace_weights(dataBlock, theEditVerts, theWeights);
    }
    return status;
}

void blurSkinDisplay::getConnectedVertices(MObject& outMesh, int nbVertices) {
    MItMeshVertex vertexIter(outMesh);
    connectedVertices.resize(nbVertices);
    connectedFaces.resize(nbVertices);
    for (int vtxTmp = 0; !vertexIter.isDone(); vertexIter.next(), ++vtxTmp) {
        MIntArray surroundingVertices, surroundingFaces;
        vertexIter.getConnectedVertices(surroundingVertices);
        connectedVertices[vtxTmp] = surroundingVertices;

        vertexIter.getConnectedFaces(surroundingFaces);
        connectedFaces[vtxTmp] = surroundingFaces;
    }
}

void blurSkinDisplay::getConnectedSkinCluster() {
    MStatus status;
    MPlug outMeshPlug(thisMObject(), blurSkinDisplay::_outMesh);
    MPlugArray connections;
    outMeshPlug.connectedTo(connections, false, true);
    // get mesh connection from outMesh
    for (unsigned int i = 0; i < connections.length(); i++) {
        MObject destConn = connections[i].node();
        if (destConn.apiType() == MFn::kSkinClusterFilter) {
            skinCluster_ = destConn;
            break;
        }
        ///////////////////////////////////////////////////////////////////////////
        // now check the connection of the plug -----------------------------------
        // we do it in maya cmds better
        /*
        MPlug weight_list_plug(thisMObject(), blurSkinDisplay::_s_skin_weights);
        MPlugArray plugs;
        weight_list_plug.connectedTo(plugs, true, false, &status);
        if (plugs.length() == 0) { // not connected to the weightList
                MFnDependencyNode skinDep(skinCluster_);
                MPlug weight_list_skin_clus = skinDep.findPlug("weightList");
                MDGModifier dg;
                //status = dg.connect( weight_list_skin_clus, weight_list_plug);

                MGlobal::displayInfo(" TRY CONNECT " + weight_list_plug.name() + " -> " +
        weight_list_skin_clus.name()); status = dg.connect( weight_list_plug,
        weight_list_skin_clus); if (MS::kSuccess != status) MGlobal::displayError("dg.connect(
        weight_list_skin_clus, weight_list_plug);"); status = dg.doIt(); if (MS::kSuccess != status)
        MGlobal::displayError("dg.doIt ;");
        }
        */
    }
}

MStatus blurSkinDisplay::refreshColors(MIntArray& theVerts, MColorArray& theEditColors) {
    MStatus status = MS::kSuccess;
    if (verbose) MGlobal::displayInfo(" refreshColors CALL ");  // beginning opening of node
    if (theEditColors.length() != theVerts.length()) {
        theEditColors.setLength(theVerts.length());
    }
    for (int i = 0; i < theVerts.length(); ++i) {
        int theVert = theVerts[i];
        MColor theColor;
        for (int j = 0; j < this->nbJoints; ++j) {  // for each joint
            theColor += jointsColors[j] * this->skinWeightList[theVert * this->nbJoints + j];
        }
        this->currColors[theVert] = theColor;
        theEditColors[i] = theColor;
    }
    return status;
}

MStatus blurSkinDisplay::fillArrayValues(bool doColors) {
    MStatus status = MS::kSuccess;

    MFnDependencyNode skinClusterDep(skinCluster_);

    MPlug weight_list_plug = skinClusterDep.findPlug("weightList");
    MPlug matrix_plug = skinClusterDep.findPlug("matrix");
    // MGlobal::displayInfo(weight_list_plug.name());
    int nbElements = weight_list_plug.numElements();
    nbJoints = matrix_plug.numElements();

    skin_weights_.resize(nbElements);
    if (doColors) {
        currColors.clear();
        currColors.setLength(nbElements);
    }
    this->skinWeightList = MDoubleArray(nbElements * nbJoints, 0.0);

    for (int i = 0; i < nbElements; ++i) {
        // weightList[i]
        MPlug ith_weights_plug = weight_list_plug.elementByPhysicalIndex(i);
        int vertexIndex = ith_weights_plug.logicalIndex();
        // MGlobal::displayInfo(ith_weights_plug.name());

        // weightList[i].weight
        MPlug plug_weights = ith_weights_plug.child(0);  // access first compound child
        int nb_weights = plug_weights.numElements();
        skin_weights_[i].resize(nb_weights);
        // MGlobal::displayInfo(plug_weights.name() + nb_weights);

        MColor theColor;
        for (int j = 0; j < nb_weights; j++) {  // for each joint
            MPlug weight_plug = plug_weights.elementByPhysicalIndex(j);
            // weightList[i].weight[j]
            int indexInfluence = weight_plug.logicalIndex();
            double theWeight = weight_plug.asDouble();

            skin_weights_[i][j] = std::make_pair(indexInfluence, theWeight);
            this->skinWeightList[vertexIndex * nbJoints + indexInfluence] = theWeight;
            if (doColors) theColor += jointsColors[indexInfluence] * theWeight;
            /*
            //get Value
            auto myPair = skin_weights_[i][j];
            int index = myPair.first;
            float  val= myPair.second;
            */
        }
        if (doColors) currColors[i] = theColor;
    }
    if (verbose) MGlobal::displayInfo(" FILLED ARRAY VALUES ");

    return status;
}

void blurSkinDisplay::set_skinning_weights(MDataBlock& block) {
    MStatus status = MS::kSuccess;
    MArrayDataHandle array_hdl = block.outputArrayValue(_s_skin_weights, &status);
    MArrayDataBuilder array_builder = array_hdl.builder(&status);

    auto nbVerts = skin_weights_.size();
    // array_builder.growArray(nbVerts);
    for (int i = 0; i < nbVerts; i++) {
        auto vertexWeight = skin_weights_[i];
        auto nbInfluences = vertexWeight.size();

        MDataHandle element_hdl = array_builder.addElement(i, &status);  // weightList[i]
        MDataHandle child = element_hdl.child(_s_per_joint_weights);     // weightList[i].weight

        MArrayDataHandle weight_list_hdl(child, &status);
        MArrayDataBuilder weight_list_builder = weight_list_hdl.builder(&status);

        for (int j = 0; j < nbInfluences; ++j) {
            auto myPair = skin_weights_[i][j];
            int index = myPair.first;
            float val = myPair.second;

            MDataHandle hdl = weight_list_builder.addElement(index, &status);
            // hdl.setDouble((double)val);
            double theWeight = this->skinWeightList[i * nbJoints + index];
            hdl.setDouble(theWeight);
        }
        weight_list_hdl.set(weight_list_builder);
    }
    array_hdl.set(array_builder);
}

void blurSkinDisplay::replace_weights(MDataBlock& block, MIntArray& theVertices,
                                      MDoubleArray& theWeights) {
    MStatus status = MS::kSuccess;
    MArrayDataHandle array_hdl = block.outputArrayValue(_s_skin_weights, &status);
    for (int i = 0; i < theVertices.length(); ++i) {
        int indexVertex = theVertices[i];
        array_hdl.jumpToArrayElement(indexVertex);
        // weightList[i]
        MDataHandle element_hdl = array_hdl.outputValue(&status);
        // weightList[i].weight
        MDataHandle child = element_hdl.child(_s_per_joint_weights);

        MArrayDataHandle weight_list_hdl(child, &status);
        MArrayDataBuilder weight_list_builder = weight_list_hdl.builder(&status);

        unsigned handle_count = weight_list_hdl.elementCount(&status);
        unsigned builder_count = weight_list_builder.elementCount(&status);
        // mayaAssert(builder_count == handle_count);

        MIntArray to_remove;
        // Scan array, update existing element, remove unsused ones
        for (unsigned j = 0; j < handle_count; ++j) {
            // weightList[i].weight[j]
            // weight_list_hdl.jumpToArrayElement(j);
            unsigned index = weight_list_hdl.elementIndex(&status);

            double weight = theWeights[i * this->nbJoints + index];
            if (weight == 0.0)
                to_remove.append(index);
            else {
                MDataHandle hdl = weight_list_hdl.outputValue(&status);
                hdl.setDouble(weight);
                theWeights[i * this->nbJoints + index] = 0.0;
            }
            weight_list_hdl.next();
        }
        for (int k = 0; k < to_remove.length(); ++k)
            weight_list_builder.removeElement(to_remove[k]);
        // add the missing
        for (unsigned j = 0; j < this->nbJoints; ++j) {
            double weight = theWeights[i * this->nbJoints + j];
            if (weight != 0.0) {
                MDataHandle hdl = weight_list_builder.addElement(j, &status);
                hdl.setDouble(weight);
            }
        }
        weight_list_hdl.set(weight_list_builder);
    }
}

MPlug blurSkinDisplay::passThroughToOne(const MPlug& plug) const {
    if (plug.attribute() == blurSkinDisplay::_inMesh) {
        return MPlug(thisMObject(), blurSkinDisplay::_outMesh);
    }

    return MPlug();
}

void* blurSkinDisplay::creator() { return (new blurSkinDisplay()); }

MStatus blurSkinDisplay::initialize() {
    MStatus status;

    MFnTypedAttribute meshAttr;
    MFnTypedAttribute tAttr;
    MFnNumericAttribute numAtt;

    blurSkinDisplay::_inMesh = meshAttr.create("inMesh", "im", MFnMeshData::kMesh, &status);
    meshAttr.setStorable(false);
    meshAttr.setConnectable(true);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_inMesh);

    // mesh output
    blurSkinDisplay::_outMesh = meshAttr.create("outMesh", "om", MFnMeshData::kMesh, &status);
    meshAttr.setStorable(false);
    meshAttr.setConnectable(true);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_outMesh);

    blurSkinDisplay::_postSetting =
        numAtt.create("postSetting", "ps", MFnNumericData::kBoolean, true, &status);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_postSetting);
    ///////////////////////////////////////////////////////////////////////////
    // paintable attributes
    ///////////////////////////////////////////////////////////////////////////
    blurSkinDisplay::_paintableAttr =
        tAttr.create("paintAttr", "pa", MFnData::kDoubleArray, &status);
    meshAttr.setStorable(true);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_paintableAttr);

    blurSkinDisplay::_clearArray =
        numAtt.create("clearArray", "ca", MFnNumericData::kBoolean, false, &status);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_clearArray);

    blurSkinDisplay::_callUndo =
        numAtt.create("callUndo", "cu", MFnNumericData::kBoolean, false, &status);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_callUndo);

    ///////////////////////////////////////////////////////////////////////////
    // creation attributes
    ///////////////////////////////////////////////////////////////////////////
    MFnEnumAttribute enumAttr;

    blurSkinDisplay::_commandAttr = enumAttr.create("command", "cmd", 0);
    CHECK_MSTATUS(enumAttr.addField("Add", 0));
    CHECK_MSTATUS(enumAttr.addField("Remove", 1));
    CHECK_MSTATUS(enumAttr.addField("AddPercent", 2));
    CHECK_MSTATUS(enumAttr.addField("Absolute", 3));
    CHECK_MSTATUS(enumAttr.addField("Smooth", 4));
    CHECK_MSTATUS(enumAttr.addField("Sharpen", 5));
    CHECK_MSTATUS(enumAttr.addField("Colors", 6));
    CHECK_MSTATUS(enumAttr.setStorable(true));
    CHECK_MSTATUS(enumAttr.setKeyable(true));
    CHECK_MSTATUS(enumAttr.setReadable(true));
    CHECK_MSTATUS(enumAttr.setWritable(true));
    CHECK_MSTATUS(enumAttr.setCached(false));

    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_commandAttr);

    blurSkinDisplay::_smoothRepeat =
        numAtt.create("smoothRepeat", "sr", MFnNumericData::kInt64, 3, &status);
    numAtt.setMin(0);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_smoothRepeat);

    blurSkinDisplay::_influenceAttr =
        numAtt.create("influenceIndex", "ii", MFnNumericData::kInt64, 0, &status);
    numAtt.setMin(0);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_influenceAttr);
    ///////////////////////////////////////////////////////////////////////////
    // Initialize skin weights multi attributes
    ///////////////////////////////////////////////////////////////////////////
    _s_per_joint_weights = numAtt.create("weights", "w", MFnNumericData::kDouble, 0.0, &status);
    numAtt.setKeyable(false);
    numAtt.setArray(true);
    numAtt.setReadable(true);
    numAtt.setWritable(false);
    numAtt.setUsesArrayDataBuilder(true);
    addAttribute(_s_per_joint_weights);

    MFnCompoundAttribute cmpAttr;
    _s_skin_weights = cmpAttr.create("weightList", "wl", &status);
    cmpAttr.setArray(true);
    cmpAttr.addChild(_s_per_joint_weights);
    cmpAttr.setKeyable(false);
    cmpAttr.setReadable(true);
    cmpAttr.setWritable(false);
    cmpAttr.setUsesArrayDataBuilder(true);
    addAttribute(_s_skin_weights);

    cmpAttr.setStorable(true);  // To be stored during file-save

    ///////////////////////////////////////////////////////////////////////////
    // attributeAffects
    ///////////////////////////////////////////////////////////////////////////
    attributeAffects(blurSkinDisplay::_inMesh, blurSkinDisplay::_outMesh);
    // attributeAffects(blurSkinDisplay::_paintableAttr, blurSkinDisplay::_fakeAttr);
    attributeAffects(blurSkinDisplay::_paintableAttr, blurSkinDisplay::_s_skin_weights);
    // attributeAffects(blurSkinDisplay::_paintableAttr, blurSkinDisplay::_outMesh);
    return status;

    MGlobal::executeCommand("makePaintable -attrType doubleArray blurSkinDisplay paintAttr");
}

MStatus blurSkinDisplay::setDependentsDirty(const MPlug& plugBeingDirtied,
                                            MPlugArray& affectedPlugs) {
    MStatus status;
    MObject thisNode = thisMObject();
    MFnDependencyNode fnThisNode(thisNode);
    this->reloadCommand =
        ((plugBeingDirtied == _commandAttr) || (plugBeingDirtied == _influenceAttr) ||
         (plugBeingDirtied == _smoothRepeat) || (plugBeingDirtied == _postSetting));
    this->clearTheArray = (plugBeingDirtied == _clearArray);
    this->callUndo = (plugBeingDirtied == _callUndo);

    if ((plugBeingDirtied == _paintableAttr) || this->reloadCommand || this->clearTheArray ||
        this->callUndo) {
        this->applyPaint = true;
        MPlug outMeshPlug(thisNode, blurSkinDisplay::_outMesh);
        affectedPlugs.append(outMeshPlug);
    }
    return (MS::kSuccess);
}
