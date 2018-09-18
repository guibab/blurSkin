#include "blurSkinEdit.h"

#include "functions.h"

MTypeId blurSkinDisplay::id(0x001226F9);

MObject blurSkinDisplay::_inMesh;
MObject blurSkinDisplay::_outMesh;
MObject blurSkinDisplay::_cpList;
MObject blurSkinDisplay::_paintableAttr;
MObject blurSkinDisplay::_clearArray;
MObject blurSkinDisplay::_callUndo;
MObject blurSkinDisplay::_postSetting;
MObject blurSkinDisplay::_commandAttr;
MObject blurSkinDisplay::_colorType;
MObject blurSkinDisplay::_soloColorType;
MObject blurSkinDisplay::_smoothRepeat;
MObject blurSkinDisplay::_smoothDepth;
MObject blurSkinDisplay::_influenceAttr;
MObject blurSkinDisplay::_influenceColor;
MObject blurSkinDisplay::_getLockWeights;
MObject blurSkinDisplay::_autoExpandAttr;
MObject blurSkinDisplay::_normalize;
MObject blurSkinDisplay::_s_per_joint_weights;
MObject blurSkinDisplay::_s_skin_weights;

blurSkinDisplay::blurSkinDisplay() {}
blurSkinDisplay::~blurSkinDisplay() {}

MStatus blurSkinDisplay::getAttributes(MDataBlock& dataBlock) {
    MStatus status;

    if (verbose) MGlobal::displayInfo(MString(" GA| reloadCommand  "));
    MDataHandle commandData = dataBlock.inputValue(_commandAttr);
    MDataHandle influenceData = dataBlock.inputValue(_influenceAttr);
    MDataHandle smoothRepeatData = dataBlock.inputValue(_smoothRepeat);
    MDataHandle smoothDepthData = dataBlock.inputValue(_smoothDepth);
    MDataHandle postSettingData = dataBlock.inputValue(_postSetting);
    MDataHandle autoExpandData = dataBlock.inputValue(_autoExpandAttr);
    MDataHandle getLockWeightsData = dataBlock.inputValue(_getLockWeights);
    MDataHandle getNormalizeData = dataBlock.inputValue(_normalize);

    MDataHandle colorTypeData = dataBlock.inputValue(_colorType);
    int prevColorCommand = this->colorCommand;
    this->colorCommand = colorTypeData.asShort();
    if (verbose) MGlobal::displayInfo(MString(" GA| colorCommand   ") + this->colorCommand + " - ");

    MDataHandle soloColorTypeData = dataBlock.inputValue(_soloColorType);
    int prevSoloColorType = this->soloColorTypeVal;
    this->soloColorTypeVal = soloColorTypeData.asShort();
    if (verbose)
        MGlobal::displayInfo(MString(" GA| soloColorType  ") + this->soloColorTypeVal + " - ");

    int prevInfluenceIndex = this->influenceIndex;
    bool prevAutoExpand = this->autoExpand;
    bool prevLockWeights = this->refreshLockWeights;

    this->refreshLockWeights = getLockWeightsData.asBool();
    this->doNormalize = getNormalizeData.asBool();
    this->influenceIndex = influenceData.asInt();
    this->commandIndex = commandData.asShort();
    this->smoothRepeat = smoothRepeatData.asInt();

    this->autoExpand = autoExpandData.asBool();
    int smoothDepthVal = smoothDepthData.asInt();

    if (smoothDepthVal != this->smoothDepth) {
        this->smoothDepth = smoothDepthVal;
        refreshVertsConnection();
    }
    this->postSetting = postSettingData.asBool();

    this->reloadCommand = false;
    this->applyPaint = false;
    this->reloadSoloColor = false;

    if (this->inputVerticesChanged) {  // the vertices are changed, time to reconnect to the
                                       // weightList of the skinCluster
        // this->inputVerticesChanged = false; //do in compute
        if (verbose) MGlobal::displayInfo(MString(" GA| enter inputVerticesChanged "));

        MDataHandle inputIDs = dataBlock.inputValue(_cpList, &status);
        MObject compList = inputIDs.data();
        MFnComponentListData compListFn(compList);
        unsigned i;
        int j;
        // MIntArray cpIds;
        this->cpIds.clear();

        MFn::Type componentType = MFn::kMeshVertComponent;
        for (i = 0; i < compListFn.length(); i++) {
            MObject comp = compListFn[i];
            if (comp.apiType() == componentType) {
                MFnSingleIndexedComponent siComp(comp);
                for (j = 0; j < siComp.elementCount(); j++) this->cpIds.append(siComp.element(j));
            }
        }
        if (this->cpIds.length() > 0) {
            // get the weights from the plug --------
            if (verbose) MGlobal::displayInfo(MString(" GA| refresh weight skinCluster "));
            MDoubleArray theWeights(this->cpIds.length() * this->nbJoints, 0.0);
            querySkinClusterValues(this->cpIds, theWeights, false);  // with colors
            // need to refreshColors refreshColors(cpIds, multiEditColors, soloEditColors);

            // set the weights in the datablock------------
            replace_weights(dataBlock, this->cpIds, theWeights);
            // reconnect ------------ done in maya
            // connectSkinClusterWL();
            // this->doConnectSkinCL = true;
            // set at zero -------- done in maya
            this->applyPaint = true;
        } else {
            this->inputVerticesChanged = false;
        }
    }
    if (verbose) MGlobal::displayInfo(MString(" GA| commandeIndex  ") + this->commandIndex + " - ");
    if (verbose)
        MGlobal::displayInfo(MString(" GA| influenceIndex ") + this->influenceIndex + " - ");
    if (verbose) MGlobal::displayInfo(MString(" GA| smoothRepeat   ") + this->smoothRepeat + " - ");
    if (verbose) MGlobal::displayInfo(MString(" GA| smoothDepth    ") + this->smoothDepth + " - ");
    if (verbose)
        MGlobal::displayInfo(MString(" GA| inputVerticesChanged  ") + this->inputVerticesChanged +
                             " - ");

    if (this->colorCommand == 1)
        this->reloadSoloColor = (prevSoloColorType != this->soloColorTypeVal) ||
                                (prevColorCommand != this->colorCommand) ||
                                (prevInfluenceIndex != this->influenceIndex);
    /*

    if (prevColorCommand != this->colorCommand) {
            MFnSkinCluster theSkinCluster(this->skinCluster_);
            MObjectArray objectsDeformed;
            theSkinCluster.getOutputGeometry(objectsDeformed);
            MFnDependencyNode deformedNameMesh(objectsDeformed[0]);
            //setAttr - type "string" Mesh_X_HeadBody_Pc_Sd1_SdDsp_Shape.currentColorSet
    "soloColorsSet"; MPlug currentColorSet = deformedNameMesh.findPlug("currentColorSet");

            if (this->colorCommand == 0)
                    currentColorSet.setValue(this->fullColorSet);
            else
                    currentColorSet.setValue(this->soloColorSet);
    }
    */

    if (verbose)
        MGlobal::displayInfo(MString(" GA| reloadSoloColor ") + this->reloadSoloColor + " - ");

    return status;
}

MStatus blurSkinDisplay::compute(const MPlug& plug, MDataBlock& dataBlock) {
    MStatus status;

    if (plug.attribute() == blurSkinDisplay::_outMesh) {
        if (verbose)
            MGlobal::displayInfo(MString("--> _outMesh CALL "));  // beginning opening of node

        MDataHandle inMeshData = dataBlock.inputValue(blurSkinDisplay::_inMesh);
        MDataHandle outMeshData = dataBlock.outputValue(blurSkinDisplay::_outMesh);

        if (this->reloadCommand) {
            if (verbose) MGlobal::displayInfo(MString("   --> GA : this->reloadCommand "));
            getAttributes(dataBlock);
        }
        if (this->changedColorInfluence != -1) {  // user changed the color of the influence
            if (verbose) MGlobal::displayInfo(MString("   --> this->changedColorInfluence "));
            if (!this->init) {
                if (verbose) MGlobal::displayInfo(MString("      --> !this->init "));
                // MGlobal::displayInfo(MString(" changedColor ") + this->changedColorInfluence);

                MArrayDataHandle influenceColorHandle =
                    dataBlock.inputValue(blurSkinDisplay::_influenceColor);
                influenceColorHandle.jumpToArrayElement(this->changedColorInfluence);
                MDataHandle theHandle = influenceColorHandle.inputValue(&status);
                float3& colorvalue = theHandle.asFloat3();
                this->jointsColors[this->changedColorInfluence] =
                    MColor(colorvalue[0], colorvalue[1], colorvalue[2]);
            } else {
                if (verbose) MGlobal::displayInfo(MString("      --> this->init "));
                this->changedColorInfluence = -1;
                dataBlock.setClean(plug);
                return status;
            }
        }
        if (skinCluster_ == MObject::kNullObj) {
            if (verbose) MGlobal::displayInfo(MString("   --> skinCluster_ == MObject::kNullObj "));
            outMeshData.copy(inMeshData);  // copy the mesh

            getConnectedSkinCluster();                              // get the skinCluster
            getListColorsJoints(skinCluster_, this->jointsColors);  // get the joints colors
            getListLockVertices(skinCluster_, this->lockVertices);

            setInfluenceColorAttr();          // set the colors in our attribute
            status = fillArrayValues(true);   // get the skin data and all the colors
            set_skinning_weights(dataBlock);  // set the skin data and all the colors
            getAttributes(dataBlock);         // for too fast a start
        }
        if (skinCluster_ != MObject::kNullObj) {
            if (verbose) MGlobal::displayInfo(MString("   --> skinCluster_  is GOOD "));
            // 1 get the colors
            MObject outMesh = outMeshData.asMesh();
            MFnMesh meshFn(outMesh);
            int nbVertices = meshFn.numVertices();
            if (this->init) {  // init of node /////////////////
                if (verbose) MGlobal::displayInfo(MString("      --> INIT"));
                this->fullVvertexList.setLength(nbVertices);
                for (int i = 0; i < nbVertices; ++i) this->fullVvertexList[i] = i;
                this->init = false;
                this->paintedValues = MDoubleArray(nbVertices, 0);
                if (verbose)
                    MGlobal::displayInfo(
                        MString("          set COLORS "));  // beginning opening of node

                // get connected vertices --------------------
                getConnectedVertices(outMesh, nbVertices);
                refreshVertsConnection();

                // get face collor assignments ----------
                MIntArray VertexCountPerPolygon, fullVvertexList;
                meshFn.getVertices(VertexCountPerPolygon, fullVvertexList);
                this->fullVertexListLength = fullVvertexList.length();
                // now setting colors		 ---
                meshFn.createColorSetDataMesh(this->fullColorSet);
                meshFn.setCurrentColorSetName(this->fullColorSet);
                meshFn.setColors(this->multiCurrentColors, &this->fullColorSet);
                meshFn.assignColors(fullVvertexList, &this->fullColorSet);
                // solo colors -----------------------
                meshFn.createColorSetDataMesh(this->soloColorSet);
                this->soloCurrentColors = MColorArray(nbVertices, MColor(0, 0, 0));
                this->soloColorsValues = MDoubleArray(nbVertices, 0.0);

                meshFn.setColors(this->soloCurrentColors, &this->soloColorSet);
                meshFn.assignColors(fullVvertexList, &this->soloColorSet);

                meshFn.createColorSetDataMesh(this->noColorSet);  // for no colors
                // locks joints // all unlock
                this->lockJoints = MIntArray(nbVertices, 0);
                getListLockJoints(skinCluster_, this->lockJoints);

            } else if (this->refreshLockWeights) {  // get list of locks
                if (verbose) MGlobal::displayInfo(MString("      --> refreshLockWeights"));
                MDataHandle callRefreshLocksData = dataBlock.inputValue(_getLockWeights);
                bool callRefreshLocksVal = callRefreshLocksData.asBool();

                if (callRefreshLocksVal) {
                    getListLockJoints(skinCluster_, this->lockJoints);
                    getListLockVertices(skinCluster_, this->lockVertices);
                    MPlug callLockWeightsPlug(thisMObject(), _getLockWeights);
                    callLockWeightsPlug.setBool(false);
                }
                this->refreshLockWeights = false;
                dataBlock.setClean(plug);
                return status;
            } else if (this->reloadSoloColor) {
                if (verbose) MGlobal::displayInfo(MString("      --> reloadSoloColor  "));
                // if (this->soloColorsValues.length() > 20266)MGlobal::displayInfo(MString("
                // previous value [20266] -  ") + this->soloColorsValues[20266]);
                editSoloColorSet(meshFn);
                // if (this->colorCommand == 0) meshFn.setCurrentColorSetName(this->fullColorSet);
                // else meshFn.setCurrentColorSetName(this->soloColorSet);
                this->reloadSoloColor = false;
                dataBlock.setClean(plug);
                return status;

                // if (this->soloColorsValues .length() > 20266) MGlobal::displayInfo(MString(" post
                // value     [20266] -  ") + this->soloColorsValues[20266]);
            } else if (this->applyPaint) {
                if (verbose) MGlobal::displayInfo(MString("      --> applyPaint  "));
                this->applyPaint = false;

                /////////////////////////////////////////////////////////////////
                // consider painted attribute --------
                /////////////////////////////////////////////////////////////////
                MColor white(1, 1, 1), multColor, soloMultColor;
                float intensity = 1.0;
                // 0 Add - 1 Remove - 2 AddPercent - 3 Absolute - 4 Smooth - 5 Sharpen - 6 Colors
                if ((commandIndex < 4) && (influenceIndex > -1)) {
                    // multColor = .7*this->jointsColors[influenceIndex] + .3*white; // paint a
                    // little whiter
                    multColor = this->jointsColors[influenceIndex];
                    if (this->soloColorTypeVal == 2)
                        soloMultColor = this->jointsColors[influenceIndex];
                    else
                        soloMultColor = white;
                } else {
                    multColor = white;
                    soloMultColor = white;
                    intensity = float(0.1);
                }
                MColorArray multiEditColors, soloEditColors;
                MIntArray editVertsIndices;
                MDoubleArray editVertsWeights;

                if (this->clearTheArray) {
                    if (verbose) MGlobal::displayInfo("         --> this->clearTheArray");
                    this->clearTheArray = false;
                    MDataHandle clearArrayData = dataBlock.inputValue(_clearArray);
                    bool clearArrayVal = clearArrayData.asBool();
                    if (verbose) {
                        MString strVal = "False";
                        if (clearArrayVal) strVal = "True";
                        MGlobal::displayInfo("         --> clearArrayVal   " + strVal);
                    }
                    if (clearArrayVal) {
                        // MGlobal::displayInfo("------------  do clear array-----------------");

                        if (this->autoExpand && this->postSetting) {  // not asking refresh of skin
                            if (verbose)
                                MGlobal::displayInfo("            --> auto Expand on postSetting ");
                            for (int k = 0; k < this->nbAutoExpand; ++k) {
                                double threshold = .6;
                                std::unordered_set<int> toFixVtx;
                                for (unsigned int i = 0; i < this->paintedValues.length(); i++) {
                                    if (this->paintedValues[i] > threshold) {
                                        MIntArray vertsAround = this->connectedVertices[i];
                                        for (int j = 0; j < vertsAround.length(); ++j) {
                                            int theVertAround = vertsAround[j];
                                            if (this->paintedValues[theVertAround] <= 0.05) {
                                                toFixVtx.insert(theVertAround);
                                            }
                                        }
                                    }
                                }
                                MIntArray vtxToFix;
                                for (int vtx : toFixVtx) vtxToFix.append(vtx);
                                for (int i = 0; i < vtxToFix.length(); ++i) {
                                    int vtx = vtxToFix[i];
                                    MIntArray vertsAround = this->connectedVertices[vtx];
                                    double newVal = 0.0;
                                    int nbAround = 0;
                                    for (int j = 0; j < vertsAround.length(); ++j) {
                                        int theVertAround = vertsAround[j];
                                        if (toFixVtx.find(theVertAround) == toFixVtx.end()) {
                                            nbAround++;
                                            newVal += this->paintedValues[theVertAround];
                                        }
                                    }
                                    newVal /= nbAround;
                                    this->paintedValues[vtx] = newVal;
                                }
                            }
                        }
                        // now set
                        for (unsigned int i = 0; i < this->paintedValues.length(); i++) {
                            if (this->paintedValues[i] != 0 &&
                                !this->lockVertices[i]) {  // not zero and not locked ----
                                // multiEditColors.append(this->multiCurrentColors[i]);
                                editVertsIndices.append(i);
                                editVertsWeights.append(this->paintedValues[i]);
                                this->paintedValues[i] = 0.0;
                            }
                        }
                        // post brushing apply values
                        // ----------------------------------------------------
                        if (this->postSetting)
                            applyCommand(dataBlock, editVertsIndices, editVertsWeights);
                        else {  // store undo
                            MDoubleArray previousWeights;

                            previousWeights.setLength(editVertsIndices.length() * this->nbJoints);
                            for (int i = 0; i < editVertsIndices.length(); ++i) {
                                int theVert = editVertsIndices[i];
                                for (int j = 0; j < this->nbJoints; ++j)
                                    previousWeights[i * this->nbJoints + j] =
                                        this->fullUndoSkinWeightList[theVert * this->nbJoints + j];
                            }
                            this->undoVertsIndices_.push_back(editVertsIndices);
                            this->undoVertsValues_.push_back(previousWeights);
                            this->postSetting_timeToStoreUndo = true;
                            this->fullUndoSkinWeightList.clear();
                        }
                        // refresh the colors with real values
                        // -------------------------------------------
                        refreshColors(editVertsIndices, multiEditColors, soloEditColors);

                        MPlug clearArrayPlug(thisMObject(), _clearArray);
                        clearArrayPlug.setBool(false);
                    }
                } else if (this->inputVerticesChanged) {
                    if (verbose)
                        MGlobal::displayInfo(MString("         --> this->inputVerticesChanged ") +
                                             this->cpIds.length());
                    this->inputVerticesChanged = false;
                    editVertsIndices.copy(this->cpIds);
                    this->cpIds.clear();
                    refreshColors(editVertsIndices, multiEditColors, soloEditColors);
                } else if (this->callUndo) {
                    if (verbose) MGlobal::displayInfo("         -- > this->callUndo");
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
                            MDoubleArray undoWeight_MArr = this->undoVertsValues_.back();
                            editVertsIndices.copy(this->undoVertsIndices_.back());

                            for (int i = 0; i < editVertsIndices.length(); ++i) {
                                int theVert = editVertsIndices[i];
                                for (int j = 0; j < this->nbJoints; ++j) {
                                    this->skinWeightList[theVert * this->nbJoints + j] =
                                        undoWeight_MArr[i * this->nbJoints + j];
                                }
                            }
                            replace_weights(dataBlock, editVertsIndices, undoWeight_MArr);

                            this->undoVertsIndices_.pop_back();
                            this->undoVertsValues_.pop_back();
                            refreshColors(editVertsIndices, multiEditColors, soloEditColors);
                        } else {
                            MGlobal::displayInfo("  NO MORE UNDOS ");
                        }
                    }
                    // now set Attr false
                    MPlug callUndoPlug(thisMObject(), _callUndo);
                    callUndoPlug.setBool(false);
                } else if (this->changedColorInfluence != -1) {
                    if (verbose)
                        MGlobal::displayInfo(MString("         --> changing Color ") +
                                             this->changedColorInfluence);
                    // now change all the colors of the multi --------------------------
                    for (int theVert = 0; theVert < this->multiCurrentColors.length(); ++theVert) {
                        double inflVal = this->skinWeightList[theVert * this->nbJoints +
                                                              this->changedColorInfluence];
                        if (inflVal != 0.0) {
                            MColor multiColor;
                            for (int j = 0; j < this->nbJoints; ++j) {  // for each joint
                                double val = this->skinWeightList[theVert * this->nbJoints + j];
                                multiColor += this->jointsColors[j] * val;
                            }
                            this->multiCurrentColors[theVert] = multiColor;
                            editVertsIndices.append(theVert);
                            multiEditColors.append(multiColor);
                        }
                    }
                    this->changedColorInfluence = -1;
                    meshFn.setSomeColors(editVertsIndices, multiEditColors, &this->fullColorSet);
                    editSoloColorSet(meshFn);

                    dataBlock.setClean(plug);
                    return status;
                } else {
                    if (verbose)
                        MGlobal::displayInfo(MString("         --> actually painting weights"));
                    // read paint values ---------------------------
                    MFnDoubleArrayData arrayData;
                    MObject dataObj = dataBlock.inputValue(_paintableAttr).data();
                    arrayData.setObject(dataObj);

                    unsigned int length = arrayData.length();
                    for (unsigned int i = 0; i < length; i++) {
                        double val = arrayData[i];
                        if ((val > 0) &&
                            (!this->lockVertices[i])) {  // superior zero and not locked
                            if (val !=
                                this->paintedValues[i]) {  // only if other zone painted ----------
                                val = std::max(0.0, std::min(val, 1.0));  // clamp
                                editVertsIndices.append(i);
                                editVertsWeights.append(val);
                                // MGlobal::displayInfo(MString(" paint value ") + i + MString(" -
                                // ") + val);
                                this->paintedValues[i] = val;  // store to not repaint

                                val *= intensity;
                                val = std::log10(val * 9 + 1);
                                multiEditColors.append(val * multColor +
                                                       (1.0 - val) * this->multiCurrentColors[i]);
                                soloEditColors.append(val * soloMultColor +
                                                      (1.0 - val) * this->soloCurrentColors[i]);
                            }
                        }
                    }
                    // during brushing apply values
                    // ---------------------------------------------------
                    if (!this->postSetting) {
                        if (this->postSetting_timeToStoreUndo) {
                            this->postSetting_timeToStoreUndo = false;
                            this->fullUndoSkinWeightList.copy(this->skinWeightList);
                        }
                        applyCommand(dataBlock, editVertsIndices, editVertsWeights, false);
                    }
                }
                meshFn.setSomeColors(editVertsIndices, multiEditColors, &this->fullColorSet);
                meshFn.setSomeColors(editVertsIndices, soloEditColors, &this->soloColorSet);
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

MStatus blurSkinDisplay::editSoloColorSet(MFnMesh& meshFn) {
    MStatus status;

    if (verbose) MGlobal::displayInfo(" editSoloColorSet CALL ");
    MColorArray colToSet;
    MIntArray vtxToSet;
    for (int theVert = 0; theVert < this->soloColorsValues.length(); ++theVert) {
        double val = this->skinWeightList[theVert * this->nbJoints + this->influenceIndex];
        // if (theVert == 20266) MGlobal::displayInfo(MString("
        // Mesh_X_HeadBody_Pc_Sd1_SdDsp_.vtx[20266] -  ") + val + MString(" - storeValue ") +
        // this->soloColorsValues[theVert]);
        bool isVtxLocked = this->lockVertices[theVert] == 1;
        if (!(this->soloColorsValues[theVert] == 0 && val == 0)) {
            MColor soloColor;
            if (this->soloColorTypeVal == 0) {
                soloColor = MColor(val, val, val);
            } else if (this->soloColorTypeVal == 1) {
                val *= 2;
                if (val > 1)
                    soloColor = MColor(val, (val - 1), 0);
                else
                    soloColor = MColor(val, 0, 0);
            } else {
                soloColor = val * this->jointsColors[this->influenceIndex];
            }
            if (isVtxLocked) soloColor = this->lockVertColor;

            this->soloCurrentColors[theVert] = soloColor;
            this->soloColorsValues[theVert] = val;
            colToSet.append(soloColor);
            vtxToSet.append(theVert);
        }
    }
    meshFn.setSomeColors(vtxToSet, colToSet, &this->soloColorSet);
    return status;
}

MStatus blurSkinDisplay::applyCommand(MDataBlock& dataBlock, MIntArray& theEditVerts,
                                      MDoubleArray& verticesWeight, bool storeUndo) {
    // 0 Add - 1 Remove - 2 AddPercent - 3 Absolute - 4 Smooth - 5 Sharpen - 6 LockVertices
    MStatus status;
    if (verbose) MGlobal::displayInfo(MString(" applyCommand ") + this->commandIndex);

    if (this->commandIndex < 6) {
        MDoubleArray previousWeights(this->nbJoints * theEditVerts.length(), 0.0);
        // std::vector< double > previousWeights;
        // std::vector< int > undoVerts;
        // undoVerts.resize(theEditVerts.length());
        // previousWeights.resize(this->nbJoints*theEditVerts.length());

        MDoubleArray theWeights(this->nbJoints * theEditVerts.length(), 0.0);
        int repeatLimit = 1;
        if (this->commandIndex == 4 || this->commandIndex == 5) repeatLimit = this->smoothRepeat;
        for (int repeat = 0; repeat < repeatLimit; ++repeat) {
            if (this->commandIndex == 4) {
                for (int i = 0; i < theEditVerts.length(); ++i) {
                    int theVert = theEditVerts[i];
                    double theVal = verticesWeight[i];

                    MIntArray vertsAround = this->allVertsAround[theVert];
                    status = setAverageWeight(vertsAround, theVert, i, this->nbJoints,
                                              this->lockJoints, this->skinWeightList, theWeights);
                }
            } else {
                if (this->lockJoints[influenceIndex] == 1 && this->commandIndex != 5)
                    return status;  //  if locked and it's not sharpen --> do nothing
                status = editArray(this->commandIndex, influenceIndex, this->nbJoints,
                                   this->lockJoints, this->skinWeightList, theEditVerts,
                                   verticesWeight, theWeights, this->doNormalize);
            }
            // now set the weights -----------------------------------------------------
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
            MIntArray undoVerts;
            undoVerts.copy(theEditVerts);
            // now store the undo ----------------
            // for (int i = 0; i < theEditVerts.length(); ++i) undoVerts[i] = theEditVerts[i];
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

void blurSkinDisplay::refreshVertsConnection() {
    this->allVertsAround.clear();
    this->allVertsAround.resize(this->connectedVertices.size());
    for (int i = 0; i < this->connectedVertices.size(); ++i) {
        MIntArray surroundingVertices = this->connectedVertices[i];
        std::unordered_set<int> setOfVerts;
        for (unsigned int itVtx = 0; itVtx < surroundingVertices.length(); itVtx++)
            setOfVerts.insert(surroundingVertices[itVtx]);
        // for the repeats
        for (int d = 1; d < this->smoothDepth; d++) {  // <= to add one more
            for (unsigned int itVtx = 0; itVtx < surroundingVertices.length(); itVtx++) {
                int vtx = surroundingVertices[itVtx];
                // for (int vtx : setOfVerts) {
                MIntArray repeatVertices = this->connectedVertices[vtx];
                for (unsigned int itVtx = 0; itVtx < repeatVertices.length(); itVtx++)
                    setOfVerts.insert(repeatVertices[itVtx]);
            }
        }
        MIntArray vertsAround;
        for (int vtx : setOfVerts) vertsAround.append(vtx);
        this->allVertsAround[i] = vertsAround;
        // this->set_vertsAround.push_back( setOfVerts );
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
    }
}

void blurSkinDisplay::connectSkinClusterWL() {
    MStatus status;
    MPlug weight_list_plug(thisMObject(), blurSkinDisplay::_s_skin_weights);
    MPlugArray plugs;
    weight_list_plug.connectedTo(plugs, true, false, &status);
    if (plugs.length() == 0) {  // not connected to the weightList
        MFnDependencyNode skinDep(skinCluster_);
        MPlug weight_list_skin_clus = skinDep.findPlug("weightList");
        MDGModifier dg;
        // status = dg.connect( weight_list_skin_clus, weight_list_plug);
        MGlobal::displayInfo(MString(" TRY CONNECT ") + weight_list_plug.name() + MString(" -> ") +
                             weight_list_skin_clus.name());
        status = dg.connect(weight_list_plug, weight_list_skin_clus);
        if (MS::kSuccess != status)
            MGlobal::displayError(
                MString("FAIL dg.connect( weight_list_skin_clus, weight_list_plug);"));
        status = dg.doIt();
        if (MS::kSuccess != status) MGlobal::displayError(MString("FAIL dg.doIt ;"));
    }

    this->doConnectSkinCL = false;
}

void blurSkinDisplay::setInfluenceColorAttr() {
    MStatus status;
    MPlug influenceColor_Plug(thisMObject(), _influenceColor);

    for (int i = 0; i < this->jointsColors.length(); ++i) {
        MPlug thecolorPlug = influenceColor_Plug.elementByLogicalIndex(i, &status);
        MColor theColor = this->jointsColors[i];
        thecolorPlug.child(0).setValue(theColor.r);
        thecolorPlug.child(1).setValue(theColor.g);
        thecolorPlug.child(2).setValue(theColor.b);
    }
    this->changedColorInfluence = -1;
}

MStatus blurSkinDisplay::refreshColors(MIntArray& editVertsIndices, MColorArray& multiEditColors,
                                       MColorArray& soloEditColors) {
    MStatus status = MS::kSuccess;
    if (verbose) MGlobal::displayInfo(MString(" refreshColors CALL "));  // beginning opening of
                                                                         // node
    if (multiEditColors.length() != editVertsIndices.length())
        multiEditColors.setLength(editVertsIndices.length());
    if (soloEditColors.length() != editVertsIndices.length())
        soloEditColors.setLength(editVertsIndices.length());

    for (int i = 0; i < editVertsIndices.length(); ++i) {
        int theVert = editVertsIndices[i];
        MColor multiColor, soloColor;
        bool isVtxLocked = this->lockVertices[theVert] == 1;
        if (!isVtxLocked) {
            for (int j = 0; j < this->nbJoints; ++j) {  // for each joint
                double val = this->skinWeightList[theVert * this->nbJoints + j];
                ;
                multiColor += jointsColors[j] * val;
                if (j == this->influenceIndex) {
                    this->soloColorsValues[i] = val;

                    if (this->soloColorTypeVal == 0) {
                        soloColor = MColor(val, val, val);
                    } else if (this->soloColorTypeVal == 1) {
                        val *= 2;
                        if (val > 1)
                            soloColor = MColor(val, (val - 1), 0);
                        else
                            soloColor = MColor(val, 0, 0);
                    } else {
                        soloColor = val * this->jointsColors[this->influenceIndex];
                    }
                }
            }
        } else {
            multiColor = this->lockVertColor;
            soloColor = this->lockVertColor;
        }
        this->multiCurrentColors[theVert] = multiColor;
        this->soloCurrentColors[theVert] = soloColor;
        multiEditColors[i] = multiColor;
        soloEditColors[i] = soloColor;
    }
    return status;
}

MStatus blurSkinDisplay::querySkinClusterValues(MIntArray& verticesIndices,
                                                MDoubleArray& theWeights, bool doColors) {
    MStatus status = MS::kSuccess;

    MFnDependencyNode skinClusterDep(skinCluster_);

    MPlug weight_list_plug = skinClusterDep.findPlug("weightList");

    for (int i = 0; i < verticesIndices.length(); ++i) {
        int vertexIndex = verticesIndices[i];
        // weightList[i]
        MPlug ith_weights_plug = weight_list_plug.elementByLogicalIndex(vertexIndex);

        // weightList[i].weight
        MPlug plug_weights = ith_weights_plug.child(0);  // access first compound child
        int nb_weights = plug_weights.numElements();

        MColor theColor;
        for (int j = 0; j < nb_weights; j++) {  // for each joint
            MPlug weight_plug = plug_weights.elementByPhysicalIndex(j);
            // weightList[i].weight[j]
            int indexInfluence = weight_plug.logicalIndex();
            double theWeight = weight_plug.asDouble();

            this->skinWeightList[vertexIndex * nbJoints + indexInfluence] = theWeight;
            theWeights[i * nbJoints + indexInfluence] = theWeight;
            if (doColors) theColor += jointsColors[indexInfluence] * theWeight;
        }
        if (doColors) this->multiCurrentColors[vertexIndex] = theColor;
    }
    if (verbose) MGlobal::displayInfo(MString(" querySkinClusterValues "));

    return status;
}

MStatus blurSkinDisplay::fillArrayValues(bool doColors) {
    MStatus status = MS::kSuccess;

    MFnDependencyNode skinClusterDep(skinCluster_);

    MPlug weight_list_plug = skinClusterDep.findPlug("weightList");
    MPlug matrix_plug = skinClusterDep.findPlug("matrix");
    // MGlobal::displayInfo(weight_list_plug.name());
    int nbElements = weight_list_plug.numElements();
    this->nbJoints = matrix_plug.numElements();

    skin_weights_.resize(nbElements);
    if (doColors) {
        this->multiCurrentColors.clear();
        this->multiCurrentColors.setLength(nbElements);
    }
    this->skinWeightList = MDoubleArray(nbElements * this->nbJoints, 0.0);

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
        bool isVtxLocked = this->lockVertices[i] == 1;
        for (int j = 0; j < nb_weights; j++) {  // for each joint
            MPlug weight_plug = plug_weights.elementByPhysicalIndex(j);
            // weightList[i].weight[j]
            int indexInfluence = weight_plug.logicalIndex();
            double theWeight = weight_plug.asDouble();

            skin_weights_[i][j] = std::make_pair(indexInfluence, theWeight);
            this->skinWeightList[vertexIndex * this->nbJoints + indexInfluence] = theWeight;
            if ((doColors) && (!isVtxLocked))  // and not locked
                theColor += this->jointsColors[indexInfluence] * theWeight;
        }
        if (doColors)
            if (isVtxLocked)
                this->multiCurrentColors[vertexIndex] = this->lockVertColor;
            else
                this->multiCurrentColors[vertexIndex] = theColor;
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

    MFnTypedAttribute attrFn;
    _cpList = attrFn.create("inputComponents", "ics", MFnComponentListData::kComponentList);
    attrFn.setStorable(false);  // To be stored during file-save
    addAttribute(_cpList);

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

    blurSkinDisplay::_getLockWeights =
        numAtt.create("getLockWeights", "glw", MFnNumericData::kBoolean, false, &status);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_getLockWeights);

    blurSkinDisplay::_normalize =
        numAtt.create("normalize", "nom", MFnNumericData::kBoolean, true, &status);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_normalize);

    ///////////////////////////////////////////////////////////////////////////
    // creation attributes
    ///////////////////////////////////////////////////////////////////////////
    MFnEnumAttribute enumAttr, enumAttr2;

    blurSkinDisplay::_commandAttr = enumAttr.create("command", "cmd", 0);
    CHECK_MSTATUS(enumAttr.addField("Add", 0));
    CHECK_MSTATUS(enumAttr.addField("Remove", 1));
    CHECK_MSTATUS(enumAttr.addField("AddPercent", 2));
    CHECK_MSTATUS(enumAttr.addField("Absolute", 3));
    CHECK_MSTATUS(enumAttr.addField("Smooth", 4));
    CHECK_MSTATUS(enumAttr.addField("Sharpen", 5));
    CHECK_MSTATUS(enumAttr.addField("Locks", 6));
    // CHECK_MSTATUS(enumAttr.addField("UpdateSkin", 7));
    CHECK_MSTATUS(enumAttr.setStorable(true));
    CHECK_MSTATUS(enumAttr.setKeyable(true));
    CHECK_MSTATUS(enumAttr.setReadable(true));
    CHECK_MSTATUS(enumAttr.setWritable(true));
    CHECK_MSTATUS(enumAttr.setCached(false));

    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_commandAttr);

    blurSkinDisplay::_autoExpandAttr =
        numAtt.create("autoExpand", "aex", MFnNumericData::kBoolean, false, &status);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_autoExpandAttr);

    blurSkinDisplay::_smoothRepeat =
        numAtt.create("smoothRepeat", "sr", MFnNumericData::kInt64, 3, &status);
    numAtt.setMin(1);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_smoothRepeat);

    blurSkinDisplay::_smoothDepth =
        numAtt.create("smoothDepth", "dpt", MFnNumericData::kInt64, 1, &status);
    numAtt.setMin(1);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_smoothDepth);

    blurSkinDisplay::_influenceAttr =
        numAtt.create("influenceIndex", "ii", MFnNumericData::kInt64, 0, &status);
    numAtt.setMin(0);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_influenceAttr);

    blurSkinDisplay::_colorType = enumAttr.create("colorType", "cty", 0);
    CHECK_MSTATUS(enumAttr.addField("Multi", 0));
    CHECK_MSTATUS(enumAttr.addField("Solo", 1));
    CHECK_MSTATUS(enumAttr.addField("None", 2));
    CHECK_MSTATUS(enumAttr.setStorable(true));
    CHECK_MSTATUS(enumAttr.setKeyable(true));
    CHECK_MSTATUS(enumAttr.setReadable(true));
    CHECK_MSTATUS(enumAttr.setWritable(true));
    CHECK_MSTATUS(enumAttr.setCached(false));

    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_colorType);

    blurSkinDisplay::_soloColorType = enumAttr.create("soloColType", "scty", 0);
    CHECK_MSTATUS(enumAttr.addField("white", 0));
    CHECK_MSTATUS(enumAttr.addField("Lava", 1));
    CHECK_MSTATUS(enumAttr.addField("Influence", 2));
    CHECK_MSTATUS(enumAttr.setStorable(true));
    CHECK_MSTATUS(enumAttr.setKeyable(true));
    CHECK_MSTATUS(enumAttr.setReadable(true));
    CHECK_MSTATUS(enumAttr.setWritable(true));
    CHECK_MSTATUS(enumAttr.setCached(false));

    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_soloColorType);

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
    // theColors
    ///////////////////////////////////////////////////////////////////////////
    MFnNumericAttribute nAttr;
    _influenceColor = nAttr.createColor("influenceColor", "ic");
    nAttr.setDefault(1.0, 1.0, 0.0);
    nAttr.setKeyable(false);
    nAttr.setUsedAsColor(true);
    nAttr.setChannelBox(false);

    nAttr.setArray(true);
    nAttr.setUsesArrayDataBuilder(true);
    addAttribute(_influenceColor);

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

MStatus blurSkinDisplay::connectionBroken(const MPlug& plug, const MPlug& otherPlug, bool asSrc) {
    if (plug == _s_skin_weights) {
        MGlobal::displayInfo(" disconnect  _s_skin_weights");
        MGlobal::displayInfo(plug.name() + " " + otherPlug.name());
    }
    return MS::kUnknownParameter;
};

MStatus blurSkinDisplay::setDependentsDirty(const MPlug& plugBeingDirtied,
                                            MPlugArray& affectedPlugs) {
    MStatus status;
    MObject thisNode = thisMObject();
    MFnDependencyNode fnThisNode(thisNode);
    this->reloadCommand = plugBeingDirtied == _commandAttr || plugBeingDirtied == _influenceAttr ||
                          plugBeingDirtied == _smoothRepeat || plugBeingDirtied == _smoothDepth ||
                          plugBeingDirtied == _postSetting || plugBeingDirtied == _colorType ||
                          plugBeingDirtied == _cpList || plugBeingDirtied == _getLockWeights ||
                          plugBeingDirtied == _soloColorType || plugBeingDirtied == _normalize ||
                          plugBeingDirtied == _autoExpandAttr;

    this->clearTheArray = (plugBeingDirtied == _clearArray);
    this->callUndo = (plugBeingDirtied == _callUndo);
    this->inputVerticesChanged = (plugBeingDirtied == _cpList);

    if (!(plugBeingDirtied == _paintableAttr || this->reloadCommand || this->clearTheArray ||
          this->callUndo)) {
        if (plugBeingDirtied.isChild()) {
            MPlug prt = plugBeingDirtied.parent();
            int ind = prt.logicalIndex();
            bool isGreenChannel =
                plugBeingDirtied.partialName() == MString("ic[") + ind + MString("].icb");
            // this->changedColor = (plugBeingDirtied == _influenceColor);
            if (isGreenChannel) {
                this->changedColorInfluence = ind;
                if (isGreenChannel && verbose)
                    MGlobal::displayInfo(" dirty dirty " + plugBeingDirtied.partialName() + " " +
                                         prt.name() + " isChild " + ind);
            }
        }
    }
    if ((plugBeingDirtied == _paintableAttr) || this->reloadCommand || this->clearTheArray ||
        this->callUndo || this->changedColorInfluence) {
        this->applyPaint = true;
        MPlug outMeshPlug(thisNode, blurSkinDisplay::_outMesh);
        affectedPlugs.append(outMeshPlug);
    }
    return (MS::kSuccess);
}

/*
MStatus blurSkinDisplay::postEvaluation(const MDGContext& context, const MEvaluationNode&
evaluationNode, PostEvaluationType evalType)
{
        MGlobal::displayInfo(" postEvaluation ");

        if (!context.isNormal())
                return MStatus::kFailure;

        MStatus status;

        //http://around-the-corner.typepad.com/adn/2017/05/improve-performance-with-mpxnodepreevaluation-mpxnodepostevaluation.html
        if (evaluationNode.dirtyPlugExists(_cpList, &status) && status)
        {
                MGlobal::displayInfo(" postEvaluation - _cpList ");

                if( this->doConnectSkinCL )connectSkinClusterWL();

        }
        return MStatus::kSuccess;
}


*/