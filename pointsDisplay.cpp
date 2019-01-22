#include "pointsDisplay.h"

MObject pointsDisplay::size;
// MObject pointsDisplay::_inMesh;
MObject pointsDisplay::_inGeometry;
MObject pointsDisplay::_inputColor;
MObject pointsDisplay::_inputAlpha;
MObject pointsDisplay::_pointWidth;
MObject pointsDisplay::_enableSmooth;
MObject pointsDisplay::_cpList;

MTypeId pointsDisplay::id(0x001226F8);
MString pointsDisplay::drawDbClassification("drawdb/geometry/pointsDisplay");
MString pointsDisplay::drawRegistrantId("PointsDisplayNodePlugin");
MObject pointsDisplay::worldS;

pointsDisplay::pointsDisplay() {}
pointsDisplay::~pointsDisplay() {}

void pointsDisplay::postConstructor() {
    MObject self = thisMObject();
    MFnDependencyNode fn_node(self);
    fn_node.setName("pointsDisplayShape#");

    _self = self;
}

MStatus pointsDisplay::compute(const MPlug& plug /*plug*/, MDataBlock& dataBlock /*data*/) {
    MStatus s;
    /*
    MDataHandle inputIDs = dataBlock.inputValue(_cpList, &s);
    MObject compList = inputIDs.data();
    MFnComponentListData compListFn(compList);
    unsigned i;
    int j;
    cpIds.clear();

    MFn::Type componentType = MFn::kMeshVertComponent;
    for (i = 0; i < compListFn.length(); i++)
    {
            MObject comp = compListFn[i];
            if (comp.apiType() == componentType)
            {
                    MFnSingleIndexedComponent siComp(comp);
                    for (j = 0; j < siComp.elementCount(); j++)
                            cpIds.append(siComp.element(j));
            }
    }
    */
    if (plug == worldS) {
        if (plug.isElement()) {
            MArrayDataHandle outputArrayHandle = dataBlock.outputArrayValue(worldS);
            outputArrayHandle.setAllClean();
        }
        dataBlock.setClean(plug);
        return MS::kSuccess;
    }

    return MS::kUnknownParameter;
    ;
}

// called by legacy default viewport
void pointsDisplay::draw(M3dView& view, const MDagPath& /*path*/, M3dView::DisplayStyle style,
                         M3dView::DisplayStatus status) {
    PointsDisplayData data;
    data.getData(_self);

    // Get the size
    //
    MObject thisNode = thisMObject();
    MPlug plug(thisNode, size);
    MDistance sizeVal;
    plug.getValue(sizeVal);

    view.beginGL();

    // draw the vertices --------------------------
    glBegin(GL_POINTS);
    /*
    if (status == M3dView::kLead) {
            glColor4f(0., 1., 0., data.color[3]);	//green
    }
    else if (status == M3dView::kActive) {
            glColor4f(1., 1., 1., data.color[3]);	//white
    }
    else {
            glColor4fv(data.color);
    }
    */
    glColor4fv(data.color);
    glPointSize(data.pointWidth);

    for (int i = 0; i < data.pointsVertices.length(); i++) {
        glVertex3f(data.pointsVertices[i][0], data.pointsVertices[i][1], data.pointsVertices[i][2]);
    }
    glEnd();

    view.endGL();
    /*
    // Draw the name of the pointsDisplay
    view.setDrawColor( MColor( 0.1f, 0.8f, 0.8f, 1.0f ) );
    view.drawText( MString("Footprint"), MPoint( 0.0, 0.0, 0.0 ), M3dView::kCenter );
    */
}

bool pointsDisplay::isBounded() const { return true; }

MBoundingBox pointsDisplay::boundingBox() const {
    // Get the size
    //
    PointsDisplayData data;
    data.getData(_self);
    return data.theBoundingBox;
}

// Called before this node is evaluated by Evaluation Manager
MStatus pointsDisplay::preEvaluation(const MDGContext& context,
                                     const MEvaluationNode& evaluationNode) {
    if (context.isNormal()) {
        MStatus status;
        if ((evaluationNode.dirtyPlugExists(_cpList, &status) && status) ||
            (evaluationNode.dirtyPlugExists(_inGeometry, &status) && status)) {
            MHWRender::MRenderer::setGeometryDrawDirty(thisMObject());
        }
    }
    return MStatus::kSuccess;
}

void* pointsDisplay::creator() { return new pointsDisplay(); }

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Viewport 2.0 override implementation
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void PointsDisplayData::getData(const MObject& node) {
    MStatus status;

    this->pointWidth = MPlug(node, pointsDisplay::_pointWidth).asFloat();
    this->enableSmooth = MPlug(node, pointsDisplay::_enableSmooth).asBool();

    MPlug inputColPlug =
        MPlug(node, pointsDisplay::_inputColor);  // fn.findPlug("overrideColorRGB");
    // MPlug inputColPlug = MPlug(node, MPxLocatorNode::overrideColorR);
    this->color[0] = inputColPlug.child(0).asFloat();
    this->color[1] = inputColPlug.child(1).asFloat();
    this->color[2] = inputColPlug.child(2).asFloat();
    this->color[3] = MPlug(node, pointsDisplay::_inputAlpha).asFloat();

    this->fColor = MColor(this->color);

    // get the geo -------------------------------------------------------------------------
    MPlug geometryPlug = MPlug(node, pointsDisplay::_inGeometry);
    MPlugArray plugs;
    geometryPlug.connectedTo(plugs, true, false, &status);
    this->pointsVertices.clear();

    if (plugs.length() == 0) {
        // MGlobal::displayError("Unable to rebind.  No geometry is connected.");
        // MGlobal::displayInfo("Unable to rebind.  No geometry is connected.");
        this->theBoundingBox = MBoundingBox();
    } else {
        MPlug cpListPlug = MPlug(node, pointsDisplay::_cpList);
        MObject compList = cpListPlug.asMObject();
        MFnComponentListData compListFn(compList);

        MObject theNode = plugs[0].node();

        /*
        MFnDagNode fn(theNode);
        MPlug matrixPlug = fn.findPlug("worldMatrix");
        matrixPlug = matrixPlug.elementByLogicalIndex(0);

        MObject matrixObject;
        matrixObject = matrixPlug.asMObject();

        MFnMatrixData worldMatrixData(matrixObject);
        MMatrix worldMatrix = worldMatrixData.matrix();
        */

        // get the transform  matrix
        MFnDagNode theShape(theNode);
        MObject prt = theShape.parent(0);

        MDagPath pth;
        status = MDagPath::getAPathTo(theNode, pth);
        MMatrix worldMatrix = pth.inclusiveMatrix();
        // MGlobal::displayInfo(MString("   --> NAME : ")+ pth.fullPathName());

        // MMatrix worldMatrix;

        if (theNode.hasFn(MFn::kMesh)) {
            MFnMesh tmpMesh(theNode);
            // get data of points and bounding box ---------------------------------------
            // const float* mayaRawPoints;
            MPointArray pointsVerticesAll;

            MFnMeshData meshData;
            MMeshSmoothOptions options;
            tmpMesh.getSmoothMeshDisplayOptions(options);
            this->theBoundingBox = tmpMesh.boundingBox(&status);
            int smoothLevel = 0;
            if (this->enableSmooth)
                smoothLevel = tmpMesh.findPlug("displaySmoothMesh", false, &status).asInt();

            if (smoothLevel > 0) {
                // options.setDivisions(smoothLevel);
                options.setDivisions(1);
                options.setSmoothUVs(false);
                // https://github.com/haggi/OpenMaya/blob/master/src/common/cpp/mayaObject.cpp

                MObject dataObject = meshData.create();
                MObject smoothedObj = tmpMesh.generateSmoothMesh(dataObject, &options, &status);
                MFnMesh smoothMesh(smoothedObj, &status);
                // mayaRawPoints = smoothMesh.getRawPoints(&status);
                smoothMesh.getPoints(pointsVerticesAll);
            } else {
                // mayaRawPoints = tmpMesh.getRawPoints(&status);
                tmpMesh.getPoints(pointsVerticesAll);
            }
            // this->pointsVertices = pointsVerticesAll;

            // get the stored vertices indices ------------------------------------------
            // MFnDependencyNode depNode(node, &status);
            // pointsDisplay* fpointsDisplay = status ?
            // dynamic_cast<pointsDisplay*>(depNode.userNode()) : NULL;

            MFn::Type componentType = MFn::kMeshVertComponent;
            for (unsigned i = 0; i < compListFn.length(); i++) {
                MObject comp = compListFn[i];
                if (comp.apiType() == componentType) {
                    MFnSingleIndexedComponent siComp(comp);
                    for (int j = 0; j < siComp.elementCount(); j++)
                        this->pointsVertices.append(pointsVerticesAll[siComp.element(j)] *
                                                    worldMatrix);  //);
                }
            }
        } else if (theNode.hasFn(MFn::kNurbsSurface)) {
            MFnNurbsSurface surfaceFn(theNode);
            this->theBoundingBox = surfaceFn.boundingBox(&status);

            MPointArray cvPoints;
            surfaceFn.getCVs(cvPoints);
            int numCVsInV = surfaceFn.numCVsInV();

            MFn::Type componentType = MFn::kSurfaceCVComponent;
            for (unsigned i = 0; i < compListFn.length(); i++) {
                MObject comp = compListFn[i];
                if (comp.apiType() == componentType) {
                    MFnDoubleIndexedComponent siComp(comp);
                    for (int j = 0; j < siComp.elementCount(); j++) {
                        int indexU, indexV;
                        // MGlobal::displayInfo(MString(" indexU, indexV")+ indexU + MString(" ")+
                        // indexV);
                        siComp.getElement(j, indexU, indexV);
                        int vertInd = numCVsInV * indexU + indexV;

                        this->pointsVertices.append(cvPoints[vertInd] * worldMatrix);
                    }
                }
            }
        } else if (theNode.hasFn(MFn::kNurbsCurve)) {
            MFnNurbsCurve curveFn(theNode);
            this->theBoundingBox = curveFn.boundingBox(&status);

            MPointArray cvPoints;
            curveFn.getCVs(cvPoints);

            // MDataHandle inputIDs = dataBlock.inputValue(_cpList, &s);
            // MObject compList = inputIDs.data();

            MFn::Type componentType = MFn::kCurveCVComponent;
            for (unsigned i = 0; i < compListFn.length(); i++) {
                MObject comp = compListFn[i];
                if (comp.apiType() == componentType) {
                    MFnSingleIndexedComponent siComp(comp);
                    for (int j = 0; j < siComp.elementCount(); j++)
                        this->pointsVertices.append(cvPoints[siComp.element(j)] * worldMatrix);
                }
            }
        } else if (theNode.hasFn(MFn::kLattice)) {  // lattice and everything else
            MFnLattice latticeFn(theNode);
            this->theBoundingBox = latticeFn.boundingBox(&status);

            MFn::Type componentType = MFn::kLatticeComponent;
            for (unsigned i = 0; i < compListFn.length(); i++) {
                MObject comp = compListFn[i];
                if (comp.apiType() == componentType) {
                    MFnTripleIndexedComponent siComp(comp);
                    for (int j = 0; j < siComp.elementCount(); j++) {
                        int s, t, u;
                        // MGlobal::displayInfo(MString(" indexU, indexV")+ indexU + MString(" ")+
                        // indexV);
                        siComp.getElement(j, s, t, u);

                        MPoint thePt = latticeFn.point(s, t, u, &status);
                        this->pointsVertices.append(thePt * worldMatrix);
                    }
                }
            }
        }
    }
}

// By setting isAlwaysDirty to false in MPxDrawOverride constructor, the
// draw override will be updated (via prepareForDraw()) only when the node
// is marked dirty via DG evaluation or dirty propagation. Additional
// callback is also added to explicitly mark the node as being dirty (via
// MRenderer::setGeometryDrawDirty()) for certain circumstances. Note that
// the draw callback in MPxDrawOverride constructor is set to NULL in order
// to achieve better performance.
PointsDisplayDrawOverride::PointsDisplayDrawOverride(const MObject& obj)
    : MHWRender::MPxDrawOverride(obj, NULL, false) {
    fModelEditorChangedCbId =
        MEventMessage::addEventCallback("modelEditorChanged", OnModelEditorChanged, this);

    MStatus status;
    MFnDependencyNode node(obj, &status);
    fPointsDisplay = status ? dynamic_cast<pointsDisplay*>(node.userNode()) : NULL;
}

PointsDisplayDrawOverride::~PointsDisplayDrawOverride() {
    fPointsDisplay = NULL;

    if (fModelEditorChangedCbId != 0) {
        MMessage::removeCallback(fModelEditorChangedCbId);
        fModelEditorChangedCbId = 0;
    }
}

void PointsDisplayDrawOverride::OnModelEditorChanged(void* clientData) {
    // Mark the node as being dirty so that it can update on display appearance
    // switch among wireframe and shaded.
    PointsDisplayDrawOverride* ovr = static_cast<PointsDisplayDrawOverride*>(clientData);
    if (ovr && ovr->fPointsDisplay) {
        MHWRender::MRenderer::setGeometryDrawDirty(ovr->fPointsDisplay->thisMObject());
    }
}

MHWRender::DrawAPI PointsDisplayDrawOverride::supportedDrawAPIs() const {
    // this plugin supports both GL and DX
    return (MHWRender::kOpenGL | MHWRender::kDirectX11 | MHWRender::kOpenGLCoreProfile);
}

float PointsDisplayDrawOverride::getMultiplier(const MDagPath& objPath) const {
    // Retrieve value of the size attribute from the node
    MStatus status;
    MObject pointsdisplayNode = objPath.node(&status);

    if (status) {
        MPlug plug(pointsdisplayNode, pointsDisplay::size);
        if (!plug.isNull()) {
            MDistance sizeVal;
            if (plug.getValue(sizeVal)) {
                return (float)sizeVal.asCentimeters();
            }
        }
    }
    return 1.0f;
}

bool PointsDisplayDrawOverride::isBounded(const MDagPath& /*objPath*/,
                                          const MDagPath& /*cameraPath*/) const {
    return true;
}

MBoundingBox PointsDisplayDrawOverride::boundingBox(const MDagPath& objPath,
                                                    const MDagPath& cameraPath) const {
    PointsDisplayData data;
    MObject node = objPath.node();
    data.getData(node);
    return data.theBoundingBox;
}

// Called by Maya each time the object needs to be drawn.
MUserData* PointsDisplayDrawOverride::prepareForDraw(const MDagPath& objPath,
                                                     const MDagPath& cameraPath,
                                                     const MHWRender::MFrameContext& frameContext,
                                                     MUserData* oldData) {
    // Any data needed from the Maya dependency graph must be retrieved and cached in this stage.
    // There is one cache data for each drawable instance, if it is not desirable to allow Maya to
    // handle data caching, simply return null in this method and ignore user data parameter in draw
    // callback method. e.g. in this sample, we compute and cache the data for usage later when we
    // create the MUIDrawManager to draw pointsdisplay in method addUIDrawables().
    PointsDisplayData* data = dynamic_cast<PointsDisplayData*>(oldData);
    if (!data) {
        data = new PointsDisplayData();
    }
    MStatus stat;
    MObject node = objPath.node(&stat);
    data->getData(node);

    // get correct color and depth priority based on the state of object, e.g. active or dormant

    // data->fColor = MHWRender::MGeometryUtilities::wireframeColor(objPath);
    switch (MHWRender::MGeometryUtilities::displayStatus(objPath)) {
        case MHWRender::kLead:
        case MHWRender::kActive:
        case MHWRender::kHilite:
        case MHWRender::kActiveComponent:
            data->fDepthPriority = MHWRender::MRenderItem::sActiveWireDepthPriority;
            break;
        default:
            data->fDepthPriority = MHWRender::MRenderItem::sDormantFilledDepthPriority;
            break;
    }
    return data;
}

// addUIDrawables() provides access to the MUIDrawManager, which can be used
// to queue up operations for drawing simple UI elements such as lines, circles and
// text. To enable addUIDrawables(), override hasUIDrawables() and make it return true.
void PointsDisplayDrawOverride::addUIDrawables(const MDagPath& objPath,
                                               MHWRender::MUIDrawManager& drawManager,
                                               const MHWRender::MFrameContext& frameContext,
                                               const MUserData* data) {
    // Get data cached by prepareForDraw() for each drawable instance, then MUIDrawManager
    // can draw simple UI by these data.
    PointsDisplayData* pLocatorData = (PointsDisplayData*)data;
    if (!pLocatorData) {
        return;
    }

    drawManager.beginDrawable(MHWRender::MUIDrawManager::kNonSelectable);

    // Draw the foot print solid/wireframe
    drawManager.setColor(pLocatorData->fColor);
    drawManager.setDepthPriority(pLocatorData->fDepthPriority);

    if (frameContext.getDisplayStyle() & MHWRender::MFrameContext::kGouraudShaded) {
        drawManager.mesh(MHWRender::MUIDrawManager::kTriangles, pLocatorData->fTriangleList);
    }

    drawManager.mesh(MHWRender::MUIDrawManager::kLines, pLocatorData->fLineList);
    // draw the vertices
    drawManager.setPointSize(pLocatorData->pointWidth);
    drawManager.mesh(MHWRender::MUIDrawManager::kPoints, pLocatorData->pointsVertices);
    // drawManager.points(pLocatorData->pointsVertices, true);

    /*
    // Draw a text "Foot"
    MPoint pos( 0.0, 0.0, 0.0 ); // Position of the text
    MColor textColor( 0.1f, 0.8f, 0.8f, 1.0f ); // Text color

    drawManager.setColor( textColor );
    drawManager.setFontSize( MHWRender::MUIDrawManager::kSmallFontSize );
    drawManager.text( pos,  MString("Footprint"), MHWRender::MUIDrawManager::kCenter );
    */

    drawManager.endDrawable();
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Plugin Registration
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

MStatus pointsDisplay::initialize() {
    MFnUnitAttribute unitFn;
    MStatus stat;

    size = unitFn.create("size", "sz", MFnUnitAttribute::kDistance);
    unitFn.setDefault(1.0);

    stat = addAttribute(size);
    if (!stat) {
        stat.perror("addAttribute");
        return stat;
    }

    MFnNumericAttribute nAttr;
    MFnEnumAttribute enumAttr;

    MFnTypedAttribute inMeshAttrFn;
    /*
    _inMesh = inMeshAttrFn.create("inMesh", "inMesh", MFnData::kMesh);
    CHECK_MSTATUS(inMeshAttrFn.setStorable(true));
    CHECK_MSTATUS(inMeshAttrFn.setKeyable(false));
    CHECK_MSTATUS(inMeshAttrFn.setReadable(true));
    CHECK_MSTATUS(inMeshAttrFn.setWritable(true));
    CHECK_MSTATUS(inMeshAttrFn.setCached(false));
    CHECK_MSTATUS(addAttribute(_inMesh));
    */
    _inGeometry = inMeshAttrFn.create("inGeometry", "inGeo", MFnGeometryData::kAny);
    inMeshAttrFn.setStorable(true);
    inMeshAttrFn.setKeyable(false);
    inMeshAttrFn.setReadable(true);
    inMeshAttrFn.setWritable(true);
    inMeshAttrFn.setCached(false);
    addAttribute(_inGeometry);

    _inputColor = nAttr.createColor("inputColor", "ico");
    nAttr.setDefault(1.0, 1.0, 0.0);
    nAttr.setKeyable(true);
    nAttr.setStorable(true);
    nAttr.setUsedAsColor(true);
    nAttr.setReadable(true);
    nAttr.setWritable(true);
    nAttr.setChannelBox(false);

    // add the color attribute to our node
    stat = addAttribute(_inputColor);

    _inputAlpha = nAttr.create("inputAlpha", "ina", MFnNumericData::kFloat, 0.5);
    nAttr.setMin(0.0);
    nAttr.setMax(1.0);
    nAttr.setKeyable(true);
    nAttr.setStorable(true);
    nAttr.setReadable(true);
    nAttr.setWritable(true);
    nAttr.setChannelBox(false);

    // add the color attribute to our node
    stat = addAttribute(_inputAlpha);

    _pointWidth = nAttr.create("pointWidth", "pw", MFnNumericData::kInt, 1);
    nAttr.setMin(1);
    nAttr.setMax(10);
    nAttr.setKeyable(true);
    nAttr.setStorable(true);
    nAttr.setReadable(true);
    nAttr.setWritable(true);
    nAttr.setChannelBox(false);

    // add the color attribute to our node
    stat = addAttribute(_pointWidth);

    MFnTypedAttribute attrFn;
    _cpList = attrFn.create("inputComponents", "ics", MFnComponentListData::kComponentList);
    attrFn.setStorable(true);  // To be stored during file-save
    addAttribute(_cpList);

    _enableSmooth = nAttr.create("enableSmooth", "es", MFnNumericData::kBoolean, false);
    nAttr.setKeyable(true);
    nAttr.setStorable(true);
    nAttr.setReadable(true);
    nAttr.setWritable(true);
    nAttr.setChannelBox(false);

    // add the color attribute to our node
    stat = addAttribute(_enableSmooth);

    worldS = unitFn.create("worldS", "ws", MFnUnitAttribute::kDistance, 1.0);
    unitFn.setWritable(true);
    unitFn.setCached(false);
    unitFn.setArray(true);
    unitFn.setUsesArrayDataBuilder(true);
    unitFn.setWorldSpace(true);

    addAttribute(worldS);
    attributeAffects(size, worldS);
    attributeAffects(_inputColor, worldS);
    attributeAffects(_cpList, worldS);

    return MS::kSuccess;
}
