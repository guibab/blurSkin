//-
// ==========================================================================
// Copyright 2015 Autodesk, Inc.  All rights reserved.
// Use of this software is subject to the terms of the Autodesk license agreement
// provided at the time of installation or download, or which otherwise
// accompanies this software in either electronic or hard copy form.
// ==========================================================================
//+

////////////////////////////////////////////////////////////////////////
// DESCRIPTION:
//
//
// This plug-in demonstrates how to draw a simple mesh like foot Print in an easy way.
//
// This easy way is supported in Viewport 2.0.
// In Viewport 2.0, MUIDrawManager can used to draw simple UI elements in method addUIDrawables().
//
// For comparison, you can reference a Maya Developer Kit sample named rawpointsDisplayNode.
// In that sample, we draw the pointsDisplay with OpenGL\DX in method
// rawPointsDisplayDrawOverride::draw().
//
// Note the method
//   pointsDisplay::draw()
// is only called in legacy default viewport to draw foot Print.
// while the methods
//   PointsDisplayDrawOverride::prepareForDraw()
//   PointsDisplayDrawOverride::addUIDrawables()
// are only called in Viewport 2.0 to prepare and draw foot Print.
//
////////////////////////////////////////////////////////////////////////

#include <maya/M3dView.h>
#include <maya/MArrayDataBuilder.h>
#include <maya/MColor.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MDistance.h>
#include <maya/MEvaluationNode.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MPlug.h>
#include <maya/MPxLocatorNode.h>
#include <maya/MString.h>
#include <maya/MTypeId.h>
#include <maya/MVector.h>

// added by me
#include <maya/MFnComponentListData.h>
#include <maya/MFnDoubleIndexedComponent.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnNurbsSurface.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MPlugArray.h>

// Viewport 2.0 includes
#include <assert.h>
#include <maya/MDrawContext.h>
#include <maya/MDrawRegistry.h>
#include <maya/MEventMessage.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MGlobal.h>
#include <maya/MHWGeometryUtilities.h>
#include <maya/MPointArray.h>
#include <maya/MPxDrawOverride.h>
#include <maya/MUserData.h>

// Foot Data
//
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Node implementation with standard viewport draw
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
class PointsDisplayData : public MUserData {
   public:
    PointsDisplayData() : MUserData(false) {}  // don't delete after draw
    ~PointsDisplayData() override {}

    MColor fColor;
    unsigned int fDepthPriority;
    MPointArray fLineList;
    MPointArray fTriangleList;

    virtual void getData(const MObject&);
    float pointWidth = 1;
    bool enableSmooth = true;
    float color[4] = {1.0f, 0.0f, 0.0f, 0.1f};
    MPointArray pointsVertices;
    MBoundingBox theBoundingBox;
};

class pointsDisplay : public MPxLocatorNode {
   public:
    pointsDisplay();
    ~pointsDisplay() override;

    MStatus compute(const MPlug& plug, MDataBlock& data) override;
    virtual void postConstructor();
    void draw(M3dView& view, const MDagPath& path, M3dView::DisplayStyle style,
              M3dView::DisplayStatus status) override;

    bool isBounded() const override;
    MBoundingBox boundingBox() const override;

    MStatus preEvaluation(const MDGContext& context,
                          const MEvaluationNode& evaluationNode) override;

    static void* creator();
    static MStatus initialize();

    static MObject size;  // The size of the foot
    // static  MObject         _inMesh;
    static MObject _inGeometry;
    static MObject _inputColor;
    static MObject _inputAlpha;
    static MObject _pointWidth;
    static MObject _cpList;
    static MObject _enableSmooth;

   public:
    static MTypeId id;
    static MString drawDbClassification;
    static MString drawRegistrantId;

    static MObject worldS;

   private:
    MObject _self;
};

class PointsDisplayDrawOverride : public MHWRender::MPxDrawOverride {
   public:
    static MHWRender::MPxDrawOverride* Creator(const MObject& obj) {
        return new PointsDisplayDrawOverride(obj);
    }

    ~PointsDisplayDrawOverride() override;

    MHWRender::DrawAPI supportedDrawAPIs() const override;

    bool isBounded(const MDagPath& objPath, const MDagPath& cameraPath) const override;

    MBoundingBox boundingBox(const MDagPath& objPath, const MDagPath& cameraPath) const override;

    MUserData* prepareForDraw(const MDagPath& objPath, const MDagPath& cameraPath,
                              const MHWRender::MFrameContext& frameContext,
                              MUserData* oldData) override;

    bool hasUIDrawables() const override { return true; }

    void addUIDrawables(const MDagPath& objPath, MHWRender::MUIDrawManager& drawManager,
                        const MHWRender::MFrameContext& frameContext,
                        const MUserData* data) override;

    bool traceCallSequence() const override {
        // Return true if internal tracing is desired.
        return false;
    }
    void handleTraceMessage(const MString& message) const override {
        MGlobal::displayInfo("pointsDisplayDrawOverride: " + message);

        // Some simple custom message formatting.
        fputs("pointsDisplayDrawOverride: ", stderr);
        fputs(message.asChar(), stderr);
        fputs("\n", stderr);
    }

   private:
    PointsDisplayDrawOverride(const MObject& obj);
    float getMultiplier(const MDagPath& objPath) const;

    static void OnModelEditorChanged(void* clientData);

    pointsDisplay* fPointsDisplay;
    MCallbackId fModelEditorChangedCbId;
};