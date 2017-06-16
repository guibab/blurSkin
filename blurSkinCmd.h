#define PI (3.141592654)
#include <math.h>
#include <maya/M3dView.h>
#include <maya/MArgDatabase.h>
#include <maya/MArgList.h>
#include <maya/MDagPath.h>
#include <maya/MDagPathArray.h>
#include <maya/MFloatArray.h>
#include <maya/MFloatMatrix.h>
#include <maya/MFloatPoint.h>
#include <maya/MFloatPointArray.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MFnDoubleIndexedComponent.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MFnMesh.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnNurbsSurface.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MFnSkinCluster.h>
#include <maya/MFnTransform.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MGlobal.h>
#include <maya/MIOStream.h>
#include <maya/MIntArray.h>
#include <maya/MItDag.h>
#include <maya/MItDependencyGraph.h>
#include <maya/MItGeometry.h>
#include <maya/MItMeshEdge.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MItMeshVertex.h>
#include <maya/MItSelectionList.h>
#include <maya/MItSurfaceCV.h>
#include <maya/MMatrix.h>
#include <maya/MPoint.h>
#include <maya/MPxCommand.h>
#include <maya/MQuaternion.h>
#include <maya/MRichSelection.h>
#include <maya/MSelectionList.h>
#include <maya/MStatus.h>
#include <maya/MString.h>
#include <maya/MSyntax.h>
#include <maya/MTransformationMatrix.h>
#include <maya/MVector.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <unordered_set>
#include <vector>

class blurSkinCmd : public MPxCommand {
   public:
    blurSkinCmd();
    virtual ~blurSkinCmd();
    enum CommandMode {
        kCommandSmooth,
        kCommandAdd,
        kCommandAbsolute,
        kCommandPercentage,
        kCommandHelp
    };

    MStatus doIt(const MArgList&);
    MStatus undoIt();
    MStatus redoIt();
    MStatus getAverageWeight(MIntArray vertices, int currentVertex);
    MStatus addWeights(int currentVertex);
    void verboseSetWeights(int currentVertex);
    void getTypeOfSurface();
    MStatus getListLockJoints();
    MStatus getAllWeights();
    MStatus useAllVertices();
    MStatus executeAction();
    MStatus printWeigth(int vertex, int u = 0, int v = 0);
    MStatus getSoftSelection();
    bool isUndoable() const;
    static void* creator();
    static MSyntax newSyntax();

    const static char* kQueryFlagShort;
    const static char* kQueryFlagLong;

    const static char* kSkinClusterNameFlagShort;
    const static char* kSkinClusterNameFlagLong;

    const static char* kMeshNameFlagShort;
    const static char* kMeshNameFlagLong;

    const static char* kIndexSkinClusterFlagShort;
    const static char* kIndexSkinClusterFlagLong;

    const static char* kPercentMovementFlagShort;
    const static char* kPercentMovementFlagLong;

    const static char* kVerboseFlagShort;
    const static char* kVerboseFlagLong;

    const static char* kListCVsIndicesFlagShort;
    const static char* kListCVsIndicesFlagLong;

    const static char* kListVerticesIndicesFlagShort;
    const static char* kListVerticesIndicesFlagLong;

    const static char* kListVerticesWeightFlagShort;
    const static char* kListVerticesWeightFlagLong;

    const static char* kListJointsFlagShort;
    const static char* kListJointsFlagLong;

    const static char* kListJointsValuesFlagShort;
    const static char* kListJointsValuesFlagLong;

    const static char* kRepeatFlagShort;
    const static char* kRepeatFlagLong;

    const static char* kDepthFlagShort;
    const static char* kDepthFlagLong;

    const static char* kRespectLocksFlagShort;
    const static char* kRespectLocksFlagLong;

    const static char* kCommandFlagShort;
    const static char* kCommandFlagLong;

    /**
    Displays help.
    */
    const static char* kHelpFlagShort;
    const static char* kHelpFlagLong;

   private:
    MStatus GatherCommandArguments(const MArgList& args);

    // MFnMesh meshFn; // the mesh
    MDagPath meshPath_;    // the mesh
    MObject skinCluster_;  // the skinCluster
    MObject component;     // the components vertices
    int nbJoints;

    MDoubleArray fullOrigWeights, weigthsForUndo, currentWeights, newWeights, weightsForSetting;
    MIntArray lockJoints, lockVertices;

    MIntArray indicesVertices_, indicesU_, indicesV_;
    MFloatArray weightVertices_;

    MStringArray listJoints_;
    MDoubleArray listJointsValues_;
    MDoubleArray perJointAddingValues_;

    MSelectionList selectionList_; /**< Selected command input nodes. */
    CommandMode command_;          // the command type
    bool verbose, respectLocks_;
    bool isNurbsSurface_, isMeshSurface_, isNurbsCurve_, isBezierCurve_;
    bool useSelection = false;
    int depth_, repeat_, indSkinCluster_;
    double percentMvt_;

    int numCVsInV_, numCVsInU_;
    int UDeg_, VDeg_;
    bool UIsPeriodic_, VIsPeriodic_;
};
