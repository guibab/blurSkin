#include <math.h>
#include <maya/MArrayDataBuilder.h>
#include <maya/MColorArray.h>
#include <maya/MDagPath.h>
#include <maya/MDagPathArray.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MDoubleArray.h>
#include <maya/MEvaluationNode.h>
#include <maya/MFnComponentListData.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MFnSkinCluster.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MIOStream.h>
#include <maya/MIntArray.h>
#include <maya/MItMeshVertex.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MPxNode.h>
#include <maya/MString.h>
#include <maya/MTypeId.h>
#include <maya/MVector.h>
#include <string.h>

#include <unordered_map>
#include <unordered_set>

class blurSkinDisplay : public MPxNode {
   private:
    bool verbose = false;
    bool init = true;
    bool doConnectSkinCL = false;
    // void displayLayerWeights(const SkinLayer &layer);
    MObject skinCluster_;
    MColorArray multiCurrentColors, jointsColors, soloCurrentColors;
    MIntArray soloColorsValues;
    MString fullColorSet = MString("multiColorsSet");
    MString soloColorSet = MString("soloColorsSet");

    MDoubleArray paintedValues;
    MIntArray fullVvertexList;
    std::vector<MIntArray> connectedVertices;  // use by MItMeshVertex getConnectedVertices
    std::vector<MIntArray> connectedFaces;     // use by MItMeshVertex getConnectedFaces
    std::vector<MIntArray> allVertsAround;     // used verts around
    int fullVertexListLength = 0;

    MIntArray lockJoints, lockVertices;
    bool applyPaint = false, clearTheArray = false, reloadCommand = true, postSetting = true;
    bool callUndo = false;
    int colorCommand = 0;
    bool reloadSoloColor = false;
    bool inputVerticesChanged = false;
    int influenceIndex = 0, commandIndex = 0, smoothRepeat = 3, smoothDepth = 1;
    int nbJoints = 0;

    // void resizeVertexIndexes(const unsigned int newSize);
    std::vector<std::vector<std::pair<int, float>>> skin_weights_;
    std::vector<MIntArray> undoVertsIndices_;
    std::vector<MDoubleArray> undoVertsValues_;
    // std::vector< std::vector< int > > undoVertsIndices_;
    // std::vector< std::vector< double > > undoVertsValues_;

    MDoubleArray skinWeightList;

    MStatus fillArrayValues(bool doColors = false);
    MStatus querySkinClusterValues(MIntArray& verticesIndices, MDoubleArray& theWeights,
                                   bool doColors = false);
    void getConnectedVertices(MObject& outMesh, int nbVertices);
    void refreshVertsConnection();
    void getConnectedSkinCluster();
    void connectSkinClusterWL();
    MStatus getAttributes(MDataBlock& dataBlock);
    MStatus applyCommand(MDataBlock& dataBlock, MIntArray& theEditVerts,
                         MDoubleArray& verticesWeight, bool storeUndo = true);
    MStatus refreshColors(MIntArray& editVertsIndices, MColorArray& multiEditColors,
                          MColorArray& soloEditColors);
    MStatus editSoloColorSet(MFnMesh& meshFn);

   public:
    blurSkinDisplay();
    virtual ~blurSkinDisplay();
    virtual MStatus compute(const MPlug& plug, MDataBlock& dataBlock);
    virtual MStatus setDependentsDirty(const MPlug& plugBeingDirtied, MPlugArray& affectedPlugs);
    MStatus connectionBroken(const MPlug& plug, const MPlug& otherPlug, bool asSrc);
    virtual MPlug passThroughToOne(const MPlug& plug) const;
    // MStatus postEvaluation(const MDGContext & 	context, const MEvaluationNode & 	evaluationNode,
    // PostEvaluationType 	evalType);
    void set_skinning_weights(MDataBlock& block);
    void replace_weights(MDataBlock& block, MIntArray& theVertices, MDoubleArray& theWeights);
    static void* creator();
    static MStatus initialize();
    static MTypeId id;

    static MObject _inMesh;
    static MObject _outMesh;
    static MObject _cpList;
    static MObject _paintableAttr;
    static MObject _clearArray;
    static MObject _callUndo;
    static MObject _postSetting;
    static MObject _commandAttr;
    static MObject _colorType;
    static MObject _smoothDepth;
    static MObject _smoothRepeat;
    static MObject _influenceAttr;

    static MObject _s_per_joint_weights;
    static MObject _s_skin_weights;
    // http://rodolphe-vaillant.fr/?e=79
};
