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
    MColorArray multiCurrentColors, jointsColors,
        soloCurrentColors;  // lock vertices color are not stored inside these arrays
    MIntArray deformersIndices;
    int nbJointsBig = 0;
    MColor lockVertColor = MColor(0.2, 0.2, 0.2);

    MDoubleArray soloColorsValues;
    MString fullColorSet = MString("multiColorsSet");
    MString soloColorSet = MString("soloColorsSet");
    MString noColorSet = MString("noColorsSet");

    MDoubleArray paintedValues;
    MIntArray fullVvertexList;
    std::vector<MIntArray> connectedVertices;  // use by MItMeshVertex getConnectedVertices
    std::vector<MIntArray> connectedFaces;     // use by MItMeshVertex getConnectedFaces
    std::vector<MIntArray> allVertsAround;     // used verts around
    int fullVertexListLength = 0;

    MIntArray lockJoints, lockVertices;
    bool applyPaint = false;
    bool clearTheArray = false;
    bool reloadCommand = true;
    bool postSetting = true;  // we apply paint as ssons as attr is changed
    bool postSetting_timeToStoreUndo =
        true;  // we store the undo for post setting after mouse release
    bool refreshLockWeights = false;
    bool callUndo = false;
    bool doNormalize = true;
    bool autoExpand = false;
    float minSoloColor = 0.0;
    float maxSoloColor = 1.0;
    int colorCommand = 0;      // multi
    int soloColorTypeVal = 1;  // 1 lava
    int nbAutoExpand = 3;      // autoExpand mode how many time we repeat it
    int changedColorInfluence = -1;
    bool reloadSoloColor = false;
    bool inputVerticesChanged = false;
    int influenceIndex = 0, commandIndex = 0, smoothRepeat = 3, smoothDepth = 1;
    int nbJoints = 0, nbVertices = 0;
    MIntArray cpIds;  // the ids of the vertices passed as to update skin for

    // mirro things -----
    bool mirrorIsActive = false;
    MIntArray mirrorInfluences, mirrorVertices;
    bool changeOfMirrorData = false;

    // void resizeVertexIndexes(const unsigned int newSize);
    std::vector<std::vector<std::pair<int, float>>> skin_weights_;
    std::vector<MIntArray> undoVertsIndices_;
    std::vector<MDoubleArray> undoVertsValues_;
    // std::vector< std::vector< int > > undoVertsIndices_;
    // std::vector< std::vector< double > > undoVertsValues_;

    MDoubleArray skinWeightList, fullUndoSkinWeightList;

    MStatus fillArrayValues(bool doColors = false);
    MStatus querySkinClusterValues(MIntArray& verticesIndices, MDoubleArray& theWeights,
                                   bool doColors = false);
    void getConnectedVertices(MObject& outMesh, int nbVertices);
    void refreshVertsConnection();
    void getConnectedSkinCluster();
    void connectSkinClusterWL();
    void setInfluenceColorAttr();
    MStatus getAttributes(MDataBlock& dataBlock);
    MStatus getMirrorInfos(MDataBlock& dataBlock);
    MStatus doStoreUndo(MIntArray& undoArray);
    MStatus applyCommand(MDataBlock& dataBlock, int influence, MIntArray& theEditVerts,
                         MDoubleArray& verticesWeight, bool storeUndo = true);
    MStatus applyCommandMirror(MDataBlock& dataBlock, MIntArray& theMirrorVerts,
                               MDoubleArray& verticesWeight);
    MStatus refreshColors(MIntArray& editVertsIndices, MColorArray& multiEditColors,
                          MColorArray& soloEditColors);
    MStatus editSoloColorSet(MFnMesh& meshFn);
    MColor getASoloColor(double val);

   public:
    blurSkinDisplay();
    virtual ~blurSkinDisplay();
    virtual MStatus compute(const MPlug& plug, MDataBlock& dataBlock);
    virtual MStatus setDependentsDirty(const MPlug& plugBeingDirtied, MPlugArray& affectedPlugs);
    MStatus connectionBroken(const MPlug& plug, const MPlug& otherPlug, bool asSrc);
    virtual MPlug passThroughToOne(const MPlug& plug) const;
    // MStatus postEvaluation(const MDGContext & 	context, const MEvaluationNode & 	evaluationNode,
    // PostEvaluationType 	evalType);
    // MStatus				shouldSave(const MPlug & plug, bool & ret);
    // virtual bool        doNotWrite() const;
    // void				beforeSave();
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
    static MObject _normalize;
    static MObject _postSetting;
    static MObject _commandAttr;
    static MObject _colorType;
    static MObject _soloColorType;
    static MObject _minSoloColor;
    static MObject _maxSoloColor;
    static MObject _smoothDepth;
    static MObject _smoothRepeat;
    static MObject _mirrorInfluenceArray;
    static MObject _mirrorActive;
    static MObject _influenceColor;
    static MObject _influenceAttr;
    static MObject _autoExpandAttr;
    static MObject _getLockWeights;

    static MObject _s_per_joint_weights;
    static MObject _s_skin_weights;
    // http://rodolphe-vaillant.fr/?e=79
};
