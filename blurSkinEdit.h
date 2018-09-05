#include <math.h>
#include <maya/MArrayDataBuilder.h>
#include <maya/MColorArray.h>
#include <maya/MDagPath.h>
#include <maya/MDagPathArray.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MDoubleArray.h>
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

#include <map>
#include <unordered_map>

class blurSkinDisplay : public MPxNode {
   private:
    bool verbose = false;
    bool init = true;

    // void displayLayerWeights(const SkinLayer &layer);
    MObject skinCluster_;
    MColorArray currColors, jointsColors;
    MDoubleArray paintedValues;
    std::vector<MIntArray> connectedVertices;  // use by MItMeshVertex getConnectedVertices
    std::vector<MIntArray> connectedFaces;     // use by MItMeshVertex getConnectedFaces
    void getConnectedVertices(MObject& outMesh, int nbVertices);

    MIntArray lockJoints, lockVertices;
    bool applyPaint = false, clearTheArray = false, reloadCommand = true, postSetting = true;
    bool callUndo = false;
    int influenceIndex = 0, commandIndex = 0;
    int nbJoints = 0;

    // void resizeVertexIndexes(const unsigned int newSize);
    std::vector<std::vector<std::pair<int, float>>> skin_weights_;

    MDoubleArray skinWeightList;

    MStatus blurSkinDisplay::fillArrayValues(bool doColors = false);
    void getConnectedSkinCluster();
    MStatus getAttributes(MDataBlock& dataBlock);
    MStatus applyCommand(MDataBlock& dataBlock, MIntArray& theEditVerts,
                         MDoubleArray& verticesWeight);
    MStatus refreshColors(MIntArray& theVerts);

   public:
    blurSkinDisplay();
    virtual ~blurSkinDisplay();
    virtual MStatus compute(const MPlug& plug, MDataBlock& dataBlock);
    virtual MStatus setDependentsDirty(const MPlug& plugBeingDirtied, MPlugArray& affectedPlugs);
    virtual MPlug passThroughToOne(const MPlug& plug) const;

    void set_skinning_weights(MDataBlock& block);
    void replace_weights(MDataBlock& block, MIntArray& theVertices, MDoubleArray& theWeights);
    static void* creator();
    static MStatus initialize();
    static MTypeId id;

    static MObject _inMesh;
    static MObject _outMesh;
    static MObject _paintableAttr;
    static MObject _clearArray;
    static MObject _callUndo;
    static MObject _postSetting;
    static MObject _commandAttr;
    static MObject _influenceAttr;

    static MObject _s_per_joint_weights;
    static MObject _s_skin_weights;
    // http://rodolphe-vaillant.fr/?e=79
};
