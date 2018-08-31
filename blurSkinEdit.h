#include <math.h>
#include <maya/MArrayDataBuilder.h>
#include <maya/MColorArray.h>
#include <maya/MDagPath.h>
#include <maya/MDagPathArray.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MFnSkinCluster.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MIOStream.h>
#include <maya/MIntArray.h>
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
    // void displayLayerWeights(const SkinLayer &layer);
    void getConnectedSkinCluster();
    MObject skinCluster_;
    MColorArray currColors, jointsColors;
    MIntArray vertexIndices;
    int init = 0;
    // void resizeVertexIndexes(const unsigned int newSize);
    std::vector<std::vector<std::pair<int, float>>> skin_weights_;

    // std::vector<  std::vector<  std::pair<  int , float  > > > skin_weights_;
    MStatus blurSkinDisplay::fillArrayValues(bool doColors = false);

   public:
    blurSkinDisplay();
    virtual ~blurSkinDisplay();
    virtual MStatus compute(const MPlug& plug, MDataBlock& dataBlock);
    // virtual MStatus     setDependentsDirty( const MPlug& plugBeingDirtied,MPlugArray
    // &affectedPlugs );
    virtual MPlug passThroughToOne(const MPlug& plug) const;

    void set_skinning_weights(MDataBlock& block);

    static void* creator();
    static MStatus initialize();
    static MTypeId id;

    static MObject _inMesh;
    static MObject _outMesh;
    static MObject _paintableAttr;
    static MObject _s_per_joint_weights;
    static MObject _s_skin_weights;
    // http://rodolphe-vaillant.fr/?e=79
};
