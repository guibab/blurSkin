#include "blurSkinEdit.h"

#include "functions.h"

MTypeId blurSkinDisplay::id(0x001226F9);

MObject blurSkinDisplay::_inMesh;
MObject blurSkinDisplay::_outMesh;
MObject blurSkinDisplay::_paintableAttr;
MObject blurSkinDisplay::_s_per_joint_weights;
MObject blurSkinDisplay::_s_skin_weights;

blurSkinDisplay::blurSkinDisplay() {}
blurSkinDisplay::~blurSkinDisplay() {}
// The compute() method does the actual work of the node using the inputs
MStatus blurSkinDisplay::compute(const MPlug& plug, MDataBlock& dataBlock) {
    MStatus status;
    if (plug.attribute() == blurSkinDisplay::_outMesh) {
        MDataHandle inMeshData = dataBlock.inputValue(blurSkinDisplay::_inMesh);
        MDataHandle outMeshData = dataBlock.outputValue(blurSkinDisplay::_outMesh);

        outMeshData.copy(inMeshData);

        MFnDagNode meshNodeFn(thisMObject());
        MDagPath meshPath;

        if (skinCluster_ == MObject::kNullObj) {
            getConnectedSkinCluster();                              // get the skinCluster
            getListColorsJoints(skinCluster_, this->jointsColors);  // get the joints colors
            status = fillArrayValues(true);  // get the skin data and all the colors
            set_skinning_weights(dataBlock);
        }

        if (skinCluster_ != MObject::kNullObj) {
            // 1 get the colors
            meshNodeFn.getPath(meshPath);
            MFnMesh meshFn(outMeshData.asMesh());
            int nbVertices = meshFn.numVertices();
            if (this->vertexIndices.length() != nbVertices) {
                this->vertexIndices.setLength(nbVertices);
                for (unsigned int i = 0; i < nbVertices; i++) {
                    this->vertexIndices[i] = i;
                }
                MGlobal::displayInfo(" set COLORS ");  // beginning opening of node
                meshFn.setVertexColors(this->currColors, this->vertexIndices);
            }
            ////////////////////////////////////////////////////////////////
            // consider painted attribute --------
            /////////////////////////////////////////////////////////////////
            MColor white(1, 1, 1);
            MColorArray theEditColors;
            theEditColors.copy(this->currColors);
            // MIntArray theEditVerts;

            MFnDoubleArrayData arrayData;
            MObject dataObj = dataBlock.inputValue(_paintableAttr).data();
            arrayData.setObject(dataObj);

            unsigned int length = arrayData.length();
            std::vector<double> val(length);

            for (unsigned int i = 0; i < length; i++) {
                double val = arrayData[i];
                if (val > 0) {
                    // MGlobal::displayInfo(MString(" paint value ") + i + MString(" - ") + val);
                    MColor newCol = val * white + (1.0 - val) * this->currColors[i];
                    theEditColors[i] = newCol;
                    // theEditColors.append(newCol);
                    // theEditVerts.append(i);
                }
            }
            meshFn.setVertexColors(theEditColors, this->vertexIndices);
        }
    }

    //}else if (plug.attribute() == blurSkinDisplay::_s_skin_weights) {

    //}
    dataBlock.setClean(plug);

    return status;
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
MStatus blurSkinDisplay::fillArrayValues(bool doColors) {
    MStatus status = MS::kSuccess;

    MFnDependencyNode skinClusterDep(skinCluster_);

    MPlug weight_list_plug = skinClusterDep.findPlug("weightList");
    // MGlobal::displayInfo(weight_list_plug.name());
    int nbElements = weight_list_plug.numElements();

    skin_weights_.resize(nbElements);
    if (doColors) {
        currColors.clear();
        currColors.setLength(nbElements);
    }
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
    MGlobal::displayInfo(" FILLED ARRAY VALUES ");

    return status;
}

void blurSkinDisplay::set_skinning_weights(MDataBlock& block) {
    MStatus status = MS::kSuccess;
    MArrayDataHandle array_hdl = block.outputArrayValue(_s_skin_weights, &status);
    MArrayDataBuilder array_builder = array_hdl.builder(&status);

    int nbVerts = skin_weights_.size();
    // array_builder.growArray(nbVerts);
    for (unsigned i = 0; i < nbVerts; i++) {
        auto vertexWeight = skin_weights_[i];
        int nbInfluences = vertexWeight.size();

        // MArrayDataHandle infLuenceWeight_hdl = array_builder.addElementArray(i);
        // MArrayDataBuilder weight_list_builder = infLuenceWeight_hdl.builder(&status);

        MDataHandle element_hdl = array_builder.addElement(i, &status);  // weightList[i]
        MDataHandle child = element_hdl.child(_s_per_joint_weights);     // weightList[i].weight

        MArrayDataHandle weight_list_hdl(child, &status);
        MArrayDataBuilder weight_list_builder = weight_list_hdl.builder(&status);

        // weight_list_builder.growArray(nbInfluences);		no need we do addElement

        // Scan array, update existing element, remove unsused ones
        for (unsigned j = 0; j < nbInfluences; ++j) {
            auto myPair = skin_weights_[i][j];
            int index = myPair.first;
            float val = myPair.second;

            MDataHandle hdl = weight_list_builder.addElement(index, &status);
            hdl.setDouble((double)val);
        }
        weight_list_hdl.set(weight_list_builder);
    }
    array_hdl.set(array_builder);
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

    blurSkinDisplay::_inMesh = meshAttr.create("inMesh", "im", MFnMeshData::kMesh, &status);
    meshAttr.setStorable(false);
    meshAttr.setConnectable(true);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_inMesh);

    // mesh output

    blurSkinDisplay::_outMesh = meshAttr.create("outMesh", "om", MFnMeshData::kMesh, &status);
    meshAttr.setStorable(false);
    meshAttr.setConnectable(true);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_outMesh);

    ///////////////////////////////////////////////////////////////////////////
    // paintable attributes
    ///////////////////////////////////////////////////////////////////////////
    blurSkinDisplay::_paintableAttr =
        tAttr.create("paintAttr", "pa", MFnData::kDoubleArray, &status);
    meshAttr.setStorable(true);
    status = blurSkinDisplay::addAttribute(blurSkinDisplay::_paintableAttr);

    ///////////////////////////////////////////////////////////////////////////
    // Initialize skin weights multi attributes
    ///////////////////////////////////////////////////////////////////////////
    MFnNumericAttribute numAtt;
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

    ///////////////////////////////////////////////////////////////////////////
    // attributeAffects
    ///////////////////////////////////////////////////////////////////////////
    attributeAffects(blurSkinDisplay::_inMesh, blurSkinDisplay::_outMesh);
    attributeAffects(blurSkinDisplay::_paintableAttr, blurSkinDisplay::_s_skin_weights);
    attributeAffects(blurSkinDisplay::_paintableAttr, blurSkinDisplay::_outMesh);
    return status;

    MGlobal::executeCommand("makePaintable -attrType doubleArray blurSkinDisplay paintAttr");
}
/*
MStatus blurSkinDisplay::setDependentsDirty( const MPlug &plugBeingDirtied,
        MPlugArray &affectedPlugs )
{
    MStatus status;
    MObject thisNode = thisMObject();
    MFnDependencyNode fnThisNode( thisNode );
    if ( plugBeingDirtied.partialName() == "A" ) {
        // "A" is dirty, so mark "B" dirty if "B" exists.
        // This implements the relationship "A blurSkinDisplay B".
        //
        fprintf(stderr,"blurSkinDisplay::setDependentsDirty, \"A\" being dirtied\n");
        MPlug pB = fnThisNode.findPlug( "B", &status );
        if ( MStatus::kSuccess == status ) {
            fprintf(stderr,"\t\t... dirtying \"B\"\n");
            CHECK_MSTATUS( affectedPlugs.append( pB ) );
        }
    }
    return( MS::kSuccess );
}
// These methods load and unload the plugin, registerNode registers the
// new node type with maya
//
*/