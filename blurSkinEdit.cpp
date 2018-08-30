#include "blurSkinEdit.h"

#include "functions.h"

MTypeId blurSkinDisplay::id(0x001226F9);

MObject blurSkinDisplay::_inMesh;
MObject blurSkinDisplay::_outMesh;
MObject blurSkinDisplay::_s_per_joint_weights;
MObject blurSkinDisplay::_s_skin_weights;

blurSkinDisplay::blurSkinDisplay() {}
blurSkinDisplay::~blurSkinDisplay() {}
// The compute() method does the actual work of the node using the inputs
MStatus blurSkinDisplay::compute(const MPlug& plug, MDataBlock& dataBlock) {
    MStatus status;
    if (plug.attribute() != blurSkinDisplay::_outMesh) {
        // ignore other outputs
        return status;
    }

    MDataHandle inMeshData = dataBlock.inputValue(blurSkinDisplay::_inMesh);
    MDataHandle outMeshData = dataBlock.outputValue(blurSkinDisplay::_outMesh);

    outMeshData.copy(inMeshData);

    MFnDagNode meshNodeFn(thisMObject());
    MDagPath meshPath;

    if (skinCluster_ == MObject::kNullObj) getConnectedSkinCluster();

    if (skinCluster_ != MObject::kNullObj) {
        // 1 get the colors
        meshNodeFn.getPath(meshPath);
        MFnMesh meshFn(outMeshData.asMesh());
        int nbVertices = meshFn.numVertices();

        if (this->currColors.length() != nbVertices) {  // beginning opening of node
            MGlobal::displayInfo(" set COLORS ");
            getListColorJoints(skinCluster_, nbVertices, this->currColors, true);
        }
        /*
        MColor yellow(1, 1, 0);
        this->currColors.setLength(nbVertices);
        for (unsigned int i = 0; i < nbVertices; i++) {
                this->currColors[i] = yellow;
        }
        */
        this->vertexIndices.setLength(nbVertices);
        for (unsigned int i = 0; i < nbVertices; i++) {
            this->vertexIndices[i] = i;
        }
        meshFn.setVertexColors(this->currColors, this->vertexIndices);
    }
    dataBlock.setClean(plug);

    return status;
}
void blurSkinDisplay::getConnectedSkinCluster() {
    MStatus status;
    MPlug outMesh(thisMObject(), blurSkinDisplay::_outMesh);
    MPlugArray connections;
    outMesh.connectedTo(connections, false, true);

    // get mesh connection from outMesh
    for (unsigned int i = 0; i < connections.length(); i++) {
        MObject mesh = connections[i].node();
        if (!mesh.hasFn(MFn::kMesh) || !mesh.hasFn(MFn::kDagNode)) continue;
        MFnDagNode meshNodeFn(mesh);
        MDagPath meshPath;
        meshNodeFn.getPath(meshPath);
        status = findSkinCluster(meshPath, skinCluster_, 0, false);

        ///////////////////////////////////////////////////////////////////////////
        // now check the connection of the plug -----------------------------------
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
void blurSkinDisplay::get_skinning_weights(MDataBlock& block) {
    MStatus status = MS::kSuccess;
    MArrayDataHandle array_hdl = block.outputArrayValue(_s_skin_weights, &status);
    /*
    for (unsigned i = 0; i < weights.size(); i++)
    {
            array_hdl.jumpToArrayElement(i);

            // weightList[i]
            MDataHandle element_hdl = array_hdl.outputValue(&status);

            // weightList[i].weight
            MDataHandle child = element_hdl.child(_s_per_joint_weights);

            MArrayDataHandle weight_list_hdl(child, &status);
            mayaCheck(status);

            MArrayDataBuilder weight_list_builder = weight_list_hdl.builder(&status);
            mayaCheck(status);

            unsigned handle_count = weight_list_hdl.elementCount(&status);
            mayaCheck(status);

            unsigned builder_count = weight_list_builder.elementCount(&status);
            mayaCheck(status);
            mayaAssert(builder_count == handle_count);

            //std::map<int  influence obj id / joint id , float> map = weights[i];
            std::map<int, float> map = weights[i];

            std::vector to_remove;
            to_remove.reserve(map.size());

            // Scan array, update existing element, remove unsused ones
            for (unsigned j = 0; j < handle_count; ++j)
            {
                    // weightList[i].weight[j]
                    mayaCheck(weight_list_hdl.jumpToArrayElement(j));
                    unsigned index = weight_list_hdl.elementIndex(&status);
                    mayaCheck(status);

                    auto elt = map.find(index);

                    if (elt != map.end())
                    {
                            MDataHandle hdl = weight_list_builder.addElement(index, &status);
                            mayaCheck(status);
                            hdl.setDouble((double)elt->second);
                            map.erase(elt);
                    }
                    else
                    {
                            to_remove.push_back(index);
                    }
            }

            for (unsigned idx : to_remove) {
                    mayaCheck(weight_list_builder.removeElement(idx));
            }

            mayaCheck(weight_list_hdl.set(weight_list_builder));
    }
    */
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
    // Initialize skin weights multi attributes
    ///////////////////////////////////////////////////////////////////////////
    MFnNumericAttribute numAtt;
    _s_per_joint_weights =
        numAtt.create("per_joint_weights", "jw", MFnNumericData::kDouble, 0.0, &status);
    numAtt.setKeyable(false);
    numAtt.setArray(true);
    numAtt.setReadable(true);
    numAtt.setWritable(false);
    numAtt.setUsesArrayDataBuilder(true);
    addAttribute(_s_per_joint_weights);

    MFnCompoundAttribute cmpAttr;
    _s_skin_weights = cmpAttr.create("skin_weight_list", "sw", &status);
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
    return status;
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