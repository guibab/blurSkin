

#ifndef _functions_h

#define _functions_h

// MAYA HEADER FILES:

#include <maya/MDagPath.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnMesh.h>
#include <maya/MFnSkinCluster.h>
#include <maya/MGlobal.h>
#include <maya/MItDependencyGraph.h>
#include <maya/MObject.h>
#include <maya/MObjectArray.h>
#include <maya/MPlug.h>

#include <vector>

// FUNCTION DECLARATION:

MStatus findSkinCluster(MDagPath MeshPath, MObject& theSkinCluster, int indSkinCluster,
                        bool verbose);
MStatus findMesh(MObject& theSkinCluster, MDagPath& theMeshPath, bool verbose);
#endif