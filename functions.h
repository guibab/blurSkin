#ifndef _functions_h

#define _functions_h

// MAYA HEADER FILES:

#include <maya/MColorArray.h>
#include <maya/MDagPath.h>
#include <maya/MDagPathArray.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MFnMesh.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MFnSkinCluster.h>
#include <maya/MGlobal.h>
#include <maya/MItDependencyGraph.h>
#include <maya/MObject.h>
#include <maya/MObjectArray.h>
#include <maya/MPlug.h>

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <unordered_set>
#include <vector>

// FUNCTION DECLARATION:
unsigned int getMIntArrayIndex(MIntArray& myArray, int searching);
void CVsAround(int storedU, int storedV, int numCVsInU, int numCVsInV, bool UIsPeriodic,
               bool VIsPeriodic, MIntArray& vertices);
MStatus findSkinCluster(MDagPath MeshPath, MObject& theSkinCluster, int indSkinCluster,
                        bool verbose);
MStatus findMesh(MObject& theSkinCluster, MDagPath& theMeshPath, bool verbose);
MStatus findOrigMesh(MObject& theSkinCluster, MObject& origMesh, bool verbose);
MStatus getListColors(MObject& skinCluster, int nbVertices, MColorArray& currColors,
                      bool useMPlug = false);
MStatus getListColorsJoints(MObject& skinCluster, MColorArray& jointsColors);
MStatus getListLockJoints(MObject& skinCluster, MIntArray& jointsLocks);
MStatus getListLockVertices(MObject& skinCluster, MIntArray& vertsLocks);
MStatus getSymetryAttributes(MObject& skinCluster, MIntArray& symetryList);
MStatus getMirrorVertices(MIntArray mirrorVertices, MIntArray& theEditVerts,
                          MIntArray& theMirrorVerts, MIntArray& editAndMirrorVerts,
                          MDoubleArray& editVertsWeights, MDoubleArray& mirrorVertsWeights,
                          MDoubleArray& editAndMirrorWeights, bool doMerge = true);
MStatus editLocks(MObject& skinCluster, MIntArray& vertsToLock, bool addToLock,
                  MIntArray& vertsLocks);
MStatus editArray(int command, int influence, int nbJoints, MIntArray& lockJoints,
                  MDoubleArray& fullWeightArray, MIntArray& vertices, MDoubleArray& verticesWeight,
                  MDoubleArray& theWeights, bool normalize = true);
MStatus setAverageWeight(MIntArray& verticesAround, int currentVertex, int indexCurrVert,
                         int nbJoints, MIntArray& lockJoints, MDoubleArray& fullWeightArray,
                         MDoubleArray& theWeights);
MStatus doPruneWeight(MDoubleArray& theWeights, int nbJoints, double pruneCutWeight);
#endif