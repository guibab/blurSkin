from maya import cmds

"""
inflColors = cmds.listConnections ("skinCluster1.influenceColor", s=True, d=False, p=True, c=True)
inflColors[1]
for jnt, inflColor in zip (inflColors[1::2],  inflColors[0::2]):    
    cmds.disconnectAttr (jnt, inflColor )
    cmds.connectAttr (jnt.replace("wireColorRGB", "objectColorRGB"), inflColor )
    #cmds.connectAttr (jnt.replace("objectColorRGB", "wireColorRGB"), inflColor )
"""

cmds.loadPlugin("blurSkin")
def setColorsOnJoints ():
	_colors = []
	for i in xrange (1,9) :
	    col = cmds.displayRGBColor( "userDefined{0}".format (i), q=True)
	    _colors.append (col)

	for jnt in cmds.ls (type = "joint") :
	    theInd = cmds.getAttr (jnt+".objectColor")
	    cmds.setAttr (jnt+".wireColorRGB",*_colors[theInd] )
	    for destConn in cmds.listConnections (jnt+".objectColorRGB", d=True, s=False, p=True, type = "skinCluster"):
	        cmds.connectAttr (jnt+".wireColorRGB", destConn , f=True)
	    


def setColorsOnSel ():
	sel = cmds.ls (sl=True, tr=True)
	msh = cmds.listRelatives (sel, type = "mesh")
	cmds.setAttr (msh [0]+".displayColors", True)
	cmds.blurSkinCmd(command="colors",meshName = msh[0], verbose = False)


def addColorNode ()	:
	sel = cmds.ls (sl=True, tr=True)
	msh = cmds.listRelatives (sel, type = "mesh")
	cmds.setAttr (msh [0]+".displayColors", True)
		
	bsd = cmds.createNode ("blurSkinDisplay")

	inConn,=cmds.listConnections (msh[0]+".inMesh", s=True, d=False, p=True)

	cmds.connectAttr (inConn, bsd+".inMesh", f=True)
	cmds.connectAttr (bsd+".outMesh", msh[0]+".inMesh",f=True)

setColorsOnJoints ()
setColorsOnSel ()
addColorNode ()	