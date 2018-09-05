"""
pth = r'C:\Users\guillaume\Documents\DEV\Maya\cpp\blurSkin\testSetColors.py'
execfile (pth, globals(),  globals())
setColorsOnJoints ()
#setColorsOnSel ()
bsd = addColorNode () 
cmds.setAttr( bsd+".influenceIndex", 10)


enterPaint (bsd)
deleteTheJobs (toSearch = "function CallAfterPaint")

"""

from maya import cmds, mel
from functools import partial
#from dcc.maya import createMelProcedure
"""
inflColors = cmds.listConnections ("skinCluster1.influenceColor", s=True, d=False, p=True, c=True)
inflColors[1]
for jnt, inflColor in zip (inflColors[1::2],  inflColors[0::2]):    
    cmds.disconnectAttr (jnt, inflColor )
    cmds.connectAttr (jnt.replace("wireColorRGB", "objectColorRGB"), inflColor )
    #cmds.connectAttr (jnt.replace("objectColorRGB", "wireColorRGB"), inflColor )
"""

mel.eval("source artAttrCreateMenuItems.mel")
def setColorsOnJoints ():
    _colors = []
    for i in xrange (1,9) :
        col = cmds.displayRGBColor( "userDefined{0}".format (i), q=True)
        _colors.append (col)
    
    for jnt in cmds.ls (type = "joint") :
        theInd = cmds.getAttr (jnt+".objectColor")
        cmds.setAttr (jnt+".wireColorRGB",*_colors[theInd] )
        for destConn in cmds.listConnections (jnt+".objectColorRGB", d=True, s=False, p=True, type = "skinCluster") or []:
            cmds.connectAttr (jnt+".wireColorRGB", destConn , f=True)
        
def setColorsOnSel ():
    cmds.loadPlugin("blurSkin")
    sel = cmds.ls (sl=True, tr=True)
    msh = cmds.listRelatives (sel, type = "mesh")
    cmds.blurSkinCmd(command="colors",meshName = msh[0], verbose = False)

def addColorNode () :
    cmds.loadPlugin("blurSkin")
    sel = cmds.ls (sl=True, tr=True)
    msh = cmds.listRelatives (sel, type = "mesh")   
    cmds.setAttr (msh [0]+".displayColors", True)
    
    hist = cmds.listHistory (sel ,lv=0,pruneDagObjects=True)       
    if hist : 
        skinClusters = cmds.ls (hist , type="skinCluster")
        if skinClusters :         
            skinCluster = skinClusters [0]
            skinConn, inConn = cmds.listConnections (skinCluster+".input[0].inputGeometry", s=True, d=False, p=True, c=True, scn=False)
            
            bsd = cmds.createNode ("blurSkinDisplay")
            
            cmds.connectAttr (inConn, bsd+".inMesh", f=True)
            cmds.connectAttr (bsd+".outMesh", skinConn, f=True)

            cmds.evalDeferred  (partial (cmds.connectAttr, bsd+".weightList", skinCluster+".weightList", f=True))
    return bsd

def enterPaint (bsd) : 
    deleteTheJobs (toSearch = "function CallAfterPaint")
    nbAtt = cmds.getAttr (bsd+".wl", size=True)
    val = [0]*nbAtt 
    cmds.setAttr (bsd+".paintAttr", val, type = "doubleArray")
    cmds.makePaintable( "blurSkinDisplay", "paintAttr")
    cmds.makePaintable( bsd, "paintAttr")
    
    msh,=cmds.ls (cmds.listHistory (bsd,af=True, f=True), type="mesh")
    prt,=cmds.listRelatives (msh, p=True, path=True)
    
    cmds.select (prt)
    mel.eval ( "artSetToolAndSelectAttr( \"artAttrCtx\", \"{0}.paintAttr\" );".format (bsd) );
    cmds.ArtPaintAttrTool ()

    #fcProc = createMelProcedure(finalPaintBrush, [('int','slot')])
    #import __main__
    #__main__.applyCallBack = True
    createScriptJob ()
    cmds.artAttrCtx( cmds.currentCtx(), edit=True, outline=True, colorfeedback = False,
    clamp = "both", clamplower = 0.0, clampupper=1.0)#, afterStrokeCmd='print "PAINT"')

def createScriptJob ():
    theJob = cmds.scriptJob( runOnce=False, attributeChange=[bsd+'.paintAttr', CallAfterPaint] )

def deleteTheJobs (toSearch = "function CallAfterPaint") : 
    res = cmds.scriptJob(listJobs=True)
    for job in res :
        if toSearch in job : 
            jobIndex = int(job.split(":")[0])
            cmds.scriptJob(kill = jobIndex)

"""
class toggleBlockPaintBrush(object):
    def __init__(self, theVariable = ""):
        self.theAttr = theVariable

    def __enter__(self):
        import __main__
        __main__.__dict__[self.theAttr] = False

    def __exit__(self, exc_type, exc_val, exc_tb):
        import __main__
        __main__.__dict__[self.theAttr] = True
"""
def CallAfterPaint ():
    #print "-- painting post --"
    currContext = cmds.currentCtx()
    if currContext == 'artAttrContext' :
        gArtAttrCurrentAttr= mel.eval('$tmp = $gArtAttrCurrentAttr')
        typeOfNode,node,attr= gArtAttrCurrentAttr.split(".")

        arrayValues = cmds.getAttr (node+'.'+attr)
        # now set the command

        doSetCommand = False
        for ind,val in enumerate (arrayValues) : 
            if val >0.0:
                doSetCommand = True
                break
        if doSetCommand : 
            zeroValues = [0]*len (arrayValues)
            cmds.setAttr ( node+'.'+attr, zeroValues, type = "doubleArray" )
            cmds.evalDeferred (partial (cmds.setAttr ,bsd+".clearArray", 1))

def finalPaintBrush(slot):    
    print "FINAL Brush"

def clearPaint (bsd):
    nbAtt = cmds.getAttr (bsd+".wl", size=True)
    val = [0]*nbAtt 
    cmds.setAttr (bsd+".paintAttr", val, type = "doubleArray")

    """
    currContext = cmds.currentCtx()
    val = cmds.artAttrCtx( currContext, query=True,value =True )        
    cmds.artAttrCtx( currContext, edit=True,value =0.0)
    cmds.artAttrCtx( currContext, edit=True, clear=True )
    cmds.artAttrCtx( currContext, edit=True,value =val)

    """
    #cmds.setAttr (bsd+".clearArray", 1)
    cmds.evalDeferred (partial (cmds.setAttr ,bsd+".clearArray", 1))

"""
bsd = "blurSkinDisplay1"




"""

"""

if not cmds.attributeQuery ("paintAttr2", node = bsd, exists = True):
    cmds.addAttr( bsd , longName="paintAttr2", dataType="doubleArray")
    cmds.makePaintable( "blurSkinDisplay", "paintAttr2")

mel.eval ( "artSetToolAndSelectAttr( \"artAttrCtx\", \"{0}.paintAttr2\" );".format (bsd) );
cmds.ArtPaintAttrTool ()


"""


"""

cmds.setAttr ("blurSkinDisplay1.weightList[7146].weights[0]",1)


conns = cmds.listConnections ("blurSkinDisplay1",p=True, c=True)
cmds.connectAttr (conns [1], conns [3], f=True)
cmds.disconnectAttr ( conns [1], conns [0])

skinConn, inConn = cmds.listConnections ("skinCluster1.input[0].inputGeometry", s=True, d=False, p=True, c=True, scn=False)
cmds.connectAttr (inConn, "blurSkinDisplay1.inMesh", f=True)
cmds.connectAttr ("blurSkinDisplay1.outMesh", skinConn, f=True)

cmds.connectAttr ("blurSkinDisplay1.weightList", "skinCluster1.weightList", f=True)





BSD_indices = cmds.getAttr ("blurSkinDisplay1.weightList", mi=True)
skin_indices = cmds.getAttr ("skinCluster1.weightList", mi=True)

theInd = 30052
BSDW_indices = cmds.getAttr ("blurSkinDisplay1.weightList[{0}].weights".format(theInd), mi=True)
weights_indices = cmds.getAttr ("skinCluster1.weightList[{0}].weights".format(theInd), mi=True)

theInfluenceIndex = BSDW_indices [6]

BSDW_value = cmds.getAttr ("blurSkinDisplay1.weightList[{0}].weights[{1}]".format(theInd, theInfluenceIndex))
weights_value = cmds.getAttr ("skinCluster1.weightList[{0}].weights[{1}]".format(theInd, theInfluenceIndex))
"""