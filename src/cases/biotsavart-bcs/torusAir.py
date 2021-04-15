#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.3.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
# import salome_notebook
# notebook = salome_notebook.NoteBook()
# sys.path.insert(0, r'/home/LNCMI-G/trophime')

import getopt
import sys
import os

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
# import SALOMEDS


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )


def Screen(r1, r2, h, z0):
  """
  Create a rectangular cross section torus 

  r1: inner radius
  r2: outer radius
  h: total height
  z0: altitude of the mid-plane

  The torus section is: [r1,r2] x [z0-h/2.,z0+h/2.]
  """

  Cylinder_1 = geompy.MakeCylinderRH(r2, h)
  Cylinder_2 = geompy.MakeCylinderRH(r1, 1.2*h)
  geompy.TranslateDXDYDZ(Cylinder_2, 0, 0, -0.1*h)
  Cut_1 = geompy.MakeCutList(Cylinder_1, [Cylinder_2], True)
  geompy.TranslateDXDYDZ(Cut_1, 0, 0, z0-h/2. )
  geompy.addToStudy( Cut_1, 'Screen_1' )
  return Cut_1

def Coil(r1, r2, h, eps, z0):
  """
  Create a coil of rectangular cross section
  with a slit 
  
  r1: inner radius
  r2: outer radius
  h: total height
  eps: thickness of the slit
  z0: altitude of the mid-plane

  The torus section is: [r1,r2] x [z0-h/2.,z0+h/2.]
  """

  # Create a rectangular cross section torus
  Cylinder_1 = geompy.MakeCylinderRH(r2, h)
  Cylinder_2 = geompy.MakeCylinderRH(r1, 1.2*h)
  geompy.TranslateDXDYDZ(Cylinder_2, 0, 0, -0.1*h)

  Cut_1 = geompy.MakeCutList(Cylinder_1, [Cylinder_2], True)
  # Create Input/Output
  InOut_1 = geompy.MakeBoxDXDYDZ(1.2*r2, 2*eps, 1.2*h)
  geompy.TranslateDXDYDZ(InOut_1, 0, -eps/2., 0)
  Coil = geompy.MakeCutList(Cut_1, [InOut_1], True)
  geompy.TranslateDXDYDZ(Coil, 0, 0, z0-h/2. )

  # geompy.addToStudy( Cylinder_1, 'Cylinder_1' )
  # geompy.addToStudy( Cylinder_2, 'Cylinder_2' )
  # geompy.addToStudy( Cut_1, 'Cut_1' )
  # geompy.addToStudy( InOut_1, 'InOut_1' )
  # geompy.addToStudy( Coil, 'Coil' )

  return Coil

# TODO:
# put data into a json?
# create Mesh
def main():
  """main"""

  print("sys.argv=", sys.argv)
  
  import argparse
  parser = argparse.ArgumentParser()

  parser.add_argument("--turn", metavar=('r1', 'r2', 'h', 'z0', 'eps'), nargs=5, action='append', help="define the coils")
  parser.add_argument("--screen", metavar=('r1', 'r2', 'h', 'z0'), nargs=4, action='append', help="define the screens")
  parser.add_argument("--output", type=str, help="specify output name", default='torusAir')
  parser.add_argument("--rinf", type=float, help="specify rinf")
  args = parser.parse_args()

  print("args=", args)
  coils_data = args.turn
  screens_data = args.screen
  output = args.output

  # extract coils data
  # coil: [r1, r2, h, z0, eps]
  print ("Create coils")
  n_coils = len(coils_data)
  Coils = []
  HPts = []
  BPts = []
  for coil in coils_data:
    print("coil=", coil)
    r1 = coil[0]
    r2 = coil[1]
    h = coil[2]
    z0 = coil[3]
    eps = coil[4]
    Coils.append(Coil(r1,r2,h,eps,z0))
    HPts.append( geompy.MakeVertex(0., (r1+r2)/2., z0-h/2.) );
    BPts.append( geompy.MakeVertex(0., (r1+r2)/2., z0+h/2.) );
  
  # extract screens data
  # screen: [r1, r2, h, z0]
  print("Create screens")
  n_screens = len(screens_data)
  Screens = []
  S_HPts=[]
  S_BPts=[]
  for screen in screens_data:
    print("screen=", screen)
    r1 = screen[0]
    r2 = screen[1]
    h = screen[2]
    z0 = screen[3]
    Screens.append(Screen(r1, r2, h, z0))
    S_HPts.append( geompy.MakeVertex(0., (r1+r2)/2., z0-h/2.) );
    S_BPts.append( geompy.MakeVertex(0., (r1+r2)/2., z0+h/2.) );  

  # create outer box
  rinf=args.rinf
  print ("Create Sphere")
  Sphere_1 = geompy.MakeSphereR(rinf)
  #geompy.addToStudy( Sphere_1, 'Sphere_1' )

  # Partition
  print ("Create Partition")
  Partition_1 = geompy.MakePartition(Coils+Screens+ [Sphere_1], [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
  geompy.addToStudy( Partition_1, 'Partition_1' )

  Solids = geompy.ExtractShapes(Partition_1, geompy.ShapeType["SOLID"], True)
  
  lGroups=[]
  print ("Create Solid and BC groups")
  for n,solid in enumerate(Solids):
    print ("Solid[%d]:"%n)
    isIn = geompy.AreCoordsInside(solid, [0,0,0])
    print ("isIn:", isIn)
    if isIn[0]:
      geompy.addToStudyInFather( Partition_1, solid, 'Air' )
      lFaces = geompy.GetShapesOnSphere(solid, geompy.ShapeType["FACE"], O, rinf-eps, GEOM.ST_ONOUT)
      if len(lFaces):
        print ("create Infty")
        infty = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"], "Infty")
        lIDFaces = geompy.GetSubShapesIDs(Partition_1, lFaces)
        geompy.UnionIDs(infty, lIDFaces)
        lGroups.append(infty)
    else:
    
      # Check if HPts in solid
      isScreen=True
      for ncoil in range(ncoils):
        Hcoord=geompy.PointCoordinates(HPts[ncoil])
        Bcoord=geompy.PointCoordinates(BPts[ncoil])
        print("coil%d:"%ncoil, (Hcoord[0]+Bcoord[0])/2.,(Hcoord[1]+Bcoord[1])/2.,(Hcoord[2]+Bcoord[2])/2.)

        isHPtsIn = geompy.AreCoordsInside(solid, Hcoord  )
        print ("ncoil=%d isHptsIn:" % ncoil, isHPtsIn)
        if isHPtsIn[0]:
          geompy.addToStudyInFather( Partition_1, solid, 'Coil%d' % ncoil )
          isScreen=False
        
          # Create face group
          lFaces = geompy.GetShapesOnCylinderWithLocation(solid, geompy.ShapeType["FACE"], OZ, O, r1+eps, GEOM.ST_ONIN)
          if len(lFaces):
            print ("create Rint")
            ri = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"], "Rint_%d"%ncoil)
            lIDFaces = geompy.GetSubShapesIDs(Partition_1, lFaces)
            geompy.UnionIDs(ri, lIDFaces)
            lGroups.append(ri)
          
          lFaces = geompy.GetShapesOnCylinderWithLocation(solid, geompy.ShapeType["FACE"], OZ, O, r2-eps, GEOM.ST_ONOUT)
          if len(lFaces):
            print ("create Rext")
            re = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"], "Rext_%d"%ncoil)
            lIDFaces = geompy.GetSubShapesIDs(Partition_1, lFaces)
            geompy.UnionIDs(re, lIDFaces)
            lGroups.append(re)
          
          lFaces = geompy.GetShapesOnPlaneWithLocation(solid, geompy.ShapeType["FACE"], OZ, HPts[ncoil], GEOM.ST_ON)
          if len(lFaces):
            print ("create HP")
            hp = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"], "HP_%d"%ncoil)
            lIDFaces = geompy.GetSubShapesIDs(Partition_1, lFaces)
            geompy.UnionIDs(hp, lIDFaces)
            lGroups.append(hp)
          
          lFaces = geompy.GetShapesOnPlaneWithLocation(solid, geompy.ShapeType["FACE"], OZ, BPts[ncoil], GEOM.ST_ON)
          if len(lFaces):
            print ("create BP")
            bp = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"], "BP_%d"%ncoil)
            lIDFaces = geompy.GetSubShapesIDs(Partition_1, lFaces)
            geompy.UnionIDs(bp, lIDFaces)
            lGroups.append(bp)

          lFaces = geompy.GetShapesOnPlane(solid, geompy.ShapeType["FACE"], OY, GEOM.ST_OUT)
          if len(lFaces):
            print ("create V0")
            v0 = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"], "V0_%d"%ncoil)
            lIDFaces = geompy.GetSubShapesIDs(Partition_1, lFaces)
            geompy.UnionIDs(v0, lIDFaces)
            lGroups.append(v0)
            
            print ("create V1")
            v1 = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"], "V1_%d"%ncoil)
            lFaces = geompy.GetShapesOnPlane(solid, geompy.ShapeType["FACE"], OY, GEOM.ST_IN)
            lIDFaces = geompy.GetSubShapesIDs(Partition_1, lFaces)
            geompy.UnionIDs(v1, lIDFaces)
            lGroups.append(v1)

      if isScreen:
        for nscreen in range(nscreens):
          Hcoord=geompy.PointCoordinates(S_HPts[nscreen])
          Bcoord=geompy.PointCoordinates(S_BPts[nscreen])
          print("coil%d:"%ncoil, (Hcoord[0]+Bcoord[0])/2.,(Hcoord[1]+Bcoord[1])/2.,(Hcoord[2]+Bcoord[2])/2.)
          
          isHPtsIn = geompy.AreCoordsInside(solid, Hcoord  )
          print ("ncoil=%d isHptsIn:" % nscreen, isHPtsIn)
          if isHPtsIn[0]:
            geompy.addToStudyInFather( Partition_1, solid, 'Screen%d' % nscreen )
        
            # TODO recover Rint, Rext
            lFaces = geompy.GetShapesOnCylinderWithLocation(solid, geompy.ShapeType["FACE"], OZ, O, r1_s+eps, GEOM.ST_ONIN)
            if len(lFaces):
              print ("create Rint")
              ri = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"], "Rint_S%d"%nscreen)
              lIDFaces = geompy.GetSubShapesIDs(Partition_1, lFaces)
              geompy.UnionIDs(ri, lIDFaces)
              lGroups.append(ri)
          
            lFaces = geompy.GetShapesOnCylinderWithLocation(solid, geompy.ShapeType["FACE"], OZ, O, r2_s-eps, GEOM.ST_ONOUT)
            if len(lFaces):
              print ("create Rext")
              re = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"], "Rext_S%d"%nscreen)
              lIDFaces = geompy.GetSubShapesIDs(Partition_1, lFaces)
              geompy.UnionIDs(re, lIDFaces)
              lGroups.append(re)
          
            lFaces = geompy.GetShapesOnPlaneWithLocation(solid, geompy.ShapeType["FACE"], OZ, S_HPts[nscreen], GEOM.ST_ON)
            if len(lFaces):
              print ("create HP")
              hp = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"], "HP_S%d"%nscreen)
              lIDFaces = geompy.GetSubShapesIDs(Partition_1, lFaces)
              geompy.UnionIDs(hp, lIDFaces)
              lGroups.append(hp)
          
            lFaces = geompy.GetShapesOnPlaneWithLocation(solid, geompy.ShapeType["FACE"], OZ, S_BPts[nscreen], GEOM.ST_ON)
            if len(lFaces):
              print ("create BP")
              bp = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"], "BP_S%d"%nscreen)
              lIDFaces = geompy.GetSubShapesIDs(Partition_1, lFaces)
              geompy.UnionIDs(bp, lIDFaces)
              lGroups.append(bp)
      

  geompy.ExportXAO(Partition_1, lGroups, [], "LNCMI - HIFIMAGNET", output + ".xao", output + ".brep")

  if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()
            
if __name__ == "__main__":
  main()
