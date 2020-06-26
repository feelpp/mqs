#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.3.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
# import salome_notebook
# notebook = salome_notebook.NoteBook()
# sys.path.insert(0, r'/home/LNCMI-G/trophime/github/mqs/src/cases/biotsavart-bcs')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )

(imported, Partition_1, [Coil0, Coil1, Screen0, Air], [Rint_0, Rext_0, HP_0, BP_0, V0_0, V1_0, Rint_1, Rext_1, HP_1, BP_1, V0_1, V1_1, Rint_S0, Rext_S0, HP_S0, BP_S0, Infty], []) = geompy.ImportXAO("torusAir.xao")
geompy.addToStudy( Partition_1, 'Partition_1' )
geompy.addToStudyInFather( Partition_1, Coil0, 'Coil0' )
geompy.addToStudyInFather( Partition_1, Coil1, 'Coil1' )
geompy.addToStudyInFather( Partition_1, Screen0, 'Screen0' )
geompy.addToStudyInFather( Partition_1, Air, 'Air' )
geompy.addToStudyInFather( Partition_1, Rint_0, 'Rint_0' )
geompy.addToStudyInFather( Partition_1, Rext_0, 'Rext_0' )
geompy.addToStudyInFather( Partition_1, HP_0, 'HP_0' )
geompy.addToStudyInFather( Partition_1, BP_0, 'BP_0' )
geompy.addToStudyInFather( Partition_1, V0_0, 'V0_0' )
geompy.addToStudyInFather( Partition_1, V1_0, 'V1_0' )
geompy.addToStudyInFather( Partition_1, Rint_1, 'Rint_1' )
geompy.addToStudyInFather( Partition_1, Rext_1, 'Rext_1' )
geompy.addToStudyInFather( Partition_1, HP_1, 'HP_1' )
geompy.addToStudyInFather( Partition_1, BP_1, 'BP_1' )
geompy.addToStudyInFather( Partition_1, V0_1, 'V0_1' )
geompy.addToStudyInFather( Partition_1, V1_1, 'V1_1' )
geompy.addToStudyInFather( Partition_1, Rint_S0, 'Rint_S0' )
geompy.addToStudyInFather( Partition_1, Rext_S0, 'Rext_S0' )
geompy.addToStudyInFather( Partition_1, HP_S0, 'HP_S0' )
geompy.addToStudyInFather( Partition_1, BP_S0, 'BP_S0' )
geompy.addToStudyInFather( Partition_1, Infty, 'Infty' )

Auto_group_for_Sub_mesh_2 = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
geompy.UnionList(Auto_group_for_Sub_mesh_2, [Rint_1, Rext_1, HP_1, BP_1, V0_1, V1_1])

Auto_group_for_Sub_mesh_3 = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
geompy.UnionList(Auto_group_for_Sub_mesh_3, [Rint_S0, Rext_S0, HP_S0, BP_S0])

Auto_group_for_Sub_mesh_1 = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
geompy.UnionList(Auto_group_for_Sub_mesh_1, [Rint_0, Rext_0, HP_0, BP_0, V0_0, V1_0])
geompy.addToStudyInFather( Partition_1, Auto_group_for_Sub_mesh_2, 'Auto_group_for_Sub-mesh_2' )
geompy.addToStudyInFather( Partition_1, Auto_group_for_Sub_mesh_3, 'Auto_group_for_Sub-mesh_3' )
geompy.addToStudyInFather( Partition_1, Auto_group_for_Sub_mesh_1, 'Auto_group_for_Sub-mesh_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Partition_1)

MG_CADSurf = Mesh_1.Triangle(algo=smeshBuilder.MG_CADSurf)
MG_CADSurf_Parameters_0 = MG_CADSurf.Parameters()


# Submesh1 (coil0)
MG_CADSurf_1 = Mesh_1.Triangle(algo=smeshBuilder.MG_CADSurf,geom=Auto_group_for_Sub_mesh_1)
MG_CADSurf_Parameters_1 = MG_CADSurf_1.Parameters()
MG_CADSurf_Parameters_1.SetPhySize( 6 )
MG_CADSurf_Parameters_1.SetMinSize( 0.69282 )
MG_CADSurf_Parameters_1.SetMaxSize( 138.564 )
MG_CADSurf_Parameters_1.SetChordalError( 34.641 )
MG_CADSurf_Parameters_1.SetAngleMesh( 4 )
MG_CADSurf_Parameters_1.SetChordalError( 15 )

# Submesh2 (coil1)
MG_CADSurf_2 = Mesh_1.Triangle(algo=smeshBuilder.MG_CADSurf,geom=Auto_group_for_Sub_mesh_2)
MG_CADSurf_Parameters_2 = MG_CADSurf_2.Parameters()
MG_CADSurf_Parameters_2.SetPhySize( 6 )
MG_CADSurf_Parameters_2.SetMinSize( 0.69282 )
MG_CADSurf_Parameters_2.SetMaxSize( 138.564 )
MG_CADSurf_Parameters_2.SetAngleMesh( 4 )
MG_CADSurf_Parameters_2.SetChordalError( 10 )

#submesh3 (screen0)
MG_CADSurf_3 = Mesh_1.Triangle(algo=smeshBuilder.MG_CADSurf,geom=Auto_group_for_Sub_mesh_3)
MG_CADSurf_Parameters_3 = MG_CADSurf_3.Parameters()
MG_CADSurf_Parameters_3.SetPhySize( 2 )
MG_CADSurf_Parameters_3.SetMinSize( 0.69282 )
MG_CADSurf_Parameters_3.SetMaxSize( 138.564 )
MG_CADSurf_Parameters_3.SetAngleMesh( 4 )
MG_CADSurf_Parameters_3.SetChordalError( 15 )

# Submesh4 (air)
MG_CADSurf_4 = Mesh_1.Triangle(algo=smeshBuilder.MG_CADSurf,geom=Infty)
MG_CADSurf_Parameters_4 = MG_CADSurf_4.Parameters()
MG_CADSurf_Parameters_4.SetPhySize( 30 )
MG_CADSurf_Parameters_4.SetMinSize( 0.69282 )
MG_CADSurf_Parameters_4.SetMaxSize( 138.564 )
MG_CADSurf_Parameters_4.SetChordalError( 15 )


# Volumic mesh
MG_Tetra = Mesh_1.Tetrahedron(algo=smeshBuilder.MG_Tetra)
MG_Tetra_Parameters_1 = MG_Tetra.Parameters()
Coil0_1 = Mesh_1.GroupOnGeom(Coil0,'Coil0',SMESH.VOLUME)
Coil1_1 = Mesh_1.GroupOnGeom(Coil1,'Coil1',SMESH.VOLUME)
Screen0_1 = Mesh_1.GroupOnGeom(Screen0,'Screen0',SMESH.VOLUME)
Air_1 = Mesh_1.GroupOnGeom(Air,'Air',SMESH.VOLUME)
Rint_0_1 = Mesh_1.GroupOnGeom(Rint_0,'Rint_0',SMESH.FACE)
Rext_0_1 = Mesh_1.GroupOnGeom(Rext_0,'Rext_0',SMESH.FACE)
HP_0_1 = Mesh_1.GroupOnGeom(HP_0,'HP_0',SMESH.FACE)
BP_0_1 = Mesh_1.GroupOnGeom(BP_0,'BP_0',SMESH.FACE)
V0_0_1 = Mesh_1.GroupOnGeom(V0_0,'V0_0',SMESH.FACE)
V1_0_1 = Mesh_1.GroupOnGeom(V1_0,'V1_0',SMESH.FACE)
Rint_1_1 = Mesh_1.GroupOnGeom(Rint_1,'Rint_1',SMESH.FACE)
Rext_1_1 = Mesh_1.GroupOnGeom(Rext_1,'Rext_1',SMESH.FACE)
HP_1_1 = Mesh_1.GroupOnGeom(HP_1,'HP_1',SMESH.FACE)
BP_1_1 = Mesh_1.GroupOnGeom(BP_1,'BP_1',SMESH.FACE)
V0_1_1 = Mesh_1.GroupOnGeom(V0_1,'V0_1',SMESH.FACE)
V1_1_1 = Mesh_1.GroupOnGeom(V1_1,'V1_1',SMESH.FACE)
Rint_S0_1 = Mesh_1.GroupOnGeom(Rint_S0,'Rint_S0',SMESH.FACE)
Rext_S0_1 = Mesh_1.GroupOnGeom(Rext_S0,'Rext_S0',SMESH.FACE)
HP_S0_1 = Mesh_1.GroupOnGeom(HP_S0,'HP_S0',SMESH.FACE)
BP_S0_1 = Mesh_1.GroupOnGeom(BP_S0,'BP_S0',SMESH.FACE)
Infty_1 = Mesh_1.GroupOnGeom(Infty,'Infty',SMESH.FACE)

isDone = Mesh_1.Compute()
[ Coil0_1, Coil1_1, Screen0_1, Air_1, Rint_0_1, Rext_0_1, HP_0_1, BP_0_1, V0_0_1, V1_0_1, Rint_1_1, Rext_1_1, HP_1_1, BP_1_1, V0_1_1, V1_1_1, Rint_S0_1, Rext_S0_1, HP_S0_1, BP_S0_1, Infty_1 ] = Mesh_1.GetGroups()
smesh.SetName(Mesh_1, 'Mesh_1')
try:
  Mesh_1.ExportMED('torusAir.med',auto_groups=0,minor=40,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')
Sub_mesh_1 = MG_CADSurf_1.GetSubMesh()
Sub_mesh_2 = MG_CADSurf_2.GetSubMesh()
Sub_mesh_3 = MG_CADSurf_3.GetSubMesh()
Sub_mesh_4 = MG_CADSurf_4.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(MG_CADSurf.GetAlgorithm(), 'MG-CADSurf')
smesh.SetName(Rext_1_1, 'Rext_1')
smesh.SetName(HP_1_1, 'HP_1')
smesh.SetName(MG_Tetra.GetAlgorithm(), 'MG-Tetra')
smesh.SetName(MG_Tetra_Parameters_1, 'MG-Tetra Parameters_1')
smesh.SetName(MG_CADSurf_Parameters_1, 'MG-CADSurf Parameters_1')
smesh.SetName(MG_CADSurf_Parameters_2, 'MG-CADSurf Parameters_2')
smesh.SetName(MG_CADSurf_Parameters_3, 'MG-CADSurf Parameters_3')
smesh.SetName(MG_CADSurf_Parameters_4, 'MG-CADSurf Parameters_4')
smesh.SetName(Rint_0_1, 'Rint_0')
smesh.SetName(Rext_0_1, 'Rext_0')
smesh.SetName(HP_0_1, 'HP_0')
smesh.SetName(BP_0_1, 'BP_0')
smesh.SetName(Sub_mesh_4, 'Sub-mesh_4')
smesh.SetName(V0_0_1, 'V0_0')
smesh.SetName(V1_0_1, 'V1_0')
smesh.SetName(Rint_1_1, 'Rint_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Sub_mesh_3, 'Sub-mesh_3')
smesh.SetName(Sub_mesh_2, 'Sub-mesh_2')
smesh.SetName(V0_1_1, 'V0_1')
smesh.SetName(BP_1_1, 'BP_1')
smesh.SetName(Rint_S0_1, 'Rint_S0')
smesh.SetName(V1_1_1, 'V1_1')
smesh.SetName(Air_1, 'Air')
smesh.SetName(HP_S0_1, 'HP_S0')
smesh.SetName(Screen0_1, 'Screen0')
smesh.SetName(Rext_S0_1, 'Rext_S0')
smesh.SetName(Coil1_1, 'Coil1')
smesh.SetName(Infty_1, 'Infty')
smesh.SetName(Coil0_1, 'Coil0')
smesh.SetName(BP_S0_1, 'BP_S0')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
