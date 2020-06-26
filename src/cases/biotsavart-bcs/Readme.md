singularity shell --nv -B /opt/DISTENE/:/opt/DISTENE:ro feelpp/research/hifimagnet/singularity/salome/salome-9.3.2-buster.sif
salome -w1 -t torusAir.py args:--turns=2,--screens=1

gmsh -0 -3 Mesh_1.med 
feelpp_mesh_partitioner --ifile Mesh_1.msh --ofile Mesh_1 --part 4 --mesh.scale=0.001

