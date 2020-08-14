singularity shell --nv -B /opt/DISTENE/:/opt/DISTENE:ro feelpp/research/hifimagnet/singularity/salome/salome-9.3.2-buster.sif
salome -w1 -t torusAir.py args:--turns=2,--screens=1

gmsh -0 -3 torusAir.med 

singularity shell hifimagnet-thermobox-P2_9.3.3-v0.108.simg
feelpp_mesh_partitioner --ifile torusAir.msh --ofile torusAir_p --part 4 --mesh.scale=0.001

