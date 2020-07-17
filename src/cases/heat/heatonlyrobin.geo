// parametres
h = 0.025;
L =1;
Lz = 5*L;

// Definition des points: stem:[(x, y, z)] et taille de l'element /
Point (1)={0,0,0,h};
Point (2)={L,0,0,h};
Point (3)={L,L,0,h};
Point (4)={0,L,0,h};
// Definition des lignes
Line (1)={1,2};
Line (2)={2,3};
Line (3)={3,4};
Line (4)={4,1};
// Definition de la surface
Curve Loop (1)={1,2,3,4} ;
Plane Surface (1)={1};

//Physical Curve ("Gamma") = {1,2,3,4};
Physical Surface ("Dirichlet", 1)={1};

out[] = Extrude {0, 0, Lz} {
  Surface{1}; 
};

Printf ("out[0]=%g", out[0]);
Printf ("out[2]=%g", out[2]);
Printf ("out[3]=%g", out[3]);
Printf ("out[4]=%g", out[4]);

Physical Volume ("Omega", out[1])={out[1]};
Physical Surface ("Robin", out[0])={out[0]};
Physical Surface ("Dirichlet", out[2])+={out[2]};
Physical Surface ("Dirichlet", out[3])+={out[3]};
Physical Surface ("Dirichlet", out[4])+={out[4]};
Physical Surface ("Dirichlet", out[5])+={out[5]};