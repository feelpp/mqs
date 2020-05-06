// parametres
h = 0.025;
L =1;
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
Line Loop (1)={1,2,3,4} ;
Plane Surface (1)={1};
Physical Surface (1)={1};