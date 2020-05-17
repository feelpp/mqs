// Gmsh project created on Sun May 17 19:22:38 2020
//+
Point(1) = {0, -0, 0, 1.0};
//+
Point(2) = {-0, 1, 0, 1.0};
//+
Point(3) = {5, 1, 0, 1.0};
//+
Point(4) = {5, 0, 0, 1.0};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Surface(1) = {1};
//+
Physical Surface("Gamma_C", 5) = {1};
//+
Extrude {0, 0, 1} {
  Surface{1}; 
}
//+
Physical Surface("gamma_I", 28) = {18};
//+
Physical Surface("gamma_O", 29) = {26};
//+
Physical Volume("Omega_C", 30) = {1};
//+
Physical Surface("Gamma_C", 31) += {27};
//+
Physical Surface("Gamma_C", 31) += {14};
//+
Physical Surface("Gamma_C", 31) += {22};