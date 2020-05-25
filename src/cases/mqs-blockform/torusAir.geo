// Define Main params
h = 0.1;
r1=1;
r2=2;
he=2.5;
eps=1.e-3;
theta1=Asin( eps/(2*r1) );
theta2=Asin( eps/(2*r2) );

r_inf = 4*r2;
h_inf = 6 *h;
IsAir=1;

// 1st quarter
Point(1) = {0, 0, -he, h};

Point(2) = {r1*Cos(theta1), eps/2., -he, h};
Point(3) = {r2*Cos(theta2), eps/2., -he, h};
Point(4) = {0, r1, -he, h};
Point(5) = {0, r2, -he, h};
Point(6) = {-r1, 0, -he, h};
Point(7) = {-r2, 0, -he, h};
Point(8) = {0, -r1, -he, h};
Point(9) = {0, -r2, -he, h};
Point(10) = {r1*Cos(-theta1), -eps/2., -he, h};
Point(11) = {r2*Cos(-theta2), -eps/2., -he, h};

Circle(1) = {2, 1, 4};
Circle(2) = {4, 1, 6};
Circle(3) = {6, 1, 8};
Circle(4) = {8, 1, 10};

Circle(5) = {3, 1, 5};
Circle(6) = {5, 1, 7};
Circle(7) = {7, 1, 9};
Circle(8) = {9, 1, 11};

Line(9) = {2, 3};
Line(10) = {10, 11};

dL=newl; Line Loop(dL) = {1:4, 10, -8, -7, -6, -5, -9};
S=news; Plane Surface(S) = {dL};

out[] = Extrude {0,0,2*he} {Surface{S};};

Physical Volume("Omega_c") = {out[1]};
Physical Surface("top") = {out[0]};
Physical Surface("bottom") = {S};
Physical Surface("Rint") = {out[2], out[3], out[4], out[5]};
Physical Surface("Rext") = {out[7], out[8], out[9], out[10]};
Physical Surface("V0") = {out[6]};
Physical Surface("V1") = {out[11]};

Bord=newl; Surface Loop(Bord) = {out[0], S, out[2], out[3], out[4], out[5],
                                            out[7], out[8], out[9], out[10],
                                            out[6], out[11]};
// Define Air (inf section)

If ( IsAir != 0 )
  P0= newp; Point(P0) = {0,0,0, h};
  P8=newp; Point(P8) = {0,     0, -r_inf, h_inf};
  P9=newp; Point(P9) = {r_inf, 0, 0, h_inf};
  P10=newp; Point(P10) = {0,   0, r_inf, h_inf};

  // L58=newl; Line(L58) = {5, P8};
  C89=newl; Circle(C89) = {P8, P0, P9};
  C910=newl; Circle(C910) = {P9, P0, P10};
  L107=newl; Line(L107) = {P8, P10};

  L_inf=newl; Line Loop(L_inf) = {C89, C910, -L107};
  // L_inf=newl; Line Loop(L_inf) = {L58, C89, C910, L107, -C67, -C56};
  S_inf=newreg; Plane Surface(S_inf) = {L_inf};

  // Build 3D geom
  quart1[] = Extrude { {0,0,1} , {0,0,0} , Pi/2. } { 
    Surface{S_inf}; //Layers{N_Layers}; //Recombine; 
  };

  quart2[] = Extrude { {0,0,1} , {0,0,0} , Pi/2 } { 
    Surface{quart1[0]}; //Layers{N_Layers/2}; //Recombine; 
  };

  quart3[] = Extrude { {0,0,1} , {0,0,0} , Pi/2 } { 
    Surface{quart2[0]}; //Layers{N_Layers/2}; //Recombine; 
  };

  quart4[] = Extrude { {0,0,1} , {0,0,0} , Pi/2 } { 
    Surface{quart3[0]}; //Layers{N_Layers/2}; //Recombine; 
  };

  // size of an array
  Printf("quart1[2]=%g", quart1[2]);
  Printf("quart1[3]=%g", quart1[3]);
  Printf("quart2[3]=%g", quart2[3]);
  Printf("quart3[3]=%g", quart3[3]);
  Printf("quart4[3]=%g", quart4[3]);
  
  Recursive Delete {
    Surface{S_inf};
    Surface{quart1[0]}; 
    Surface{quart2[0]}; 
    Surface{quart3[0]}; 
    Surface{quart4[0]}; 
  }

  Infty=newl; Surface Loop(Infty) = {quart1[2], quart2[2], quart3[2], quart4[2],
                                     quart1[3], quart2[3], quart3[3], quart4[3]};
  Air = newv; Volume(Air) = {Infty, -Bord};
  Physical Volume("air") = {Air};
  Physical Surface("Infty") = {quart1[2], quart2[2], quart3[2], quart4[2],
                               quart1[3], quart2[3], quart3[3], quart4[3]};

  Recursive Delete {
    Volume{quart1[1]};
    Volume{quart2[1]};
    Volume{quart3[1]};
    Volume{quart4[1]};
  }

  Coherence;
EndIf
