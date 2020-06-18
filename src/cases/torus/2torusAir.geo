// Define Main params
Unit=1.e-3;

h = 5*Unit;

r1=75.*Unit;
r2=100.2*Unit;
he=50/2.*Unit;
eps=1*Unit;
theta1=Asin( eps/(2*r1) );
theta2=Asin( eps/(2*r2) );

r_inf = 1.2*10*r2; // need to be adapted to NbTurns and z0
h_inf = 10*h;
IsAir=1;

// Define Torus with a slit
Macro torus
  P1= newp; Point(P1) = {0, 0, z0-he, h};

  P2= newp; Point(P2) = {r1*Cos(theta1), eps/2., z0-he, h};
  P3= newp; Point(P3) = {r2*Cos(theta2), eps/2., z0-he, h};
  P4= newp; Point(P4) = {0, r1, z0-he, h};
  P5= newp; Point(P5) = {0, r2, z0-he, h};
  P6= newp; Point(P6) = {-r1, 0, z0-he, h};
  P7= newp; Point(P7) = {-r2, 0, z0-he, h};
  P8= newp; Point(P8) = {0, -r1, z0-he, h};
  P9= newp; Point(P9) = {0, -r2, z0-he, h};
  P10= newp; Point(P10) = {r1*Cos(-theta1), -eps/2., z0-he, h};
  P11= newp; Point(P11) = {r2*Cos(-theta2), -eps/2., z0-he, h};

  C1= newl; Circle(C1) = {P2, P1, P4};
  C2= newl; Circle(C2) = {P4, P1, P6};
  C3= newl; Circle(C3) = {P6, P1, P8};
  C4= newl; Circle(C4) = {P8, P1, P10};

  C5=newl; Circle(C5) = {P3, P1, P5};
  C6=newl; Circle(C6) = {P5, P1, P7};
  C7=newl; Circle(C7) = {P7, P1, P9};
  C8=newl; Circle(C8) = {P9, P1, P11};

  L9=newl; Line(L9) = {P2, P3};
  L10=newl; Line(L10) = {P10, P11};

  dL=newl; Line Loop(dL) = {C1:C4, L10, -C8, -C7, -C6, -C5, -L9};
  S=news; Plane Surface(S) = {dL};

  out[] = Extrude {0,0,2*he} {Surface{S};};

  Physical Volume(Sprintf("coil%g", t)) = {out[1]};
  Physical Surface(Sprintf("BP_%g", t)) = {out[0]};
  Physical Surface(Sprintf("HP_%g", t)) = {S};
  Physical Surface(Sprintf("Rint_%g", t)) = {out[2], out[3], out[4], out[5]};
  Physical Surface(Sprintf("Rext_%g", t)) = {out[7], out[8], out[9], out[10]};
  Physical Surface(Sprintf("V0_%g", t)) = {out[6]};
  Physical Surface(Sprintf("V1_%g", t)) = {out[11]};

  Bord[t]=news; Surface Loop(Bord[t]) = {out[0], S, out[2], out[3], out[4], out[5],
                                            out[7], out[8], out[9], out[10],
                                            out[6], out[11]};
Return

// Create a stack of Torus
Nturn=2;
epsz=1*Unit;
z0=-(3*(2*he)+2*epsz)/2.;
dz=2*Fabs(z0);
For t In {0:Nturn-1}
  Call torus;
  z0+=dz;
EndFor
  
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
