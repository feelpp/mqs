// -*- vim:filetype=gmsh
Mesh.OptimizeNetgen=0; // otherwise it crashes when IsFull=1

IsFull=0;
Unit=1.e-3;
h = 5*Unit;
h_inf = 10*h;
h_ext = 10*h;

r1=75.*Unit;
r2=100.2*Unit;
L=50/2.*Unit;

r_ext=10*r2;
r_inf=1.2*r_ext;

// Define Center
P0= newp; Point(P0) = {0,0,0, h};

Section[]={};
dSection[]={};

// Define torus section
Macro turn
  P1= newp; Point(P1) = {r1,0,z0-L, h};
  P2= newp; Point(P2) = {r2,0,z0-L, h};
  P3= newp; Point(P3) = {r2,0,z0+L, h};
  P4= newp; Point(P4) = {r1,0,z0+L, h};

  L12=newl; Line(L12) = {P1, P2};
  L23=newl; Line(L23) = {P2, P3};
  L34=newl; Line(L34) = {P3, P4};
  L41=newl; Line(L41) = {P4, P1};

  L0=newl; Line Loop(L0) = {L12, L23, L34, L41};
  S=newreg; Plane Surface(S) = {L0};
  V=Extrude { {0,0,1} , {0,0,0} , Pi/2. } {Surface{S};};
  
  dSection[t]=L0;
  Section[t]=S;
  coil[t]=V[1];
  V1[t]=V[0];

  Ri[t]=V[5];
  Re[t]=V[3];
  HP[t]=V[4];
  BP[t]=V[2];
  
Return

Nturn=2;
eps=1*Unit;
z0=-(3*(2*L)+2*eps)/2.;
dz=2*z0;
For t In {0:Nturn-1}
  Call turn;
  z0+=dz;
EndFor

// Define ext section
P5=newp; Point(P5) = {0,0, -r_ext, h_ext};
P6=newp; Point(P6) = {r_ext, 0, 0, h_ext};
P7=newp; Point(P7) = {0, 0, r_ext, h_ext};

L05=newl; Line(L05) = {P0, P5};
C56=newl; Circle(C56) = {P5, P0, P6};
C67=newl; Circle(C67) = {P6, P0, P7};
L70=newl; Line(L70) = {P7, P0};

L_ext=newl; Line Loop(L_ext) = {L05, C56, C67, L70};
S_ext=newreg; Plane Surface(S_ext) = {L_ext, -dSection[]};

// Define inf section
P8=newp; Point(P8) = {0,     0, -r_inf, h_inf};
P9=newp; Point(P9) = {r_inf, 0, 0, h_inf};
P10=newp; Point(P10) = {0,   0, r_inf, h_inf};

L58=newl; Line(L58) = {P5, P8};
C89=newl; Circle(C89) = {P8, P0, P9};
C910=newl; Circle(C910) = {P9, P0, P10};
L107=newl; Line(L107) = {P10, P7};

L_inf=newl; Line Loop(L_inf) = {L58, C89, C910, L107, -C67, -C56};
S_inf=newreg; Plane Surface(S_inf) = {L_inf};

// Build 3D geom

out0[] = Extrude { {0,0,1} , {0,0,0} , Pi/2. } { 
  Surface{S_ext}; //Layers{N_Layers}; //Recombine; 
};
out1[] = Extrude { {0,0,1} , {0,0,0} , Pi/2. } { 
  Surface{S_inf}; //Layers{N_Layers}; //Recombine; 
};


//  Define Physical

Printf("coil: %g", #coil[]-1);
Printf("out0: %g", #out0[]-1);
Printf("out1: %g", #out1[]-1);

For t In {0:Nturn-1}
  Physical Volume(Sprintf("coil%g", t)) = {coil[t]}; // Tore
  Physical Surface(Sprintf("V0_%g", t)) = {Section[t]};
  Physical Surface(Sprintf("V1_%g", t)) = {V1[t]};
  Physical Surface(Sprintf("Rint_%g", t)) = {Ri[t]};
  Physical Surface(Sprintf("Rext_%g", t)) = {Re[t]};
  Physical Surface(Sprintf("HP_%g", t)) = {HP[t]};
  Physical Surface(Sprintf("BP_%g", t)) = {BP[t]};  
EndFor
Physical Volume("air") = {out0[1], out1[1]}; // Infini

Physical Surface("OXOZ") = {S_ext, S_inf};
Physical Surface("OYOZ") = {out0[0], out1[0]};

n=0;
BorderIds[]={};
// For id In {2 : out1[]-1}
For id In {2 : 3} 
    BorderIds[n]=out1[id];
    Printf("Border[%g]=%g", n, out1[id]);
    n++;
EndFor

Physical Surface("Border") = {BorderIds[]};

//  Physical Surface("Border") = {85, 88, 163, 166, 222, 225, 241, 244, 300, 303 };  // Inf

