#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/checker.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>



int main(int argc, char**argv )
{
  using namespace Feel;
  po::options_description options( "MQS options" );
  options.add_options()
    ( "sigma", po::value<std::string>()->default_value( "1" ), "electrical conductivity" )
    ( "mu_mag", po::value<std::string>()->default_value( "1" ), "relative magnetic permeability" )
    ( "A0", po::value<std::string>()->default_value( "{0,0,0}" ), "initial A" )
    ( "V0", po::value<std::string>()->default_value( "0" ), "initial V" )
    ( "Aexact", po::value<std::string>()->default_value( "{0,0,0}" ), "exact A" )
    ( "Vexact", po::value<std::string>()->default_value( "0" ), "exact V" )
    ( "Ad", po::value<std::string>()->default_value( "{0,0,0}" ), "Ad" )
    ( "v0", po::value<std::string>()->default_value( "0" ), "v0" )
    ( "v1", po::value<std::string>()->default_value( "0" ), "v1" );

  Environment env( _argc=argc, _argv=argv,_desc=options,
		   _about=about(_name="mqs",
				_author="Feel++ Consortium",
				_email="feelpp-devel@feelpp.org"));

  // Dirichlet for Magnetic potential
  auto Ad = expr<3,1>(soption(_name="Ad"));
  Feel::cout << "Ad=" << Ad << std::endl;

  // Dirichlet for electric potential
  auto v1 = expr(soption(_name="v1"));
  Feel::cout << "v1=" << v1 << std::endl;

  auto v0 = expr(soption(_name="v0"));
  Feel::cout << "vO=" << v0 << std::endl;

  //Recuperer time frame

  double dt = doption(_name = "ts.time-step");
  std::cout << "time-step=" << dt << std::endl;

  double tmax = doption(_name = "ts.time-final");
  std::cout << "time-final=" << tmax << std::endl;

  // Init solution for Magnetic Potential
  auto A0 = expr<3,1>(soption(_name="A0"));
  Feel::cout << "A0=" << A0 << std::endl;

  // Init solution for Potential
  auto V0 = expr(soption(_name="V0"));
  Feel::cout << "V0=" << V0 << std::endl;
 
  // Define sigma and mu
  auto sigma = expr(soption(_name = "sigma"));
  Feel::cout << "sigma=" << sigma << std::endl;

  auto mur = expr(soption(_name = "mu_mag"));
  Feel::cout << "mur=" << mur << std::endl;

  // Enforce a solution
  auto Aexact_g = expr<3, 1>(soption(_name = "Aexact"));
  Feel::cout << "Aexact=" << Aexact_g << std::endl;

  auto Vexact_g = expr(soption(_name = "Vexact"));
  Feel::cout << "Vexact=" << Vexact_g << std::endl;

  // Load Mesh and define Product space
  
  auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);
  auto cond_mesh = createSubmesh(mesh,markedelements(mesh,"Omega_C"));

  auto Ah = Pchv<1>( mesh );
  auto Vh = Pch<1>( cond_mesh );

  auto A = Ah->element(A0); // how to init A to A0?;
  auto V = Vh->element(V0);

  auto Aexact = Ah->element();
  auto Vexact = Vh->element();

  auto phi = Ah->element();
  auto psi = Vh->element();
  
  auto Zh = product(Ah,Vh);
  auto U = Zh.element();
#if 0
  auto cAh = Zh[0_c];
  auto cVh = Zh[1_c];
#endif  
  

  auto rhs = blockform1( Zh );
  auto lhs = blockform2( Zh );

  double t = 0;

  auto e = exporter( _mesh=mesh );

  Aexact_g.setParameterValues({{"t", t}});
  Aexact = project(_space = Ah, _expr = Aexact_g);
  
  Vexact_g.setParameterValues({{"t", t}});
  Vexact = project(_space = Vh, _expr = Vexact_g);
  
  e->step(t)->add("A", A0);
  e->step(t)->add("Aexact", Aexact);
  e->step(t)->add("V", V0);
  e->step(t)->add("Vexact", Vexact);
  e->save();


  double H1Aerror = 0;
  double L2Aerror = 0;
  double H1Verror = 0;
  double L2Verror = 0;
  
  auto mu0 = 4.e-7 * M_PI ; // SI Unit : H/m = m.kg/s2/A2

  for (t = dt; t <= tmax; t += dt){
    Aexact_g.setParameterValues({{"t", t}});
    Aexact = project(_space = Ah, _expr = Aexact_g);
    Vexact_g.setParameterValues({{"t", t}});
    Vexact = project(_space = Vh, _expr = Vexact_g);
    v0.setParameterValues({{"t", t}});
    v1.setParameterValues({{"t", t}});
    Ad.setParameterValues({{"t", t}});
    
    lhs.zero();
    rhs.zero();
    tic();
    // Ampere law: sigma dA/dt + rot(1/(mu-r*mu_0) rotA) + sigma grad(V) = Js
    lhs(0_c, 0_c) = integrate( _range=elements(mesh),
		     _expr = dt * inner(curl(phi) , curlt(A)) );
    lhs(0_c, 0_c) += integrate( _range=elements(cond_mesh),
		     _expr = mur * mu0 * sigma * inner(id(phi) , idt(A) ));

    lhs(0_c, 1_c) += integrate(_range=elements(cond_mesh),
         _expr = dt * mu0 * mur * sigma*inner(id(phi),trans(gradt(V))) );

    rhs(0_c) += integrate(_range=elements(cond_mesh),
                        _expr = mu0 * mur * sigma * inner(id(phi) , idv(A)));

    // Current conservation
    lhs(1_c, 0_c) += integrate( _range=elements(cond_mesh),
			       _expr = sigma * inner(idt(A), trans(grad(psi))) );
      
    lhs(1_c, 1_c) += integrate( _range=elements(cond_mesh),
			       _expr = sigma * dt * inner(gradt(V), grad(psi)) );

    rhs(1_c) += integrate(_range=elements(cond_mesh),
                        _expr = sigma * inner(idv(A), trans(grad(psi))) );

    /* Add Boundary conditions */
#if 0
    lhs(0_c, 0_c) += on(_range=markedfaces(mesh,"Infty"), _rhs=rhs(0_c), _element=phi, _expr= Ad);
#endif

#if 1
    lhs(0_c, 0_c) += integrate( _range=boundaryfaces(mesh),
                       		     _expr = dt * inner(cross(id(phi),N()) , curlt(A)) );
    lhs(0_c, 0_c) += integrate( _range=boundaryfaces(mesh),
                       		     _expr = dt * inner(cross(idt(phi),N()) , curl(A)) );
    lhs(0_c, 0_c) += integrate( _range=boundaryfaces(mesh),
                       		     _expr = inner(cross(id(phi),N()) , cross(idt(A),N()))/hFace() );  
    rhs(0_c) += integrate(_range=boundaryfaces(mesh),
                        _expr = dt * inner(cross(idt(phi),N()) , curl(A)) );                                    
    rhs(0_c) += integrate(_range=boundaryfaces(mesh),
                        _expr = inner(cross(id(phi),N()) , cross(Ad,N()))/hFace());                                                                                              
#else
#if 1
    /* 1/4th of a torus + Air */
    lhs(0_c, 0_c) += on(_range=markedfaces(mesh,"V0"), _rhs=rhs(0_c), _element=phi, _expr= Ad);
    lhs(0_c, 0_c) += on(_range=markedfaces(mesh,"V1"), _rhs=rhs(0_c), _element=phi, _expr= Ad);
    lhs(0_c, 0_c) += on(_range=markedfaces(mesh,"Gamma_C"), _rhs=rhs(0_c), _element=phi, _expr= Ad);
#endif
      
    lhs(1_c, 1_c) += on(_range=markedfaces(cond_mesh,"V0"), _rhs=rhs(1_c), _element=psi, _expr= v0);
    lhs(1_c, 1_c) += on(_range=markedfaces(cond_mesh,"V1"), _rhs=rhs(1_c), _element=psi, _expr= v1);
    lhs(1_c, 1_c) += on(_range=markedfaces(cond_mesh,"Gamma_C"), _rhs=rhs(1_c), _element=psi, _expr= Vexact_g);
#endif    
    toc("assembling", true);

    /* Solve */
    tic();
    lhs.solve(_rhs=rhs,_solution=U);
    toc("solve", true);

    tic();
    e->step(t)->add( "A", U(0_c));
    e->step(t)->add( "V", U(1_c));
    e->step(t)->add( "Aexact", Aexact);
    e->step(t)->add( "Vexact", Vexact);
    e->save();
    toc("export", true);

    double L2Aexact = normL2(_range = elements(mesh), _expr = Aexact_g);
    L2Aerror = normL2(_range=elements(mesh), _expr=(idv(U(0_c)) - idv(Aexact)));
    H1Aerror = normH1(_range=elements(mesh), _expr = (idv(U(0_c)) - idv(Aexact)), _grad_expr = (gradv(U(0_c)) - gradv(Aexact)));

    Feel::cout << "A error: " << "t="<< t << " " 
               << L2Aexact << " " 
               << L2Aerror << " " << L2Aerror / L2Aexact << " " << H1Aerror << std::endl;

    double L2Vexact = normL2(_range = elements(mesh), _expr = Vexact_g);
    L2Verror = normL2(_range=elements(mesh), _expr=(idv(U(1_c)) - idv(Vexact)));
    H1Verror = normH1(_range=elements(mesh), _expr = (idv(U(1_c)) - idv(Vexact)), _grad_expr = (gradv(U(1_c)) - gradv(Vexact)));

    Feel::cout << "V error: " << "t="<< t << " " 
               << L2Vexact << " " 
               << L2Verror << " " << L2Verror / L2Vexact << " " << H1Verror << std::endl;

    #if 0
    A = U(0_c); 
    V = U(1_c);
    #endif
  }
}
