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
    ( "mu_mag", po::value<std::string>()->default_value( "1" ), "relative magnetic permeability" );
  Environment env( _argc=argc, _argv=argv,_desc=options,
		   _about=about(_name="Maxwell Quasi-Static",
				_author="Feel++ Consortium",
				_email="feelpp-devel@feelpp.org"));

  // Dirichlet for Magnetic potential
  auto Ad = expr(soption(_name="functions.Ad"));
  Feel::cout << "Ad=" << Ad << std::endl;

  // Dirichlet for electric potential
  auto v1 = expr(soption(_name="functions.v1"));
  Feel::cout << "v1=" << v1 << std::endl;

  auto v0 = expr(soption(_name="functions.vO"));
  Feel::cout << "vO=" << v0 << std::endl;

  //Recuperer time frame

  double dt = doption(_name = "ts.time-step");
  std::cout << "time-step=" << dt << std::endl;

  double tmax = doption(_name = "ts.time-final");
  std::cout << "time-final=" << tmax << std::endl;

  // Init solution for Magnetic Potential
  auto A0 = expr<3,1>(soption(_name="functions.A0"));
  Feel::cout << "A0=" << A0 << std::endl;
     
  // Define sigma and mu
  auto sigma = doption(_name = "sigma");
  auto mur = doption(_name = "mu_mag");

  // Load Mesh and define Product space
  
  auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);
  auto cond_mesh = createSubmesh(mesh,markedelements(mesh,"Omega_c"));

  auto Ah = Pchv<1>( mesh );
  auto Vh = Pch<1>( cond_mesh );

  auto A = Ah->element(A0); // how to init A to A0?;
  auto V = Vh->element();

  auto phi = Ah->element();
  auto psi = Vh->element();
  
  auto Zh = product(Ah,Vh);
  auto U = Zh->element();
#if 1
  auto cAh = Zh.functionSpace(0_c);
  auto cVh = Zh.functionSpace(1_c);
#endif  
  

  auto rhs = blockform1( Zh );
  auto lhs = blockform2( Zh );

  double t = dt;

  auto e = exporter( _mesh=mesh );

  auto mu0 = 4.e-7 * M_PI ; // SI Unit : H/m = m.kg/s2/A2

  while(t < tmax){

#if 1
    tic();
    // Ampere law
    lhs(0_c, 0_c) += integrate( _range=elements(mesh),
		   _expr = dt * inner(curl(phi) , curlt(A)) );
    lhs(0_c, 0_c) += integrate( _range=elements(cond_mesh),
		     _expr = mur * mu0 * sigma * inner(id(phi) , idt(A) ));

    lhs(0_c, 1_c) += integrate(_range=elements(cond_mesh),_expr = dt * mu0 * mur * sigma*inner(trans(grad(V)),id(phi)));

    // Current conservation
    lhs(1_c, 0_c) += integrate( _range=elements(cond_mesh),
			       _expr = sigma * inner(idt(A), grad(psi)) );
      
    lhs(1_c, 1_c) = integrate( _range=elements(cond_mesh),
			       _expr = sigma * dt * inner(idt(V), grad(psi)) );

    /* Add Boundary conditions */
    lhs(0_c, 0_c) += on(_range=markedfaces(mesh,"Infty"), _rhs=rhs(0_c), _element=phi, _expr= Ad);

#if 1
    /* 1/4th of a torus + Air */
    lhs(0_c, 0_c) += on(_range=markedfaces(mesh,"V0"), _rhs=rhs(0_c), _element=phi, _expr= Ad);
    lhs(0_c, 0_c) += on(_range=markedfaces(mesh,"V1"), _rhs=, _element=phi, _expr= Ad);
#endif
      
    lhs(1_c, 1_c) += on(_range=markedfaces(cond_mesh,"VO"), _rhs=rhs(1_c), _element=psi, _expr= v0);
    lhs(1_c, 1_c) += on(_range=markedfaces(cond_mesh,"V1"), _rhs=rhs(1_c), _element=psi, _expr= v1);
    toc("assembling", true);

    /* Solve */
    tic();
    lhs.solve(_rhs=rhs,_solution=U);
    toc("solve", true);

    tic();
    e->step(t)->add( "A", U(0_c));
    e->step(t)->add( "V", U(1_c));
    e->save();
    toc("export", true);
#endif
    t += dt;
  }
}
