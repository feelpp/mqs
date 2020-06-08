
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/checker.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>

#include <feel/feelmodels/modelproperties.hpp>

int main(int argc, char**argv )
{
  using namespace Feel;
  po::options_description options( "MQS options" );
  options.add_options()
    ("model-file", Feel::po::value<std::string>()->default_value( "" ), "file describing model properties")
    ( "weakdir", po::value<bool>()->default_value( "false" ), "use Dirichlet weak formulation" )
    ( "penalty-coeff", po::value<double>()->default_value( 1.e+3 ), "penalty coefficient for weak Dirichlet" )
    ( "A0", po::value<std::string>()->default_value( "{0,0,0}" ), "initial A" )
    ( "V0", po::value<std::string>()->default_value( "0" ), "initial V" )
    ( "Aexact", po::value<std::string>()->default_value( "" ), "exact A" )
    ( "Vexact", po::value<std::string>()->default_value( "" ), "exact V" );

  Environment env( _argc=argc, _argv=argv,_desc=options.add(Feel::backend_options("mqs")),
		   _about=about(_name="mqs",
				_author="Feel++ Consortium",
				_email="feelpp-devel@feelpp.org"));

  //Recuperer time frame
  double dt = doption(_name = "ts.time-step");
  std::cout << "time-step=" << dt << std::endl;

  double tmax = doption(_name = "ts.time-final");
  std::cout << "time-final=" << tmax << std::endl;

  // Eventually get a solution
  bool Uexact = false;
  auto Aexact_g = expr<3, 1>("{0,0,0}");
  auto Vexact_g = expr("0");

  std::string Aexact_s = soption(_name = "Aexact");
  Feel::cout << "Aexact=" << Aexact_s << std::endl;
  std::string Vexact_s = soption(_name = "Vexact");
  Feel::cout << "Vexact=" << Vexact_s << std::endl;

  if ( !Aexact_s.empty() )
    {
      Uexact = true;
      Aexact_g = expr<3, 1>(Aexact_s);
      Vexact_g = expr(Vexact_s);
      Feel::cout << "* Aexact=" << Aexact_s << std::endl;
      Feel::cout << "* Vexact=" << Vexact_s << std::endl;
    }

  // Load Mesh
  auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);

  // Load json model file
  std::shared_ptr<ModelProperties> M_modelProps;

  std::string modelPropFilename = Environment::expand( soption( _name="model-file") );
  if ( !modelPropFilename.empty() )
    M_modelProps = std::make_shared<ModelProperties>( modelPropFilename );
  else
    throw std::logic_error( "model-file: " + soption(_name="model-file") + " no such file" );

  auto M_materials = M_modelProps->materials().materialWithPhysic(std::vector<std::string>({"electric","thermo-electric"}));
  std::vector<std::string> range;
  for( auto const& mp : M_materials )
    for (auto const& marker : mp.second.meshMarkers() )
      range.push_back(marker);
  Feel::cout << "Electric Materials markers: " << range << std::endl;
  
  // Define SpaceFunctions
  auto Ah = Pchv<1>( mesh );
  auto Vh = Pch<1>( mesh, markedelements(mesh, range) );
  auto cond_mesh = Vh->mesh();
  
  // Space Product
  auto Zh = product(Ah,Vh);
  auto U = Zh.element();

  tic();
  solve::strategy strategy = solve::strategy::monolithic; // if it enough if I want use fieldsplit?
  auto rhs = blockform1( Zh, strategy ,backend(_name="mqs") );
  auto lhs = blockform2( Zh, strategy ,backend(_name="mqs") );
  toc("create blockforms", true);

  double t = 0;

  auto e = exporter( _mesh=mesh );

  tic();
  auto Aexact = Ah->element();
  auto Vexact = Vh->element();

  auto A = Ah->element(); //Ah->element(A0); // how to init A to A0?;
  auto V = Vh->element(); //Vh->element(V0);
  if ( Uexact )
    {
      Aexact_g.setParameterValues({{"t", t}});
      Aexact = project(_space = Ah, _expr = Aexact_g);
      Feel::cout << "Define Aexact" << std::endl;

      Vexact_g.setParameterValues({{"t", t}});
      Vexact = project(_space = Vh, _expr = Vexact_g);
      Feel::cout << "Define Vexact" << std::endl;
  
      A = Aexact;
      V = Vexact;
    }
  toc("init exact solution", true);

  auto phi = Ah->element();
  auto psi = Vh->element();
  
  
#if 1
  node_type pt(3);
  pt[0] = 0.5; pt[1] = 0.5; pt[2] = 2.5;
  auto val = A(pt);
  auto Ax = val(0,0,0); // evaluation de Ay
  auto Ay = val(1,0,0); // evaluation de Az
  auto Az = val(2,0,0); // evaluation de Az
  auto Vval = V(pt);
  Feel::cout << "t=" << t << ", ";
  Feel::cout << "A(" << pt[0] << "," << pt[1] << "," << pt[2] << ") = {" << Ax << "," << Ay << "," << Az << "}, ";
  Feel::cout << "V(" << pt[0] << "," << pt[1] << "," << pt[2] << ")=" << Vval(0,0,0) << std::endl;
#endif

  tic();
  e->step(t)->add("A", A);
  e->step(t)->add("V", V);
  if ( Uexact )
    {
      e->step(t)->add("Aexact", Aexact);
      e->step(t)->add("Vexact", Vexact);
    }
  e->save();
  toc("export init solution", true);
  
  double L2Aexact, H1Aerror, L2Aerror;
  double L2Vexact, H1Verror, L2Verror;
  
  if ( Uexact )
    {
      L2Aexact = normL2(_range = elements(mesh), _expr = Aexact_g);
      H1Aerror = 0;
      L2Aerror = 0;
      L2Vexact = normL2(_range = elements(cond_mesh), _expr = Vexact_g);
      H1Verror = 0;
      L2Verror = 0;
    }
  
  auto mu0 = 4.e-7 * M_PI ; // SI Unit : H/m = m.kg/s2/A2

  
  for (t = dt; t < tmax; t += dt)
    {

      tic();
      for( auto const& pairMat : M_modelProps->materials() )
	{
	  auto name = pairMat.first;
	  auto material = pairMat.second;

	  auto mur = material.getScalar("mu_mag");

	  // Ampere law: sigma dA/dt + rot(1/(mu_r*mu_0) rotA) + sigma grad(V) = Js
	  lhs(0_c, 0_c) += integrate( _range=markedelements(mesh, material.meshMarkers()),
				      _expr = dt * 1/mur * inner(curl(phi) , curlt(A)) );
	  Feel::cout << "create lhs(0,0)" << std::endl;

	}

      for( auto const& pairMat : M_materials )
	{
	  auto name = pairMat.first;
	  auto material = pairMat.second;

	  auto sigma = material.getScalar("sigma");

	  // Ampere law: sigma dA/dt + rot(1/(mu_r*mu_0) rotA) + sigma grad(V) = Js
	  lhs(0_c, 0_c) += integrate( _range=markedelements(cond_mesh, material.meshMarkers()),
				      _expr = mu0 * sigma * inner(id(phi) , idt(A) ));
	  Feel::cout << "create lhs(0,0)" << std::endl;

	  lhs(0_c, 1_c) += integrate(_range=markedelements(cond_mesh, material.meshMarkers()),
				     _expr = dt * mu0 * sigma * inner(id(phi),trans(gradt(V))) );
	  Feel::cout << "create lhs(0,1)" << std::endl;

	  rhs(0_c) += integrate(_range=markedelements(cond_mesh, material.meshMarkers()),
				_expr = mu0 * sigma * inner(id(phi) , idv(A)));
	  Feel::cout << "create row(0)" << std::endl;

	  // Current conservation: div( -sigma grad(V) -sigma*dA/dt) = Qs
	  lhs(1_c, 1_c) += integrate( _range=markedelements(cond_mesh, material.meshMarkers()),
				      _expr = sigma * dt * inner(gradt(V), grad(psi)) );
	  Feel::cout << "create lhs(1,1)" << std::endl;
	  lhs(1_c, 0_c) += integrate( _range=markedelements(cond_mesh, material.meshMarkers()),
				      _expr = sigma * inner(idt(A), trans(grad(psi))) );
	  Feel::cout << "create lhs(1,0)" << std::endl;

	  rhs(1_c) += integrate( _range=markedelements(cond_mesh, material.meshMarkers()),
				 _expr = sigma * inner(idv(A), trans(grad(psi))) );
	  Feel::cout << "create rhs(1)" << std::endl;

	}
      toc("assembling", true);
     
      tic();
      // define M_weakdir
      auto weakdir = boption("weakdir");
      auto itField = M_modelProps->boundaryConditions().find( "magnetic-potential");
      if ( itField != M_modelProps->boundaryConditions().end() )
	{
	  auto mapField = (*itField).second;
	  auto itType = mapField.find( "Dirichlet" );
	  if ( itType != mapField.end() )
	    {
	      for ( auto const& exAtMarker : (*itType).second )
		{
		  std::string marker = exAtMarker.marker();
		  auto g = expr<3,1>(exAtMarker.expression());
		  g.setParameterValues({{"t", t}});
		  Feel::cout << "A Dirichlet[" << marker << "] : " << g << " (weak=" << weakdir << ")" << std::endl;

		  if (! weakdir )
		    {
		      lhs(0_c, 0_c) += on(_range=markedfaces(mesh,marker), _rhs=rhs(0_c), _element=phi, _expr= g);
		    }
		  else
		    {
		      auto gamma = doption("penalty-coeff");
		      lhs(0_c, 0_c) += integrate( _range=markedfaces(mesh,marker),
						  _expr = mu0 * dt * inner(cross(id(phi),N()) , curlt(A)) );
		      lhs(0_c, 0_c) += integrate( _range=markedfaces(mesh,marker),
						  _expr = mu0 * dt * inner(cross(idt(A),N()) , curl(phi)) );
		      lhs(0_c, 0_c) += integrate( _range=markedfaces(mesh,marker),
						  _expr = mu0 * gamma * inner(cross(id(phi),N()) , cross(idt(A),N()))/hFace() );  
		      rhs(0_c) += integrate(_range=markedfaces(mesh,marker),
					    _expr = mu0 * dt * inner(g , curl(phi)));
		      rhs(0_c) += integrate(_range=markedfaces(mesh,marker),
					    _expr = mu0 * gamma * inner(cross(id(phi),N()) , g)/hFace());                                                                      }
		  Feel::cout << "block(0,0) on " << marker << std::endl;
		}
	    }
	  // 	ItType = mapField.find( "Neumann" );
	  // 	if ( itType != mapField.end() )
	  // 	  {
	  // 	    for ( auto const& exAtMarker : (*itType).second )
	  // 	      {
	  // 		std::string marker = exAtMarker.marker();
	  // 		auto g = expr<3,1>(exAtMarker.expression());
	  //            g.setParameterValues({{"t", t}});
	  // 		Feel::cout << "Neuman[" << marker << "] : " << exAtMarker.expression() << std::endl;
	  //            lhs(0_c, 0_c) += integrate(_range=markedfaces(mesh,marker), ....);
	  //            Feel::cout << "block(0,0) on " << marker << std::endl;
	  // 	      }
	  // 	  }
	}   
      itField = M_modelProps->boundaryConditions().find( "electric-potential");
      if ( itField != M_modelProps->boundaryConditions().end() )
	{
	  auto mapField = (*itField).second;
	  auto itType = mapField.find( "Dirichlet" );
	  if ( itType != mapField.end() )
	    {
	      for ( auto const& exAtMarker : (*itType).second )
		{
		  std::string marker = exAtMarker.marker();
		  auto g = expr(exAtMarker.expression());
		  g.setParameterValues({{"t", t}});
		  Feel::cout << "V Dirichlet[" << marker << "] : " << g << " (weak=" << weakdir << ")" << std::endl;
		  if (! weakdir )
		    {
		      lhs(1_c, 1_c) += on(_range=markedfaces(cond_mesh,marker), _rhs=rhs(1_c), _element=psi, _expr= g);
		    }
		  else
		    {
		      // see: http://docs.feelpp.org/math/fem/laplacian/nitsche/#_weak_treatment_of_dirichlet_boundary_conditions
		      // shall be multiply by sigma but howto get sigma
		      
		      auto gamma = doption("penalty-coeff");
		      lhs(1_c, 1_c) += integrate( _range=markedfaces(cond_mesh, marker),
						  _expr = dt *  -(gradt(V)*N())*id(psi) );
		      lhs(1_c, 1_c) += integrate( _range=markedfaces(cond_mesh, marker),
						  _expr = dt * -(grad(psi)*N())*idt(V) );
		      lhs(1_c, 1_c) += integrate( _range=markedfaces(cond_mesh, marker),
						  _expr = dt * gamma*id(psi)*idt(V)/hFace() );  
		      rhs(1_c) += integrate(_range=markedfaces(cond_mesh, marker),
					    _expr = dt * -(grad(psi)*N())*g );
		      rhs(1_c) += integrate(_range=markedfaces(cond_mesh, marker),
					    _expr = dt * gamma*id(psi)*g/hFace());                                                                                      }
		  Feel::cout << "block(1,1) on " << marker << std::endl;
		}
	    }
	}       
      toc("boundary conditions", true);
    
      /* Solve */
      tic();
      auto result = lhs.solve(_rhs=rhs,_solution=U);
      toc("solve", true);

#if 1
      A = U(0_c); 
      V = U(1_c);

      val = A(pt);
      Ax = val(0,0,0); // evaluation de Ay
      Ay = val(1,0,0); // evaluation de Az
      Az = val(2,0,0); // evaluation de Az
      Vval = V(pt);
      Feel::cout << "t=" << t << ", ";
      Feel::cout << "A(" << pt[0] << "," << pt[1] << "," << pt[2] << ") = {" << Ax << "," << Ay << "," << Az << "}, ";
      Feel::cout << "V(" << pt[0] << "," << pt[1] << "," << pt[2] << ")=" << Vval(0,0,0) << std::endl;
#endif

      tic();
      e->step(t)->add( "A", U(0_c));
      e->step(t)->add( "V", U(1_c));
    
      if ( Uexact )
	{
	  Aexact_g.setParameterValues({{"t", t}});
	  Aexact = project(_space = Ah, _expr = Aexact_g);
	  Vexact_g.setParameterValues({{"t", t}});
	  Vexact = project(_space = Vh, _expr = Vexact_g);

	  e->step(t)->add( "Aexact", Aexact);
	  e->step(t)->add( "Vexact", Vexact);

	  L2Aexact = normL2(_range = elements(mesh), _expr = Aexact_g);
	  L2Aerror = normL2(elements(mesh), (idv(U(0_c)) - idv(Aexact)));
	  H1Aerror = normH1(elements(mesh), _expr = (idv(U(0_c)) - idv(Aexact)), _grad_expr = (gradv(U(0_c)) - gradv(Aexact)));

	  Feel::cout << "error: " << "t="<< t;
	  Feel::cout << " A:" << L2Aerror << " " << L2Aerror / L2Aexact << " " << H1Aerror;

	  L2Vexact = normL2(_range = elements(cond_mesh), _expr = Vexact_g);
	  L2Verror = normL2(elements(cond_mesh), (idv(U(1_c)) - idv(Vexact)));
	  H1Verror = normH1(elements(cond_mesh), _expr = (idv(U(1_c)) - idv(Vexact)), _grad_expr = (gradv(U(1_c)) - gradv(Vexact)));

	  Feel::cout << " V: " << L2Verror << " " << L2Verror / L2Vexact << " " << H1Verror << std::endl;
	}
      e->save();
      toc("export", true);

#if 1
      A = U(0_c); 
      V = U(1_c);
#endif
      
      /* reinit  */
      lhs.zero();
      rhs.zero();
    }
}
