
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
  if (Environment::worldComm().isMasterRank())
    {
      std::cout << "mesh->numGlobalElements() "<< mesh->numGlobalElements() << std::endl;
      std::cout << "cond_mesh->numGlobalElements() "<< cond_mesh->numGlobalElements() << std::endl;
      std::cout << "Ah->nDof() "<<Ah->nDof() << std::endl;
      std::cout << "Vh->nDof() "<<Vh->nDof() << std::endl;
    }

  auto A = Ah->elementPtr(); //Ah->element(A0); // how to init A to A0?;
  auto V = Vh->elementPtr(); //Vh->element(V0);

  // Vincent way
  tic();
  BlocksBaseGraphCSR myblockGraph(2,2);
  myblockGraph(0,0) = stencil(_test=Ah,_trial=Ah, _diag_is_nonzero=false, _close=false)->graph();
  myblockGraph(0,1) = stencil(_test=Ah,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
  myblockGraph(1,0) = stencil(_test=Vh,_trial=Ah, _diag_is_nonzero=false, _close=false)->graph();
  myblockGraph(1,1) = stencil(_test=Vh,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
  auto M = backend()->newBlockMatrix(_block=myblockGraph);

  BlocksBaseVector<double> myblockVec(2);
  myblockVec(0,0) = backend()->newVector( Ah );
  myblockVec(1,0) = backend()->newVector( Vh );
  auto F = backend()->newBlockVector(_block=myblockVec, _copy_values=false);

  BlocksBaseVector<double> myblockVecSol(2);
  myblockVecSol(0,0) = A;
  myblockVecSol(1,0) = V;
  auto U = backend()->newBlockVector(_block=myblockVecSol, _copy_values=false);
  toc("create blockforms", true);

  auto mybackend = backend(_name="mqs");

  double t = 0;

  auto e = exporter( _mesh=mesh );

  tic();
  auto Aexact = Ah->element();
  auto Vexact = Vh->element();

  if ( Uexact )
    {
      Aexact_g.setParameterValues({{"t", t}});
      Aexact = project(_space = Ah, _expr = Aexact_g);
      Feel::cout << "Define Aexact" << std::endl;

      Vexact_g.setParameterValues({{"t", t}});
      Vexact = project(_space = Vh, _expr = Vexact_g);
      Feel::cout << "Define Vexact" << std::endl;
  
      (*A) = Aexact;
      (*V) = Vexact;
    }
  toc("init exact solution", true);

  
#if 1
  node_type pt(3);
  pt[0] = 0.5; pt[1] = 0.5; pt[2] = 2.5;
  auto val = (*A)(pt);
  auto Ax = val(0,0,0); // evaluation de Ay
  auto Ay = val(1,0,0); // evaluation de Az
  auto Az = val(2,0,0); // evaluation de Az
  auto Vval = (*V)(pt);
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
  
  auto mu0 = 1; // 4.e-7 * M_PI ; // SI Unit : H/m = m.kg/s2/A2

  
  for (t = dt; t < tmax; t += dt)
    {

      tic();
      for( auto const& pairMat : M_modelProps->materials() )
	{
	  auto name = pairMat.first;
	  auto material = pairMat.second;

	  auto mur = material.getScalar("mu_mag");

	  // Ampere law: sigma dA/dt + rot(1/(mu_r*mu_0) rotA) + sigma grad(V) = Js
	  form2( _trial=Ah, _test=Ah ,_matrix=M )
	    += integrate( _range=markedelements(mesh, material.meshMarkers()),
			  _expr = dt * 1/mur * inner(curl(A) , curlt(A)) );
	  Feel::cout << "create lhs(0,0)" << std::endl;

	}

      for( auto const& pairMat : M_materials )
	{
	  auto name = pairMat.first;
	  auto material = pairMat.second;

	  auto sigma = material.getScalar("sigma");

	  // Ampere law: sigma dA/dt + rot(1/(mu_r*mu_0) rotA) + sigma grad(V) = Js
	  form2( _trial=Ah, _test=Ah ,_matrix=M )
	    += integrate( _range=markedelements(cond_mesh, material.meshMarkers()),
			  _expr = mu0 * sigma * inner(id(A) , idt(A) ));
	  Feel::cout << "create lhs(0,0)" << std::endl;

	  form2( _trial=Vh, _test=Ah ,_matrix=M, _rowstart=0, _colstart=1 )
	    += integrate(_range=markedelements(cond_mesh, material.meshMarkers()),
			 _expr = dt * mu0 * sigma * inner(id(A),trans(gradt(V))) );
	  Feel::cout << "create lhs(0,1)" << std::endl;

	  form1( _test=Ah, _vector=F )
	    += integrate(_range=markedelements(cond_mesh, material.meshMarkers()),
			 _expr = mu0 * sigma * inner(id(A) , idv(A)));
	  Feel::cout << "create rhs(0)" << std::endl;

	  // auto Js = ;
	  // form1( _test=Ah, _vector=F )
	  //   += integrate(_range=markedelements(cond_mesh, material.meshMarkers()),
	  // 		 _expr = dt * mu0 * inner(id(A) , Js));
	  // Feel::cout << "create rhs(0)" << std::endl;

	  // Current conservation: div( -sigma grad(V) -sigma*dA/dt) = Qs
	  form2( _trial=Vh, _test=Vh ,_matrix=M, _rowstart=1, _colstart=1 )
	    += integrate( _range=markedelements(cond_mesh, material.meshMarkers()),
			  _expr = sigma * dt * inner(gradt(V), grad(V)) );
	  Feel::cout << "create lhs(1,1)" << std::endl;

	  form2( _trial=Ah, _test=Vh ,_matrix=M, _rowstart=1, _colstart=0 )
	    += integrate( _range=markedelements(cond_mesh, material.meshMarkers()),
			  _expr = sigma * inner(idt(A), trans(grad(V))) );
	  Feel::cout << "create lhs(1,0)" << std::endl;

	  form1( _test=Vh ,_vector=F, _rowstart=1 )
	    += integrate( _range=markedelements(cond_mesh, material.meshMarkers()),
			  _expr = sigma * inner(idt(A), trans(grad(V))) );
	  Feel::cout << "create rhs(1)" << std::endl;

	  // auto Qs = ...;
	  // form1( _test=Vh, _vector=F, _rowstart=1 )
	  //   += integrate(_range=markedelements(cond_mesh, material.meshMarkers()),
	  // 		 _expr = dt * Qs * id(V);
	  // Feel::cout << "create row(1)" << std::endl;
	}
      toc("assembling", true);
     
      tic();
      // Implement Dirichlet fort
      // TODO: define M_weakdir
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
		  Feel::cout << "A Dirichlet[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
		  form2( _trial=Ah, _test=Ah ,_matrix=M, _rowstart=0, _colstart=0 )
		    += on(_range=markedfaces(mesh,marker), _rhs=F, _element=*A, _expr= g);
		}
	    }
	  // 	ItType = mapField.find( "Neumann" );
	  // 	if ( itType != mapField.end() )
	  // 	  {
	  // 	    for ( auto const& exAtMarker : (*itType).second )
	  // 	      {
	  // 		std::string marker = exAtMarker.marker();
	  // 		auto g = expr<3,1>(exAtMarker.expression());
	  //          g.setParameterValues({{"t", t}});
	  // 		Feel::cout << "Neuman[" << marker << "] : " << exAtMarker.expression() << std::endl;
	  //          lhs(0_c, 0_c) += integrate(_range=markedfaces(mesh,marker), ....);
	  //          Feel::cout << "block(0,0) on " << marker << std::endl;
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
		  Feel::cout << "V Dirichlet[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
		  form2( _trial=Vh, _test=Vh ,_matrix=M, _rowstart=1, _colstart=1 )
		    += on(_range=markedfaces(cond_mesh,marker), _rhs=F, _element=*V, _expr= g);
		}
	    }
	}       
      toc("boundary conditions", true);
    
      /* Solve */
      tic();
      auto result = mybackend->solve( _matrix=M, _rhs=F, _solution=U);
      toc("solve", true);

      myblockVecSol.localize(U);
#if 1
      val = (*A)(pt);
      Ax = val(0,0,0); // evaluation de Ay
      Ay = val(1,0,0); // evaluation de Az
      Az = val(2,0,0); // evaluation de Az
      Vval = (*V)(pt);
      Feel::cout << "t=" << t << ", ";
      Feel::cout << "A(" << pt[0] << "," << pt[1] << "," << pt[2] << ") = {" << Ax << "," << Ay << "," << Az << "}, ";
      Feel::cout << "V(" << pt[0] << "," << pt[1] << "," << pt[2] << ")=" << Vval(0,0,0) << std::endl;
#endif

      tic();
      e->step(t)->add( "A", A);
      e->step(t)->add( "V", V);
    
      if ( Uexact )
	{
	  Aexact_g.setParameterValues({{"t", t}});
	  Aexact = project(_space = Ah, _expr = Aexact_g);
	  Vexact_g.setParameterValues({{"t", t}});
	  Vexact = project(_space = Vh, _expr = Vexact_g);

	  e->step(t)->add( "Aexact", Aexact);
	  e->step(t)->add( "Vexact", Vexact);
	}
      e->save();
      toc("export", true);

      // Compute error
      if ( Uexact )
	{
	  L2Aexact = normL2(_range = elements(mesh), _expr = Aexact_g);
	  L2Aerror = normL2(elements(mesh), (idv(A) - idv(Aexact)));
	  H1Aerror = normH1(elements(mesh), _expr = (idv(A) - idv(Aexact)), _grad_expr = (gradv(A) - gradv(Aexact)));
#if 0
	  L2Aerror = normL2(elements(mesh), (idv(U(0_c)) - idv(Aexact)));
	  H1Aerror = normH1(elements(mesh), _expr = (idv(U(0_c)) - idv(Aexact)), _grad_expr = (gradv(U(0_c)) - gradv(Aexact)));
#endif
	  Feel::cout << "error: " << "t="<< t;
	  Feel::cout << " A: " << L2Aerror << " " << L2Aerror / L2Aexact << " " << H1Aerror << " ";

	  L2Vexact = normL2(_range = elements(cond_mesh), _expr = Vexact_g);
	  L2Verror = normL2(elements(cond_mesh), (idv(V) - idv(Vexact)));
	  H1Verror = normH1(elements(cond_mesh), _expr = (idv(V) - idv(Vexact)), _grad_expr = (gradv(V) - gradv(Vexact)));
#if 0
	  L2Verror = normL2(elements(cond_mesh), (idv(U(1_c)) - idv(Vexact)));
	  H1Verror = normH1(elements(cond_mesh), _expr = (idv(U(1_c)) - idv(Vexact)), _grad_expr = (gradv(U(1_c)) - gradv(Vexact)));
#endif
	  Feel::cout << " V: " << L2Verror << " " << L2Verror / L2Vexact << " " << H1Verror << std::endl;
	}
    
      /* reinit  */
      M->zero();
      F->zero();
#if 0
      lhs.zero();
      rhs.zero();
#endif
    }
}
