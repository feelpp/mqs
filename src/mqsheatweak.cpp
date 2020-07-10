#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/checker.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelalg/vectorblock.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>

#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelts/bdf.hpp>

int main(int argc, char**argv )
{

  using namespace Feel;
  po::options_description options( "MQS options" );
  options.add_options()
    ( "model-file", Feel::po::value<std::string>()->default_value( "" ), "file describing model properties")
    ( "verbosity", po::value<int>()->default_value( 0 ), "set verbosisity level" )
    ( "weakdir", po::value<bool>()->default_value( "false" ), "use Dirichlet weak formulation" )
    ( "penalty-coeff", po::value<double>()->default_value( 1.e+3 ), "penalty coefficient for weak Dirichlet" )
    ( "A0", po::value<std::string>()->default_value( "{0,0,0}" ), "initial A" )
    ( "V0", po::value<std::string>()->default_value( "0" ), "initial V" )
    ( "T0", po::value<std::string>()->default_value( "0" ), "initial T" )
    ( "Aexact", po::value<std::string>()->default_value( "" ), "exact A" )
    ( "Vexact", po::value<std::string>()->default_value( "" ), "exact V" )
    ( "Texact", po::value<std::string>()->default_value( "" ), "exact T" );

  Environment env( _argc=argc, _argv=argv,_desc=options.add(Feel::backend_options("mqsheat")),
		   _about=about(_name="mqsheat",
				_author="Feel++ Consortium",
				_email="feelpp-devel@feelpp.org"));

  int M_verbose = ioption(_name="verbosity");
  
  //Recuperer time frame
  double dt = doption(_name = "ts.time-step");
  Feel::cout << "time-step=" << dt << std::endl;

  double tmax = doption(_name = "ts.time-final");
  Feel::cout << "time-final=" << tmax << std::endl;

  // Eventually get a solution
  bool Uexact = false;
  bool UexactT = false;

  std::string Aexact_s = soption(_name = "Aexact");
  std::string Vexact_s = soption(_name = "Vexact");
  std::string Texact_s = soption(_name = "Texact");

  if ( !Aexact_s.empty() && !Vexact_s.empty() )
    {
      Uexact = true;
      Feel::cout << "* Aexact=" << Aexact_s << std::endl;
      Feel::cout << "* Vexact=" << Vexact_s << std::endl;
    }

  if ( !Texact_s.empty() )
    {
      UexactT = true;
      Feel::cout << "* Texact=" << Texact_s << std::endl;
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

  auto M_materials = M_modelProps->materials().materialWithPhysic(std::vector<std::string>({"electric","thermo-electric","heat"}));
  std::vector<std::string> range;
  for( auto const& mp : M_materials )
    for (auto const& marker : mp.second.meshMarkers() )
      range.push_back(marker);
  Feel::cout << "Electric Materials markers: " << range << std::endl;
  
  // Define SpaceFunctions
  tic();
  auto Ah = Pchv<1>( mesh );
  auto Vh = Pch<1>( mesh, markedelements(mesh, range) );
  auto Th = Pch<1>( mesh );

  auto cond_mesh = Vh->mesh();
  if (Environment::worldComm().isMasterRank())
    {
      std::cout << "mesh->numGlobalElements() "<< mesh->numGlobalElements() << std::endl;
      std::cout << "cond_mesh->numGlobalElements() "<< cond_mesh->numGlobalElements() << std::endl;
      std::cout << "Ah->nDof() "<<Ah->nDof() << std::endl;
      std::cout << "Vh->nDof() "<<Vh->nDof() << std::endl;
    }
#if 1
  auto Jh = Pdhv<0>( mesh, markedelements(mesh, range) );
#endif
  auto Bh = Pdhv<0>( mesh );

  toc("define space functions", (M_verbose > 0));

  // init solutions
  tic();
  auto A0 = expr<3, 1>(soption(_name="A0"));
  auto V0 = expr(soption(_name="V0"));
  auto T0 = expr(soption(_name="T0"));

  auto A = Ah->elementPtr(); //Ah->element(A0); // how to init A to A0?;
  auto V = Vh->elementPtr(); //Vh->element(V0);
  auto T = Th->element(T0);

  auto A0e = Ah->element();
  auto V0e = Vh->element();
  auto T0e = Th->element();

  (*A) = project(_space = Ah, _expr = A0);
  (*V) = project(_space = Vh, _expr = V0);

  auto Aold = (*A);
  auto Vold = (*V);
  toc("init solutions", (M_verbose > 0));

    
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
  toc("create Algebric blockforms", (M_verbose > 0));

  auto mybackend = backend(_name="mqsheat");

  //double t = 0;

  double L2Aexact, H1Aerror, L2Aerror;
  double L2Vexact, H1Verror, L2Verror;
  double L2Texact, H1Terror, L2Terror;

  auto Aexact = Ah->element();
  auto Vexact = Vh->element();
  auto Texact = Th->element();

  auto Aexact_g = expr<3, 1>("{0,0,0}");
  auto Vexact_g = expr("0");
  auto Texact_g = expr("0");

  auto mybdfA = bdf(_space = Ah, _name="mybdfA");
  auto mybdfV = bdf(_space = Vh, _name="mybdfV");
  auto mybdfT = bdf(_space = Th, _name="mybdfT");

  for (auto time : mybdfA -> priorTimes() )
  {
    if (Environment::worldComm().isMasterRank())
    {
      std::cout << "Initialize prior times (from timeInitial()) : " << time.second << "s index: " << time.first << "\n";
    }
    A0.setParameterValues({{"t",time.second}});
    A0e = project(_space=Ah, _expr=A0);
    mybdfA->setUnknown(time.first,A0e);
    //mybdfA->setUnknown(time.first,*A);
  }
  for (auto time : mybdfV -> priorTimes() )
  {
    if (Environment::worldComm().isMasterRank())
    {
      std::cout << "Initialize prior times (from timeInitial()) : " << time.second << "s index: " << time.first << "\n";
    }
    V0.setParameterValues({{"t",time.second}});
    V0e = project(_space=Vh, _expr=V0);
    mybdfV->setUnknown(time.first,V0e);
  }
  for (auto time : mybdfT -> priorTimes() )
  {
    if (Environment::worldComm().isMasterRank())
    {
      std::cout << "Initialize prior times (from timeInitial()) : " << time.second << "s index: " << time.first << "\n";
    }
    T0.setParameterValues({{"t",time.second}});
    T0e = project(_space=Th, _expr=T0);
    mybdfT->setUnknown(time.first,T0e);
  }

  if ( Uexact )
  {
    tic();
    Aexact_g = expr<3, 1>(Aexact_s);
    Aexact_g.setParameterValues({{"t", mybdfA->timeInitial()}});
    Aexact = project(_space = Ah, _expr = Aexact_g);
    Feel::cout << "Define Aexact" << std::endl;
    (*A) = Aexact;
    
    Vexact_g = expr(Vexact_s);
    Vexact_g.setParameterValues({{"t", mybdfV->timeInitial()}});
    Vexact = project(_space = Vh, _expr = Vexact_g);
    Feel::cout << "Define Vexact" << std::endl;
    (*V) = Vexact;

    L2Aexact = normL2(_range = elements(mesh), _expr = Aexact_g);
    H1Aerror = 0;
    L2Aerror = 0;
    L2Vexact = normL2(_range = elements(cond_mesh), _expr = Vexact_g);
    H1Verror = 0;
    L2Verror = 0;
    toc("init A and V exact solution", (M_verbose > 0));
  }

  if(UexactT)
  {
    tic();
    Texact_g = expr(Texact_s);
    Texact_g.setParameterValues({{"t", mybdfT->timeInitial()}});
    Texact = project(_space = Th, _expr = Texact_g);
    Feel::cout << "Define Texact" << std::endl;

    L2Texact = normL2(_range = elements(cond_mesh), _expr = Texact_g);
    H1Verror = 0;
    L2Verror = 0;
    toc("init T exact solution", (M_verbose > 0));
  }
  
  // Compute Magnetic Field
  tic();
  node_type pt(3);
  pt[0] = 0.; pt[1] = 0.; pt[2] = 0.;
  auto M_B = Bh->element();
  M_B = vf::project(_space=Bh, _range=elements(mesh), _expr=curlv(A));
  auto val = M_B(pt);
  auto Bx = val(0,0,0); // evaluation de Bx
  auto By = val(1,0,0); // evaluation de By
  auto Bz = val(2,0,0); // evaluation de Bz
#if 0
  node_type vpt(3);
  vpt[0] = 0.; vpt[1] = 87.5e-3; vpt[2] = 0.;
  auto Vval = (*V)(vpt);
#endif

//test T(0,0,0)
  node_type Tpt(3);
  Tpt[0] = 0.; Tpt[1] = 87.5e-3; Tpt[2] = 0.;
  node_type Tpt0(3);
  Tpt0[0] = 0.; Tpt0[1] = 0; Tpt0[2] = 0.;
  auto Tval = T(Tpt); 
  auto Tval0 = T(Tpt0);
//
  Feel::cout << "t=" << mybdfA->timeInitial() << ", ";
  Feel::cout << "B(" << pt[0] << "," << pt[1] << "," << pt[2] << ") = {" << Bx << "," << By << "," << Bz << "}, ";
  //Feel::cout << "V(" << pt[0] << "," << pt[1] << "," << pt[2] << ")=" << Vval(0,0,0);
  Feel::cout << std::endl;
  toc("compute induction field", (M_verbose > 0));

  tic();
  auto e = exporter( _mesh=mesh );

  e->step(0)->add("A", A);
  e->step(0)->add("V", V);
  e->step(0)->add("T", T0);
  e->step(0)->add("B", M_B);

  if ( Uexact )
  {
    e->step(0)->add("Aexact", Aexact);
    e->step(0)->add("Vexact", Vexact);
  }

  if ( UexactT )
  {
    e->step(0)->add("Texact", Texact);
  }

  Aold = (*A);
  Vold = (*V);

auto bdfA_poly = mybdfA->polyDeriv();
#if 1
  // Feel::cout << "Compute Current density" << std::endl;
  auto J_cond = Jh->element();
  auto J_induct = Jh->element();
  for( auto const& pairMat : M_materials )
  {
    auto name = pairMat.first;
    auto material = pairMat.second;

    auto sigma = material.getScalar("sigma");
    // Feel::cout << "Material:" << material.meshMarkers() << " ";
	  
    J_cond += vf::project( _space=Jh, _range=markedelements(cond_mesh, material.meshMarkers()),
  			     _expr=-sigma * trans(gradv(V)) );
    // Feel::cout << "J_cond:" << material.meshMarkers() << " ";
	  
    J_induct += vf::project( _space=Jh, _range=markedelements(cond_mesh, material.meshMarkers()),
  			       _expr=-sigma * (mybdfA->polyDerivCoefficient(0) * idv(A)-idv(bdfA_poly)) );
    // Feel::cout << "J_induct:" << material.meshMarkers() << std::endl;
  }
  e->step(0)->add( "Jcond", J_cond );
  e->step(0)->add( "Jinduct", J_induct );
  e->step(0)->add( "J", idv(J_cond)+idv(J_induct) );
#endif
  e->save();
  toc("export init solution", (M_verbose > 0));
  
  auto mu0 = 4.e-7 * M_PI ; // SI Unit : H/m = m.kg/s2/A2


  std::string Vname[8];
  std::string Vfirst[8];
  std::string Bname;
  int ii = 0;
  int firstStep = 0; 


  std::ofstream ofile;
  ofile.open("data.csv");

  mybdfA->start();
  mybdfV->start();
  mybdfT->start();

  auto a1 = form2(_trial = Th, _test = Th);
  auto l1 = form1(_test=Th);

  for (double t = dt; mybdfA->isFinished() == false; )
  {
    auto bdfA_poly = mybdfA->polyDeriv();
    tic();
    auto M00 = form2( _trial=Ah, _test=Ah ,_matrix=M, _rowstart=0, _colstart=0 ); 
    for( auto const& pairMat : M_modelProps->materials() )
	  {
	    auto name = pairMat.first;
	    auto material = pairMat.second;

	    auto mur = material.getScalar("mu_mag");

	    // Ampere law: sigma dA/dt + rot(1/(mu_r*mu_0) rotA) + sigma grad(V) = Js
	    // M00 += integrate( _range=markedelements(mesh, material.meshMarkers()),
	    // 		  _expr = dt * 1/mur * inner(curl(A) , curlt(A)) );
	  
	    M00 += integrate( _range=markedelements(mesh, material.meshMarkers()),
			      _expr = 1/mur * trace(trans(gradt(A))*grad(A)) );
	    //Feel::cout << "create lhs(0,0):" << material.meshMarkers() << std::endl;

	  }

    auto M01 = form2( _trial=Vh, _test=Ah ,_matrix=M, _rowstart=0, _colstart=1 );
    auto F0 = form1( _test=Ah, _vector=F, _rowstart=0 );

    auto M11 = form2( _trial=Vh, _test=Vh ,_matrix=M, _rowstart=1, _colstart=1 );
    auto M10 = form2( _trial=Ah, _test=Vh ,_matrix=M, _rowstart=1, _colstart=0 );
    auto F1 = form1( _test=Vh ,_vector=F, _rowstart=1 );
      
    for( auto const& pairMat : M_materials )
	  {
	    auto name = pairMat.first;
	    auto material = pairMat.second;

	    auto sigma = material.getScalar("sigma");

	    // Ampere law: sigma dA/dt + rot(1/(mu_r*mu_0) rotA) + sigma grad(V) = Js
	    M00 += integrate( _range=markedelements(mesh, material.meshMarkers()),
			      _expr = mu0 * mybdfA->polyDerivCoefficient(0) * sigma * inner(id(A) , idt(A) ));
	    //Feel::cout << "create lhs(0,0):" << material.meshMarkers() << std::endl;

	    M01  += integrate(_range=markedelements(mesh, material.meshMarkers()),
			      _expr = mu0 * sigma * inner(id(A),trans(gradt(V))) );
	    //Feel::cout << "create lhs(0,1)" << std::endl;

	    F0 += integrate(_range=markedelements(mesh, material.meshMarkers()),
			    _expr = mu0 * sigma * inner(id(A) , idv(bdfA_poly)));
	    //Feel::cout << "create rhs(0)" << std::endl;

	    // auto Js = ;
	    // F0 += integrate(_range=markedelements(cond_mesh, material.meshMarkers()),
	    // 		 _expr = dt * mu0 * inner(id(A) , Js));
	    // Feel::cout << "create rhs(0)" << std::endl;

	    // Current conservation: div( -sigma grad(V) -sigma*dA/dt) = Qs
	  
	    M11  += integrate( _range=markedelements(cond_mesh, material.meshMarkers()),
			      _expr = mu0 * sigma *inner(gradt(V), grad(V)) );
	    //Feel::cout << "create lhs(1,1)" << std::endl;

	  
	    M10  += integrate( _range=markedelements(cond_mesh, material.meshMarkers()),
			      _expr = mu0 * sigma * mybdfV->polyDerivCoefficient(0) * inner(idt(A), trans(grad(V))) );
	    //Feel::cout << "create lhs(1,0)" << std::endl;

	    F1 += integrate( _range=markedelements(cond_mesh, material.meshMarkers()),
			    _expr = mu0 * sigma * inner(idv(bdfA_poly), trans(grad(V))) );
	    //Feel::cout << "create rhs(1)" << std::endl;

	    // auto Qs = ...;
	    // F1 += integrate(_range=markedelements(cond_mesh, material.meshMarkers()),
	    // 		 _expr = dt * Qs * id(V);
	    // Feel::cout << "create row(1)" << std::endl;
	  }
    toc("assembling", (M_verbose > 0));
     
    tic();
    // Implement Dirichlet fort
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
		      g.setParameterValues({{"t", mybdfA->time()}});
		      //Feel::cout << "A Dirichlet[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
		      M00 += on(_range=markedfaces(mesh,marker), _rhs=F, _element=*A, _expr= g);
		    }
	    }
	    itType = mapField.find( "DirichletX" );
	    if ( itType != mapField.end() )
	    {
	      for ( auto const& exAtMarker : (*itType).second )
		    {
		      std::string marker = exAtMarker.marker();
		      auto g = expr(exAtMarker.expression());
		      g.setParameterValues({{"t", mybdfA->time()}});
		      //Feel::cout << "A DirichletX[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
		      M00 += on(_range=markedfaces(mesh,marker), _rhs=F, _element=(*A)[Component::X], _expr= g);
		    }
	    }
	    itType = mapField.find( "DirichletY" );
	    if ( itType != mapField.end() )
	    {
	      for ( auto const& exAtMarker : (*itType).second )
		    {
		      std::string marker = exAtMarker.marker();
		      auto g = expr(exAtMarker.expression());
		      g.setParameterValues({{"t",mybdfA->time()}});
		      //Feel::cout << "A DirichletY[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
		      M00 += on(_range=markedfaces(mesh,marker), _rhs=F, _element=(*A)[Component::Y], _expr= g);
		    }
	    }
	    itType = mapField.find( "DirichletZ" );
	    if ( itType != mapField.end() )
	    {
	      for ( auto const& exAtMarker : (*itType).second )
		    {
		      std::string marker = exAtMarker.marker();
		      auto g = expr(exAtMarker.expression());
		      g.setParameterValues({{"t", mybdfA->time()}});
		      //Feel::cout << "A DirichletZ[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
		      M00 += on(_range=markedfaces(mesh,marker), _rhs=F, _element=(*A)[Component::Z], _expr= g);
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
		      g.setParameterValues({{"t", mybdfV->time()}});
		      // Feel::cout << "V[" << marker << "]=" << g.evaluate()(0,0) << ", ";
		      M11 += on(_range=markedfaces(cond_mesh,marker), _rhs=F, _element=*V, _expr= g);
		    }
	    }
	  }       
    toc("boundary conditions", (M_verbose > 0));
    
    /* Solve */
    tic();
    auto result = mybackend->solve( _matrix=M, _rhs=F, _solution=U, _rebuild=true);
    std::string msg = (boost::format("[MQS %2%] t=%1% NbIter=%3% Residual=%4%") % mybdfA->time()
			 % soption("mqsheat.pc-type")
			 % result.nIterations()
			 % result.residual()).str();
    if (result.isConverged())
	  {
	    Feel::cout << tc::green << msg << tc::reset << " "; // << std::endl;
	  }
    else
	  {
	    std::string errmsg = msg + " Failed to converge";
	    throw std::logic_error( errmsg );
	  }
      
    toc("solve", (M_verbose > 0));

    // update A and V pointers from U
    myblockVecSol.localize(U);

    // Display Magnetic Field
    M_B = vf::project(_space=Bh, _range=elements(mesh), _expr=curlv(A));
    val = M_B(pt);
    Bx = val(0,0,0); // evaluation de Bx
    By = val(1,0,0); // evaluation de By
    Bz = val(2,0,0); // evaluation de Bz
#if 0
    Vval = (*V)(vpt);
    Feel::cout << "V(" << pt[0] << "," << pt[1] << "," << pt[2] << ")=" << Vval(0,0,0);
    Feel::cout << std::endl;
#endif
    tic();
    e->step(mybdfA->time())->add( "A", A);
    e->step(mybdfV->time())->add( "V", V);
      
    e->step(mybdfA->time())->add( "B", M_B );
    // M_gradV = vf::project(_space=Jh, _range=elements(cond_mesh), _expr=trans(gradv(V))); // breaks in // why?
    // e->step(t)->add( "E", M_gradV );

    // Update current densities
#if 1
    J_cond = vf::project(_space=Jh, _range=elements(cond_mesh), _expr=expr<3, 1>("{0,0,0}")); //Jh->element();
    J_induct = vf::project(_space=Jh, _range=elements(cond_mesh), _expr=expr<3, 1>("{0,0,0}")); //Jh->element();
    for( auto const& pairMat : M_materials )
	  {
	    auto name = pairMat.first;
	    auto material = pairMat.second;

	    auto sigma = material.getScalar("sigma");
	    J_cond += vf::project( _space=Jh, _range=markedelements(cond_mesh, material.meshMarkers()),
				        _expr=-sigma * trans(gradv(V)) );
	    J_induct += vf::project( _space=Jh, _range=markedelements(cond_mesh, material.meshMarkers()),
				          _expr=-sigma * (mybdfA->polyDerivCoefficient(0)*idv(A)-idv(bdfA_poly)) );
	  }
    e->step(mybdfA->time())->add( "Jcond", J_cond );
    e->step(mybdfA->time())->add( "Jinduct", J_induct );
    e->step(mybdfA->time())->add( "J", idv(J_cond)+idv(J_induct) );
#endif
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
		      g.setParameterValues({{"t", mybdfV->time()}});
		      Feel::cout << "V[" << marker << "]=" << g.evaluate()(0,0) << ", ";
          Vname[ii] = std::to_string(g.evaluate()(0,0));
          ii ++;
#if 1
		      double I = integrate( markedfaces( cond_mesh, marker ), inner(idv(J_induct),N()) + inner(idv(J_cond),N()) ).evaluate()(0,0);
#endif		      
          Feel::cout << "I[" << marker << "]=" << I << ", ";
          Vname[ii] = std::to_string(I);
          ii ++;
          if (firstStep == 0)
          {
            Vfirst[ii-2] = "V[" + marker + "]";
            Vfirst[ii-1] = "I[" + marker + "]";;
          }
		    }
	    }
	  }

    Feel::cout << " B(" << pt[0] << "," << pt[1] << "," << pt[2] << ") = {" << Bx << "," << By << "," << Bz << "}";
    Feel::cout << std::endl;

    if(firstStep == 0)
    {
      if (ii == 4)
      { 
        if (ofile)
        {
          ofile << std::setprecision(10) << "t,NbIter,Residual," << Vfirst[0] << "," 
                << Vfirst[1] << "," << Vfirst[2] << "," << Vfirst[3] << "," << "Bz," << "T" << std::endl;
        }
      }
      else
      {
        if (ofile)
        {
          ofile << std::setprecision(10) << "t,NbIter,Residual," << Vfirst[0] << "," 
                 << Vfirst[1] << "," << Vfirst[2] << "," << Vfirst[3] << "," << Vfirst[4] 
                 << "," << Vfirst[5] << "," << Vfirst[6] << "," << Vfirst[7] << "," << "Bz,T" << std::endl;
        }
      }
      firstStep = 1;
    }


    //Bname = "{"+std::to_string(Bx) + "," + std::to_string(By) + "," + std::to_string(Bz) + "}";
    Bname = std::to_string(Bz);
    if (ii == 4)
    {
      if (ofile)
      {
        ofile << std::setprecision(10) << mybdfA->time() << "," << result.nIterations() << "," << result.residual() 
             << "," << Vname[0] << "," << Vname[1] << "," << Vname[2] << "," << Vname[3] 
             << "," << Bname;
      }
    }
    else
    {
      if (ofile)
      {
        ofile << std::setprecision(10) << mybdfA->time() << "," << result.nIterations() << "," << result.residual() 
             << "," << Vname[0] << "," << Vname[1] << "," << Vname[2] << "," << Vname[3] << ","
             << Vname[4] << "," << Vname[5] << "," << Vname[6] << "," << Vname[7] << "," 
             << Bname;
      }
    }
    ii = 0;

    if ( Uexact )
	  {
	    Aexact_g.setParameterValues({{"t", mybdfA->time()}});
	    Aexact = project(_space = Ah, _expr = Aexact_g);
	    Vexact_g.setParameterValues({{"t", mybdfV->time()}});
	    Vexact = project(_space = Vh, _expr = Vexact_g);

	    e->step(mybdfA->time())->add( "Aexact", Aexact);
	    e->step(mybdfV->time())->add( "Vexact", Vexact);
	  }
    
    if ( UexactT )
	  {
	    Texact_g.setParameterValues({{"t", mybdfT->time()}});
	    Texact = project(_space = Th, _expr = Texact_g);

	    e->step(mybdfT->time())->add( "Texact", Texact);
	  }


//********************** T calcul *********************
    auto bdfT_poly = mybdfT->polyDeriv();

    for( auto const& pairMat : M_modelProps->materials() )
	  {
	    auto name = pairMat.first;
	    auto material = pairMat.second;

      auto k = material.getScalar("k");
	    auto rho = material.getScalar("rho");
      auto Cp = material.getScalar("Cp");
        
      //heat 
	    a1 += integrate( _range=markedelements(mesh, material.meshMarkers()),
			      _expr = k * inner( gradt(T),grad(T) ) );

	    // heat
	    a1 += integrate( _range=markedelements(mesh, material.meshMarkers()),
			      _expr = Cp * rho * mybdfT->polyDerivCoefficient(0) * id(T) * idt(T) );

	    l1 += integrate(_range=markedelements(mesh, material.meshMarkers()),
			    _expr = Cp * rho * id(T) * idv(bdfT_poly) );
	  }

    for( auto const& pairMat : M_materials )
	  {
	    auto name = pairMat.first;
	    auto material = pairMat.second;

      auto sigma = material.getScalar("sigma");

      //source term from mqs sigma * ||E||^2 = 1/sigma ||J||^2
	    l1 += integrate(_range=markedelements(mesh, material.meshMarkers()),
	                    _expr = (1/sigma) * id(T) * inner(idv(J_cond)+idv(J_induct),idv(J_cond)+idv(J_induct)) );
	  }

    itField = M_modelProps->boundaryConditions().find( "temperature");
    if ( itField != M_modelProps->boundaryConditions().end() )
	  {
	    auto mapField = (*itField).second;

      auto itType = mapField.find( "Neumann" );
	    if ( itType != mapField.end() )
	    {
	      for ( auto const& exAtMarker : (*itType).second )
	     	{
	     		std::string marker = exAtMarker.marker();
	     		auto g = expr(exAtMarker.expression());
	        g.setParameterValues({{"t", mybdfT->time()}});
	       	Feel::cout << "Neuman[" << marker << "] : " << exAtMarker.expression() << std::endl;
	        l1 += integrate(_range=markedfaces(mesh,marker), 
                           _expr=  - g * id(T) );
        }
      }

      itType = mapField.find( "Robin" );
	    if ( itType != mapField.end() )
	    {
	      for ( auto const& exAtMarker : (*itType).second )
	     	{
	     		std::string marker = exAtMarker.marker();
          auto h = expr(exAtMarker.expression1());
	     		auto Tw = expr(exAtMarker.expression2());
	        Tw.setParameterValues({{"t", mybdfT->time()}});
          h.setParameterValues({{"t", mybdfT->time()}});
	       	//Feel::cout << "Robin[" << marker << "] : " << exAtMarker.expression1() << std::endl;
          //Feel::cout << "Robin[" << marker << "] : " << exAtMarker.expression2() << std::endl;
	        a1 += integrate(_range=markedfaces(mesh,marker), 
                           _expr= h * idt(T) * id(T) );
          l1 += integrate(_range=markedfaces(mesh,marker), 
                          _expr= h * Tw * id(T) );
        }
      }
	    itType = mapField.find( "Dirichlet" );
	    if ( itType != mapField.end() )
	    {
	      for ( auto const& exAtMarker : (*itType).second )
		    {
		      std::string marker = exAtMarker.marker();
		      auto g = expr(exAtMarker.expression());
		      g.setParameterValues({{"t", mybdfT->time()}});
		      Feel::cout << "T Dirichlet[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
		      a1 += on(_range=markedfaces(mesh,marker), _rhs=F, _element=T, _expr= g);
		    }
	    }      
	  }

    a1.solve(_rhs = l1, _solution = T);
    e->step(mybdfT->time())->add("T",T);

    Tval = T(Tpt);
    Feel::cout << "t = " << mybdfT->time() << std::endl;
    Feel::cout << "T(" << Tpt[0] << "," << Tpt[1] << "," << Tpt[2] << ") = " << Tval(0,0,0) << std::endl;
    Feel::cout << "T(" << Tpt0[0] << "," << Tpt0[1] << "," << Tpt0[2] << ") = " << Tval0(0,0,0) << std::endl;

    if (ofile)
    {
      ofile << "," << Tval(0,0,0) << std::endl;
    }
//*****************************************************

    e->save();
    toc("export", (M_verbose > 0));

    // Compute error
    if ( Uexact )
	  {
	    L2Aexact = normL2(_range = elements(mesh), _expr = Aexact_g);
	    L2Aerror = normL2(elements(mesh), (idv(A) - idv(Aexact)));
	    H1Aerror = normH1(elements(mesh), _expr = (idv(A) - idv(Aexact)), _grad_expr = (gradv(A) - gradv(Aexact)));
	    Feel::cout << "error: " << "t="<< mybdfA->time();
	    Feel::cout << " A: " << L2Aerror << " " << L2Aerror / L2Aexact << " " << H1Aerror << " ";

  	  L2Vexact = normL2(_range = elements(cond_mesh), _expr = Vexact_g);
  	  L2Verror = normL2(elements(cond_mesh), (idv(V) - idv(Vexact)));
	    H1Verror = normH1(elements(cond_mesh), _expr = (idv(V) - idv(Vexact)), _grad_expr = (gradv(V) - gradv(Vexact)));
	    Feel::cout << " V: " << L2Verror << " " << L2Verror / L2Vexact << " " << H1Verror << std::endl;
	  }
    
    if(UexactT)
    {
  	  L2Texact = normL2(_range = elements(mesh), _expr = Texact_g);
  	  L2Terror = normL2(elements(mesh), (idv(T) - idv(Texact)));
	    H1Terror = normH1(elements(mesh), _expr = (idv(T) - idv(Texact)), _grad_expr = (gradv(T) - gradv(Texact)));
	    Feel::cout << " T: " << L2Terror << " " << L2Terror / L2Texact << " " << H1Terror << std::endl;
	  }
    mybdfA->next(*A);
    mybdfV->next(*V);
    mybdfT->next(T);

    /* reinit  */
    M->zero();
    F->zero();

    l1.zero();
    a1.zero();

    Aold = (*A);
    Vold = (*V);

  }
  //std::cout << mqs << std::endl;
  ofile.close();
}
