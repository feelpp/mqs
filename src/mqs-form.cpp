#include <tabulate/table.hpp>
#include <tabulate/markdown_exporter.hpp>
using namespace tabulate;

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

#include <feel/feelmodels/maxwell/biotsavart.hpp>

int main(int argc, char**argv )
{

  using namespace Feel;
  po::options_description options( "MQS options" );
  options.add_options()
    ( "model-file", Feel::po::value<std::string>()->default_value( "" ), "file describing model properties")
    ( "adaptive", po::value<bool>()->default_value( false ), "activate dt apdative scheme" )
    ( "dttol", po::value<double>()->default_value( 0. ), "dt tolerance" )
    ( "forced-sequence", po::value< std::vector<double> >()->default_value(std::vector<double>()), "list of forced times" )
    ( "verbosity", po::value<int>()->default_value( 0 ), "set verbosisity level" )
    ( "weakdir", po::value<bool>()->default_value( "false" ), "use Dirichlet weak formulation" )
    ( "penalty-coeff", po::value<double>()->default_value( 1.e+3 ), "penalty coefficient for weak Dirichlet" )
    ( "A0", po::value<std::string>()->default_value( "{0,0,0}" ), "initial A" )
    ( "V0", po::value<std::string>()->default_value( "0" ), "initial V" )
    ( "Aexact", po::value<std::string>()->default_value( "" ), "exact A" )
    ( "Vexact", po::value<std::string>()->default_value( "" ), "exact V" );

  Environment env( _argc=argc, _argv=argv,_desc=options.add(Feel::backend_options("mqs")).add(Feel::biotsavart_options()),
		   _about=about(_name="mqs",
				_author="Feel++ Consortium",
				_email="feelpp-devel@feelpp.org"));

  int M_verbose = ioption(_name="verbosity");
  
  //Recuperer time frame
  double dt = doption(_name = "ts.time-step");
  Feel::cout << "time-step=" << dt << std::endl;

  double tmax = doption(_name = "ts.time-final");
  Feel::cout << "time-final=" << tmax << std::endl;

  double dttol = doption(_name="dttol");
  if ( boption("adaptive") && dttol == 0)
    dttol = dt/100.;
  Feel::cout << "time-dttol=" << dttol << std::endl;
  double dtprev = dt;

  double dt_min = dttol/100.;
  double dt_max = dt*100.;
  Feel::cout << "time-dtmin=" << dt_min << std::endl;
  Feel::cout << "time-dtmax=" << dt_max << std::endl;
  
  // Eventually get a solution
  bool Uexact = false;

  std::string Aexact_s = soption(_name = "Aexact");
  std::string Vexact_s = soption(_name = "Vexact");

  if ( !Aexact_s.empty() && !Vexact_s.empty() )
    {
      Uexact = true;
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
  std::vector<std::string> range_conductors;
  for( auto const& mp : M_materials )
    for (auto const& marker : mp.second.meshMarkers() )
      range_conductors.push_back(marker);
  // Feel::cout << "Electric Materials markers: " << range_conductors << std::endl;
  std::set<std::string> conductors(std::begin(range_conductors), std::end(range_conductors));
  Feel::cout << "Electric Materials markers (set): " << conductors << std::endl;
  
  // Define SpaceFunctions
  tic();
  auto Ah = Pchv<1>( mesh );
  auto Vh = Pch<1>( mesh, markedelements(mesh, range_conductors) );

  auto cond_mesh = Vh->mesh();
  if (Environment::worldComm().isMasterRank())
    {
      std::cout << "mesh->numGlobalElements() "<< mesh->numGlobalElements() << std::endl;
      std::cout << "cond_mesh->numGlobalElements() "<< cond_mesh->numGlobalElements() << std::endl;
      std::cout << "Ah->nDof() "<<Ah->nDof() << std::endl;
      std::cout << "Vh->nDof() "<<Vh->nDof() << std::endl;
    }

  auto Jh = Pdhv<0>( mesh, markedelements(mesh, range_conductors) );
  auto Bh = Pdhv<0>( mesh );

  auto A = Ah->elementPtr(); //Ah->element(A0); // how to init A to A0?;
  auto V = Vh->elementPtr(); //Vh->element(V0);
  toc("define space functions", (M_verbose > 0));

  // init solutions
  tic();
  auto A0 = expr<3, 1>(soption(_name="A0"));
  auto V0 = expr(soption(_name="V0"));
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

  auto mybackend = backend(_name="mqs");

  double t = 0;
  double epsNL = 1.e-3;
  double errorNL, normA;
  double Residual;
  int nIterations;
  
  double L2Aexact, H1Aerror, L2Aerror;
  double L2Vexact, H1Verror, L2Verror;

  auto Aexact = Ah->element();
  auto Vexact = Vh->element();

  auto Aexact_g = expr<3, 1>("{0,0,0}");
  auto Vexact_g = expr("0");
  if ( Uexact )
    {
      tic();
      Aexact_g = expr<3, 1>(Aexact_s);
      Aexact_g.setParameterValues({{"t", t}});
      Aexact = project(_space = Ah, _expr = Aexact_g);
      Feel::cout << "Define Aexact" << std::endl;
      (*A) = Aexact;
      
      Vexact_g = expr(Vexact_s);
      Vexact_g.setParameterValues({{"t", t}});
      Vexact = project(_space = Vh, _expr = Vexact_g);
      Feel::cout << "Define Vexact" << std::endl;
      (*V) = Vexact;

      L2Aexact = normL2(_range = elements(mesh), _expr = Aexact_g);
      H1Aerror = 0;
      L2Aerror = 0;
      L2Vexact = normL2(_range = elements(cond_mesh), _expr = Vexact_g);
      H1Verror = 0;
      L2Verror = 0;
      toc("init exact solution", (M_verbose > 0));
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
  Feel::cout << "t=" << t << ", ";
  Feel::cout << "B(" << pt[0] << "," << pt[1] << "," << pt[2] << ") = {" << Bx << "," << By << "," << Bz << "}, ";
  //Feel::cout << "V(" << pt[0] << "," << pt[1] << "," << pt[2] << ")=" << Vval(0,0,0);
  Feel::cout << std::endl;
  toc("compute induction field", (M_verbose > 0));

  tic();
  auto e = exporter( _mesh=mesh );

  e->step(t)->add("A", A);
  e->step(t)->add("V", V);
  e->step(t)->add("B", M_B);
  
  // Feel::cout << "Compute Electric Field" << std::endl;
  // auto M_gradV = Jh->element(); 
  // M_gradV = vf::project(_space=Jh, _range=elements(cond_mesh), _expr=trans(gradv(V))); // breaks in // why?
  // e->step(t)->add("E", M_gradV);

  if ( Uexact )
    {
      e->step(t)->add("Aexact", Aexact);
      e->step(t)->add("Vexact", Vexact);
    }

  Aold = (*A);
  Vold = (*V);

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
  			       _expr=-sigma * (idv(A)-idv(Aold))/dt );
      // Feel::cout << "J_induct:" << material.meshMarkers() << std::endl;
    }
  e->step(t)->add( "Jcond", J_cond );
  e->step(t)->add( "Jinduct", J_induct );
  e->step(t)->add( "J", idv(J_cond)+idv(J_induct) );

  e->save();
  toc("export init solution", (M_verbose > 0));
  
  auto mu0 = 4.e-7 * M_PI ; // SI Unit : H/m = m.kg/s2/A2

  Table mqs;
  std::string Vname[8];
  std::string Vfirst[8];
  std::string Bname;
  int ii = 0;
  int firstStep = 0;

  int iterNL = 0;
  int maxiterNL = 10;
  double initResidual;

  // define sequence of forced time steps
  double epstime = 1.e-3;
  bool reached = false;
  std::vector<double> forced_times;
  forced_times.push_back(0.1);
  forced_times.push_back(0.5);
  forced_times.push_back(tmax);
  Feel::cout << "Forced sequence:" << forced_times << std::endl;

  int n_forced = 0;
  double forced_t = forced_times[n_forced];
  
  for(t = dt; t <= tmax+1e-10; )
    {
      // Update current densities
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
				   _expr=-sigma * (idv(A)-idv(Aold))/dt );
	}

      // if BiotSavart Bc shall loop until ||A-Aold||<eps
      bool nonlinear = false;

      tic();
      do {
	tic();

	auto Anl = (*A);
	auto Vnl = (*V);

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
			      _expr = dt * 1/mur * trace(trans(gradt(A))*grad(A)) );
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
			      _expr = mu0 * sigma * inner(id(A) , idt(A) ));
	    //Feel::cout << "create lhs(0,0):" << material.meshMarkers() << std::endl;

	    M01  += integrate(_range=markedelements(mesh, material.meshMarkers()),
			      _expr = dt * mu0 * sigma * inner(id(A),trans(gradt(V))) );
	    //Feel::cout << "create lhs(0,1)" << std::endl;

	    F0 += integrate(_range=markedelements(mesh, material.meshMarkers()),
			    _expr = mu0 * sigma * inner(id(A) , idv(Aold)));
	    //Feel::cout << "create rhs(0)" << std::endl;

	    // auto Js = ;
	    // F0 += integrate(_range=markedelements(cond_mesh, material.meshMarkers()),
	    // 		 _expr = dt * mu0 * inner(id(A) , Js));
	    // Feel::cout << "create rhs(0)" << std::endl;

	    // Current conservation: div( -sigma grad(V) -sigma*dA/dt) = Qs
	  
	    M11  += integrate( _range=markedelements(cond_mesh, material.meshMarkers()),
			       _expr = mu0 * sigma * dt * inner(gradt(V), grad(V)) );
	    //Feel::cout << "create lhs(1,1)" << std::endl;

	  
	    M10  += integrate( _range=markedelements(cond_mesh, material.meshMarkers()),
			       _expr = mu0 * sigma * inner(idt(A), trans(grad(V))) );
	    //Feel::cout << "create lhs(1,0)" << std::endl;

	    F1 += integrate( _range=markedelements(cond_mesh, material.meshMarkers()),
			     _expr = mu0 * sigma * inner(idv(Aold), trans(grad(V))) );
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
		    g.setParameterValues({{"t", t}});
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
		    g.setParameterValues({{"t", t}});
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
		    g.setParameterValues({{"t", t}});
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
		    g.setParameterValues({{"t", t}});
		    //Feel::cout << "A DirichletZ[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
		    M00 += on(_range=markedfaces(mesh,marker), _rhs=F, _element=(*A)[Component::Z], _expr= g);
		  }
	      }
	    itType = mapField.find( "BiotSavart" );
	    if ( itType != mapField.end() )
	      {
		nonlinear = true;

		for ( auto const& exAtMarker : (*itType).second )
		  {
		    std::string marker = exAtMarker.marker(); // make a list instead
		    std::set<std::string> markers;
		    markers.insert(marker);
		    
		    auto As = BiotSavart<3>(mesh, markers);
		    auto jEx = idv(J_cond)+idv(J_induct);

		    As.compute(jEx, false, true, conductors);
		    
		    //Feel::cout << "A BiotSavart[" << marker << "] : " << std::endl;
		    auto Abc = As.magneticPotential();
		    M00 += on(_range=markedfaces(mesh,marker), _rhs=F, _element=(*A), _expr= idv(Abc) );
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
		    // Feel::cout << "V[" << marker << "]=" << g.evaluate()(0,0) << ", ";
		    M11 += on(_range=markedfaces(cond_mesh,marker), _rhs=F, _element=*V, _expr= g);
		  }
	      }
	  }       
	toc("boundary conditions", (M_verbose > 0));
    
	/* Solve */
	tic();
	auto result = mybackend->solve( _matrix=M, _rhs=F, _solution=U, _rebuild=true);
	std::string msg = (boost::format("[MQS %2%] t=%1% NbIter=%3% Residual=%4%") % t
			   % soption("mqs.pc-type")
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

	if ( iterNL == 0 )
	  initResidual = result.residual();
	
	Residual =  result.residual();
	nIterations = result.nIterations();
	toc("solve", (M_verbose > 0));

	// update A and V pointers from U
	myblockVecSol.localize(U);

	// Update current densities
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
				     _expr=-sigma * (idv(A)-idv(Aold))/dt );
	  }

	// compute errorNL (see V. Chabannes comment for more precise handling)
	if ( nonlinear==true)
	  {
	    errorNL = normL2(_range = elements(mesh), _expr = (idv(A)-idv(Anl)) );
	    normA = normL2(_range = elements(mesh), _expr = idv(A) );
	    Feel::cout << "iterNL=" << iterNL << " ,";
	    Feel::cout << "errorNL=" << errorNL << " ,";
	    Feel::cout << "nomrA=" << normA << ", ";
	    Feel::cout << "epsNL*normA=" << epsNL*normA << ", ";
	    // Feel::cout << std::endl;

	    iterNL++;
	  }
	
      } while ( (nonlinear==true) && (errorNL > epsNL*normA) && (iterNL < maxiterNL) );
      toc("non-linear step", ( (M_verbose > 0) && nonlinear) );

      // reset NL counter
      iterNL = 0;
    
      bool do_export=true;
      
      /* Solve */
      if ( boption( "adaptive") )
        {
	  tic();
	  std::string adapt_msg;
	  // time filtering , get order 2
	  auto filter = [&dt, &dtprev]( auto const& in, auto const& inprev, auto& out ) { 
			  double nu = dt*(dt+dtprev)/(dtprev*(2*dt+dtprev));
			  double c1 = 2*dtprev/(dt+dtprev);
			  double c2 = 2*dt/(dt+dtprev);
			  Feel::cout << "  adaptive time stepping nu=" << nu << " c1=" << c1 << " c2=" << c2 << "; "; //<< std::endl;
			  out.on( _range=elements(out.mesh()), _expr=idv(in)-(nu/2)*(c1*idv(inprev) - 2*idv(in) + c2*idv(inprev) )); 
			};
	  auto Apost = (*A);
	  auto Vpost = (*V);
	  filter( A, Aold, Apost );
	  filter( V, Vold, Vpost );
	  auto estA = normL2( _range=elements(mesh), _expr=idv(A)-idv(Apost));
	  auto estV = normL2( _range=elements(cond_mesh), _expr=idv(V)-idv(Vpost));
	  auto est = std::max( estA, estV );
	  Feel::cout << "est : " << std::scientific << std::setprecision(3) << est << " estA : " << estA << " estV : " << estV << " (dttol=" << dttol << "); ";// << std::endl;

	  Feel::cout << "forced_time=" << forced_t << ", ";
	  Feel::cout << "t=" << t << ", ";
	  Feel::cout << "allmost=" << fabs(1-forced_t/t) << " (" << (fabs(1-forced_t/t) <= epstime) << ") ";
	  Feel::cout << "reached" << reached << std::endl;
	  if ( est > dttol )
	    {  
	      //Feel::cout << "reject (>dttol): dt estimate: " << 0.7 * dt * sqrt(dttol/est);
	      t -= dt;
	      dt/=2.;
	      adapt_msg = "refining(/2) the time step";
	      if ( dt < dt_min )
		{
		  dt = dt_min;
		  adapt_msg = "refining the time step to dt_min";
		}
	      // no export
	      do_export=false;

	      // time rejected
	      if ( reached )
		{
		  reached = false;
		  Feel::cout << "***";
		}
	    }
	  else //if ( est < dttol )
	    {
	      //Feel::cout << "accepted (<ddtol): dt estimate: " << 0.9 * dt * sqrt(dttol/est);
	      
	      dtprev=dt;

	      Aold = (*A);
	      Vold = (*V);

	      // export
	      do_export=true; //false;

	      bool allmost = ( fabs(1-forced_t/t) <= epstime );
	      if ( !allmost )
		{
		  if ( est <= dttol/8. )
		    {
		      dt*=2;
		      adapt_msg = "increasing(x2) the time step";
		      if ( dt > dt_max )
			{
			  dt = dt_max;
			  adapt_msg = "increasing the time step to dt_max";
			}
		    }
		  else
		    {
		      adapt_msg = "keeping the time step";
		    }
		}

	      // time accepted
	      if ( reached || allmost )
		if ( n_forced < forced_times.size()-1 )
		  {
		    reached = false;
		    n_forced++;
		    forced_t = forced_times[n_forced] ;
		    dt = doption(_name = "ts.time-step");
		    dtprev=dt;
		    Feel::cout << "go to next sequence" << std::endl;
		  }
	    }
	   
	  std::string msg = (boost::format("[adapt dt=%1%] ") % dt).str();
	  msg += adapt_msg;
	  Feel:cout << msg << std::endl;
	  toc( msg, (M_verbose > 0));
        }
      else
	{
	  if ( reached )
	    {
	      Feel::cout << "reached[" << forced_t << "]: ";
	      Feel::cout << ( n_forced < forced_times.size()-1 ) << std::endl;
	      if ( n_forced < forced_times.size()-1 )
		{
		  n_forced++;
		  forced_t = forced_times[n_forced];
		  reached = false;
		  dt = doption(_name = "ts.time-step");
		  Feel::cout << "go to next sequence" << std::endl;
		}
	    }
	  Aold = (*A);
	  Vold = (*V);
	}
      
      
      if ( do_export)
	{
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
	  e->step(t)->add( "A", A);
	  e->step(t)->add( "V", V);
      
	  e->step(t)->add( "B", M_B );
	  // M_gradV = vf::project(_space=Jh, _range=elements(cond_mesh), _expr=trans(gradv(V))); // breaks in // why?
	  // e->step(t)->add( "E", M_gradV );

	  e->step(t)->add( "Jcond", J_cond );
	  e->step(t)->add( "Jinduct", J_induct );
	  e->step(t)->add( "J", idv(J_cond)+idv(J_induct) );

	  auto itField = M_modelProps->boundaryConditions().find( "electric-potential");
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
		      Feel::cout << "V[" << marker << "]=" << g.evaluate()(0,0) << ", ";
		      Vname[ii] = std::to_string(g.evaluate()(0,0));
		      ii ++;
		      
		      double I = integrate( markedfaces( cond_mesh, marker ), inner(idv(J_induct),N()) + inner(idv(J_cond),N()) ).evaluate()(0,0);
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
	      Bname = "Bz("+std::to_string(pt[0]) + "," + std::to_string(pt[1]) + "," + std::to_string(pt[2]) + ")";
	      if (ii == 4)
		{
		  mqs.add_row({"t","NbIter","Residual",Vfirst[0],Vfirst[1],Vfirst[2],Vfirst[3],Bname});
		}
	      else
		{
		  mqs.add_row({"t","NbIter","Residual",Vfirst[0],Vfirst[1],Vfirst[2],Vfirst[3]
			       ,Vfirst[4],Vfirst[5],Vfirst[6],Vfirst[7],Bname});
		}
	      firstStep = 1;
	    }


	  //Bname = "{"+std::to_string(Bx) + "," + std::to_string(By) + "," + std::to_string(Bz) + "}";
	  Bname = std::to_string(Bz);
	  if (ii == 4)
	    {
	      mqs.add_row({std::to_string(t),std::to_string(nIterations),
			   std::to_string(Residual),
			   Vname[0],Vname[1],Vname[2],Vname[3],Bname});
	    }
	  else
	    {
	      mqs.add_row({std::to_string(t),std::to_string(nIterations),
			   std::to_string(Residual),
			   Vname[0],Vname[1],Vname[2],Vname[3],Vname[4],Vname[5],Vname[6],Vname[7],Bname});
	    }
	  ii = 0;

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
	  toc("export", (M_verbose > 0));

	  // Compute error
	  if ( Uexact )
	    {
	      L2Aexact = normL2(_range = elements(mesh), _expr = Aexact_g);
	      L2Aerror = normL2(elements(mesh), (idv(A) - idv(Aexact)));
	      H1Aerror = normH1(elements(mesh), _expr = (idv(A) - idv(Aexact)), _grad_expr = (gradv(A) - gradv(Aexact)));
	      Feel::cout << "error: " << "t="<< t;
	      Feel::cout << " A: " << L2Aerror << " " << L2Aerror / L2Aexact << " " << H1Aerror << " ";

	      L2Vexact = normL2(_range = elements(cond_mesh), _expr = Vexact_g);
	      L2Verror = normL2(elements(cond_mesh), (idv(V) - idv(Vexact)));
	      H1Verror = normH1(elements(cond_mesh), _expr = (idv(V) - idv(Vexact)), _grad_expr = (gradv(V) - gradv(Vexact)));
	      Feel::cout << " V: " << L2Verror << " " << L2Verror / L2Vexact << " " << H1Verror << std::endl;
	    }
	}

      /* reinit  */
      M->zero();
      F->zero();

      t += dt;

      // force time step
      if ( !forced_times.empty() && !reached )
	{
	  forced_t = forced_times[n_forced];
	  if ( t > forced_t )
	    {
	      reached = true;
	      if ( fabs(1-forced_t/t) > epstime )
		{
		  dt -= t-forced_t;
		  t = forced_t;
		  Feel::cout << "forced_time=" << forced_t << ", ";
		  Feel::cout << "forced_dt=" << dt << ", ";
		  Feel::cout << "t=" << t << std::endl;
		}
	    }
	}
    }

  // export as markdow table
  MarkdownExporter exporter;
  auto markdown = exporter.dump(mqs);
  Feel::cout << markdown << std::endl;
}
