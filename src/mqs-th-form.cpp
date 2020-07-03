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
  po::options_description options( "MQS Th options" );
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
    ( "T0", po::value<std::string>()->default_value( "273.15" ), "initial T" )
    ( "Aexact", po::value<std::string>()->default_value( "" ), "exact A" )
    ( "Vexact", po::value<std::string>()->default_value( "" ), "exact V" );

  Environment env( _argc=argc, _argv=argv,_desc=options.add(Feel::backend_options("mqs")).add(Feel::biotsavart_options()),
		   _about=about(_name="mqs-th",
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

  // Materials in heat ONLY
  auto M_th_materials = M_modelProps->materials().materialWithPhysic(std::vector<std::string>({"heat"}));
  std::vector<std::string> range_th;
  for( auto const& mp : M_th_materials )
    for (auto const& marker : mp.second.meshMarkers() )
      range_th.push_back(marker);
  // Feel::cout << "Electric Materials markers: " << range_conductors << std::endl;
  std::set<std::string> thdomains(std::begin(range_th), std::end(range_th));
  
  // Define SpaceFunctions
  tic();
  auto Ah = Pchv<1>( mesh );
  auto Vh = Pch<1>( mesh, markedelements(mesh, conductors) );
  auto Th = Pch<1>( mesh, markedelements(mesh, thdomains) );
  
  auto cond_mesh = Vh->mesh();
  auto th_mesh = Th->mesh();
  if (Environment::worldComm().isMasterRank())
    {
      std::cout << "mesh->numGlobalElements() "<< mesh->numGlobalElements() << std::endl;
      std::cout << "cond_mesh->numGlobalElements() "<< cond_mesh->numGlobalElements() << std::endl;
      std::cout << "Ah->nDof() "<<Ah->nDof() << std::endl;
      std::cout << "Vh->nDof() "<<Vh->nDof() << std::endl;
    }

  auto Jh = Pdhv<0>( mesh, markedelements(mesh, conductors) );
  auto Bh = Pdhv<0>( mesh );

  auto A = Ah->elementPtr(); //Ah->element(A0); // how to init A to A0?;
  auto V = Vh->elementPtr(); //Vh->element(V0);
  auto Heat_T = Th->elementPtr(); //Vh->element(V0);
  toc("define space functions", (M_verbose > 0));

  // init solutions
  tic();
  auto A0 = expr<3, 1>(soption(_name="A0"));
  auto V0 = expr(soption(_name="V0"));
  auto T0 = expr(soption(_name="T0"));
  A->on(_range=elements(mesh), _expr=A0); //(*A) = project(_space = Ah, _expr = A0);
  V->on(_range=elements(support(Vh)), _expr=V0); //(*V) = project(_space = Vh, _expr = V0);
  Heat_T->on(_range=elements(support(Vh)), _expr=T0); //(*Heat_T) = project(_space = Th, _expr = T0);

  auto Aold = (*A);
  auto Vold = (*V);
  auto Heat_Told = (*Heat_T);
  toc("init solutions", (M_verbose > 0));
    
  // Vincent way
  tic();
  BlocksBaseGraphCSR myblockGraph(3,3);
  myblockGraph(0,0) = stencil(_test=Ah,_trial=Ah, _diag_is_nonzero=false, _close=false)->graph();
  myblockGraph(0,1) = stencil(_test=Ah,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
  myblockGraph(1,0) = stencil(_test=Vh,_trial=Ah, _diag_is_nonzero=false, _close=false)->graph();
  myblockGraph(1,1) = stencil(_test=Vh,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
  myblockGraph(2,2) = stencil(_test=Th,_trial=Th, _diag_is_nonzero=false, _close=false)->graph();
  auto M = backend()->newBlockMatrix(_block=myblockGraph);

  BlocksBaseVector<double> myblockVec(3);
  myblockVec(0,0) = backend()->newVector( Ah );
  myblockVec(1,0) = backend()->newVector( Vh );
  myblockVec(2,0) = backend()->newVector( Th );
  auto F = backend()->newBlockVector(_block=myblockVec, _copy_values=false);

  BlocksBaseVector<double> myblockVecSol(3);
  myblockVecSol(0,0) = A;
  myblockVecSol(1,0) = V;
  myblockVecSol(2,0) = Heat_T;
  auto U = backend()->newBlockVector(_block=myblockVecSol, _copy_values=false);
  toc("create Algebric blockforms", (M_verbose > 0));

  auto mybackend = backend(_name="mqs");

  double t = 0;
  double epsNL = 1.e-3, epsNL_th = 1.e-5;

  double errorNL, normA;
  double errorNL_th, normT;
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
  toc("compute induction field", (M_verbose > 0));

  auto Tm = minmax( _range=elements(support(Th)), _pset=_Q<2>(), _expr=idv(Heat_T));
  double Heat_Tmax = Tm.max();
  double Heat_Tmin = Tm.min();
  Feel::cout << "Tmin=" << Heat_Tmin << ", ";
  Feel::cout << "Tmax=" << Heat_Tmax << ", ";
  
  auto Tmean = mean( _range=elements(support(Th)), _expr=idv(Heat_T))(0,0);
  double measure = integrate( elements(support(Th)), cst(1.0) ).evaluate()(0,0);
  Feel::cout << "Tmean=" << Tmean << ", ";
  Feel::cout << "measure=" << measure << ", ";

  double Tstd_dev = normL2( elements(support(Th)), (idv(Heat_T)-cst(Tmean)) );
  Tstd_dev = math::sqrt(Tstd_dev / measure);
  Feel::cout << "Tstd_dev=" << Tstd_dev << ", ";

  Feel::cout << std::endl;
  
  tic();
  auto e = exporter( _mesh=mesh );

  e->step(t)->add("A", A);
  e->step(t)->add("V", V);
  e->step(t)->add("T", Heat_T);
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
  Heat_Told = (*Heat_T);

  // if BiotSavart Bc shall loop until ||A-Aold||<eps
  bool nonlinear = false;
  bool Tdepend = false;

  for( auto const& pairMat : M_materials )
    {
      auto name = pairMat.first;
      auto material = pairMat.second;


      auto sigma = material.getScalar("sigma");
      if ( sigma.expression().hasSymbol( "heat_T" ) )
	{
	  nonlinear = true;
	  Tdepend = true;
	  break;
	}
    }

  for( auto const& pairMat : M_th_materials )
    {
      auto name = pairMat.first;
      auto material = pairMat.second;

      auto k = material.getScalar("k");
      if ( k.expression().hasSymbol( "heat_T" ) )
	{
	  nonlinear = true;
	  Tdepend = true;
	  break;
	}
      auto rho = material.getScalar("rho");
      if ( rho.expression().hasSymbol( "heat_T" ) )
	{
	  nonlinear = true;
	  Tdepend = true;
	  break;
	}
      auto Cp = material.getScalar("Cp");
      if ( Cp.expression().hasSymbol( "heat_T" ) )
	{
	  nonlinear = true;
	  Tdepend = true;
	  break;
	}
    }

  // Feel::cout << "Compute Current density" << std::endl;
  auto J_cond = Jh->element();
  auto J_induct = Jh->element();
  for( auto const& pairMat : M_materials )
    {
      auto name = pairMat.first;
      auto material = pairMat.second;

      auto sigma = material.getScalar("sigma", "heat_T", idv(Heat_T) );
      // if ( sigma.expression().hasSymbol( "heat_T" ) )
      // 	{
      // 	  auto sExpr = material.getScalar("sigma", "heat_T", idv(Heat_T) );
      // 	  M_sigma.on(_range=markedelements(cond_mesh, material.meshMarkers()),_expr=sExpr );
      // 	}
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

  auto M_sigma = Th->element();
  auto M_k = Th->element();
  auto M_rho = Th->element();
  auto M_Cp = Th->element();

  if ( Tdepend )
    {
      auto M_sigma = Th->element();
      for( auto const& pairMat : M_materials )
	{
	  auto name = pairMat.first;
	  auto material = pairMat.second;

	  auto sigma = material.getScalar("sigma", "heat_T", idv(Heat_T) );
	  M_sigma.on(_range=markedelements(cond_mesh, material.meshMarkers()),_expr=sigma );
	}
      e->step(t)->add( "sigma", M_sigma );
      
      for( auto const& pairMat : M_th_materials )
	{
	  auto name = pairMat.first;
	  auto material = pairMat.second;

	  auto k = material.getScalar("k", "heat_T", idv(Heat_T) );
	  M_k.on(_range=markedelements(th_mesh, material.meshMarkers()),_expr=k );

	  auto rho = material.getScalar("rho", "heat_T", idv(Heat_T) );
	  M_rho.on(_range=markedelements(th_mesh, material.meshMarkers()),_expr=rho );
	  
	  auto Cp = material.getScalar("Cp", "heat_T", idv(Heat_T) );
	  M_Cp.on(_range=markedelements(th_mesh, material.meshMarkers()),_expr=Cp );
	}
      e->step(t)->add( "k", M_k );
      e->step(t)->add( "rho", M_rho );
      e->step(t)->add( "Cp", M_Cp );

    }
  e->save();
  toc("export init solution", (M_verbose > 0));
  
  auto mu0 = 4.e-7 * M_PI ; // SI Unit : H/m = m.kg/s2/A2

  Table mqs;
  std::string Vname[12];
  std::string Vfirst[12];
  std::string Bname;
  int ii = 0;
  int firstStep = 0;

  int iterNL = 0;
  int maxiterNL = 10;
  double initResidual;

  // define sequence of forced time steps
  double epstime = 1.e-3;
  bool reached = false;
  std::vector<double> forced_times = vdoption("forced-sequence");;
  forced_times.push_back(tmax);
  Feel::cout << "Forced sequence:" << forced_times << std::endl;

  int n_forced = 0;
  double forced_t = forced_times[n_forced];
  
  for(t = dt; t <= tmax+1e-10; )
    {

      tic();
      do {
	auto Anl = (*A);
	auto Vnl = (*V);
	auto Tnl = (*Heat_T);

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
			      _expr = dt * 1/mur * trace(trans(gradt(A))*grad(A)) );
	    //Feel::cout << "create lhs(0,0):" << material.meshMarkers() << std::endl;

	  }

	auto M01 = form2( _trial=Vh, _test=Ah ,_matrix=M, _rowstart=0, _colstart=1 );
	auto F0 = form1( _test=Ah, _vector=F, _rowstart=0 );

	auto M11 = form2( _trial=Vh, _test=Vh ,_matrix=M, _rowstart=1, _colstart=1 );
	auto M10 = form2( _trial=Ah, _test=Vh ,_matrix=M, _rowstart=1, _colstart=0 );
	auto F1 = form1( _test=Vh ,_vector=F, _rowstart=1 );

	auto M22 = form2( _trial=Th, _test=Th ,_matrix=M, _rowstart=2, _colstart=2 );
	auto F2 = form1( _test=Th ,_vector=F, _rowstart=2 );
      
	for( auto const& pairMat : M_materials )
	  {
	    auto name = pairMat.first;
	    auto material = pairMat.second;

	    auto sigma = material.getScalar("sigma", "heat_T", idv(Tnl) );

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

	    // heat equation (only contribution for Joule losses)
	    F2 += integrate( _range=markedelements(th_mesh, material.meshMarkers()),
			     _expr = mu0 * dt * 1/sigma * inner((idv(J_cond)+idv(J_induct)), (idv(J_cond)+idv(J_induct))) * id(Heat_T) );
	  
	  }

	// heat equation
	for( auto const& pairMat : M_th_materials )
	  {
	    auto name = pairMat.first;
	    auto material = pairMat.second;

	    auto k = material.getScalar("k", "heat_T", idv(Tnl));
	    auto rho = material.getScalar("rho", "heat_T", idv(Tnl));
	    auto Cp = material.getScalar("Cp", "heat_T", idv(Tnl));

	    // heat equation
	    M22 += integrate( _range=markedelements(th_mesh, material.meshMarkers()),
			      _expr = mu0 * rho * Cp * (id(Heat_T) * idt(Heat_T) ));
	    M22  += integrate( _range=markedelements(th_mesh, material.meshMarkers()),
			       _expr = mu0 * dt * k * inner(gradt(Heat_T), grad(Heat_T)) );
	    F2 += integrate( _range=markedelements(th_mesh, material.meshMarkers()),
			     _expr = mu0 *rho * Cp * idv(Heat_Told) * id(Heat_T) );
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
	itField = M_modelProps->boundaryConditions().find( "temperature");
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
		    //Feel::cout << "T Dirichlet[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
		    M22 += on(_range=markedfaces(mesh,marker), _rhs=F, _element=*Heat_T, _expr= g);
		  }
	      }
	    itType = mapField.find( "Neumann" );
	    if ( itType != mapField.end() )
	      {
		for ( auto const& exAtMarker : (*itType).second )
		  {
		    std::string marker = exAtMarker.marker();
		    auto g = expr(exAtMarker.expression());
		    g.setParameterValues({{"t", t}});
		    //Feel::cout << "T Neumann[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
		    F2 += integrate( markedfaces(th_mesh, marker), mu0*dt*g*id(Heat_T) );
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

		    h.setParameterValues({{"t", t}});
		    Tw.setParameterValues({{"t", t}});

		    //Feel::cout << "T Robin[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
		    M22 += integrate( markedfaces(th_mesh, marker), mu0*dt*h*idt(Heat_T)*id(Heat_T) );
		    F2 += integrate( markedfaces(th_mesh, marker), mu0*dt*h*Tw*id(Heat_T) );
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

	    auto sigma = material.getScalar("sigma", "heat_T", idv(Heat_T));
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

	    errorNL_th = normL2(_range = elements(th_mesh), _expr = (idv(Heat_T)-idv(Tnl)) );
	    normT = normL2(_range = elements(th_mesh), _expr = idv(Heat_T) );
	    Feel::cout << "errorNL_th=" << errorNL_th << " ,";
	    Feel::cout << "nomrT=" << normT << ", ";
	    Feel::cout << "epsNL_th*normT=" << epsNL_th*normT << ", ";
	    Feel::cout << std::endl;

	    iterNL++;
	  }
	
      } while ( (nonlinear==true) && (errorNL > epsNL*normA) &&  (errorNL_th > epsNL_th*normT) && (iterNL < maxiterNL) );
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
	  auto Heat_Tpost = (*Heat_T);
	  filter( A, Aold, Apost );
	  filter( V, Vold, Vpost );
	  filter( Heat_T, Heat_Told, Heat_Tpost );
	  auto estA = normL2( _range=elements(mesh), _expr=idv(A)-idv(Apost));
	  auto estV = normL2( _range=elements(cond_mesh), _expr=idv(V)-idv(Vpost));
	  auto estT = normL2( _range=elements(th_mesh), _expr=idv(Heat_T)-idv(Heat_Tpost));
	  auto est = std::max( estA, estV );
	  est = std::max( est, estT );
	  Feel::cout << "est : " << std::scientific << std::setprecision(3) << est << " ";
	  Feel::cout << "estA : " << std::scientific << std::setprecision(3) << estA << " ";
	  Feel::cout << "estV : "  << std::scientific << std::setprecision(3) << estV << " ";
	  Feel::cout << "estT : "  << std::scientific << std::setprecision(3) << estT << " ";
	  Feel::cout << "(dttol=" << std::scientific << std::setprecision(3) << dttol << "); ";// << std::endl;

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
	      Heat_Told = (*Heat_T);

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
	  bool allmost = ( fabs(1-forced_t/t) <= epstime );
	  if ( reached || allmost )
	    {
	      Feel::cout << "reached[" << forced_t << "]: ";
	      Feel::cout << "allmost=" << allmost << ", ";
	      Feel::cout << "reached=" << reached << ", ";
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
	  Heat_Told = (*Heat_T);
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

	  Tm = minmax( _range=elements(support(Th)), _pset=_Q<2>(), _expr=idv(Heat_T));
	  Heat_Tmax = Tm.max();
	  Heat_Tmin = Tm.min();
	  Feel::cout << "Tmin=" << Heat_Tmin << ", ";
	  Feel::cout << "Tmax=" << Heat_Tmax << ", ";
  
	  Tmean = mean( _range=elements(support(Th)), _expr=idv(Heat_T))(0,0);
	  Feel::cout << "Tmean=" << Tmean << ", ";

	  double Tstd_dev = normL2( elements(support(Th)), (idv(Heat_T)-cst(Tmean)) );
	  Tstd_dev = math::sqrt(Tstd_dev / measure);

	  e->step(t)->add( "T", Heat_T);
      
	  e->step(t)->add( "B", M_B );
	  // M_gradV = vf::project(_space=Jh, _range=elements(cond_mesh), _expr=trans(gradv(V))); // breaks in // why?
	  // e->step(t)->add( "E", M_gradV );

	  e->step(t)->add( "Jcond", J_cond );
	  e->step(t)->add( "Jinduct", J_induct );
	  e->step(t)->add( "J", idv(J_cond)+idv(J_induct) );

	  if ( Tdepend )
	    {
	      for( auto const& pairMat : M_materials )
		{
		  auto name = pairMat.first;
		  auto material = pairMat.second;

		  auto sigma = material.getScalar("sigma", "heat_T", idv(Heat_T) );
		  M_sigma.on(_range=markedelements(cond_mesh, material.meshMarkers()),_expr=sigma );

		  auto k = material.getScalar("k", "heat_T", idv(Heat_T) );
		  M_k.on(_range=markedelements(th_mesh, material.meshMarkers()),_expr=k );

		  auto rho = material.getScalar("rho", "heat_T", idv(Heat_T) );
		  M_rho.on(_range=markedelements(th_mesh, material.meshMarkers()),_expr=rho );
	  
		  auto Cp = material.getScalar("Cp", "heat_T", idv(Heat_T) );
		  M_Cp.on(_range=markedelements(th_mesh, material.meshMarkers()),_expr=Cp );
		}
	      e->step(t)->add( "sigma", M_sigma );
      
	      for( auto const& pairMat : M_th_materials )
		{
		  auto name = pairMat.first;
		  auto material = pairMat.second;

		  auto k = material.getScalar("k", "heat_T", idv(Heat_T) );
		  M_k.on(_range=markedelements(th_mesh, material.meshMarkers()),_expr=k );

		  auto rho = material.getScalar("rho", "heat_T", idv(Heat_T) );
		  M_rho.on(_range=markedelements(th_mesh, material.meshMarkers()),_expr=rho );
	  
		  auto Cp = material.getScalar("Cp", "heat_T", idv(Heat_T) );
		  M_Cp.on(_range=markedelements(th_mesh, material.meshMarkers()),_expr=Cp );
		}
	      e->step(t)->add( "k", M_k );
	      e->step(t)->add( "rho", M_rho );
	      e->step(t)->add( "Cp", M_Cp );

	    }

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
		  mqs.add_row({"t","NbIter","Residual",Vfirst[0],Vfirst[1],Vfirst[2],Vfirst[3],Bname, "Tmin", "Tmean", "Tmax", "Tstd_dev"});
		}
	      else
		{
		  mqs.add_row({"t","NbIter","Residual",Vfirst[0],Vfirst[1],Vfirst[2],Vfirst[3]
			       ,Vfirst[4],Vfirst[5],Vfirst[6],Vfirst[7],Bname, "Tmin", "Tmean", "Tmax", "Tstd_dev"});
		}
	      firstStep = 1;
	    }


	  //Bname = "{"+std::to_string(Bx) + "," + std::to_string(By) + "," + std::to_string(Bz) + "}";
	  Bname = std::to_string(Bz);
	  if (ii == 4)
	    {
	      mqs.add_row({std::to_string(t),std::to_string(nIterations),
			   std::to_string(Residual),
			   Vname[0],Vname[1],Vname[2],Vname[3],Bname, std::to_string(Heat_Tmin), std::to_string(Tmean), std::to_string(Heat_Tmax), std::to_string(Tstd_dev)});
	    }
	  else
	    {
	      mqs.add_row({std::to_string(t),std::to_string(nIterations),
			   std::to_string(Residual),
			   Vname[0],Vname[1],Vname[2],Vname[3],Vname[4],Vname[5],Vname[6],Vname[7],Bname, std::to_string(Heat_Tmin), std::to_string(Tmean), std::to_string(Heat_Tmax), std::to_string(Tstd_dev)});
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
	  if ( t >= forced_t )
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
