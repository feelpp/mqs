#include <tabulate/table.hpp>
#include <tabulate/markdown_exporter.hpp>
using namespace tabulate;

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/checker.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>
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
    ( "relax", po::value<double>()->default_value( 0 ), "relaxation parameter for picard" )
    ( "epsNL", po::value<double>()->default_value( 1.e-5 ), "eps for A picard" )
    ( "epsNL-V", po::value<double>()->default_value( 1.e-5 ), "eps for V picard" )
    ( "epsNL-Th", po::value<double>()->default_value( 1.e-5 ), "eps for Heat equation picard" )
    ( "itermaxNL", po::value<int>()->default_value( 10 ), "max iteration number for picard" )
    ( "epstime", po::value<double>()->default_value( 1.e-10 ), "eps for force time step detection" )
    ( "adaptive", po::value<bool>()->default_value( false ), "activate dt apdative scheme" )
    ( "dttol", po::value<double>()->default_value( 0. ), "dt tolerance" )
    ( "dt_min", po::value<double>()->default_value( 0. ), "dt min" )
    ( "dt_max", po::value<double>()->default_value( 0. ), "dt max" )
    ( "forced-sequence", po::value< std::vector<double> >()->default_value(std::vector<double>()), "list of forced times" )
    ( "verbosity", po::value<int>()->default_value( 0 ), "set verbosisity level" )
    ( "weakdir", po::value<bool>()->default_value( "false" ), "use Dirichlet weak formulation" )
    ( "penalty-coeff", po::value<double>()->default_value( 1.e+3 ), "penalty coefficient for weak Dirichlet" )
    ( "A0", po::value<std::string>()->default_value( "{0,0,0}" ), "initial A" )
    ( "V0", po::value<std::string>()->default_value( "0" ), "initial V" )
    ( "T0", po::value<std::string>()->default_value( "293.15" ), "initial T" )
    ( "Aexact", po::value<std::string>()->default_value( "" ), "exact A" )
    ( "Vexact", po::value<std::string>()->default_value( "" ), "exact V" );

#ifdef STRONGCOUPLING
  Environment env( _argc=argc, _argv=argv,_desc=options.add(Feel::backend_options("mqs")).add(Feel::biotsavart_options()),
		   _about=about(_name="mqs-th",
				_author="Feel++ Consortium",
				_email="feelpp-devel@feelpp.org"));
#else
  Environment env( _argc=argc, _argv=argv,_desc=options.add(Feel::backend_options("mqs")).add(Feel::backend_options("heat")).add(Feel::biotsavart_options()),
		   _about=about(_name="mqs-th",
				_author="Feel++ Consortium",
				_email="feelpp-devel@feelpp.org"));
#endif
  
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

  double dt_min = doption("dt_min");
  if ( boption("adaptive") && dt_min == 0)
    dt_min = std::max(dt/100., dttol);
  Feel::cout << "time-dtmin=" << dt_min << std::endl;
  
  double dt_max = std::min(dt*1000., tmax/10.);
  if ( boption("adaptive") && dt_max == 0)
    dt_max = std::min(dt*1000., tmax/10.);
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
  std::set<std::string> conductors(std::begin(range_conductors), std::end(range_conductors));
  Feel::cout << "Electric Materials markers (set): " << conductors << std::endl;

  // Materials in heat ONLY
  auto M_th_materials = M_modelProps->materials().materialWithPhysic(std::vector<std::string>({"heat"}));
  std::vector<std::string> range_th;
  for( auto const& mp : M_th_materials )
    for (auto const& marker : mp.second.meshMarkers() )
      range_th.push_back(marker);
  std::set<std::string> thdomains(std::begin(range_th), std::end(range_th));
  Feel::cout << "Thermic Materials markers (set): " << thdomains << std::endl;
  
  // Define SpaceFunctions
  tic();
  auto Ah = Pchv<1>( mesh );
  auto Vh = Pch<1>( mesh, markedelements(mesh, conductors) );
  auto Th = Pch<1>( mesh, markedelements(mesh, thdomains) );

  auto Jh = Pchv<1>( mesh, markedelements(mesh, conductors) );
  auto Eh = Pdhv<0>( mesh, markedelements(mesh, conductors) );
  auto Bh = Pdhv<0>( mesh );
  auto Qh = Pdh<0>( mesh, markedelements(mesh, conductors) );
  toc("define space functions", (M_verbose > 0));

  if (Environment::worldComm().isMasterRank())
    {
      std::cout << "mesh->numGlobalElements() "<< mesh->numGlobalElements() << std::endl;
      std::cout << "Ah->nDof() "<< Ah->nDof() << std::endl;
      std::cout << "Vh->nDof() "<< Vh->nDof() << std::endl;
      std::cout << "Th->nDof() "<< Th->nDof() << std::endl;
    }

  tic();
  auto A = Ah->elementPtr(); //Ah->element(A0); // how to init A to A0?;
  auto V = Vh->elementPtr(); //Vh->element(V0);
#ifdef STRONGCOUPLING
  auto Heat_T = Th->elementPtr(); //Vh->element(V0);
#else
  auto Heat_T = Th->element(); //Vh->element(V0);
#endif
  toc("define fields pointers", (M_verbose > 0));

  // init solutions
  tic();
  auto A0 = expr<3, 1>(soption(_name="A0"));
  auto V0 = expr(soption(_name="V0"));
  auto T0 = expr(soption(_name="T0"));
  A->on(_range=elements(mesh), _expr=A0); //(*A) = project(_space = Ah, _expr = A0);
  V->on(_range=elements(support(Vh)), _expr=V0); //(*V) = project(_space = Vh, _expr = V0);
#ifdef STRONGCOUPLING
  Heat_T->on(_range=elements(support(Th)), _expr=T0); //(*Heat_T) = project(_space = Th, _expr = T0);
#else
  Heat_T.on(_range=elements(support(Th)), _expr=T0); //(*Heat_T) = project(_space = Th, _expr = T0);
#endif

  auto Aold = (*A);
  auto Vold = (*V);
#ifdef STRONGCOUPLING
  auto Heat_Told = (*Heat_T);
#else
  auto Heat_Told = Heat_T;
#endif
  toc("init solutions", (M_verbose > 0));
    
  // Vincent way
  tic();
#ifdef STRONGCOUPLING
  BlocksBaseGraphCSR myblockGraph(3,3);
#else
  BlocksBaseGraphCSR myblockGraph(2,2);
#endif
  myblockGraph(0,0) = stencil(_test=Ah,_trial=Ah, _diag_is_nonzero=false, _close=false)->graph();
  myblockGraph(0,1) = stencil(_test=Ah,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
  myblockGraph(1,0) = stencil(_test=Vh,_trial=Ah, _diag_is_nonzero=false, _close=false)->graph();
  myblockGraph(1,1) = stencil(_test=Vh,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
#ifdef STRONGCOUPLING
  myblockGraph(0,2) = stencil(_test=Ah,_trial=Th, _diag_is_nonzero=false, _close=false)->graph();
  myblockGraph(1,2) = stencil(_test=Vh,_trial=Th, _diag_is_nonzero=false, _close=false)->graph();
  myblockGraph(2,0) = stencil(_test=Th,_trial=Ah, _diag_is_nonzero=false, _close=false)->graph();
  myblockGraph(2,1) = stencil(_test=Th,_trial=Vh, _diag_is_nonzero=false, _close=false)->graph();
  myblockGraph(2,2) = stencil(_test=Th,_trial=Th, _diag_is_nonzero=false, _close=false)->graph();
#else
  auto Mth = backend()->newMatrix(_test=Th,_trial=Th);
#endif
  auto M = backend()->newBlockMatrix(_block=myblockGraph);

#ifdef STRONGCOUPLING
  BlocksBaseVector<double> myblockVec(3);
#else
  BlocksBaseVector<double> myblockVec(2);
#endif
  myblockVec(0,0) = backend()->newVector( Ah );
  myblockVec(1,0) = backend()->newVector( Vh );
#ifdef STRONGCOUPLING
  myblockVec(2,0) = backend()->newVector( Th );
#else
  auto Fth = backend()->newVector( Th );
#endif
  auto F = backend()->newBlockVector(_block=myblockVec, _copy_values=false);

#ifdef STRONGCOUPLING
  BlocksBaseVector<double> myblockVecSol(3);
#else
  BlocksBaseVector<double> myblockVecSol(2);
#endif
  myblockVecSol(0,0) = A;
  myblockVecSol(1,0) = V;
#ifdef STRONGCOUPLING
  myblockVecSol(2,0) = Heat_T;
#endif
  auto U = backend()->newBlockVector(_block=myblockVecSol, _copy_values=false);
  toc("create Algebric blockforms", (M_verbose > 0));

  auto mqsbackend = backend(_name="mqs");
#ifndef STRONGCOUPLING
  auto heatbackend = backend(_name="heat");
#endif
  
  double t = 0;
  double Residual, Residual_Th;
  int nIterations, nIterations_Th ;
  
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
      (*A) = Aexact;
      
      Vexact_g = expr(Vexact_s);
      Vexact_g.setParameterValues({{"t", t}});
      Vexact = project(_space = Vh, _expr = Vexact_g);
      (*V) = Vexact;

      L2Aexact = normL2(_range = elements(mesh), _expr = Aexact_g);
      H1Aerror = 0;
      L2Aerror = 0;
      L2Vexact = normL2(_range = elements(support(Vh)), _expr = Vexact_g);
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

  Feel::cout << "t=" << t << ", ";
  Feel::cout << "B(" << pt[0] << "," << pt[1] << "," << pt[2] << ") = {" << Bx << "," << By << "," << Bz << "}, ";

  auto Tm = minmax( _range=elements(support(Th)), _pset=_Q<2>(), _expr=idv(Heat_T));
  double Heat_Tmax = Tm.max();
  double Heat_Tmin = Tm.min();
  Feel::cout << "Tmin=" << Heat_Tmin << ", ";
  Feel::cout << "Tmax=" << Heat_Tmax << ", ";
  
  auto Tmean = mean( _range=elements(support(Th)), _expr=idv(Heat_T))(0,0);
  double thmeasure = integrate( elements(support(Th)), cst(1.0) ).evaluate()(0,0);
  double vhmeasure = integrate( elements(support(Vh)), cst(1.0) ).evaluate()(0,0);
  double ahmeasure = integrate( elements(support(Ah)), cst(1.0) ).evaluate()(0,0);
  Feel::cout << "Tmean=" << Tmean << ", ";
  Feel::cout << "thmeasure=" << thmeasure << ", ";

  double Tstd_dev = normL2( elements(support(Th)), (idv(Heat_T)-cst(Tmean)) );
  Tstd_dev = math::sqrt(Tstd_dev / thmeasure);
  Feel::cout << "Tstd_dev=" << Tstd_dev << ", ";

  Feel::cout << std::endl;
  toc("some stats", (M_verbose > 0));
  
  tic();
  auto e = exporter( _mesh=mesh );

  e->step(t)->add("A", A);
  e->step(t)->add("V", V);
  e->step(t)->add("T", Heat_T);
  e->step(t)->add("B", M_B);
  
  if ( Uexact )
    {
      e->step(t)->add("Aexact", Aexact);
      e->step(t)->add("Vexact", Vexact);
    }

  Aold = (*A);
  Vold = (*V);
#ifdef STRONGCOUPLING
  Heat_Told = (*Heat_T);
#else
  Heat_Told = Heat_T;
#endif
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
	  Tdepend = true;
	  break;
	}
    }

  for( auto const& pairMat : M_th_materials )
    {
      auto name = pairMat.first;
      auto material = pairMat.second;

      for (auto const& prop : {"k", "rho", "Cp"})
	{
	  auto exp_prop = material.getScalar(prop);
	  if ( exp_prop.expression().hasSymbol( "heat_T" ) )
	    {
	      Tdepend = true;
	      break;
	    }
	}
    }
  //Feel::cout << "nonlinear=" << nonlinear << ", Tdepend=" << Tdepend << std::endl;
  if ( !nonlinear ) nonlinear = Tdepend;
  
  double relax = doption("relax");
  double epsNL = doption("epsNL");
  double epsNL_V = doption("epsNL-V");
  double epsNL_th = doption("epsNL-Th");
  int maxiterNL = ioption("itermaxNL");
  if ( nonlinear )
    {
      Feel::cout << "*** NonLinear problem detected ***" << std::endl;
      Feel::cout << "maxiterNL=" << maxiterNL << std::endl;
      Feel::cout << "epsNL=" << epsNL << std::endl;
      Feel::cout << "epsNL-V=" << epsNL_V << std::endl;
      Feel::cout << "epsNL-Th=" << epsNL_th << std::endl;
      Feel::cout << "relax=" << relax << std::endl;
      Feel::cout << "**********************************" << std::endl;
    }
  
  double rtol = doption("mqs.ksp-rtol");
  double atol = doption("mqs.ksp-atol");
  double stol = doption("mqs.snes-stol");
  Feel::cout << "rtol=" << rtol << std::endl;
  Feel::cout << "atol=" << atol << std::endl;
  Feel::cout << "stol=" << stol << std::endl;
  
  double errorNL, normA;
  double errorNL_V, normV;
  double errorNL_th, normT;

  // Feel::cout << "Compute Current density" << std::endl;
  auto Efield = Eh->element();
  auto J_cond = Eh->element();
  auto dAdt = Jh->element();
  auto J_induct = Jh->element();
  //auto J = Eh->element();
  auto Qth = Qh->element();
  Efield.zero();
  J_cond.zero();
  dAdt.zero();
  J_induct.zero();
  //J.zero();
  Qth.zero();
  for( auto const& pairMat : M_materials )
    {
      auto name = pairMat.first;
      auto material = pairMat.second;

      Efield.on(_range=markedelements(mesh, material.meshMarkers()), _expr= -trans(gradv(V)) );
      dAdt.on(_range=markedelements(mesh, material.meshMarkers()), _expr= -(idv(A)-idv(Aold))/dt );

      auto sigma = material.getScalar("sigma", "heat_T", idv(Heat_T) );
      J_cond.on(_range=markedelements(mesh, material.meshMarkers()), _expr=-sigma * trans(gradv(V)) );
      J_induct.on(_range=markedelements(mesh, material.meshMarkers()), _expr=-sigma * (idv(A)-idv(Aold))/dt );
	    
      Qth.on( _range=markedelements(mesh, material.meshMarkers()), _expr= inner(idv(J_cond)+idv(J_induct),idv(Efield)+idv(dAdt)) );
    }
  e->step(t)->add( "Jcond", J_cond );
  e->step(t)->add( "Jinduct", J_induct );
  //e->step(t)->add( "J", J );
  e->step(t)->add( "Qth", Qth );

  // Feel::cout << "Create Fields for Physical properties" << std::endl;
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
	  M_sigma.on(_range=markedelements(mesh, material.meshMarkers()),_expr=sigma );
	}
      e->step(t)->add( "sigma", M_sigma );

      for( auto const& pairMat : M_th_materials )
	{
	  auto name = pairMat.first;
	  auto material = pairMat.second;

	  auto k = material.getScalar("k", "heat_T", idv(Heat_T) );
	  M_k.on(_range=markedelements(mesh, material.meshMarkers()),_expr=k );

	  auto rho = material.getScalar("rho", "heat_T", idv(Heat_T) );
	  M_rho.on(_range=markedelements(mesh, material.meshMarkers()),_expr=rho );
	  
	  auto Cp = material.getScalar("Cp", "heat_T", idv(Heat_T) );
	  M_Cp.on(_range=markedelements(mesh, material.meshMarkers()),_expr=Cp );
	}
      e->step(t)->add( "k", M_k );
      e->step(t)->add( "rho", M_rho );
      e->step(t)->add( "Cp", M_Cp );
      
    }
  e->save();
  toc("export init solution", (M_verbose > 0));
  
  auto mu0 = 4.e-7 * M_PI ; // SI Unit : H/m = m.kg/s2/A2

  Table mqs;
  std::vector<std::variant<std::string, Table>> Vfirst; //[12];

#ifdef STRONGCOUPLING
  for (auto const& field : {"t","NbIter[MQS+Th]","Residual[MQS+Th]"})
    {
      Vfirst.push_back(field);
    }
#else
  for (auto const& field : {"t","NbIter[MQS]","Residual[MQS]","NbIter[Th]","Residual[Th]"})
    {
      Vfirst.push_back(field);
    }
#endif  
  if ( nonlinear ) Vfirst.push_back("NLIter");
  
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
	      Vfirst.push_back("V[" + marker + "]");
	      Vfirst.push_back("I[" + marker + "]");
	    }
	}
    }
  Vfirst.push_back("Bz("+std::to_string(pt[0]) + "," + std::to_string(pt[1]) + "," + std::to_string(pt[2]) + ")");

  for( auto const& mp : M_th_materials )
    for (auto const& marker : mp.second.meshMarkers() )
      for (auto const& field : {"Tmin", "Tmean", "Tmax", "Tstd_dev"})
	{
	  std::string sfield = field;
	  sfield += "[" + marker + "]";
	  Vfirst.push_back(sfield);
	}
  mqs.add_row(Vfirst);
  //Feel::cout << "MQS table: " << Vfirst.size() << std::endl;
  
  int iterNL = 0;

  // define sequence of forced time steps
  double epstime = doption("epstime");
  bool reached = false;
  bool allmost = false;
  std::vector<double> forced_times = vdoption("forced-sequence");;
  forced_times.push_back(tmax);
  std::sort(forced_times.begin(), forced_times.end());
  Feel::cout << "Forced sequence:" << forced_times << std::endl;
  Feel::cout << "epstime:" << epstime << std::endl;

  int n_forced = 0;
  double forced_t = forced_times[n_forced];

  bool converged = false;
  for(t = dt; t <= tmax; )
    {
      // force time step
      if ( !forced_times.empty() && !reached )
	{
	  forced_t = forced_times[n_forced];
	  allmost = ( fabs(1-forced_t/t) <= epstime );
	  if ( t >= forced_t )
	    {
	      reached = true;
	      if ( !allmost )
		{
		  dt -= t-forced_t;
		  t = forced_t;
		  Feel::cout << "forced_dt=" << dt <<  std::endl;
		}
	    }
	}
      Feel::cout << "** forced_time=" << forced_t << ", ";
      Feel::cout << "t=" << t << ", ";
      Feel::cout << "allmost=" << fabs(1-forced_t/t) << " (" << (fabs(1-forced_t/t) <= epstime) << ") ";
      Feel::cout << "n_forced=" << n_forced << ",";
      Feel::cout << "forced_times.size()=" << forced_times.size() << ",";
      Feel::cout << "reached=" << reached << std::endl;

      tic();
      do {

	auto Anl = (*A); // use deepCopy??
	auto Vnl = (*V); // use deepCopy??	
#ifdef STRONGCOUPLING
	auto Tnl = (*Heat_T); // use deepCopy??
#else
	auto Tnl = Heat_T;
#endif

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
			      _expr = 1/(mu0*mur) * trace(trans(gradt(A))*grad(A)) );
	    //Feel::cout << "create lhs(0,0):" << material.meshMarkers() << std::endl;

	  }

	auto M01 = form2( _trial=Vh, _test=Ah ,_matrix=M, _rowstart=0, _colstart=1 );
	auto F0 = form1( _test=Ah, _vector=F, _rowstart=0 );

	auto M11 = form2( _trial=Vh, _test=Vh ,_matrix=M, _rowstart=1, _colstart=1 );
	auto M10 = form2( _trial=Ah, _test=Vh ,_matrix=M, _rowstart=1, _colstart=0 );
	auto F1 = form1( _test=Vh ,_vector=F, _rowstart=1 );

#ifdef STRONGCOUPLING
	//Feel::cout << "strong coupling declare blocks M02 and M12" << std::endl;
	auto M02 = form2( _trial=Th, _test=Ah ,_matrix=M, _rowstart=0, _colstart=2 );	
	auto M12 = form2( _trial=Th, _test=Vh ,_matrix=M, _rowstart=1, _colstart=2 );

	//Feel::cout << "strong coupling declare blocks M20 and M21" << std::endl;
	auto M20 = form2( _trial=Ah, _test=Th ,_matrix=M, _rowstart=2, _colstart=0 );	
	auto M21 = form2( _trial=Vh, _test=Th ,_matrix=M, _rowstart=2, _colstart=1 );

	//Feel::cout << "strong coupling declare blocks M20 and M21" << std::endl;
	auto M22 = form2( _trial=Th, _test=Th ,_matrix=M, _rowstart=2, _colstart=2 );
	auto F2 = form1( _test=Th ,_vector=F, _rowstart=2 );
#else
	auto M22 = form2( _trial=Th, _test=Th, _matrix=Mth);
	auto F2 = form1( _test=Th, _vector=Fth);
#endif
	
	for( auto const& pairMat : M_materials )
	  {
	    auto name = pairMat.first;
	    auto material = pairMat.second;

	    auto sigma = material.getScalar("sigma", "heat_T", idv(Tnl) );

	    // Ampere law: sigma dA/dt + rot(1/(mu_r*mu_0) rotA) + sigma grad(V) = Js
	    M00 += integrate( _range=markedelements(mesh, material.meshMarkers()),
			      _expr = sigma * inner(id(A) , idt(A) )/dt);
	    //Feel::cout << "create lhs(0,0):" << material.meshMarkers() << std::endl;

	    M01  += integrate(_range=markedelements(mesh, material.meshMarkers()),
			      _expr = sigma * inner(id(A),trans(gradt(V))) );
	    //Feel::cout << "create lhs(0,1)" << std::endl;

	    F0 += integrate(_range=markedelements(mesh, material.meshMarkers()),
			    _expr = sigma * inner(id(A) , idv(Aold))/dt);
	    //Feel::cout << "create rhs(0)" << std::endl;

	    // auto Js = ;
	    // F0 += integrate(_range=markedelements(mesh, material.meshMarkers()),
	    // 		 _expr = dt * mu0 * inner(id(A) , Js));
	    // Feel::cout << "create rhs(0)" << std::endl;

	    // Current conservation: div( -sigma grad(V) -sigma*dA/dt) = Qs
	  
	    M11  += integrate( _range=markedelements(mesh, material.meshMarkers()),
			       _expr = sigma * inner(gradt(V), grad(V)) );
	    //Feel::cout << "create lhs(1,1)" << std::endl;

	    M10  += integrate( _range=markedelements(mesh, material.meshMarkers()),
			       _expr = sigma * inner(idt(A), trans(grad(V)))/dt );
	    //Feel::cout << "create lhs(1,0)" << std::endl;

	    F1 += integrate( _range=markedelements(mesh, material.meshMarkers()),
			     _expr = sigma * inner(idv(Aold), trans(grad(V)))/dt );
	    //Feel::cout << "create rhs(1)" << std::endl;

	    // auto Qs = ...;
	    // F1 += integrate(_range=markedelements(mesh, material.meshMarkers()),
	    // 		 _expr = dt * Qs * id(V);
	    // Feel::cout << "create row(1)" << std::endl;

	    // heat equation (only contribution for Joule losses)
#ifndef STRONGCOUPLING
	    //Feel::cout << "weak coupling F2 with Joules(n-1)" << std::endl;
	    F2 += integrate( _range=markedelements(mesh, material.meshMarkers()),
	    		     _expr = inner(idv(J_cond),idv(Efield)) * id(Heat_T) );
	    F2 += integrate( _range=markedelements(mesh, material.meshMarkers()),
	    		     _expr = inner(idv(J_induct),idv(dAdt)) * id(Heat_T) );
	    F2 += integrate( _range=markedelements(mesh, material.meshMarkers()),
	    		     _expr = inner(idv(J_cond),idv(dAdt)) * id(Heat_T) );
	    F2 += integrate( _range=markedelements(mesh, material.meshMarkers()),
	    		     _expr = inner(idv(J_induct),idv(Efield)) * id(Heat_T) );
#else
#ifdef FORM1
	    // Feel::cout << "strong coupling M20" << std::endl;
	    M20 += integrate( _range=markedelements(mesh, material.meshMarkers()),
	     		      _expr = sigma * inner(id(A), idv(Efield))/dt * idt(Heat_T) );
	    M20 += integrate( _range=markedelements(mesh, material.meshMarkers()),
	     		      _expr = sigma * inner(id(A), idv(dAdt))/dt * idt(Heat_T) );
	    // Feel::cout << "strong coupling F2" << std::endl;
	    F2 += integrate( _range=markedelements(mesh, material.meshMarkers()),
	     		     _expr = sigma * inner(idv(Aold), idv(Efield))/dt * id(Heat_T) );
	    F2 += integrate( _range=markedelements(mesh, material.meshMarkers()),
	     		     _expr = sigma * inner(idv(Aold), idv(dAdt))/dt * id(Heat_T) );
	    // Feel::cout << "strong coupling M21" << std::endl;
	    M21 += integrate( _range=markedelements(mesh, material.meshMarkers()),
	    		      _expr = sigma * inner(trans(grad(V)), idv(Efield)) * idt(Heat_T) );
	    M21 += integrate( _range=markedelements(mesh, material.meshMarkers()),
	    		      _expr = sigma * inner(trans(grad(V)), idv(dAdt)) * idt(Heat_T) );
#else
	    F2 += integrate( _range=markedelements(mesh, material.meshMarkers()),
	    		     _expr = inner(idv(J_cond),idv(Efield)) * id(Heat_T) );
	    F2 += integrate( _range=markedelements(mesh, material.meshMarkers()),
	    		     _expr = inner(idv(J_induct),idv(dAdt)) * id(Heat_T) );
	    F2 += integrate( _range=markedelements(mesh, material.meshMarkers()),
	    		     _expr = inner(idv(J_cond),idv(dAdt)) * id(Heat_T) );
	    F2 += integrate( _range=markedelements(mesh, material.meshMarkers()),
	    		     _expr = inner(idv(J_induct),idv(Efield)) * id(Heat_T) );
#endif
#endif
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
	    M22 += integrate( _range=markedelements(mesh, material.meshMarkers()),
			      _expr = rho * Cp * (id(Heat_T) * idt(Heat_T) )/dt);
	    M22  += integrate( _range=markedelements(mesh, material.meshMarkers()),
			       _expr = k * inner(gradt(Heat_T), grad(Heat_T)) );
	    F2 += integrate( _range=markedelements(mesh, material.meshMarkers()),
			     _expr = rho * Cp * idv(Heat_Told) * id(Heat_T) / dt );
	  }
	toc("assembling", (M_verbose > 0));
     
	tic();
	// Implement Dirichlet fort
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
		    LOG(INFO) << "T Dirichlet[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
#ifdef STRONGCOUPLING
		    M22 += on(_range=markedfaces(mesh, marker), _rhs=F, _element=*Heat_T, _expr= g);
#else
		    M22 += on(_range=markedfaces(mesh, marker), _rhs=Fth, _element=Heat_T, _expr= g);
#endif
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
		    LOG(INFO) << "T Neumann[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
		    F2 += integrate( markedfaces(mesh, marker), g*id(Heat_T) );
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

		    LOG(INFO) << "T Robin[" << marker << "] : ";
		    LOG(INFO) << "h=" << exAtMarker.expression1() << ", h=" << h << ", ";
		    LOG(INFO) << "Tw=" << exAtMarker.expression2() << ", Tw=" << Tw << ", "<< std::endl;
		    M22 += integrate( markedfaces(mesh, marker), h*idt(Heat_T)*id(Heat_T) );
		    F2 += integrate( markedfaces(mesh, marker), h*Tw*id(Heat_T) );
		  }
	      }
	  }
	itField = M_modelProps->boundaryConditions().find( "magnetic-potential");
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
		    LOG(INFO) << "A Dirichlet[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
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
		    LOG(INFO) << "A DirichletX[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
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
		    LOG(INFO) << "A DirichletY[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
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
		    LOG(INFO) << "A DirichletZ[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
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
		    auto jEx = idv(J_cond) + idv(J_induct);

		    As.compute(jEx, false, true, conductors);
		    
		    LOG(INFO) << "A BiotSavart[" << marker << "] : " << std::endl;
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
		    LOG(INFO) << "V[" << marker << "]=" << g.evaluate()(0,0) << std::endl;
		    M11 += on(_range=markedfaces(mesh,marker), _rhs=F, _element=*V, _expr= g);
		  }
	      }
	  }
	toc("boundary conditions", (M_verbose > 0));
    
	/* Solve */
	tic();
	auto result = mqsbackend->solve( _matrix=M, _rhs=F, _solution=U, _rebuild=true);
#ifdef STRONGCOUPLING
	std::string msg = (boost::format("[MQS+Heat %2%] t=%1% NbIter=%3% Residual=%4%") % t
			   % soption("mqs.pc-type")
			   % result.nIterations()
			   % result.residual()).str();
#else
	std::string msg = (boost::format("[MQS %2%] t=%1% NbIter=%3% Residual=%4%") % t
			   % soption("mqs.pc-type")
			   % result.nIterations()
			   % result.residual()).str();
#endif
	if (result.isConverged())
	  {
	    Feel::cout << tc::green << msg << tc::reset << " "; // << std::endl;
	  }
	else
	  {
	    std::string errmsg = msg + " Failed to converge";
	    throw std::logic_error( errmsg );
	  }
#ifndef STRONGCOUPLING
	Feel::cout << "Solve Th (weak coupling)" << std::endl;
	auto thresult = heatbackend->solve( _matrix=Mth, _rhs=Fth, _solution=Heat_T, _rebuild=true);
	std::string thmsg = (boost::format("[Heat %2%] t=%1% NbIter=%3% Residual=%4%") % t
			     % soption("heat.pc-type")
			     % thresult.nIterations()
			     % thresult.residual()).str();
	if (thresult.isConverged())
	  {
	    Feel::cout << tc::green << thmsg << tc::reset << " "; // << std::endl;
	  }
	else
	  {
	    std::string errmsg = thmsg + " Failed to converge";
	    throw std::logic_error( errmsg );
	  }
#else
	Feel::cout << "Solve MQS+Th (strong coupling)";
#ifdef FORM1
	Feel::cout << " -- form1";
#endif
	Feel::cout << std::endl;
#endif
	
	Residual =  result.residual();
	nIterations = result.nIterations();
#ifndef STRONGCOUPLING
	Residual_Th =  thresult.residual();
	nIterations_Th = thresult.nIterations();
#endif
	toc("solve", (M_verbose > 0));

	// update A and V pointers from U
	myblockVecSol.localize(U);

	// compute errorNL (see V. Chabannes comment for more precise handling)
	if ( nonlinear )
	  {
	    errorNL = normL2(_range = elements(mesh), _expr = (idv(A)-idv(Anl)) );
	    normA = normL2(_range = elements(mesh), _expr = idv(A) );
	    bool cvgA =(errorNL <= std::max(epsNL*normA, atol));
	    Feel::cout << "iterNL=" << iterNL << " ,";
	    Feel::cout << "errorNL=" << errorNL << (cvgA? "*":"") << " ,";

	    errorNL_V = normL2(_range = elements(support(Vh)), _expr = (idv(V)-idv(Vnl)) );
	    normV = normL2(_range = elements(support(Vh)), _expr = idv(V) );
	    bool cvgV= (errorNL_V <= std::max(epsNL_V*normV, atol));
	    Feel::cout << "errorNL_V=" << errorNL_V << (cvgV? "*":"") << " ,";

	    // 	    // relaxation
	    // 	    Trelax.zero();
	    // #ifdef STRONGCOUPLING
	    // 	    Trelax.add( 1-relax, (*Heat_T) ); // use deepCopy??
	    // #else
	    // 	    Trelax.add( 1-relax, Heat_T );
	    // #endif
	    // 	    Trelax.add( relax, Tnl );
    
	    errorNL_th = normL2(_range = elements(support(Th)), _expr = (idv(Heat_T)-idv(Tnl)) );
	    normT = normL2(_range = elements(support(Th)), _expr = idv(Heat_T) );
	    bool cvgT = (errorNL_th <= std::max(epsNL_th*normT, atol));
	    Feel::cout << "errorNL_th=" << errorNL_th << (cvgT? "*":"") << " ,";
	    Feel::cout << std::endl;

	    iterNL++;

	    converged = cvgA && cvgV && cvgT;
	  }
	
	/* reinit  */
	M->zero();
	F->zero();

#ifndef STRONGCOUPLING
	Mth->zero();
	Fth->zero();
#endif

      } while ( (nonlinear==true) && !converged	&& (iterNL < maxiterNL) );
      toc("non-linear step", ( (M_verbose > 0) && nonlinear) );

      if ( nonlinear && !converged )
	{
	  if ( iterNL >=  maxiterNL )
	    throw std::logic_error( "NL: max picard iteration reached" );
	  if ( errorNL > std::max(epsNL*normA, atol) )
	    throw std::logic_error( "NL: picard on A failed" );
	  if ( errorNL_V > std::max(epsNL_V*normV, atol) )
	    throw std::logic_error( "NL: picard on V failed" );
	  if ( errorNL_th > std::max(epsNL_th*normT, atol) )
	    throw std::logic_error( "NL: picard on T failed" );
	}
	
      bool do_export=true;
      
      /* Adapt time step */
      if ( boption( "adaptive") )
        {
	  tic();
	  //Feel::cout << "adaptive time stepping t=" << t << std::endl;
	  std::string adapt_msg;
	  // time filtering , get order 2
	  auto filter = [&dt, &dtprev]( auto const& in, auto const& inprev, auto& out ) { 
			  double nu = dt*(dt+dtprev)/(dtprev*(2*dt+dtprev));
			  double c1 = 2*dtprev/(dt+dtprev);
			  double c2 = 2*dt/(dt+dtprev);
			  // Feel::cout << "  adaptive time stepping nu=" << nu << " c1=" << c1 << " c2=" << c2 << "; "; //<< std::endl;
			  out.on( _range=elements(out.mesh()), _expr=idv(in)-(nu/2)*(c1*idv(inprev) - 2*idv(in) + c2*idv(inprev) )); 
			};
	  auto Apost = (*A);
	  auto Vpost = (*V);
#ifdef STRONGCOUPLING
	  auto Heat_Tpost = (*Heat_T);
#else
	  auto Heat_Tpost = Heat_T;
#endif
	  filter( A, Aold, Apost );
	  filter( V, Vold, Vpost );
	  filter( Heat_T, Heat_Told, Heat_Tpost );
	  auto estA = normL2( _range=elements(mesh), _expr=idv(A)-idv(Apost));
	  auto estV = normL2( _range=elements(support(Vh)), _expr=idv(V)-idv(Vpost));
	  auto estT = normL2( _range=elements(support(Th)), _expr=idv(Heat_T)-idv(Heat_Tpost));
	  auto est = std::max( estA, vhmeasure/ahmeasure*estV );
	  est = std::max( est, thmeasure/ahmeasure*estT );
	  Feel::cout << "est : " << std::scientific << std::setprecision(3) << est << " ";
	  Feel::cout << "estA : " << std::scientific << std::setprecision(3) << estA << " ";
	  Feel::cout << "estV : "  << std::scientific << std::setprecision(3) << vhmeasure/ahmeasure*estV << " ";
	  Feel::cout << "estT : "  << std::scientific << std::setprecision(3) << thmeasure/ahmeasure*estT << " ";
	  Feel::cout << "(dttol=" << std::scientific << std::setprecision(3) << dttol << "); ";// << std::endl;

	  //Feel::cout << "forced_time=" << forced_t << ", ";
	  //Feel::cout << "t=" << t << ", ";
	  //Feel::cout << "allmost=" << fabs(1-forced_t/t) << " (" << (fabs(1-forced_t/t) <= epstime) << ") " << "[" << allmost << "]";
	  allmost = (fabs(1-forced_t/t) <= epstime);
	  //Feel::cout << "nallmost=" << allmost << ",";
	  //Feel::cout << "reached=" << reached;
	  Feel::cout << std::endl;
	  if ( est > dttol )
	    {  
	      if ( dt == dt_min )
		{
		  // Update current densities
		  //Feel::cout << "dt_min: update Qth" << std::endl;
		  Efield.zero();
		  J_cond.zero();
		  dAdt.zero();
		  J_induct.zero();
		  //J.zero();
		  Qth.zero();
		  for( auto const& pairMat : M_materials )
		    {
		      auto name = pairMat.first;
		      auto material = pairMat.second;

		      Efield.on(_range=markedelements(mesh, material.meshMarkers()), _expr= -trans(gradv(V)) );
		      dAdt.on(_range=markedelements(mesh, material.meshMarkers()), _expr= -(idv(A)-idv(Aold))/dt );

		      auto sigma = material.getScalar("sigma", "heat_T", idv(Heat_T));
		      J_cond.on(_range=markedelements(mesh, material.meshMarkers()), _expr=-sigma * trans(gradv(V)) );
		      J_induct.on(_range=markedelements(mesh, material.meshMarkers()), _expr=-sigma * (idv(A)-idv(Aold))/dt );
	    
		      Qth.on( _range=markedelements(mesh, material.meshMarkers()), _expr= inner(idv(J_cond)+idv(J_induct),idv(Efield)+idv(dAdt)) );
		    }

		  dtprev=dt;

		  Aold = (*A);
		  Vold = (*V);
#ifdef STRONGCOUPLING
		  Heat_Told = (*Heat_T);
#else
		  Heat_Told = Heat_T;
#endif
		  // export
		  do_export=true; //false;

		  adapt_msg = "keeping the time step (dt_min)";
		}
	      else
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
		      //Feel::cout << "***";
		    }
		}
	    }
	  else //if ( est < dttol )
	    {
	      // Update current densities
	      //Feel::cout << "est<dttol: update Qth" << std::endl;
	      Efield.zero();
	      J_cond.zero();
	      J_induct.zero();
	      dAdt.zero();
	      Qth.zero();
	      for( auto const& pairMat : M_materials )
		{
		  auto name = pairMat.first;
		  auto material = pairMat.second;

		  Efield.on(_range=markedelements(mesh, material.meshMarkers()), _expr= -trans(gradv(V)) );
		  dAdt.on(_range=markedelements(mesh, material.meshMarkers()), _expr= -(idv(A)-idv(Aold))/dt );

		  auto sigma = material.getScalar("sigma", "heat_T", idv(Heat_T));
		  J_cond.on(_range=markedelements(mesh, material.meshMarkers()), _expr=-sigma * trans(gradv(V)) );
		  J_induct.on(_range=markedelements(mesh, material.meshMarkers()), _expr=-sigma * (idv(A)-idv(Aold))/dt );
	    
		  Qth.on( _range=markedelements(mesh, material.meshMarkers()), _expr= inner(idv(J_cond)+idv(J_induct),idv(Efield)+idv(dAdt)) );
	      
		}

	      dtprev=dt;

	      Aold = (*A);
	      Vold = (*V);
#ifdef STRONGCOUPLING
	      Heat_Told = (*Heat_T);
#else
	      Heat_Told = Heat_T;
#endif
	      // export
	      do_export=true; //false;

	      if ( est <= dttol/8. && dt != dt_max )
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
		  if ( dt == dt_max ) adapt_msg += " (dt_max)";
		}

	    }
	   
	  std::string msg = (boost::format("[adapt dt=%1%] ") % dt).str();
	  msg += adapt_msg;

	  // time accepted
	  //Feel::cout << " xx allmost (" << allmost << ") xx " << std::endl;
	  if ( reached && allmost )
	    if ( n_forced < forced_times.size()-1 )
	      {
		reached = false;
		n_forced++;
		forced_t = forced_times[n_forced] ;
		//dt = doption(_name = "ts.time-step");
		std::string strdt = (boost::format(" - go to next sequence t=%1%") % forced_t).str();
		msg += strdt;
	      }
	      
	  Feel:cout << msg << std::endl;
	  toc( msg, (M_verbose > 0));
        }
      else
	{

	  // Update current densities
	  //Feel::cout << "\nno adapt: update Qth" << std::endl;
	  Efield.zero();
	  J_cond.zero();
	  J_induct.zero();
	  dAdt.zero();
	  Qth.zero();
	  for( auto const& pairMat : M_materials )
	    {
	      auto name = pairMat.first;
	      auto material = pairMat.second;

	      Efield.on(_range=markedelements(mesh, material.meshMarkers()), _expr= -trans(gradv(V)) );
	      dAdt.on(_range=markedelements(mesh, material.meshMarkers()), _expr= -(idv(A)-idv(Aold))/dt );

	      auto sigma = material.getScalar("sigma", "heat_T", idv(Heat_T));
	      J_cond.on(_range=markedelements(mesh, material.meshMarkers()), _expr=-sigma * trans(gradv(V)) );
	      J_induct.on(_range=markedelements(mesh, material.meshMarkers()), _expr=-sigma * (idv(A)-idv(Aold))/dt );
	    
	      Qth.on( _range=markedelements(mesh, material.meshMarkers()), _expr= inner(idv(J_cond)+idv(J_induct),idv(Efield)+idv(dAdt)) );
	    }

	  //Feel::cout << " ** allmost (" << allmost << ") ** " << std::endl;
	  if ( allmost )
	    {
	      if ( n_forced < forced_times.size()-1 )
		{
		  n_forced++;
		  forced_t = forced_times[n_forced];
		  reached = false;
		  //dt = doption(_name = "ts.time-step");
		  std::string strdt = (boost::format(" - go to next sequence t=%1%") % forced_t).str();
		  Feel::cout << strdt << std::endl;
		}
	    }

	  Aold = (*A);
	  Vold = (*V);
#ifdef STRONGCOUPLING
	  Heat_Told = (*Heat_T);
#else
	  Heat_Told = Heat_T;
#endif
	}
      
      
      if ( do_export)
	{
	  // Display Magnetic Field
	  M_B = vf::project(_space=Bh, _range=elements(mesh), _expr=curlv(A));
	  val = M_B(pt);
	  Bx = val(0,0,0); // evaluation de Bx
	  By = val(1,0,0); // evaluation de By
	  Bz = val(2,0,0); // evaluation de Bz

	  tic();
	  e->step(t)->add( "A", A);
	  e->step(t)->add( "V", V);

	  double P = integrate(_range=elements(support(Qh)), _expr=idv(Qth)).evaluate()(0,0);
	  Feel::cout << "P[]=" << P << ", ";
	  
	  Tm = minmax( _range=elements(support(Th)), _pset=_Q<2>(), _expr=idv(Heat_T));
	  Heat_Tmax = Tm.max();
	  Heat_Tmin = Tm.min();
	  Feel::cout << "Tmin=" << Heat_Tmin << ", ";
	  Feel::cout << "Tmax=" << Heat_Tmax << ", ";
  
	  Tmean = mean( _range=elements(support(Th)), _expr=idv(Heat_T))(0,0);
	  Feel::cout << "Tmean=" << Tmean << ", ";

	  double Tstd_dev = normL2( elements(support(Th)), (idv(Heat_T)-cst(Tmean)) );
	  Tstd_dev = math::sqrt(Tstd_dev / thmeasure);

	  e->step(t)->add( "T", Heat_T);
      
	  e->step(t)->add( "B", M_B );

	  e->step(t)->add( "Jcond", J_cond );
	  e->step(t)->add( "Jinduct", J_induct );
	  //e->step(t)->add( "J", J );
	  e->step(t)->add( "Qth", Qth );

	  if ( Tdepend )
	    {
	      for( auto const& pairMat : M_materials )
		{
		  auto name = pairMat.first;
		  auto material = pairMat.second;

		  auto sigma = material.getScalar("sigma", "heat_T", idv(Heat_T) );
		  M_sigma.on(_range=markedelements(mesh, material.meshMarkers()),_expr=sigma );
		}
	      e->step(t)->add( "sigma", M_sigma );
      
	      for( auto const& pairMat : M_th_materials )
		{
		  auto name = pairMat.first;
		  auto material = pairMat.second;

		  auto k = material.getScalar("k", "heat_T", idv(Heat_T) );
		  M_k.on(_range=markedelements(mesh, material.meshMarkers()),_expr=k );

		  auto rho = material.getScalar("rho", "heat_T", idv(Heat_T) );
		  M_rho.on(_range=markedelements(mesh, material.meshMarkers()),_expr=rho );
	  
		  auto Cp = material.getScalar("Cp", "heat_T", idv(Heat_T) );
		  M_Cp.on(_range=markedelements(mesh, material.meshMarkers()),_expr=Cp );
		}
	      e->step(t)->add( "k", M_k );
	      e->step(t)->add( "rho", M_rho );
	      e->step(t)->add( "Cp", M_Cp );

	    }

	  // export data to table
	  std::vector<std::variant<std::string, Table>> exported_data;
	  exported_data.push_back(std::to_string(t));
	  exported_data.push_back(std::to_string(nIterations));
	  exported_data.push_back((boost::format("%1%") % Residual).str());
#ifndef STRONGCOUPLING
	  exported_data.push_back(std::to_string(nIterations_Th));
	  exported_data.push_back((boost::format("%1%") % Residual_Th).str());
#endif
	  if ( nonlinear) exported_data.push_back(std::to_string(iterNL));
	  
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
		      Feel::cout << "V[" << marker << "]=" << g.evaluate()(0,0) << ", ";
		      exported_data.push_back( std::to_string(g.evaluate()(0,0)) );
		      
		      double I = integrate( markedfaces( mesh, marker ), inner(idv(J_cond),N()) ).evaluate()(0,0);
		      I += integrate( markedfaces( mesh, marker ), inner(idv(J_induct),N()) ).evaluate()(0,0);
		      exported_data.push_back( std::to_string(I) );
		      Feel::cout << "I[" << marker << "]=" << I << ", ";
		    }
		}
	    }

	  Feel::cout << " B(" << pt[0] << "," << pt[1] << "," << pt[2] << ") = {" << Bx << "," << By << "," << Bz << "}";
	  Feel::cout << std::endl;

	  //Bname = "{"+std::to_string(Bx) + "," + std::to_string(By) + "," + std::to_string(Bz) + "}";
	  exported_data.push_back( std::to_string(Bz) );
	  
	  for( auto const& mp : M_th_materials )
	    for (auto const& marker : mp.second.meshMarkers() )
	      {
		Tm = minmax( _range=markedelements(support(Th), marker), _pset=_Q<2>(), _expr=idv(Heat_T));
		Heat_Tmax = Tm.max();
		Heat_Tmin = Tm.min();
  
		Tmean = mean( _range=markedelements(support(Th), marker), _expr=idv(Heat_T))(0,0);
		double dthmeasure = integrate( markedelements(support(Th),marker), cst(1.0) ).evaluate()(0,0);
		  
		Tstd_dev = normL2( elements(support(Th)), (idv(Heat_T)-cst(Tmean)) );
		Tstd_dev = math::sqrt(Tstd_dev / dthmeasure);
		  
		exported_data.push_back( std::to_string(Heat_Tmin) );
		exported_data.push_back( std::to_string(Tmean) );
		exported_data.push_back( std::to_string(Heat_Tmax) );
		exported_data.push_back( std::to_string(Tstd_dev) );
	      }
	  //Feel::cout << "MQS table data: " << exported_data.size() << std::endl;
	  mqs.add_row(exported_data);

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

	      L2Vexact = normL2(_range = elements(support(Vh)), _expr = Vexact_g);
	      L2Verror = normL2(elements(support(Vh)), (idv(V) - idv(Vexact)));
	      H1Verror = normH1(elements(support(Vh)), _expr = (idv(V) - idv(Vexact)), _grad_expr = (gradv(V) - gradv(Vexact)));
	      Feel::cout << " V: " << L2Verror << " " << L2Verror / L2Vexact << " " << H1Verror << std::endl;
	    }
	}
      
      // reset NL counter
      iterNL = 0;
    
      t += dt;

      // force final time if necessary
      Feel::cout << "next t=" << t << " (reached=" << reached << ")" << std::endl;
      if ( t>tmax && !reached)
	{
	  dt = t - tmax;
	  t = tmax;
	  Feel::cout <<  (boost::format(" ** force t to %1% and dt to %2%") % t % dt).str() << std::endl;
	  reached = true;
	}
    }

  // export as markdow table
  std::ofstream tablefile;
  tablefile.open("1Dmodel.adoc");
  if (!tablefile)
    throw std::logic_error( "1Dmodel.adoc: cannot create file" );

  MarkdownExporter exporter;
  auto markdown = exporter.dump(mqs);
  tablefile << markdown << std::endl;
  tablefile.close();
}
