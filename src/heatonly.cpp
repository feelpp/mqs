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
    ( "T0", po::value<std::string>()->default_value( "{0,0,0}" ), "initial T" )
    ( "Texact", po::value<std::string>()->default_value( "" ), "exact T" )

  Environment env( _argc=argc, _argv=argv,_desc=options.add(Feel::backend_options("mqs")),
		   _about=about(_name="heatmqs",
				_author="Feel++ Consortium",
				_email="feelpp-devel@feelpp.org"));

  int M_verbose = ioption(_name="verbosity");

    //Recuperer time frame
  double dt = doption(_name = "ts.time-step");
  Feel::cout << "time-step=" << dt << std::endl;

  double tmax = doption(_name = "ts.time-final");
  Feel::cout << "time-final=" << tmax << std::endl;

  // Eventually get a solution
  bool Texact = false;

  std::string Texact_s = soption(_name = "Texact");

  if ( !Aexact_s.empty() && !Vexact_s.empty() )
    {
      Uexact = true;
      Feel::cout << "* Texact=" << Texact_s << std::endl;
    }

  // Load Mesh
  auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);

#if 0 //json load, adapt it to heat
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
#endif

  // Define SpaceFunctions
  tic();
  auto Th = Pchv<1>( mesh );

#if 0 // ??????
auto cond_mesh = Vh->mesh();
  if (Environment::worldComm().isMasterRank())
  {
    std::cout << "mesh->numGlobalElements() "<< mesh->numGlobalElements() << std::endl;
    std::cout << "cond_mesh->numGlobalElements() "<< cond_mesh->numGlobalElements() << std::endl;
    std::cout << "Ah->nDof() "<<Ah->nDof() << std::endl;
    std::cout << "Vh->nDof() "<<Vh->nDof() << std::endl;
  }
#endif 

  auto T = Th->elementPtr();  

  toc("define space functions", (M_verbose > 0));

  // init solutions
  tic();
  auto T0 = expr<3, 1>(soption(_name="T0"));
  auto T0e = Th->element();
  (*T) = project(_space = Th, _expr = T0);

  auto Told = (*T);
  toc("init solutions", (M_verbose > 0));
    
  // Vincent way
  tic();
  BlocksBaseGraphCSR myblockGraph(2,2);
  myblockGraph(0,0) = stencil(_test=Th,_trial=Th, _diag_is_nonzero=false, _close=false)->graph();
  myblockGraph(0,1) = stencil(_test=Th,_trial=Th, _diag_is_nonzero=false, _close=false)->graph();
  auto M = backend()->newBlockMatrix(_block=myblockGraph);

  BlocksBaseVector<double> myblockVec(1);
  myblockVec(0,0) = backend()->newVector( Th );
  auto F = backend()->newBlockVector(_block=myblockVec, _copy_values=false);

  BlocksBaseVector<double> myblockVecSol(1);
  myblockVecSol(0,0) = T;
  auto U = backend()->newBlockVector(_block=myblockVecSol, _copy_values=false);
  toc("create Algebric blockforms", (M_verbose > 0));

  auto mybackend = backend(_name="mqs");

  //double t = 0;

  double L2Texact, H1Terror, L2Terror;

  auto Texact = Th->element();

  auto Texact_g = expr<3, 1>("{0,0,0}");

  auto mybdfT = bdf(_space = Th, _name="mybdfT");

  for (auto time : mybdfT -> priorTimes() )
  {
    if (Environment::worldComm().isMasterRank())
    {
      std::cout << "Initialize prior times (from timeInitial()) : " << time.second << "s index: " << time.first << "\n";
    }
    T0.setParameterValues({{"t",time.second}});
    T0e = project(_space=Th, _expr=T0);
    mybdfA->setUnknown(time.first,T0e);
  }

  if ( Texact )
  {
    tic();
    Texact_g = expr<3, 1>(Texact_s);
    Texact_g.setParameterValues({{"t", mybdfT->timeInitial()}});
    Texact = project(_space = Th, _expr = Texact_g);
    Feel::cout << "Define Texact" << std::endl;
    (*T) = Texact;

    L2Texact = normL2(_range = elements(mesh), _expr = Texact_g);
    H1Terror = 0;
    L2Terror = 0;
    toc("init exact solution", (M_verbose > 0));
  }

  tic();
  auto e = exporter( _mesh=mesh );

  e->step(0)->add("T", T);

  if ( Texact )
  {
    e->step(0)->add("Texact", Texact);
  }

  Told = (*T);

  e->save();
  toc("export init solution", (M_verbose > 0));
  
  auto mu0 = 4.e-7 * M_PI ; // SI Unit : H/m = m.kg/s2/A2

  mybdfT->start();

  for (double t = dt; mybdfT->isFinished() == false; )
  {
    auto bdfT_poly = mybdfT->polyDeriv();
    tic();
    auto M00 = form2( _trial=Th, _test=Th ,_matrix=M, _rowstart=0, _colstart=0 ); 
    for( auto const& pairMat : M_modelProps->materials() )
	  {
	    auto name = pairMat.first;
	    auto material = pairMat.second;

	    auto mur = material.getScalar("mu_mag");
  
      //heat
	    M00 += integrate( _range=markedelements(mesh, material.meshMarkers()),
			      _expr = trace(trans(gradt(T))*grad(T)) );

	  }

    auto M01 = form2( _trial=Th, _test=Th ,_matrix=M, _rowstart=0, _colstart=1 );
    auto F0 = form1( _test=Th, _vector=F, _rowstart=0 );
      
    for( auto const& pairMat : M_materials )
	  {
	    auto name = pairMat.first;
	    auto material = pairMat.second;

	    auto sigma = material.getScalar("sigma");

	    // heat
	    M00 += integrate( _range=markedelements(mesh, material.meshMarkers()),
			      _expr = mybdfT->polyDerivCoefficient(0) * inner(id(T) , idt() ));
T
      //heat
	    M01  += integrate(_range=markedelements(mesh, material.meshMarkers()),
			      _expr = 0 );

      //heat
	    F0 += integrate(_range=markedelements(mesh, material.meshMarkers()),
			    _expr = inner(id(T) , idv(bdfT_poly)));
	  }
    toc("assembling", (M_verbose > 0));

    tic();
    // Implement Dirichlet fort
    auto itField = M_modelProps->boundaryConditions().find( "temperature");
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
		      g.setParameterValues({{"t", mybdfT->time()}});
		      //Feel::cout << "A Dirichlet[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
		      M00 += on(_range=markedfaces(mesh,marker), _rhs=F, _element=*T, _expr= g);
		    }
	    }
#if 0      
	    itType = mapField.find( "DirichletX" );
	    if ( itType != mapField.end() )
	    {
	      for ( auto const& exAtMarker : (*itType).second )
		    {
		      std::string marker = exAtMarker.marker();
		      auto g = expr(exAtMarker.expression());
		      g.setParameterValues({{"t", mybdfT->time()}});
		      //Feel::cout << "A DirichletX[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
		      M00 += on(_range=markedfaces(mesh,marker), _rhs=F, _element=(*T)[Component::X], _expr= g);
		    }
	    }
	    itType = mapField.find( "DirichletY" );
	    if ( itType != mapField.end() )
	    {
	      for ( auto const& exAtMarker : (*itType).second )
		    {
		      std::string marker = exAtMarker.marker();
		      auto g = expr(exAtMarker.expression());
		      g.setParameterValues({{"t",mybdfT->time()}});
		      //Feel::cout << "A DirichletY[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
		      M00 += on(_range=markedfaces(mesh,marker), _rhs=F, _element=(*T)[Component::Y], _expr= g);
		    }
	    }
	    itType = mapField.find( "DirichletZ" );
	    if ( itType != mapField.end() )
	    {
	      for ( auto const& exAtMarker : (*itType).second )
		    {
		      std::string marker = exAtMarker.marker();
		      auto g = expr(exAtMarker.expression());
		      g.setParameterValues({{"t", mybdfT->time()}});
		      //Feel::cout << "A DirichletZ[" << marker << "] : " << exAtMarker.expression() << ", g=" << g << std::endl;
		      M00 += on(_range=markedfaces(mesh,marker), _rhs=F, _element=(*T)[Component::Z], _expr= g);
		    }
	    }
#endif      
	  }
               
    toc("boundary conditions", (M_verbose > 0));
    
    /* Solve */
    tic();
    auto result = mybackend->solve( _matrix=M, _rhs=F, _solution=U, _rebuild=true);
    std::string msg = (boost::format("[MQS %2%] t=%1% NbIter=%3% Residual=%4%") % mybdfT->time()
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
      
    toc("solve", (M_verbose > 0));

    // update A and V pointers from U
    myblockVecSol.localize(U);

    tic();
    e->step(mybdfT->time())->add( "T", T);

    if ( Texact )
	  {
	    Texact_g.setParameterValues({{"t", mybdfT->time()}});
	    Texact = project(_space = Th, _expr = Texact_g);

	    e->step(mybdfT->time())->add( "Aexact", Texact);
	  }
    e->save();
    toc("export", (M_verbose > 0));

    // Compute error
    if ( Texact )
	  {
	    L2Texact = normL2(_range = elements(mesh), _expr = Texact_g);
	    L2Terror = normL2(elements(mesh), (idv(T) - idv(Texact)));
	    H1Terror = normH1(elements(mesh), _expr = (idv(T) - idv(Texact)), _grad_expr = (gradv(T) - gradv(Texact)));
	    Feel::cout << "error: " << "t="<< mybdfT->time();
	    Feel::cout << " T: " << L2Terror << " " << L2Terror / L2Texact << " " << H1Terror << " ";
	  }
    
    mybdfA->next(*T);

    /* reinit  */
    M->zero();
    F->zero();

    Told = (*T);
  }
}