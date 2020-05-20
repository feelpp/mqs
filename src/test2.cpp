#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>


int main(int argc, char**argv )
{
     using namespace Feel;
    po::options_description options( "Laplacian options" );
    options.add_options()
      ( "hc", po::value<std::string>()->default_value( "1" ), "convective heat transfer coefficient" )
      ( "Tinf", po::value<std::string>()->default_value( "1" ), "Temperature far from boundary" );
     Environment env( _argc=argc, _argv=argv,_desc=options,
                      _about=about(_name="heat",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org"));
     
    auto gI = expr<3,1>(soption(_name="functions.i"));
    Feel::cout << "gI=" << gI << std::endl;

    auto gO = expr<3,1>(soption(_name="functions.o"));
    Feel::cout << "gO=" << gO << std::endl;

    auto dA = expr<3,1>(soption(_name="functions.a"));
    Feel::cout << "A=" << A << std::endl;

    //Recuperer sigma,

    double sigma = 58000;

    double dt = doption(_name = "ts.time-step");
    std::cout << "time-step=" << dt << std::endl;

    double tmax = doption(_name = "ts.time-final");
    std::cout << "time-final=" << tmax << std::endl;

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);
    auto cond_mesh = createSubmesh(mesh,markedelements(mesh,"Omega_C"));

    auto Ah = Pchv<3>( mesh );
    auto Vh = Pch<3>( cond_mesh );

    auto V = Vh->element();

    auto psi = Vh->element();

    auto l2 = form1( _test=Vh );

    double t = dt;

    auto e = exporter( _mesh=mesh );
    auto a2 = form2( _trial=Vh, _test=Ah);

    while(t < tmax){

        l2 = integrate(_range=elements(cond_mesh),
                    _expr = sigma * inner( dA, grad(psi) );
        
        a2 += integrate(_range=elements(cond_mesh),
                    _expr = sigma * dt * inner(idt(V), grad(psi) );

        a2 += on(_range=markedfaces(cond_mesh,"Gamma_I"), _rhs=l2, _element=psi, 
                _expr= gI );
        a2 += on(_range=markedfaces(cond_mesh,"Gamma_O"), _rhs=l2, _element=psi, 
                _expr= gO );

        a2.solve(_rhs=l2,_solution=V);
       
        e->step(t)->add( "V", V);
        e->save();
        t += dt;
    }
}