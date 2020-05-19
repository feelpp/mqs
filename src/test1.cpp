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

    auto A0 = expr<3,1>(soption(_name="functions.a"));
    Feel::cout << "A0=" << A0 << std::endl;
     
    auto gI = expr<3,1>(soption(_name="functions.i"));
    Feel::cout << "gI=" << gI << std::endl;

    auto gO = expr<3,1>(soption(_name="functions.o"));
    Feel::cout << "gO=" << gO << std::endl;

#if 0
    auto gC = expr<3,1>(soption(_name="functions.C"));
    Feel::cout << "gC=" << gC << std::endl;
#endif

    auto V = expr(soption(_name="functions.v"));
    Feel::cout << "V=" << V << std::endl;

    //Recuperer mu,sigma,

    double mu = 1;

    double sigma = 58000;

    double dt = doption(_name = "ts.time-step");
    std::cout << "time-step=" << dt << std::endl;

    double tmax = doption(_name = "ts.time-final");
    std::cout << "time-final=" << tmax << std::endl;

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);
    auto cond_mesh = createSubmesh(mesh,markedelements(mesh,"Omega_C"));

    auto Ah = Pchv<3>( mesh );

    auto A = Ah->element(A0);

    auto phi = Ah->element();

    auto l1 = form1( _test=Ah );

    double t = dt;

    auto e = exporter( _mesh=mesh );
    auto a1 = form2( _trial=Ah, _test=Ah);

    while(t < tmax){

        l1 = integrate(_range=elements(cond_mesh),
                        _expr = sigma * inner(id(phi) , idv(A) - dt*trans(grad<3>(V))) );
        
        a1 = integrate(_range=elements(mesh),
                    _expr = (dt/mu) * inner(curl(phi) , curlt(A)) );
        a1 += integrate(_range=elements(cond_mesh),
                    _expr = sigma * inner(id(phi) , idt(A) ));
        a1 += on(_range=markedfaces(mesh,"Gamma_I"), _rhs=l1, _element=phi, 
                _expr= gI );
        a1 += on(_range=markedfaces(mesh,"Gamma_O"), _rhs=l1, _element=phi, 
                _expr= gO);
#if 0
        a1 += on(_range=markedfaces(mesh,"Gamma_C"), _rhs=l1, _element=phi, 
                _expr= gC);
#endif
        a1.solve(_rhs=l1,_solution=A);
       
        e->step(t)->add( "A", A);
        e->save();
        t += dt;
    }
}