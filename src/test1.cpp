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

    /*
    auto Ad = expr(soption(_name="functions.Ad"));
    Feel::cout << "Ad=" << Ad << std::endl;
    */

    auto A0 = expr(soption(_name="functions.A"));
    Feel::cout << "A0=" << A0 << std::endl;
     
    auto gI = expr(soption(_name="functions.gI"));
    Feel::cout << "gI=" << gI << std::endl;

    auto gO = expr(soption(_name="functions.gO"));
    Feel::cout << "gO=" << gO << std::endl;

    auto gC = expr(soption(_name="functions.gC"));
    Feel::cout << "gC=" << gC << std::endl;

    auto V = expr(soption(_name="functions.V"));
    Feel::cout << "V=" << V << std::endl;

    //Recuperer mu,sigma,

    double mu = 1

    double sigma = 58000

    double dt = doption(_name = "ts.time-step");
    std::cout << "time-step=" << u0 << std::endl;

    double tmax = doption(_name = "ts.time-final");
    std::cout << "time-final=" << u0 << std::endl;

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
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
                        _expr = sigma * inner(id(phi) , idv(A) - dt*grad(V)) );
        
        a1 = integrate(_range=elements(mesh),
                    _expr = (dt/mu) * inner(curl(phi) , curlt(A)) );
        a1 += integrate( range=elements(cond_mesh),
                    _expr = sigma * inner(id(phi) , idt(A) ));
        a1 += on(_range=markedfaces(mesh,"Gamma_I"), _rhs=l1, _element=phi, 
                _expr= gI );
        a1 += on(_range=markedfaces(mesh,"Gamma_O"), _rhs=l1, _element=phi, 
                _expr= gO);
        a1 += on(_range=markedfaces(mesh,"Gamma_C"), _rhs=l1, _element=phi, 
                _expr= gC);
        /*
        a1 += on(_range=markedfaces(mesh,"Gamma_D"), _rhs=l1, _element=phi, 
                _expr= Ad );
        */

        a1.solve(_rhs=l1,_solution=a1);
       
        e->step(t)->add( "a1", a1);
        e->save();
        t += dt;
    }
}