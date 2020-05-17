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


    auto g = expr(soption(_name="functions.g"));
    Feel::cout << "g=" << g << std::endl;

    auto f = expr(soption(_name="functions.f"));
    Feel::cout << "f=" << f << std::endl;

    auto h = expr(soption(_name="functions.h"));
    Feel::cout << "h=" << h << std::endl;

    //auto Ad = expr(soption(_name="functions.Ad"));
    //Feel::cout << "h=" << h << std::endl;

    //Recuperer mu,sigma,

    //auto u0 = expr(soption(_name = "functions.u"));
    //Feel::cout << "u0=" << u0 << std::endl;

    double dt = doption(_name = "ts.time-step");
    std::cout << "time-step=" << u0 << std::endl;

    double tmax = doption(_name = "ts.time-final");
    std::cout << "time-final=" << u0 << std::endl;

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto cond_mesh = createSubmesh(mesh,markedelements(mesh,"Omega_c"));

    auto Ah = Pchv<3>( mesh );
    auto Vh = Pch<3>( cond_mesh );

    auto A = Ah->element();
    auto V = Vh->element();

    auto phi = Ah->element();
    auto psi = Vh->element();

    auto l1 = form1( _test=Ah );
    auto l2 = form1( _test=Vh );

    double t = dt;

    auto e = exporter( _mesh=mesh );
    auto a1 = form2( _trial=Ah, _test=Ah);
    auto a2 = form2( _trial=Vh, _test=Ah);

    while(t < tmax){

        l1 = integrate(_range=elements(cond_mesh),
                        _expr = sigma * inner(id(phi) , idv(A) - dt*grad(V)) );
        //l1 = integrate(_range=elements(cond_mesh),_expr = sigma*grad(V)*id(phi) + trans(id(phi))*idv(A)/dt));
        
        a1 = integrate(_range=elements(mesh),
                    _expr = (dt/mu) * inner(curl(phi) , curlt(A)) );
        a1 += integrate( range=elements(cond_mesh),
                    _expr = sigma * inner(id(phi) , idt(A) ));
        a1 += on(_range=markedfaces(mesh,"Gamma_D"), _rhs=l1, _element=phi, 
                _expr= Ad );

        a1.solve(_rhs=l1,_solution=a1);

        /*l2 = integrate(_range=elements(cond_mesh),
                    _expr = sigma * inner( idv(A) - idt(A), grad(phi) );
        
        a2 += integrate( range=elements(cond_mesh),
                    _expr = sigma * dt * inner(idt(V), grad(phi) );
        */           
        e->step(t)->add( "a1", a1);
        e->save();
        t += dt;
    }
}