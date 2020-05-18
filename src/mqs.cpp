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

    auto Ad = expr(soption(_name="functions.Ad"));
    Feel::cout << "Ad=" << Ad << std::endl;

    auto A0 = expr(soption(_name="functions.A"));
    Feel::cout << "A0=" << A0 << std::endl;
     
    auto gI = expr(soption(_name="functions.gI"));
    Feel::cout << "gI=" << gI << std::endl;

    auto g0 = expr(soption(_name="functions.gO"));
    Feel::cout << "gO=" << gO << std::endl;

    //Recuperer mu,sigma,

    double dt = doption(_name = "ts.time-step");
    std::cout << "time-step=" << u0 << std::endl;

    double tmax = doption(_name = "ts.time-final");
    std::cout << "time-final=" << u0 << std::endl;

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto cond_mesh = createSubmesh(mesh,markedelements(mesh,"Omega_c"));

    auto Ah = Pchv<3>( mesh );
    auto Vh = Pch<3>( cond_mesh );

    auto A = Ah->element(A0);
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

        /*
        l2 = integrate(_range=elements(cond_mesh),
                    _expr = sigma * inner( idv(A) - idt(A), grad(psi) );
        
        a2 += integrate( range=elements(cond_mesh),
                    _expr = sigma * dt * inner(idt(V), grad(psi) );

        a2 += on(_range=markedfaces(cond_mesh,"Omega_I"), _rhs=l2, _element=psi, 
                _expr= gI );
        a2 += on(_range=markedfaces(cond_mesh,"Omega_O"), _rhs=l2, _element=psi, 
                _expr= gO );
        */              
        e->step(t)->add( "a1", a1);
        e->save();
        t += dt;
    }
}