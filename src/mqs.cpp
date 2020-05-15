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

    auto Vh = Pchv<3>( mesh );
    auto Ah = Pch<3>( cond_mesh );

    auto V = Vh->element();
    auto A = Ah->element();

    auto phi = Vh->element();
    //auto psi = Vh->element();

    auto l1 = form1( _test=Vh );
    //auto l2 = form1( _test=Vh );

    double t = dt;

    auto e = exporter( _mesh=mesh );
    auto a = form2( _trial=Vh, _test=Vh);
    while(t < tmax){

        l1 = integrate(_range=elements(cond_mesh),_expr = -sigma*inner(id(phi),trans(grad(V)) - idv(A)/dt));
        //l1 = integrate(_range=elements(cond_mesh),_expr = sigma*grad(V)*id(phi) + trans(id(phi))*idv(A)/dt));
        
        a1 = integrate(_range=elements(mesh),
                    _expr = (1/mu)*inner(curl(phi),curlt(A)));
        a1 += on(_range=markedfaces(mesh,"Omega_D"), _rhs=l1, _element=phi, _expr= Ad*inner(curlt(A)) );
        a1 += on( range=elements(cond_mesh),_expr = sigma*inner(id(phi),trans(grad(V)) + idt(A)/dt));

        a1.solve(_rhs=l1,_solution=a1);

        /*l2 = integrate(_range=elements(mesh),_expr = 0);
        auto a = form2( _trial=Vh, _test=Vh);
        a2 = integrate(_range=elements(mesh),
                    _expr = );*/

        e->step(t)->add( "a1", a1);
        e->save();
        t += dt;
    }
}