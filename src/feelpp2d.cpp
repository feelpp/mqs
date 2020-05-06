#include <feel/feel.hpp>

int main(int argc, char**argv )
{
     using namespace Feel;

     Environment env( _argc=argc, _argv=argv,
                      _about=about(_name="qs_laplacian",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto Vh = Pch<2>( mesh );

    auto g = expr(soption(_name="functions.g"));
    Feel::cout << "g=" << g << std::endl;

    auto f = expr<2,1>(soption(_name="functions.f"));
    Feel::cout << "f=" << f << std::endl;

    double dt = doption(_name = "a");
    Feel::cout << "dt=" << dt << std::endl;

    double u0 = expr(soption(_name = "functions.u");
    Feel::cout << "u0=" << u0 << std::endl;

    double tmax = doption(_name = "tmax");
    Feel::cout << "tmax=" << tmax << std::endl;

    auto u = Vh->element(u0);
    auto un = Vh->element();
    auto v = Vh->element();

    auto l = form1( _test=Vh );

    double t = 0;

    while(t < tmax){
        l = integrate(_range=elements(mesh),
                    _expr =  (id(u) + dt*id(f)) * trans(id(v)) );

        auto a = form2( _trial=Vh, _test=Vh);
        a = integrate(_range=elements(mesh),
                    _expr = id(un)*trans(id(v)) + dt*grad(un)*trans(grad(v)) );
        a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=g );
        a.solve(_rhs=l,_solution=u);
        t += dt;
    }

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->save();
}