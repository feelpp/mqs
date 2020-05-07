#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>


int main(int argc, char**argv )
{
     using namespace Feel;

     Environment env( _argc=argc, _argv=argv,
                      _about=about(_name="heat",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org"));


    auto g = expr(soption(_name="functions.g"));
    Feel::cout << "g=" << g << std::endl;

    auto f = expr(soption(_name="functions.f"));
    Feel::cout << "f=" << f << std::endl;

    auto u0 = expr(soption(_name = "functions.u"));
    Feel::cout << "u0=" << u0 << std::endl;

    double tmax = doption(_name = "ts.tmax");

    double dt = doption(_name = "ts.dt");

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto Vh = Pch<2>( mesh );
    auto u = Vh->element(u0);
    auto v = Vh->element();

    auto l = form1( _test=Vh );

    double t = dt;

    auto e = exporter( _mesh=mesh );

    while(t < tmax){
        l = integrate(_range=elements(mesh),
                    _expr =  (idv(u) + dt*f) * id(v));

        auto a = form2( _trial=Vh, _test=Vh);
        a = integrate(_range=elements(mesh),
                    _expr = idt(u)*id(v) + dt*gradt(u)*trans(grad(v)) );
        a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=g );
        a.solve(_rhs=l,_solution=u);
        e->step(t)->add( "u", u);
        e->save();
        t += dt;
    }
}