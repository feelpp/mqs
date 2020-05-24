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
                      _about=about(_name="test2",
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org"));

    auto V0 = expr(soption(_name = "functions.v"));
    Feel::cout << "V0=" << V0 << std::endl;
    
    auto gI = expr(soption(_name="functions.i"));
    Feel::cout << "gI=" << gI << std::endl;

    auto gO = expr(soption(_name="functions.o"));
    Feel::cout << "gO=" << gO << std::endl;

    auto Ad = expr(soption(_name = "functions.d"));
    Feel::cout << "Ad=" << Ad << std::endl;

    auto dA = expr<3,1>(soption(_name="functions.a"));
    Feel::cout << "dA=" << dA << std::endl;

    auto sigma = expr(soption(_name="functions.s"));
    Feel::cout << "sigma=" << sigma << std::endl;

    auto Vexact_g = expr(soption(_name = "functions.e"));
    Feel::cout << "Vexact=" << Vexact_g << std::endl; 

    double dt = doption(_name = "ts.time-step");
    std::cout << "time-step=" << dt << std::endl;

    double tmax = doption(_name = "ts.time-final");
    std::cout << "time-final=" << tmax << std::endl;

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);
    auto cond_mesh = createSubmesh(mesh,markedelements(mesh,"Omega_C"));

    auto Ah = Pchv<1>( mesh );
    auto Vh = Pch<1>( cond_mesh );

    auto V = Vh->element(V0);

    auto psi = Vh->element();

    auto Vexact = Vh->element();

    auto l2 = form1( _test=Vh );

    double t = 0;

    auto e = exporter( _mesh=mesh );
    auto a2 = form2( _trial=Vh, _test=Vh);

    Vexact_g.setParameterValues({{"t", t}});
    Vexact = project(_space = Vh, _expr = Vexact_g);

    dA.setParameterValues({{"t", t}});
    e->step(t)->add("V", V0);
    e->step(t)->add("Vexact", Vexact);
    e->save();
    
    double L2Vexact = normL2(_range = elements(mesh), _expr = Vexact_g);
    double H1error = 0;
    double L2error = 0;
    Feel::cout << "H1 error at t = " << t << ": " << H1error << std::endl;
    Feel::cout << "L2 error at t = " << t << ": " << L2error << std::endl;

    for (t = dt; t < tmax; t += dt){
        Vexact_g.setParameterValues({{"t", t}});
        Vexact = project(_space = Vh, _expr = Vexact_g);
        dA.setParameterValues({{"t", t}});
        gO.setParameterValues({{"t", t}});
        gI.setParameterValues({{"t", t}});
        Ad.setParameterValues({{"t", t}});

        l2.zero();
        a2.zero();

        l2 = integrate(_range=elements(cond_mesh),
                    _expr = sigma * inner( -dA, trans(grad(psi)) ));
        
        a2 = integrate(_range=elements(cond_mesh),
                    _expr = sigma * inner(gradt(V), grad(psi) ));

        a2 += on(_range=markedfaces(cond_mesh,"Gamma_I"), _rhs=l2, _element=psi, 
                _expr= gI );
        a2 += on(_range=markedfaces(cond_mesh,"Gamma_O"), _rhs=l2, _element=psi, 
                _expr= gO );
        a2 += on(_range=markedfaces(cond_mesh,"Gamma_C"), _rhs=l2, _element=psi, 
                _expr= Ad );

        a2.solve(_rhs=l2,_solution=V);
        e->step(t)->add( "V", V);
        e->step(t)->add("Vexact", Vexact_g);
        e->save();

        L2Vexact = normL2(_range = elements(mesh), _expr = Vexact_g);
        L2error = normL2(elements(mesh), (idv(V) - idv(Vexact)));
        H1error = normH1(elements(mesh), _expr = (idv(V) - idv(Vexact)), _grad_expr = (gradv(V) - gradv(Vexact)));

        Feel::cout << t << "" << L2error << " " << L2error / L2Vexact << " " << H1error << std::endl;

    }
}