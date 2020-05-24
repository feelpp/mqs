#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

int main(int argc, char **argv)
{
        using namespace Feel;
        po::options_description options("Laplacian options");
        options.add_options()("hc", po::value<std::string>()->default_value("1"), "convective heat transfer coefficient")("Tinf", po::value<std::string>()->default_value("1"), "Temperature far from boundary");
        Environment env(_argc = argc, _argv = argv, _desc = options,
                        _about = about(_name = "test1",
                                       _author = "Feel++ Consortium",
                                       _email = "feelpp-devel@feelpp.org"));

        auto A0 = expr<3, 1>(soption(_name = "functions.a"));
        Feel::cout << "A0=" << A0 << std::endl;

        auto gI = expr<3, 1>(soption(_name = "functions.i"));
        Feel::cout << "gI=" << gI << std::endl;

        auto gO = expr<3, 1>(soption(_name = "functions.o"));
        Feel::cout << "gO=" << gO << std::endl;
#if 0
    auto V = expr(soption(_name="functions.v"));
    Feel::cout << "V=" << V << std::endl;
#endif
        auto gradV = expr<3, 1>(soption(_name = "functions.v"));
        Feel::cout << "gradV=" << gradV << std::endl;

        auto Ad = expr<3, 1>(soption(_name = "functions.d"));
        Feel::cout << "Ad=" << Ad << std::endl;

        auto mu = expr(soption(_name = "functions.m"));
        Feel::cout << "mu=" << mu << std::endl;

        auto sigma = expr(soption(_name = "functions.s"));
        Feel::cout << "sigma=" << sigma << std::endl;

        auto Aexact_g = expr<3, 1>(soption(_name = "functions.e"));
        Feel::cout << "Aexact=" << Aexact_g << std::endl;

        double dt = doption(_name = "ts.time-step");
        std::cout << "time-step=" << dt << std::endl;

        double tmax = doption(_name = "ts.time-final");
        std::cout << "time-final=" << tmax << std::endl;

        auto mesh = loadMesh(_mesh = new Mesh<Simplex<3>>);
        auto cond_mesh = createSubmesh(mesh, markedelements(mesh, "Omega_C"));

        auto Ah = Pchv<1>(mesh);

        auto A = Ah->element(A0);

        auto phi = Ah->element();

        auto Aexact = Ah->element();

        auto l1 = form1(_test = Ah);

        double t = 0;

        auto e = exporter(_mesh = mesh);
        auto a1 = form2(_trial = Ah, _test = Ah);

        Aexact_g.setParameterValues({{"t", t}});
        Aexact = project(_space = Ah, _expr = Aexact_g);

        gradV.setParameterValues({{"t", t}});
        e->step(t)->add("A", A0);
        e->step(t)->add("Aexact", Aexact);
        e->save();

        double L2Aexact = normL2(_range = elements(mesh), _expr = Aexact_g);
        double H1error = 0;
        double L2error = 0;
        Feel::cout << "H1 error at t = " << t << ": " << H1error << std::endl;
        Feel::cout << "L2 error at t = " << t << ": " << L2error << std::endl;

        for (t = dt; t < max; t += dt)
#if 0
        l1 = integrate(_range=elements(cond_mesh),
                        _expr = sigma * inner(id(phi) , idv(A) - dt*trans(grad<3>(V))) );
#endif
        Aexact_g.setParameterValues({{"t", t}});
        Aexact = project(_space = Ah, _expr = Aexact_g);
        gradV.setParameterValues({{"t", t}});
        gO.setParameterValues({{"t", t}});
        gI.setParameterValues({{"t", t}});
        Ad.setParameterValues({{"t", t}});

        l1.zero();
        a1.zero();
        l1 = integrate(_range = elements(cond_mesh),
                       _expr = sigma * inner(id(phi), idv(A) - dt * gradV));
        a1 = integrate(_range = elements(mesh),
                       _expr = (dt / mu) * inner(curl(phi), curlt(A)));
        a1 += integrate(_range = elements(cond_mesh),
                        _expr = sigma * inner(id(phi), idt(A)));
        a1 += on(_range = markedfaces(mesh, "Gamma_I"), _rhs = l1, _element = phi,
                 _expr = gI);
        a1 += on(_range = markedfaces(mesh, "Gamma_O"), _rhs = l1, _element = phi,
                 _expr = gO);
        a1 += on(_range = markedfaces(mesh, "Gamma_C"), _rhs = l1, _element = phi,
                 _expr = Ad);

        a1.solve(_rhs = l1, _solution = A);

        e->step(t)->add("A", A);
        e->step(t)->add("Aexact", Aexact_g);
        e->save();

        L2Aexact = normL2(_range = elements(mesh), _expr = Aexact_g);
        L2error = normL2(elements(mesh), (idv(A) - idv(Aexact)));
        H1error = normH1(elements(mesh), _expr = (idv(A) - idv(Aexact)), _grad_expr = (gradv(A) - gradv(Aexact)));

        Feel::cout << t << "" << L2error << " " << L2error / L2Aexact << " " << H1error << std::endl;
}
}