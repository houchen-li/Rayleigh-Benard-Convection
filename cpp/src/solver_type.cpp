#include "../include/solver_type.h"

RBCSystem::SolverType::SolverType(const Real &Prandtl_number, const Real &Rayleigh_number, const Real &wave_number,
        const StateType &start_state, const Real &start_time, const Real &end_time):
        l(0), Pr(Prandtl_number), Ra(Rayleigh_number), a(wave_number), f(), t(), observer(l, Pr, Ra, a)
{
    f[0] = start_state;
    t[0] = start_time;
    t[1] = end_time;

    return;
}

RBCSystem::SolverType::SolverType(const SolverType &rhs):
        l(rhs.l), Pr(rhs.Pr), Ra(rhs.Ra), a(rhs.a), f(), t(), observer(l, Pr, Ra, a)
{
    f[0] = rhs.f[0];
    f[1] = rhs.f[1];
    t[0] = rhs.t[0];
    t[1] = rhs.t[1];

    return;
}

RBCSystem::SolverType & RBCSystem::SolverType::operator=(const SolverType &rhs)
{
    l = rhs.l;
    Pr = rhs.Pr;
    Ra = rhs.Ra;
    a = rhs.a;
    f[0] = rhs.f[0];
    f[1] = rhs.f[1];
    t[0] = rhs.t[0];
    t[1] = rhs.t[1];

    return *this;
}

RBCSystem::SolverType::~SolverType(void) { return; }

void RBCSystem::SolverType::evaluate(void)
{
    using namespace boost::numeric::odeint;

    f[1] = f[0];
    integrate_adaptive(make_controlled<runge_kutta_dopri5<StateType>>(1e-10, 1e-6), boost::ref(*this), f[1], t[0], t[1], 0.01, observer);

	return;
}

void RBCSystem::SolverType::saveFile(H5::Group &h5group, const String &subgroup_name) const
{
    Uint i;
    String dset_name;
    H5::Group subgrp, subsubgrp;
    H5::DataSet dset;

    subgrp = h5group.createGroup(subgroup_name);
    addAttr(subgrp, "Nx", Nx);
    addAttr(subgrp, "Ny", Ny);
    addAttr(subgrp, "length", l);
    addAttr(subgrp, "Prandtl_number", Pr);
    addAttr(subgrp, "Rayleigh_number", Ra);
    addAttr(subgrp, "wave_number", a);
    subsubgrp = subgrp.createGroup("states");
    for (i = 0; i < 2; i++) {
        dset_name = std::to_string(i);
        f[i].saveFile(subsubgrp, dset_name);
        dset = subsubgrp.openDataSet(dset_name);
        addAttr(dset, "t", t[i]);
        addAttr(dset, "Nu", Nuof(f[i]));
        dset.close();
    }
    subsubgrp.close();
    subgrp.close();

	return;
}
