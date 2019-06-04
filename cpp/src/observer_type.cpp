//
// Created by houchen_li on 5/1/19.
//

#include "../include/observer_type.h"

RBCSystem::ObserverType::ObserverType(Uint &l, Real &Pr, Real &Ra, Real &a) : l(l), Pr(Pr), Ra(Ra), a(a) { return; }

RBCSystem::ObserverType::~ObserverType(void) { return; }

void RBCSystem::ObserverType::startObservation(void)
{
    FILE *fp = fopen("Nu.dat", "w");

    fprintf(fp, "************************************\n");
    fprintf(fp, "%16s\t%16.1lf\n", "Prandtl_number", this -> Pr);
    fprintf(fp, "%16s\t%16.1lf\n", "Rayleigh_number", this -> Ra);
    fprintf(fp, "%16s\t%16.1lf\n", "wave_number", this -> a);
    fprintf(fp, "************************************\n");
    fprintf(fp, "               t\t              Nu\n");
    fprintf(fp, "************************************\n");
    fclose(fp);

    return;
}

void RBCSystem::ObserverType::endObservation(void)
{
    FILE *fp = fopen("Nu.dat", "a");
    fprintf(fp, "************************************\n");
    fclose(fp);

    return;
}

void RBCSystem::ObserverType::operator()(const RBCSystem::StateType &f, const Real &t)
{
    Real Nu;
    FILE *fp = fopen("Nu.dat", "a");
    H5::Group grp, subgrp;
    H5::DataSet dset;
#ifdef OBSERVE_STATES
    H5::H5File outf(std::to_string(l)+".h5", H5F_ACC_TRUNC);
#else
    H5::H5File outf("current_state.h5", H5F_ACC_TRUNC);
#endif

    Nu = Nuof(f);
    fprintf(fp, "%16.12lf\t%16.12lf\n", t, Nu);
    fclose(fp);
    grp = outf.createGroup("group");
    addAttr(grp, "Nx", Nx);
    addAttr(grp, "Ny", Ny);
    addAttr(grp, "length", l);
    addAttr(grp, "Prandtl_number", Pr);
    addAttr(grp, "Rayleigh_number", Ra);
    addAttr(grp, "wave_number", a);
    subgrp = grp.createGroup("states");
    f.saveFile(subgrp, "1");
    dset = subgrp.openDataSet("1");
    addAttr(dset, "t", t);
    addAttr(dset, "Nu", Nu);
    dset.close();
    subgrp.close();
    grp.close();
    outf.close();
    l++;

    return;
}
