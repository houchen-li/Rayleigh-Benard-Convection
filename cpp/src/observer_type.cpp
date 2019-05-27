//
// Created by houchen_li on 5/1/19.
//

#include "../include/observer_type.h"

RBCSystem::ObserverType::ObserverType(Uint &l, Real &Pr, Real &Ra, Real &a) : l(l), Pr(Pr), Ra(Ra), a(a)
{
    FILE *fp = fopen("Nu.dat", "w");

    fprintf(fp, "************************\n");
    fprintf(fp, "%s\t%lf\n", "Prandtl_number", this -> Pr);
    fprintf(fp, "%s\t%lf\n", "Rayleigh_number", this -> Ra);
    fprintf(fp, "%s\t%lf\n", "wave_number", this -> a);
    fprintf(fp, "************************\n");
    fprintf(fp, "t\t\tNu\n");
    fprintf(fp, "************************\n");
    fclose(fp);

    return;
}

RBCSystem::ObserverType::~ObserverType(void) { return; }

void RBCSystem::ObserverType::operator()(const RBCSystem::StateType &f, const Real &t)
{
    Uint j;
    Real Nu;
    FILE *fp = fopen("Nu.dat", "a");
    H5::DataSet dset;
#ifdef OBSERVE_STATES
    H5::H5File outf(std::to_string(l)+".h5", H5F_ACC_TRUNC);
#else
    H5::H5File outf("current_state.h5", H5F_ACC_TRUNC);
#endif

    Nu = Nuof(f);
    fprintf(fp, "%8lf\t%8lf\n", t, Nu);
    fclose(fp);
    f.saveFile(outf, "data");
    dset = outf.openDataSet("data");
    addAttr(dset, "t", t);
    addAttr(dset, "Nu", Nu);
    dset.close();
    outf.close();
    l++;

    return;
}


