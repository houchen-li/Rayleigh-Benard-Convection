//
// Created by houchen_li on 5/1/19.
//

#include "../include/state_type.h"

RBCSystem::StateType::StateType(void): boost::numeric::ublas::fixed_vector<Real, (2 * Nx + 1) * Ny>() { return; }

RBCSystem::StateType::StateType(const StateType &rhs): boost::numeric::ublas::fixed_vector<Real, (2 * Nx + 1) * Ny>(rhs) { return; }

RBCSystem::StateType & RBCSystem::StateType::operator=(const StateType &rhs) { boost::numeric::ublas::fixed_vector<Real, (2 * Nx + 1) * Ny>::operator=(rhs); return *this; }

RBCSystem::StateType::~StateType(void) { return; }

void RBCSystem::StateType::loadFile(const H5::Group &h5group, const String &dset_name)
{
    Uint nx, ny;
    H5::DataSet dset;
    H5::Attribute attr;

    dset = h5group.openDataSet(dset_name);
    attr = dset.openAttribute("Nx");
    attr.read(H5::PredType::NATIVE_UINT64, &nx);
    attr.close();
    attr = dset.openAttribute("Ny");
    attr.read(H5::PredType::NATIVE_UINT64, &ny);
    attr.close();
    if (nx != Nx || ny != Ny)
        throw H5::FileIException("RBCSystem::StateType::loadFile", "the size of the reading data set is not consistent with the structure in used.");
    dset.read(data().data(), H5::PredType::NATIVE_DOUBLE);
    dset.close();

    return;
}

void RBCSystem::StateType::saveFile(H5::Group &h5group, const String &dset_name) const
{
    hsize_t rank, dims[1];
    H5::DataSpace dspace;
    H5::DataSet dset;

    rank = 1;
    dims[0] = (2*Nx+1)*Ny;
    dspace = H5::DataSpace(rank, dims);
    dset = h5group.createDataSet(dset_name, H5::PredType::NATIVE_DOUBLE, dspace);
    dset.write(data().data(), H5::PredType::NATIVE_DOUBLE);
    addAttr(dset, "Nx", Nx);
    addAttr(dset, "Ny", Ny);
    dset.close();

    return;
}

RBCSystem::StateType RBCSystem::StateType::TrivialState(void)
{
    Uint i;
    StateType res;

    for (i = 0; i < (2*Nx+1)*Ny; i++) {
        res[i] = 0.0;
    }
    res.psi(1, 1) = 1.0;

    return res;
}

Real RBCSystem::Nuof(const StateType &f)
{
    Uint j;
    Real res;

    res = 1.0;
    for (j = 2; j < 2*Ny+2; j+=2)
        res -= j*PI*f.phi(j);

    return res;
}
