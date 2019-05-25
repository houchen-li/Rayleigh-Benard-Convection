#include "../include/utils.h"

Real Power(const Real &val, const Uint &p)
{
    Uint i;
    Real res;

    res = 0;
    for (i = 0; i < p; i++) {
        res *= val;
    }

    return res;
}

void addAttr(H5::H5Object &h5obj, const String &attr_name, const Uint &val)
{
    hsize_t rank, dims[1];
    H5::DataSpace attr_space(H5S_SCALAR);
    H5::Attribute attr;

    rank = 1;
    dims[0] = 1;
    attr = h5obj.createAttribute(attr_name, H5::PredType::NATIVE_UINT64, attr_space);
    attr.write(H5::PredType::NATIVE_UINT64, &val);
    attr.close();

    return;
}

void addAttr(H5::H5Object &h5obj, const String &attr_name, const Real &val)
{
    hsize_t rank, dims[1];
    H5::DataSpace attr_space(H5S_SCALAR);
    H5::Attribute attr;

    rank = 1;
    dims[0] = 1;
    attr = h5obj.createAttribute(attr_name, H5::PredType::NATIVE_DOUBLE, attr_space);
    attr.write(H5::PredType::NATIVE_DOUBLE, &val);
    attr.close();

    return;
}

void addArray(H5::Group &h5group, const String &array_name, const Real *a, const Uint &n)
{
	hsize_t rank, dims[1];
	H5::DataSpace dspace;
	H5::DataSet dset;

	rank = 1;
	dims[0] = n;
	dspace = H5::DataSpace(rank, dims);
	dset = h5group.createDataSet(array_name, H5::PredType::NATIVE_DOUBLE, dspace);
	addAttr(dset, "length", n);
	dset.write(a, H5::PredType::NATIVE_DOUBLE);
	dset.close();

	return;
}
