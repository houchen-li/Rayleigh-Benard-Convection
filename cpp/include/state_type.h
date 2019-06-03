//
// Created by houchen_li on 5/1/19.
//

#pragma once
#ifndef STATE_TYPE_H
#define STATE_TYPE_H

#include <exception>
#include <boost/numeric/ublas/vector.hpp>
#include <H5Cpp.h>
#include "def.h"
#include "utils.h"

namespace RBCSystem {
    constexpr Uint Nx = 10;
    constexpr Uint Ny = 20;

    class StateType : public boost::numeric::ublas::fixed_vector<Real, (2 * Nx + 1) * Ny> {
    public:
        StateType(void);

        StateType(const StateType &rhs);

        StateType &operator=(const StateType &rhs);

        ~StateType();

        inline Real &phi(const Uint &j) { return operator()(j / 2 - 1); }

        inline const Real &phi(const Uint &j) const { return operator()(j / 2 - 1); }

        inline Real &theta(const Uint &i, const Uint &j) { return operator()(Ny * i + (j + i % 2) / 2 - 1); }

        inline const Real &theta(const Uint &i, const Uint &j) const { return operator()(Ny * i + (j + i % 2) / 2 - 1); }

        inline Real &psi(const Uint &i, const Uint &j) { return operator()(Ny * (Nx + i) + (j + i % 2) / 2 - 1); }

        inline const Real &psi(const Uint &i, const Uint &j) const { return operator()(Ny * (Nx + i) + (j + i % 2) / 2 - 1); }

        void loadFile(const H5::Group &h5group, const String &dset_name);

        inline void loadFile(const H5::Group &h5group, const char *dset_name) { loadFile(h5group, static_cast<String>(dset_name)); return; }

        void saveFile(H5::Group &h5group, const String &dset_name) const;

        inline void saveFile(H5::Group &h5group, const char *dset_name) const { saveFile(h5group, static_cast<String>(dset_name)); return; }

        static StateType TrivialState(void);
    };

    Real Nuof(const StateType &f);
}
#endif //STATE_TYPE_H
