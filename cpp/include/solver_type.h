//
// Created by houchen_li on 5/1/19.
//

#pragma once
#ifndef SOLVER_TYPE_H
#define SOLVER_TYPE_H

#include <boost/numeric/odeint.hpp>
#include <boost/ref.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/constants/constants.hpp>
#include <H5Cpp.h>
#include "def.h"
#include "utils.h"
#include "state_type.h"
#include "observer_type.h"

namespace RBCSystem {

    class SolverType {
    public:
        explicit SolverType(const Real &Prandtl_number = 1.0, const Real &Rayleigh_number = 1000.0,
                const Real &wave_number = 1.0, const StateType &start_state = StateType::TrivialState(),
                const Real &start_time = 0.0, const Real &end_time = 1.0);

        SolverType(const SolverType &rhs);

        SolverType &operator=(const SolverType &rhs);

        ~SolverType(void);

        inline const Uint &length(void) const { return l; }

        inline const Real &Prandtl_number(void) const { return Pr; }

        inline void Prandtl_number(const Real &Prandtl_number) { Pr = Prandtl_number; return; }

        inline const Real &Rayleigh_number(void) const { return Ra; }

        inline void Rayleigh_number(const Real &Rayleigh_number) { Ra = Rayleigh_number; return; }

        inline const Real &wave_number(void) const { return a; }

        inline void wave_number(const Real &wave_number) { a = wave_number; return; }

        inline const RBCSystem::StateType &start_state(void) const { return f[0]; }

        inline void start_state(const RBCSystem::StateType &start_state) { f[0] = start_state; return; }

        inline const RBCSystem::StateType &end_state(void) const { return f[1]; }

        inline const Real &start_time() const { return t[0]; }

        inline void start_time(const Real &start_time) { t[0] = start_time; return; }

        inline const Real &end_time() const { return t[1]; }

        inline void end_time(const Real &end_time) { t[1] = end_time; return; }

        void evaluate();

        void saveFile(H5::Group &h5group, const String &subgroup_name) const;

        inline void saveFile(H5::Group &h5group, const char *subgroup_name) const { saveFile(h5group, static_cast<String>(subgroup_name)); return; }

        virtual void operator()(const StateType &f, StateType &dfdt, const Real &t) const = 0;

    private:
        Uint l;
        Real Pr, Ra, a;
        StateType f[2];
        Real t[2];
        ObserverType observer;
    };

};

#endif /* end of include guard: SOLVER_TYPE_H */
