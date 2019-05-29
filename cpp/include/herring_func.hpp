//
// Created by houchen_li on 5/1/19.
//

#pragma once
#ifndef HERRING_FUNC_H
#define HERRING_FUNC_H

#include "def.h"
#include "utils.h"
#include "state_type.h"

void func(const RBCSystem::StateType &f, RBCSystem::StateType &dfdt, const Real &t, const Real &Pr, const Real &Ra, const Real &a) {
    Uint i, j, k;
    Real sum;

    for (j = 2; j < 2 * RBCSystem::Ny + 2; j += 2) {
        dfdt.phi(j) = -(j * j * PI_2 * f.phi(j));
    }
    for (i = 1; i < RBCSystem::Nx + 1; i++) {
        for (j = 2; j < 2 * RBCSystem::Ny + 2; j += 2) {
            sum = 0.0;
            for (k = 2 - i % 2; k < 2 * RBCSystem::Ny + 2 - i % 2 - j; k += 2)
                sum += f.psi(i, k) * f.theta(i, j + k);
            for (k = 2 - i % 2; k < j + i % 2; k += 2)
                sum -= f.psi(i, k) * f.theta(i, j - k);
            for (k = j + 2 - i % 2; k < 2 * RBCSystem::Ny + 2 - i % 2; k += 2)
                sum += f.psi(i, k) * f.theta(i, k - j);
            dfdt.phi(j) += (j / 2) * PI_2 * a * i * sum;
        }
    }

    for (i = 1; i < RBCSystem::Nx + 1; i++) {
        for (j = 2 - i % 2; j < 2 * RBCSystem::Ny + 2 - i % 2; j += 2) {
            dfdt.theta(i, j) = 2.0 * PI * a * i * f.psi(i, j) - (4 * a * a * i * i + j * j) * PI_2 * f.theta(i, j);
            sum = 0.0;
            for (k = 2; k < 2 * RBCSystem::Ny + 2 - i % 2 - j; k += 2)
                sum += k * f.phi(k) * f.psi(i, j + k);
            for (k = 2; k < j + i % 2; k += 2)
                sum += k * f.phi(k) * f.psi(i, j - k);
            for (k = j + 2 - i % 2; k < 2 * RBCSystem::Ny + 2; k += 2)
                sum -= k * f.phi(k) * f.psi(i, k - j);
            dfdt.theta(i, j) -= PI_2 * a * i * sum;
        }
    }

    for (i = 1; i < RBCSystem::Nx + 1; i++)
        for (j = 2 - i % 2; j < 2 * RBCSystem::Ny + 2 - i % 2; j += 2)
            dfdt.psi(i, j) = 2.0 * PI * a * i * Pr * Ra / (4 * a * a * i * i + j * j) * PI_2 * f.theta(i, j) -
                             (4 * a * a * i * i + j * j) * PI_2 * Pr * f.psi(i, j);
    /*for (i = 1; i < RBCSystem::Nx; i++) {
        for (j = 2 - i % 2; j < 2 * RBCSystem::Ny + i % 2; j += 2)
            dfdt.psi(i, j) -= PI_2 * a / 2 * f.psi(1, 1) * f.psi(i + 1, j + 1) *
                              (4 * a * a * ((i + 1) * (i + 1) - 1) + (j + 1) * (j + 1) - 1) * ((Real) i - (Real) j) /
                              (4 * a * a * i * i + j * j);
        for (j = 2 + i % 2; j < 2 * RBCSystem::Ny + 2 - i % 2; j += 2)
            dfdt.psi(i, j) -= PI_2 * a / 2 * f.psi(1, 1) * f.psi(i + 1, j - 1) *
                              (4 * a * a * ((i + 1) * (i + 1) - 1) + (j - 1) * (j - 1) - 1) * (i + j) /
                              (4 * a * a * i * i + j * j);
    }
    for (i = 2; i < RBCSystem::Nx + 1; i++) {
        for (j = 2 + i % 2; j < 2 * RBCSystem::Ny + 2 - i % 2; j += 2)
            dfdt.psi(i, j) += PI_2 * a / 2 * f.psi(1, 1) * f.psi(i - 1, j - 1) *
                              (4 * a * a * ((i - 1) * (i - 1) - 1) + (j - 1) * (j - 1) - 1) * ((Real) i - (Real) j) /
                              (4 * a * a * i * i + j * j);
        for (j = 2 - i % 2; j < 2 * RBCSystem::Ny + i % 2; j += 2)
            dfdt.psi(i, j) += PI_2 * a / 2 * f.psi(1, 1) * f.psi(i - 1, j + 1) *
                              (4 * a * a * ((i - 1) * (i - 1) - 1) + (j + 1) * (j + 1) - 1) * (i + j) /
                              (4 * a * a * i * i + j * j);
    }*/

    return;
};

#endif /* ifndef SYMBOL */
