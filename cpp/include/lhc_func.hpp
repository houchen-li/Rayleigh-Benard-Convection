//
// Created by houchen_li on 5/1/19.
//

#pragma once
#ifndef LHC_FUNC_H
#define LHC_FUNC_H

#include "def.h"
#include "utils.h"
#include "state_type.h"

void func(const RBCSystem::StateType &f, RBCSystem::StateType &dfdt, const Real &t, const Real &Pr, const Real &Ra, const Real &a)
{
    const Real a_2 = a * a;
    Uint i, j, k, i_2, j_2;
    Real sum;

    for (j = 2; j < 2 * RBCSystem::Ny + 2; j += 2) {
        j_2 = j * j;
        dfdt.phi(j) = -(j_2 * PI_2 * f.phi(j));
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
        i_2 = i * i;
        for (j = 2 - i % 2; j < 2 * RBCSystem::Ny + 2 - i % 2; j += 2) {
            j_2 = j * j;
            dfdt.theta(i, j) = 2.0 * PI * a * i * f.psi(i, j) - (4 * a_2 * i_2 + j_2) * PI_2 * f.theta(i, j);
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
    for (i = 1; i < RBCSystem::Nx; i++) {
        i_2 = i * i;
        for (j = 2 - i % 2; j < 2 * RBCSystem::Ny + i % 2; j += 2) {
            j_2 = j * j;
        }
        for (j = 2 + i % 2; j < 2 * RBCSystem::Ny + 2 - i % 2; j += 2) {
            j_2 = j * j;
        }
    }
    for (i = 2; i < RBCSystem::Nx + 1; i++) {
        i_2 = i * i;
        for (j = 2 + i % 2; j < 2 * RBCSystem::Ny + 2 - i % 2; j += 2) {
            j_2 = j * j;
        }
        for (j = 2 - i % 2; j < 2 * RBCSystem::Ny + i % 2; j += 2) {
            j_2 = j * j;
        }
    }

    for (i = 1; i < RBCSystem::Nx + 1; i++) {
        i_2 = i * i;
        for (j = 2 - i % 2; j < 2 * RBCSystem::Ny + 2 - i % 2; j += 2) {
            j_2 = j * j;
            dfdt.psi(i, j) = 2.0 * a * i * Pr * Ra / ((4.0 * a_2 * i_2 + j_2) * PI) * f.theta(i, j) -
                             (4.0 * a_2 * i_2 + j_2) * PI_2 * Pr * f.psi(i, j);
        }
    }
    for (i = 1; i < RBCSystem::Nx; i++) {
        i_2 = i * i;
        for (j = 2 - i % 2; j < 2 * RBCSystem::Ny + i % 2; j += 2) {
            j_2 = j * j;
            dfdt.psi(i, j) -= PI_2 * a / 2.0 * f.psi(1, 1) * f.psi(i + 1, j + 1) *
                              (1 + (8.0 * a_2 * i + 2.0 * j) / (4.0 * a_2 * i_2 + j_2)) * (static_cast<Real>(i) - static_cast<Real>(j));
        }
        for (j = 2 + i % 2; j < 2 * RBCSystem::Ny + 2 - i % 2; j += 2) {
            j_2 = j * j;
            dfdt.psi(i, j) -= PI_2 * a / 2.0 * f.psi(1, 1) * f.psi(i + 1, j - 1) *
                              (1 + (8.0 * a_2 * i - 2.0 * j) / (4.0 * a_2 * i_2 + j_2)) * (i + j);
        }
    }
    for (i = 2; i < RBCSystem::Nx + 1; i++) {
        i_2 = i * i;
        for (j = 2 + i % 2; j < 2 * RBCSystem::Ny + 2 - i % 2; j += 2) {
            j_2 = j * j;
            dfdt.psi(i, j) += PI_2 * a / 2.0 * f.psi(1, 1) * f.psi(i - 1, j - 1) *
                              (1 - (8.0 * a_2 * i + 2.0 * j) / (4.0 * a_2 * i_2 + j_2)) * (static_cast<Real>(i) - static_cast<Real>(j));
        }
        for (j = 2 - i % 2; j < 2 * RBCSystem::Ny + i % 2; j += 2) {
            j_2 = j * j;
            dfdt.psi(i, j) += PI_2 * a / 2.0 * f.psi(1, 1) * f.psi(i - 1, j + 1) *
                              (1 - (8.0 * a_2 * i - 2.0 * j) / (4.0 * a_2 * i_2 + j_2)) * (i + j);
        }
    }

    return;
};

#endif /* ifndef SYMBOL */
