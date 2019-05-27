//
// Created by houchen_li on 5/1/19.
//

#pragma once
#ifndef HERRING_FUNC_H
#define HERRING_FUNC_H

#include "def.h"
#include "utils.h"
#include "state_type.h"

void func(const RBCSystem::StateType &f, RBCSystem::StateType &dfdt, const Real &t, const Real &Pr, const Real &Ra, const Real &a)
{
    Uint i, j, k;
    Real laplacian, sum;

    for (j = 2; j < 2*RBCSystem::Ny+2; j+=2) {
        dfdt.phi(j) = -(j * j * PI_2 * f.phi(j));
    }
    for (i = 1; i < RBCSystem::Nx+1; i++) {
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
    for (i = 1; i < RBCSystem::Nx+1; i++) {
        for (j = 2 - i % 2; j < 2 * RBCSystem::Ny + 2 - i % 2; j += 2) {
            laplacian = (4 * a * a * i * i + j * j) * PI_2;
            dfdt.theta(i, j) = 2.0 * PI * a * i * f.psi(i, j) - laplacian * f.theta(i, j);
            sum = 0.0;
            for (k = 2; k < 2 * RBCSystem::Ny + 2 - i % 2 - j; k += 2)
                sum += k * f.phi(k) * f.psi(i, j + k);
            for (k = 2; k < j + i % 2; k += 2)
                sum += k * f.phi(k) * f.psi(i, j - k);
            for (k = j + 2 - i % 2; k < 2 * RBCSystem::Ny + 2; k += 2)
                sum -= k * f.phi(k) * f.psi(i, k - j);
            dfdt.theta(i, j) -= PI_2 * a * i * sum;
            dfdt.psi(i, j) = 2.0 * PI * a * i * Pr * Ra / laplacian * f.theta(i, j) - laplacian * Pr * f.psi(i, j);
        }
    }
    for (j = 3; j < 2*RBCSystem::Ny + 1; j += 2) {
        dfdt.psi(1, j) =
                PI_4 * a / 2.0 * f.psi(1, 1) * (-(16.0 * a * a + (j + 1) * (j + 1)) * (j - 1) * f.psi(2, j + 1) -
                                                (16.0 * a * a + (j - 1) * (j - 1)) * (j + 1) * f.psi(2, j - 1) -
                                                (4.0 * a * a + 1) *
                                                ((j - 1) * f.psi(2, j + 1) + (j + 1) * f.psi(2, j - 1)));
    }

    return;
};

#endif /* ifndef SYMBOL */
