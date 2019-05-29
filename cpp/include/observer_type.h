//
// Created by houchen_li on 5/1/19.
//

#pragma once
#ifndef OBSERVER_TYPE_H
#define OBSERVER_TYPE_H

#include <fstream>
#include <H5Cpp.h>
#include "def.h"
#include "state_type.h"

//#define OBSERVE_STATES

namespace RBCSystem {

    class ObserverType {
    public:
        explicit ObserverType(Uint &l, Real &Pr, Real &Ra, Real &a);

        ~ObserverType(void);

        void startObservation(void);

        void endObservation(void);

        void operator()(const StateType &f, const Real &t);

    private:
        Uint &l;
        Real &Pr, &Ra, &a;
    };

}

#endif //OBSERVER_TYPE_H
