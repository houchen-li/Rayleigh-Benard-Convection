//
// Created by houchen_li on 5/1/19.
//

#pragma once
#ifndef UTILS_H
#define UTILS_H

#include <boost/math/constants/constants.hpp>
#include <H5Cpp.h>
#include "def.h"

constexpr Real PI = boost::math::constants::pi<Real>();

constexpr Real PI_2 = PI*PI;

constexpr Real PI_4 = PI*PI;

Real Power(const Real &val, const Uint &p);

void addAttr(H5::H5Object &h5obj, const String &attr_name, const Uint &val);

inline void addAttr(H5::H5Object &h5obj, const char *attr_name, const Uint &val) { addAttr(h5obj, static_cast<String>(attr_name), val); return; }

void addAttr(H5::H5Object &h5obj, const String &attr_name, const Real &val);

inline void addAttr(H5::H5Object &h5obj, const char *attr_name, const Real &val) { addAttr(h5obj, static_cast<String>(attr_name), val); return; }

void addArray(H5::Group &h5group, const String &array_name, const Real *a, const Uint &n);

inline void addArray(H5::Group &h5group, const char *array_name, const Real *a, const Uint &n) { addArray(h5group, static_cast<String>(array_name), a, n); return; }

#endif
