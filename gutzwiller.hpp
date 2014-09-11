/*
 * gutzwiller.hpp
 *
 *  Created on: Sep 11, 2014
 *      Author: Abuenameh
 */

#ifndef GUTZWILLER_HPP_
#define GUTZWILLER_HPP_

#include "complex.hpp"
#include "configuration.h"

#define L 1
#define nmax 2
#define dim (nmax+1)

__host__ __device__ inline int mod(int i) {
	return (i + L) % L;
}

__host__ __device__ inline real g(int n, int m) {
	return sqrt(1.0*(n + 1) * m);
}

__host__ __device__ inline double eps(real* U, int i, int j, int n, int m) {
	return n * U[i] - (m - 1) * U[j];
}

struct parameters {
	real* U;
	real mu;
	real* J;
    real theta;
};




#endif /* GUTZWILLER_HPP_ */
