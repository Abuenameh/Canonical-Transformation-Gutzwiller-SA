/*
 * energy.hpp
 *
 *  Created on: Sep 11, 2014
 *      Author: Abuenameh
 */

#ifndef ENERGY_HPP_
#define ENERGY_HPP_

#include "gutzwiller.hpp"

template<class T>
__host__ __device__ static inline T Sqr(T x) {
	return x * x;
}

template<class T>
class Energy {
public:
	__host__ __device__ T operator()(const T *f, unsigned int n,
			void *f_data) const {
		parameters* parms = (parameters*) f_data;
		real* U = parms->U;
		real* J = parms->J;
		real mu = parms->mu;
		real theta = parms->theta;

		return
//		return (-240*Sqr(f[2])*Sqr(f[3])*Sqr(J[1]) + 240*Sqr(f[0])*Sqr(f[4])*Sqr(J[1]) + 240*Sqr(f[1])*Sqr(f[4])*Sqr(J[1]) - 90*Sqr(f[2])*Sqr(f[4])*Sqr(J[1]) - 90*Sqr(f[3])*Sqr(f[4])*Sqr(J[1]) + 240*Sqr(f[0])*Sqr(f[5])*Sqr(J[1]) + 240*Sqr(f[1])*Sqr(f[5])*Sqr(J[1]) - 90*Sqr(f[2])*Sqr(f[5])*Sqr(J[1]) - 90*Sqr(f[3])*Sqr(f[5])*Sqr(J[1]) - 720*Sqr(f[4])*Sqr(f[5])*Sqr(J[1]) + 90*Sqr(f[0])*Sqr(f[6])*Sqr(J[1]) + 90*Sqr(f[1])*Sqr(f[6])*Sqr(J[1]) + 640*Sqr(f[2])*Sqr(f[6])*Sqr(J[1]) + 640*Sqr(f[3])*Sqr(f[6])*Sqr(J[1]) - 240*Sqr(f[4])*Sqr(f[6])*Sqr(J[1]) - 240*Sqr(f[5])*Sqr(f[6])*Sqr(J[1]) + 90*Sqr(f[0])*Sqr(f[7])*Sqr(J[1]) + 90*Sqr(f[1])*Sqr(f[7])*Sqr(J[1]) + 640*Sqr(f[2])*Sqr(f[7])*Sqr(J[1]) + 640*Sqr(f[3])*Sqr(f[7])*Sqr(J[1]) - 240*Sqr(f[4])*Sqr(f[7])*Sqr(J[1]) - 240*Sqr(f[5])*Sqr(f[7])*Sqr(J[1]) - 1440*Sqr(f[6])*Sqr(f[7])*Sqr(J[1]) + 80*Sqr(f[0])*Sqr(f[8])*Sqr(J[1]) + 80*Sqr(f[1])*Sqr(f[8])*Sqr(J[1]) + 165*Sqr(f[2])*Sqr(f[8])*Sqr(J[1]) + 165*Sqr(f[3])*Sqr(f[8])*Sqr(J[1]) + 1240*Sqr(f[4])*Sqr(f[8])*Sqr(J[1]) + 1240*Sqr(f[5])*Sqr(f[8])*Sqr(J[1]) - 450*Sqr(f[6])*Sqr(f[8])*Sqr(J[1]) - 450*Sqr(f[7])*Sqr(f[8])*Sqr(J[1]) + 80*Sqr(f[0])*Sqr(f[9])*Sqr(J[1]) + 80*Sqr(f[1])*Sqr(f[9])*Sqr(J[1]) + 165*Sqr(f[2])*Sqr(f[9])*Sqr(J[1]) + 165*Sqr(f[3])*Sqr(f[9])*Sqr(J[1]) + 1240*Sqr(f[4])*Sqr(f[9])*Sqr(J[1]) + 1240*Sqr(f[5])*Sqr(f[9])*Sqr(J[1]) - 450*Sqr(f[6])*Sqr(f[9])*Sqr(J[1]) - 450*Sqr(f[7])*Sqr(f[9])*Sqr(J[1]) - 2400*Sqr(f[8])*Sqr(f[9])*Sqr(J[1]) + 75*Sqr(f[0])*Sqr(f[10])*Sqr(J[1]) + 75*Sqr(f[1])*Sqr(f[10])*Sqr(J[1]) + 128*Sqr(f[2])*Sqr(f[10])*Sqr(J[1]) + 128*Sqr(f[3])*Sqr(f[10])*Sqr(J[1]) + 270*Sqr(f[4])*Sqr(f[10])*Sqr(J[1]) + 270*Sqr(f[5])*Sqr(f[10])*Sqr(J[1]) + 2040*Sqr(f[6])*Sqr(f[10])*Sqr(J[1]) + 2040*Sqr(f[7])*Sqr(f[10])*Sqr(J[1]) - 720*Sqr(f[8])*Sqr(f[10])*Sqr(J[1]) - 720*Sqr(f[9])*Sqr(f[10])*Sqr(J[1]) + 75*Sqr(f[0])*Sqr(f[11])*Sqr(J[1]) + 75*Sqr(f[1])*Sqr(f[11])*Sqr(J[1]) + 128*Sqr(f[2])*Sqr(f[11])*Sqr(J[1]) + 128*Sqr(f[3])*Sqr(f[11])*Sqr(J[1]) + 270*Sqr(f[4])*Sqr(f[11])*Sqr(J[1]) + 270*Sqr(f[5])*Sqr(f[11])*Sqr(J[1]) + 2040*Sqr(f[6])*Sqr(f[11])*Sqr(J[1]) + 2040*Sqr(f[7])*Sqr(f[11])*Sqr(J[1]) - 720*Sqr(f[8])*Sqr(f[11])*Sqr(J[1]) - 720*Sqr(f[9])*Sqr(f[11])*Sqr(J[1]) - 3600*Sqr(f[10])*Sqr(f[11])*Sqr(J[1]) + 72*Sqr(f[0])*Sqr(f[12])*Sqr(J[1]) + 72*Sqr(f[1])*Sqr(f[12])*Sqr(J[1]) + 110*Sqr(f[2])*Sqr(f[12])*Sqr(J[1]) + 110*Sqr(f[3])*Sqr(f[12])*Sqr(J[1]) + 192*Sqr(f[4])*Sqr(f[12])*Sqr(J[1]) + 192*Sqr(f[5])*Sqr(f[12])*Sqr(J[1]) + 405*Sqr(f[6])*Sqr(f[12])*Sqr(J[1]) + 405*Sqr(f[7])*Sqr(f[12])*Sqr(J[1]) + 3040*Sqr(f[8])*Sqr(f[12])*Sqr(J[1]) + 3040*Sqr(f[9])*Sqr(f[12])*Sqr(J[1]) - 1050*Sqr(f[10])*Sqr(f[12])*Sqr(J[1]) - 1050*Sqr(f[11])*Sqr(f[12])*Sqr(J[1]) + 72*Sqr(f[0])*Sqr(f[13])*Sqr(J[1]) + 72*Sqr(f[1])*Sqr(f[13])*Sqr(J[1]) + 110*Sqr(f[2])*Sqr(f[13])*Sqr(J[1]) + 110*Sqr(f[3])*Sqr(f[13])*Sqr(J[1]) + 192*Sqr(f[4])*Sqr(f[13])*Sqr(J[1]) + 192*Sqr(f[5])*Sqr(f[13])*Sqr(J[1]) + 405*Sqr(f[6])*Sqr(f[13])*Sqr(J[1]) + 405*Sqr(f[7])*Sqr(f[13])*Sqr(J[1]) + 3040*Sqr(f[8])*Sqr(f[13])*Sqr(J[1]) + 3040*Sqr(f[9])*Sqr(f[13])*Sqr(J[1]) - 1050*Sqr(f[10])*Sqr(f[13])*Sqr(J[1]) - 1050*Sqr(f[11])*Sqr(f[13])*Sqr(J[1]) - 5040*Sqr(f[12])*Sqr(f[13])*Sqr(J[1]) + 70*Sqr(f[0])*Sqr(f[14])*Sqr(J[1]) + 70*Sqr(f[1])*Sqr(f[14])*Sqr(J[1]) + 168*Sqr(f[2])*Sqr(f[14])*Sqr(J[1]) + 168*Sqr(f[3])*Sqr(f[14])*Sqr(J[1]) + 315*Sqr(f[4])*Sqr(f[14])*Sqr(J[1]) + 315*Sqr(f[5])*Sqr(f[14])*Sqr(J[1]) + 560*Sqr(f[6])*Sqr(f[14])*Sqr(J[1]) + 560*Sqr(f[7])*Sqr(f[14])*Sqr(J[1]) + 1050*Sqr(f[8])*Sqr(f[14])*Sqr(J[1]) + 1050*Sqr(f[9])*Sqr(f[14])*Sqr(J[1]) + 5040*Sqr(f[10])*Sqr(f[14])*Sqr(J[1]) + 5040*Sqr(f[11])*Sqr(f[14])*Sqr(J[1]) + 70*Sqr(f[0])*Sqr(f[15])*Sqr(J[1]) + 70*Sqr(f[1])*Sqr(f[15])*Sqr(J[1]) + 168*Sqr(f[2])*Sqr(f[15])*Sqr(J[1]) + 168*Sqr(f[3])*Sqr(f[15])*Sqr(J[1]) + 315*Sqr(f[4])*Sqr(f[15])*Sqr(J[1]) + 315*Sqr(f[5])*Sqr(f[15])*Sqr(J[1]) + 560*Sqr(f[6])*Sqr(f[15])*Sqr(J[1]) + 560*Sqr(f[7])*Sqr(f[15])*Sqr(J[1]) + 1050*Sqr(f[8])*Sqr(f[15])*Sqr(J[1]) + 1050*Sqr(f[9])*Sqr(f[15])*Sqr(J[1]) + 5040*Sqr(f[10])*Sqr(f[15])*Sqr(J[1]) + 5040*Sqr(f[11])*Sqr(f[15])*Sqr(J[1]) - 120*Sqr(J[1])*Sqr(Sqr(f[2])) - 120*Sqr(J[1])*Sqr(Sqr(f[3])) - 360*Sqr(J[1])*Sqr(Sqr(f[4])) - 360*Sqr(J[1])*Sqr(Sqr(f[5])) - 720*Sqr(J[1])*Sqr(Sqr(f[6])) - 720*Sqr(J[1])*Sqr(Sqr(f[7])) - 1200*Sqr(J[1])*Sqr(Sqr(f[8])) - 1200*Sqr(J[1])*Sqr(Sqr(f[9])) - 1800*Sqr(J[1])*Sqr(Sqr(f[10])) - 1800*Sqr(J[1])*Sqr(Sqr(f[11])) - 2520*Sqr(J[1])*Sqr(Sqr(f[12])) - 2520*Sqr(J[1])*Sqr(Sqr(f[13])) + 60*(Sqr(f[0]) + Sqr(f[1]) + Sqr(f[2]) + Sqr(f[3]) + Sqr(f[4]) + Sqr(f[5]) + Sqr(f[6]) + Sqr(f[7]) + Sqr(f[8]) + Sqr(f[9]) + Sqr(f[10]) + Sqr(f[11]) + Sqr(f[12]) + Sqr(f[13]) + Sqr(f[14]) + Sqr(f[15]))*(Sqr(f[4]) + Sqr(f[5]) + 3*Sqr(f[6]) + 3*Sqr(f[7]) + 6*Sqr(f[8]) + 6*Sqr(f[9]) + 10*Sqr(f[10]) + 10*Sqr(f[11]) + 3*(5*Sqr(f[12]) + 5*Sqr(f[13]) + 7*(Sqr(f[14]) + Sqr(f[15]))))*Sqr(U[1]) - 60*(J[1]*(Sqr(f[0])*(Sqr(f[2]) + Sqr(f[3])) + Sqr(f[1])*(Sqr(f[2]) + Sqr(f[3])) + 2*Sqr(f[2])*Sqr(f[4]) + 2*Sqr(f[3])*Sqr(f[4]) + 2*Sqr(f[2])*Sqr(f[5]) + 2*Sqr(f[3])*Sqr(f[5]) + 3*Sqr(f[4])*Sqr(f[6]) + 3*Sqr(f[5])*Sqr(f[6]) + 3*Sqr(f[4])*Sqr(f[7]) + 3*Sqr(f[5])*Sqr(f[7]) + 4*Sqr(f[6])*Sqr(f[8]) + 4*Sqr(f[7])*Sqr(f[8]) + 4*Sqr(f[6])*Sqr(f[9]) + 4*Sqr(f[7])*Sqr(f[9]) + 5*Sqr(f[8])*Sqr(f[10]) + 5*Sqr(f[9])*Sqr(f[10]) + 5*Sqr(f[8])*Sqr(f[11]) + 5*Sqr(f[9])*Sqr(f[11]) + 6*Sqr(f[10])*Sqr(f[12]) + 6*Sqr(f[11])*Sqr(f[12]) + 6*Sqr(f[10])*Sqr(f[13]) + 6*Sqr(f[11])*Sqr(f[13]) + 7*Sqr(f[12])*Sqr(f[14]) + 7*Sqr(f[13])*Sqr(f[14]) + 7*(Sqr(f[12]) + Sqr(f[13]))*Sqr(f[15])) + mu*(Sqr(f[0]) + Sqr(f[1]) + Sqr(f[2]) + Sqr(f[3]) + Sqr(f[4]) + Sqr(f[5]) + Sqr(f[6]) + Sqr(f[7]) + Sqr(f[8]) + Sqr(f[9]) + Sqr(f[10]) + Sqr(f[11]) + Sqr(f[12]) + Sqr(f[13]) + Sqr(f[14]) + Sqr(f[15]))*(Sqr(f[2]) + Sqr(f[3]) + 2*Sqr(f[4]) + 2*Sqr(f[5]) + 3*Sqr(f[6]) + 3*Sqr(f[7]) + 4*Sqr(f[8]) + 4*Sqr(f[9]) + 5*Sqr(f[10]) + 5*Sqr(f[11]) + 6*Sqr(f[12]) + 6*Sqr(f[13]) + 7*(Sqr(f[14]) + Sqr(f[15]))))*U[1])/(60.*Sqr(Sqr(f[0]) + Sqr(f[1]) + Sqr(f[2]) + Sqr(f[3]) + Sqr(f[4]) + Sqr(f[5]) + Sqr(f[6]) + Sqr(f[7]) + Sqr(f[8]) + Sqr(f[9]) + Sqr(f[10]) + Sqr(f[11]) + Sqr(f[12]) + Sqr(f[13]) + Sqr(f[14]) + Sqr(f[15]))*U[1]);
	}

	__host__ __device__ T operator()(const T *x, unsigned int n,
			void *f_data, int qwe) const {
		typedef typename complextype<T>::type complex_t;

		parameters* parms = (parameters*) f_data;
		real* U = parms->U;
		real* J = parms->J;
		real mu = parms->mu;
		real theta = parms->theta;

		complex_t expth = complex_t::make_cudacomplex(cos(theta), sin(theta)); //complextype<T>::make_complex(cos(theta), sin(theta));
		complex_t expmth = ~expth;
		complex_t exp2th = expth * expth;
		complex_t expm2th = ~exp2th;

		complex_t Ec = complex_t::zero();

		const complex_t* f[L];
		T norm[L];
		for (int i = 0; i < L; i++) {
			f[i] = reinterpret_cast<const complex_t*>(&x[2 * i * dim]);
			norm[i] = 0;
			for (int n = 0; n <= nmax; n++) {
				norm[i] += f[i][n].norm();
			}
		}

		for (int i = 0; i < L; i++) {

			int k1 = mod(i - 2);
			int j1 = mod(i - 1);
			int j2 = mod(i + 1);
			int k2 = mod(i + 2);

//			complex_t E0 = complex_t::zero();
//			complex_t E1j1 = complex_t::zero();
//			complex_t E1j2 = complex_t::zero();
//			complex_t E2j1 = complex_t::zero();
//			complex_t E2j2 = complex_t::zero();
//			complex_t E3j1 = complex_t::zero();
//			complex_t E3j2 = complex_t::zero();
//			complex_t E4j1j2 = complex_t::zero();
//			complex_t E4j1k1 = complex_t::zero();
//			complex_t E4j2k2 = complex_t::zero();
//			complex_t E5j1j2 = complex_t::zero();
//			complex_t E5j1k1 = complex_t::zero();
//			complex_t E5j2k2 = complex_t::zero();

			complex_t Ei = complex_t::zero();
			complex_t Ej1 = complex_t::zero();
			complex_t Ej2 = complex_t::zero();
			complex_t Ej1j2 = complex_t::zero();
			complex_t Ej1k1 = complex_t::zero();
			complex_t Ej2k2 = complex_t::zero();

			for (int n = 0; n <= nmax; n++) {
				Ei += (0.5 * U[i] * n * (n - 1) - mu * n) * ~f[i][n] * f[i][n];

				/*if (n < nmax) {
					Ej1 += -J[j1] * expth * g(n, n + 1) * ~f[i][n + 1]
							* ~f[j1][n] * f[i][n] * f[j1][n + 1];
					Ej2 += -J[i] * expmth * g(n, n + 1) * ~f[i][n + 1]
							* ~f[j2][n] * f[i][n] * f[j2][n + 1];

					if (n > 0) {
						Ej1 += 0.5 * J[j1] * J[j1] * exp2th * g(n, n)
								* g(n - 1, n + 1) * ~f[i][n + 1] * ~f[j1][n - 1]
								* f[i][n - 1] * f[j1][n + 1]
								* (1 / eps(U, i, j1, n, n)
										- 1 / eps(U, i, j1, n - 1, n + 1));
						Ej2 += 0.5 * J[i] * J[i] * expm2th * g(n, n)
								* g(n - 1, n + 1) * ~f[i][n + 1] * ~f[j2][n - 1]
								* f[i][n - 1] * f[j2][n + 1]
								* (1 / eps(U, i, j2, n, n)
										- 1 / eps(U, i, j2, n - 1, n + 1));
					}

					for (int m = 1; m <= nmax; m++) {
						if (n != m - 1) {
							Ej1 += 0.5 * (J[j1] * J[j1] / eps(U, i, j1, n, m))
									* g(n, m) * g(m - 1, n + 1)
									* (~f[i][n + 1] * ~f[j1][m - 1]
											* f[i][n + 1] * f[j1][m - 1]
											- ~f[i][n] * ~f[j1][m] * f[i][n]
													* f[j1][m]);
							Ej2 += 0.5 * (J[i] * J[i] / eps(U, i, j2, n, m))
									* g(n, m) * g(m - 1, n + 1)
									* (~f[i][n + 1] * ~f[j2][m - 1]
											* f[i][n + 1] * f[j2][m - 1]
											- ~f[i][n] * ~f[j2][m] * f[i][n]
													* f[j2][m]);
						}
					}

					if (n > 0) {
						Ej1j2 += 0.5 * (J[j1] * J[i] / eps(U, i, j1, n, n))
								* g(n, n) * g(n - 1, n + 1) * ~f[i][n + 1]
								* ~f[j1][n - 1] * ~f[j2][n] * f[i][n - 1]
								* f[j1][n] * f[j2][n + 1];
						Ej1j2 += 0.5 * (J[i] * J[j1] / eps(U, i, j2, n, n))
								* g(n, n) * g(n - 1, n + 1) * ~f[i][n + 1]
								* ~f[j2][n - 1] * ~f[j1][n] * f[i][n - 1]
								* f[j2][n] * f[j1][n + 1];
						Ej1k1 += 0.5 * (J[j1] * J[k1] / eps(U, i, j1, n, n))
								* g(n, n) * g(n - 1, n + 1) * ~f[i][n + 1]
								* ~f[j1][n - 1] * ~f[k1][n] * f[i][n]
								* f[j1][n + 1] * f[k1][n - 1];
						Ej2k2 += 0.5 * (J[i] * J[j2] / eps(U, i, j2, n, n))
								* g(n, n) * g(n - 1, n + 1) * ~f[i][n + 1]
								* ~f[j2][n - 1] * ~f[k2][n] * f[i][n]
								* f[j2][n + 1] * f[k2][n - 1];
						Ej1j2 -= 0.5
								* (J[j1] * J[i] / eps(U, i, j1, n - 1, n + 1))
								* g(n, n) * g(n - 1, n + 1) * ~f[i][n + 1]
								* ~f[j1][n] * ~f[j2][n - 1] * f[i][n - 1]
								* f[j1][n + 1] * f[j2][n];
						Ej1j2 -= 0.5
								* (J[i] * J[j1] / eps(U, i, j2, n - 1, n + 1))
								* g(n, n) * g(n - 1, n + 1) * ~f[i][n + 1]
								* ~f[j2][n] * ~f[j1][n - 1] * f[i][n - 1]
								* f[j2][n + 1] * f[j1][n];
						Ej1k1 -= 0.5
								* (J[j1] * J[k1] / eps(U, i, j1, n - 1, n + 1))
								* g(n, n) * g(n - 1, n + 1) * ~f[i][n]
								* ~f[j1][n - 1] * ~f[k1][n + 1] * f[i][n - 1]
								* f[j1][n + 1] * f[k1][n];
						Ej2k2 -= 0.5
								* (J[i] * J[j2] / eps(U, i, j2, n - 1, n + 1))
								* g(n, n) * g(n - 1, n + 1) * ~f[i][n]
								* ~f[j2][n - 1] * ~f[k2][n + 1] * f[i][n - 1]
								* f[j2][n + 1] * f[k2][n];
					}

					for (int m = 1; m <= nmax; m++) {
						if (n != m - 1 && n < nmax) {
							Ej1j2 += 0.5
									* (J[j1] * J[i] * exp2th
											/ eps(U, i, j1, n, m)) * g(n, m)
									* g(m - 1, n + 1) * ~f[i][n + 1]
									* ~f[j1][m - 1] * ~f[j2][m] * f[i][n + 1]
									* f[j1][m] * f[j2][m - 1];
							Ej1j2 += 0.5
									* (J[i] * J[j1] * expm2th
											/ eps(U, i, j2, n, m)) * g(n, m)
									* g(m - 1, n + 1) * ~f[i][n + 1]
									* ~f[j2][m - 1] * ~f[j1][m] * f[i][n + 1]
									* f[j2][m] * f[j1][m - 1];

							Ej1k1 += 0.5
									* (J[j1] * J[k1] * exp2th
											/ eps(U, i, j1, n, m)) * g(n, m)
									* g(m - 1, n + 1) * ~f[i][n + 1]
									* ~f[j1][m - 1] * ~f[k1][n] * f[i][n]
									* f[j1][m - 1] * f[k1][n + 1];
							Ej2k2 += 0.5
									* (J[i] * J[j2] * expm2th
											/ eps(U, i, j2, n, m)) * g(n, m)
									* g(m - 1, n + 1) * ~f[i][n + 1]
									* ~f[j2][m - 1] * ~f[k2][n] * f[i][n]
									* f[j2][m - 1] * f[k2][n + 1];

							Ej1j2 -= 0.5
									* (J[j1] * J[i] * exp2th
											/ eps(U, i, j1, n, m)) * g(n, m)
									* g(m - 1, n + 1) * ~f[i][n] * ~f[j1][m - 1]
									* ~f[j2][m] * f[i][n] * f[j1][m]
									* f[j2][m - 1];
							Ej1j2 -= 0.5
									* (J[i] * J[j1] * expm2th
											/ eps(U, i, j2, n, m)) * g(n, m)
									* g(m - 1, n + 1) * ~f[i][n] * ~f[j2][m - 1]
									* ~f[j1][m] * f[i][n] * f[j2][m]
									* f[j1][m - 1];

							Ej1k1 -= 0.5
									* (J[j1] * J[k1] * exp2th
											/ eps(U, i, j1, n, m)) * g(n, m)
									* g(m - 1, n + 1) * ~f[i][n + 1] * ~f[j1][m]
									* ~f[k1][n] * f[i][n] * f[j1][m]
									* f[k1][n + 1];
							Ej2k2 -= 0.5
									* (J[i] * J[j2] * expm2th
											/ eps(U, i, j2, n, m)) * g(n, m)
									* g(m - 1, n + 1) * ~f[i][n + 1] * ~f[j2][m]
									* ~f[k2][n] * f[i][n] * f[j2][m]
									* f[k2][n + 1];
						}
					}
				}*/
			}

			Ec += Ei / norm[i];
			Ec += Ej1 / (norm[i] * norm[j1]);
			Ec += Ej2 / (norm[i] * norm[j2]);
			Ec += Ej1j2 / (norm[i] * norm[j1] * norm[j2]);
			Ec += Ej1k1 / (norm[i] * norm[j1] * norm[k1]);
			Ec += Ej2k2 / (norm[i] * norm[j2] * norm[k2]);

//			Ec += E0 / norm[i];
//
//			Ec += E1j1 / (norm[i] * norm[j1]);
//			Ec += E1j2 / (norm[i] * norm[j2]);
//
//			Ec += E2j1 / (norm[i] * norm[j1]);
//			Ec += E2j2 / (norm[i] * norm[j2]);
//
//			Ec += E3j1 / (norm[i] * norm[j1]);
//			Ec += E3j2 / (norm[i] * norm[j2]);
//
//			Ec += E4j1j2 / (norm[i] * norm[j1] * norm[j2]);
//			Ec += E4j1k1 / (norm[i] * norm[j1] * norm[k1]);
//			Ec += E4j2k2 / (norm[i] * norm[j2] * norm[k2]);
//
//			Ec += E5j1j2 / (norm[i] * norm[j1] * norm[j2]);
//			Ec += E5j1k1 / (norm[i] * norm[j1] * norm[k1]);
//			Ec += E5j2k2 / (norm[i] * norm[j2] * norm[k2]);*/
		}

		return Ec.real();
	}
};

double f_nelderMead(unsigned int n, const double *x, double *grad,
		void *f_data) {
	return Energy<double>()(x, n, f_data);
}

#endif /* ENERGY_HPP_ */
