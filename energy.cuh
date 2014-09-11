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
class Energy {
public:
	__host__ __device__ T operator()(const T *x, unsigned int n,
		void *f_data) const {
		typedef typename complextype<T>::type complex_t;

		parameters* parms = (parameters*)f_data;
		real* U = parms->U;
		real* J = parms->J;
		real mu = parms->mu;
		real theta = parms->theta;

		complex_t expth = complex_t::make_cudacomplex(cos(theta), sin(theta));//complextype<T>::make_complex(cos(theta), sin(theta));
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

			complex_t E0 = complex_t::zero();
			complex_t E1j1 = complex_t::zero();
			complex_t E1j2 = complex_t::zero();
			complex_t E2j1 = complex_t::zero();
			complex_t E2j2 = complex_t::zero();
			complex_t E3j1 = complex_t::zero();
			complex_t E3j2 = complex_t::zero();
			complex_t E4j1j2 = complex_t::zero();
			complex_t E4j1k1 = complex_t::zero();
			complex_t E4j2k2 = complex_t::zero();
			complex_t E5j1j2 = complex_t::zero();
			complex_t E5j1k1 = complex_t::zero();
			complex_t E5j2k2 = complex_t::zero();

			for (int n = 0; n <= nmax; n++) {
				E0 += (0.5 * U[i] * n * (n - 1) - mu * n) * ~f[i][n]
					* f[i][n];

				if (n < nmax) {
					E1j1 += -J[j1] * expth * g(n, n + 1) * ~f[i][n + 1]
						* ~f[j1][n] * f[i][n] * f[j1][n + 1];
					E1j2 += -J[i] * expmth * g(n, n + 1) * ~f[i][n + 1]
						* ~f[j2][n] * f[i][n] * f[j2][n + 1];

					if (n > 0) {
						E2j1 += 0.5 * J[j1] * J[j1] * exp2th * g(n, n)
							* g(n - 1, n + 1) * ~f[i][n + 1] * ~f[j1][n - 1]
							* f[i][n - 1] * f[j1][n + 1]
							* (1 / eps(U, i, j1, n, n)
								- 1 / eps(U, i, j1, n - 1, n + 1));
						E2j2 += 0.5 * J[i] * J[i] * expm2th * g(n, n)
							* g(n - 1, n + 1) * ~f[i][n + 1] * ~f[j2][n - 1]
							* f[i][n - 1] * f[j2][n + 1]
							* (1 / eps(U, i, j2, n, n)
								- 1 / eps(U, i, j2, n - 1, n + 1));
					}

					for (int m = 1; m <= nmax; m++) {
						if (n != m - 1) {
							E3j1 +=
								0.5 * (J[j1] * J[j1] / eps(U, i, j1, n, m))
									* g(n, m) * g(m - 1, n + 1)
									* (~f[i][n + 1] * ~f[j1][m - 1]
										* f[i][n + 1] * f[j1][m - 1]
										- ~f[i][n] * ~f[j1][m] * f[i][n]
											* f[j1][m]);
							E3j2 +=
								0.5 * (J[i] * J[i] / eps(U, i, j2, n, m))
									* g(n, m) * g(m - 1, n + 1)
									* (~f[i][n + 1] * ~f[j2][m - 1]
										* f[i][n + 1] * f[j2][m - 1]
										- ~f[i][n] * ~f[j2][m] * f[i][n]
											* f[j2][m]);
						}
					}

					if (n > 0) {
						E4j1j2 += 0.5 * (J[j1] * J[i] / eps(U, i, j1, n, n))
							* g(n, n) * g(n - 1, n + 1) * ~f[i][n + 1]
							* ~f[j1][n - 1] * ~f[j2][n] * f[i][n - 1] * f[j1][n]
							* f[j2][n + 1];
						E4j1j2 += 0.5 * (J[i] * J[j1] / eps(U, i, j2, n, n))
							* g(n, n) * g(n - 1, n + 1) * ~f[i][n + 1]
							* ~f[j2][n - 1] * ~f[j1][n] * f[i][n - 1] * f[j2][n]
							* f[j1][n + 1];
						E4j1k1 += 0.5 * (J[j1] * J[k1] / eps(U, i, j1, n, n))
							* g(n, n) * g(n - 1, n + 1) * ~f[i][n + 1]
							* ~f[j1][n - 1] * ~f[k1][n] * f[i][n] * f[j1][n + 1]
							* f[k1][n - 1];
						E4j2k2 += 0.5 * (J[i] * J[j2] / eps(U, i, j2, n, n))
							* g(n, n) * g(n - 1, n + 1) * ~f[i][n + 1]
							* ~f[j2][n - 1] * ~f[k2][n] * f[i][n] * f[j2][n + 1]
							* f[k2][n - 1];
						E4j1j2 -= 0.5
							* (J[j1] * J[i] / eps(U, i, j1, n - 1, n + 1))
							* g(n, n) * g(n - 1, n + 1) * ~f[i][n + 1]
							* ~f[j1][n] * ~f[j2][n - 1] * f[i][n - 1]
							* f[j1][n + 1] * f[j2][n];
						E4j1j2 -= 0.5
							* (J[i] * J[j1] / eps(U, i, j2, n - 1, n + 1))
							* g(n, n) * g(n - 1, n + 1) * ~f[i][n + 1]
							* ~f[j2][n] * ~f[j1][n - 1] * f[i][n - 1]
							* f[j2][n + 1] * f[j1][n];
						E4j1k1 -= 0.5
							* (J[j1] * J[k1] / eps(U, i, j1, n - 1, n + 1))
							* g(n, n) * g(n - 1, n + 1) * ~f[i][n]
							* ~f[j1][n - 1] * ~f[k1][n + 1] * f[i][n - 1]
							* f[j1][n + 1] * f[k1][n];
						E4j2k2 -= 0.5
							* (J[i] * J[j2] / eps(U, i, j2, n - 1, n + 1))
							* g(n, n) * g(n - 1, n + 1) * ~f[i][n]
							* ~f[j2][n - 1] * ~f[k2][n + 1] * f[i][n - 1]
							* f[j2][n + 1] * f[k2][n];
					}

					for (int m = 1; m <= nmax; m++) {
						if (n != m - 1 && n < nmax) {
							E5j1j2 +=0.5
								* (J[j1] * J[i] * exp2th / eps(U, i, j1, n, m))
								* g(n, m) * g(m - 1, n + 1) * ~f[i][n + 1]
								* ~f[j1][m - 1] * ~f[j2][m] * f[i][n + 1]
								* f[j1][m] * f[j2][m - 1];
							E5j1j2 +=0.5
								* (J[i] * J[j1] * expm2th / eps(U, i, j2, n, m))
								* g(n, m) * g(m - 1, n + 1) * ~f[i][n + 1]
								* ~f[j2][m - 1] * ~f[j1][m] * f[i][n + 1]
								* f[j2][m] * f[j1][m - 1];

							E5j1k1 +=0.5
								* (J[j1] * J[k1] * exp2th / eps(U, i, j1, n, m))
								* g(n, m) * g(m - 1, n + 1) * ~f[i][n + 1]
								* ~f[j1][m - 1] * ~f[k1][n] * f[i][n]
								* f[j1][m - 1] * f[k1][n + 1];
							E5j2k2 +=0.5
								* (J[i] * J[j2] * expm2th / eps(U, i, j2, n, m))
								* g(n, m) * g(m - 1, n + 1) * ~f[i][n + 1]
								* ~f[j2][m - 1] * ~f[k2][n] * f[i][n]
								* f[j2][m - 1] * f[k2][n + 1];

							E5j1j2 -=0.5
								* (J[j1] * J[i] * exp2th / eps(U, i, j1, n, m))
								* g(n, m) * g(m - 1, n + 1) * ~f[i][n]
								* ~f[j1][m - 1] * ~f[j2][m] * f[i][n] * f[j1][m]
								* f[j2][m - 1];
							E5j1j2 -=0.5
								* (J[i] * J[j1] * expm2th / eps(U, i, j2, n, m))
								* g(n, m) * g(m - 1, n + 1) * ~f[i][n]
								* ~f[j2][m - 1] * ~f[j1][m] * f[i][n] * f[j2][m]
								* f[j1][m - 1];

							E5j1k1 -= 0.5
								* (J[j1] * J[k1] * exp2th / eps(U, i, j1, n, m))
								* g(n, m) * g(m - 1, n + 1) * ~f[i][n + 1]
								* ~f[j1][m] * ~f[k1][n] * f[i][n] * f[j1][m]
								* f[k1][n + 1];
							E5j2k2 -= 0.5
								* (J[i] * J[j2] * expm2th / eps(U, i, j2, n, m))
								* g(n, m) * g(m - 1, n + 1) * ~f[i][n + 1]
								* ~f[j2][m] * ~f[k2][n] * f[i][n] * f[j2][m]
								* f[k2][n + 1];
						}
					}
				}
			}

			Ec += E0 / norm[i];

			Ec += E1j1 / (norm[i] * norm[j1]);
			Ec += E1j2 / (norm[i] * norm[j2]);

			Ec += E2j1 / (norm[i] * norm[j1]);
			Ec += E2j2 / (norm[i] * norm[j2]);

			Ec += E3j1 / (norm[i] * norm[j1]);
			Ec += E3j2 / (norm[i] * norm[j2]);

			Ec += E4j1j2 / (norm[i] * norm[j1] * norm[j2]);
			Ec += E4j1k1 / (norm[i] * norm[j1] * norm[k1]);
			Ec += E4j2k2 / (norm[i] * norm[j2] * norm[k2]);

			Ec += E5j1j2 / (norm[i] * norm[j1] * norm[j2]);
			Ec += E5j1k1 / (norm[i] * norm[j1] * norm[k1]);
			Ec += E5j2k2 / (norm[i] * norm[j2] * norm[k2]);
		}

		return Ec.real();
	}
};

#endif /* ENERGY_HPP_ */
