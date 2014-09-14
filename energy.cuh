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
__host__ __device__ static inline T Sqr2(T x) {
	return x * x * x * x;
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

#include "vars.cuh"
#include "return.cuh"

//		real f_0_2 = Sqr(f[0]);
//		real f_1_2 = Sqr(f[1]);
//		real f_2_2 = Sqr(f[2]);
//		real f_3_2 = Sqr(f[3]);
//		real f_4_2 = Sqr(f[4]);
//		real f_5_2 = Sqr(f[5]);
//		real f_6_2 = Sqr(f[6]);
//		real f_7_2 = Sqr(f[7]);
//		real f_8_2 = Sqr(f[8]);
//		real f_9_2 = Sqr(f[9]);
//		real f_10_2 = Sqr(f[10]);
//		real f_11_2 = Sqr(f[11]);
//		real f_12_2 = Sqr(f[12]);
//		real f_13_2 = Sqr(f[13]);
//		real f_14_2 = Sqr(f[14]);
//		real f_15_2 = Sqr(f[15]);
//		real f_16_2 = Sqr(f[16]);
//		real f_17_2 = Sqr(f[17]);
//		real f_18_2 = Sqr(f[18]);
//		real f_19_2 = Sqr(f[19]);
//		real f_20_2 = Sqr(f[20]);
//		real f_21_2 = Sqr(f[21]);
//		real f_22_2 = Sqr(f[22]);
//		real f_23_2 = Sqr(f[23]);
//		real f_24_2 = Sqr(f[24]);
//		real f_25_2 = Sqr(f[25]);
//		real f_26_2 = Sqr(f[26]);
//		real f_27_2 = Sqr(f[27]);
//		real f_28_2 = Sqr(f[28]);
//		real f_29_2 = Sqr(f[29]);
//		real f_30_2 = Sqr(f[30]);
//		real f_31_2 = Sqr(f[31]);
//		real f_32_2 = Sqr(f[32]);
//		real f_33_2 = Sqr(f[33]);
//		real f_34_2 = Sqr(f[34]);
//		real f_35_2 = Sqr(f[35]);
//		real f_36_2 = Sqr(f[36]);
//		real f_37_2 = Sqr(f[37]);
//		real f_38_2 = Sqr(f[38]);
//		real f_39_2 = Sqr(f[39]);
//		real f_40_2 = Sqr(f[40]);
//		real f_41_2 = Sqr(f[41]);
//		real f_42_2 = Sqr(f[42]);
//		real f_43_2 = Sqr(f[43]);
//		real f_44_2 = Sqr(f[44]);
//		real f_45_2 = Sqr(f[45]);
//		real f_46_2 = Sqr(f[46]);
//		real f_47_2 = Sqr(f[47]);
//		real f_48_2 = Sqr(f[48]);
//		real f_49_2 = Sqr(f[49]);
//		real f_50_2 = Sqr(f[50]);
//		real f_51_2 = Sqr(f[51]);
//		real f_52_2 = Sqr(f[52]);
//		real f_53_2 = Sqr(f[53]);
//		real f_54_2 = Sqr(f[54]);
//		real f_55_2 = Sqr(f[55]);
//		real f_56_2 = Sqr(f[56]);
//		real f_57_2 = Sqr(f[57]);
//		real f_58_2 = Sqr(f[58]);
//		real f_59_2 = Sqr(f[59]);
//		real f_60_2 = Sqr(f[60]);
//		real f_61_2 = Sqr(f[61]);
//		real f_62_2 = Sqr(f[62]);
//		real f_63_2 = Sqr(f[63]);
//		real f_64_2 = Sqr(f[64]);
//		real f_65_2 = Sqr(f[65]);
//		real f_66_2 = Sqr(f[66]);
//		real f_67_2 = Sqr(f[67]);
//		real f_68_2 = Sqr(f[68]);
//		real f_69_2 = Sqr(f[69]);
//		real f_70_2 = Sqr(f[70]);
//		real f_71_2 = Sqr(f[71]);
//		real f_72_2 = Sqr(f[72]);
//		real f_73_2 = Sqr(f[73]);
//		real f_74_2 = Sqr(f[74]);
//		real f_75_2 = Sqr(f[75]);
//		real f_76_2 = Sqr(f[76]);
//		real f_77_2 = Sqr(f[77]);
//		real f_78_2 = Sqr(f[78]);
//		real f_79_2 = Sqr(f[79]);
//
////		real f_0_4 = Sqr2(f[0]);
////		real f_1_4 = Sqr2(f[1]);
////		real f_2_4 = Sqr2(f[2]);
////		real f_3_4 = Sqr2(f[3]);
////		real f_4_4 = Sqr2(f[4]);
////		real f_5_4 = Sqr2(f[5]);
////		real f_6_4 = Sqr2(f[6]);
////		real f_7_4 = Sqr2(f[7]);
////		real f_8_4 = Sqr2(f[8]);
////		real f_9_4 = Sqr2(f[9]);
////		real f_10_4 = Sqr2(f[10]);
////		real f_11_4 = Sqr2(f[11]);
////		real f_12_4 = Sqr2(f[12]);
////		real f_13_4 = Sqr2(f[13]);
////		real f_14_4 = Sqr2(f[14]);
////		real f_15_4 = Sqr2(f[15]);
//
//		real J_0_2 = Sqr(J[0]);
//		real J_1_2 = Sqr(J[1]);
//		real J_2_2 = Sqr(J[2]);
//		real J_3_2 = Sqr(J[3]);
//		real J_4_2 = Sqr(J[4]);
//
//		real U_0_2 = Sqr(U[0]);
//		real U_1_2 = Sqr(U[1]);
//		real U_2_2 = Sqr(U[2]);
//		real U_3_2 = Sqr(U[3]);
//		real U_4_2 = Sqr(U[4]);
//
//
//#include "energyret.cuh"
////		return (75*f_0_2*f_10_2*J_0_2 - 1800*f_10_4*J_0_2 + 75*f_0_2*f_11_2*J_0_2 - 3600*f_10_2*f_11_2*J_0_2 - 1800*f_11_4*J_0_2 + 72*f_0_2*f_12_2*J_0_2 - 1050*f_10_2*f_12_2*J_0_2 - 1050*f_11_2*f_12_2*J_0_2 - 2520*f_12_4*J_0_2 + 72*f_0_2*f_13_2*J_0_2 - 1050*f_10_2*f_13_2*J_0_2 - 1050*f_11_2*f_13_2*J_0_2 - 5040*f_12_2*f_13_2*J_0_2 - 2520*f_13_4*J_0_2 + 70*f_0_2*f_14_2*J_0_2 + 5040*f_10_2*f_14_2*J_0_2 + 5040*f_11_2*f_14_2*J_0_2 + 70*f_0_2*f_15_2*J_0_2 + 5040*f_10_2*f_15_2*J_0_2 + 5040*f_11_2*f_15_2*J_0_2 + 75*f_10_2*f_1_2*J_0_2 + 75*f_11_2*f_1_2*J_0_2 + 72*f_12_2*f_1_2*J_0_2 + 72*f_13_2*f_1_2*J_0_2 + 70*f_14_2*f_1_2*J_0_2 + 70*f_15_2*f_1_2*J_0_2 + 128*f_10_2*f_2_2*J_0_2 + 128*f_11_2*f_2_2*J_0_2 + 110*f_12_2*f_2_2*J_0_2 + 110*f_13_2*f_2_2*J_0_2 + 168*f_14_2*f_2_2*J_0_2 + 168*f_15_2*f_2_2*J_0_2 - 120*f_2_4*J_0_2 + 128*f_10_2*f_3_2*J_0_2 + 128*f_11_2*f_3_2*J_0_2 + 110*f_12_2*f_3_2*J_0_2 + 110*f_13_2*f_3_2*J_0_2 + 168*f_14_2*f_3_2*J_0_2 + 168*f_15_2*f_3_2*J_0_2 - 240*f_2_2*f_3_2*J_0_2 - 120*f_3_4*J_0_2 + 240*f_0_2*f_4_2*J_0_2 + 270*f_10_2*f_4_2*J_0_2 + 270*f_11_2*f_4_2*J_0_2 + 192*f_12_2*f_4_2*J_0_2 + 192*f_13_2*f_4_2*J_0_2 + 315*f_14_2*f_4_2*J_0_2 + 315*f_15_2*f_4_2*J_0_2 + 240*f_1_2*f_4_2*J_0_2 - 90*f_2_2*f_4_2*J_0_2 - 90*f_3_2*f_4_2*J_0_2 - 360*f_4_4*J_0_2 + 240*f_0_2*f_5_2*J_0_2 + 270*f_10_2*f_5_2*J_0_2 + 270*f_11_2*f_5_2*J_0_2 + 192*f_12_2*f_5_2*J_0_2 + 192*f_13_2*f_5_2*J_0_2 + 315*f_14_2*f_5_2*J_0_2 + 315*f_15_2*f_5_2*J_0_2 + 240*f_1_2*f_5_2*J_0_2 - 90*f_2_2*f_5_2*J_0_2 - 90*f_3_2*f_5_2*J_0_2 - 720*f_4_2*f_5_2*J_0_2 - 360*f_5_4*J_0_2 + 90*f_0_2*f_6_2*J_0_2 + 2040*f_10_2*f_6_2*J_0_2 + 2040*f_11_2*f_6_2*J_0_2 + 405*f_12_2*f_6_2*J_0_2 + 405*f_13_2*f_6_2*J_0_2 + 560*f_14_2*f_6_2*J_0_2 + 560*f_15_2*f_6_2*J_0_2 + 90*f_1_2*f_6_2*J_0_2 + 640*f_2_2*f_6_2*J_0_2 + 640*f_3_2*f_6_2*J_0_2 - 240*f_4_2*f_6_2*J_0_2 - 240*f_5_2*f_6_2*J_0_2 - 720*f_6_4*J_0_2 + 90*f_0_2*f_7_2*J_0_2 + 2040*f_10_2*f_7_2*J_0_2 + 2040*f_11_2*f_7_2*J_0_2 + 405*f_12_2*f_7_2*J_0_2 + 405*f_13_2*f_7_2*J_0_2 + 560*f_14_2*f_7_2*J_0_2 + 560*f_15_2*f_7_2*J_0_2 + 90*f_1_2*f_7_2*J_0_2 + 640*f_2_2*f_7_2*J_0_2 + 640*f_3_2*f_7_2*J_0_2 - 240*f_4_2*f_7_2*J_0_2 - 240*f_5_2*f_7_2*J_0_2 - 1440*f_6_2*f_7_2*J_0_2 - 720*f_7_4*J_0_2 + 80*f_0_2*f_8_2*J_0_2 - 720*f_10_2*f_8_2*J_0_2 - 720*f_11_2*f_8_2*J_0_2 + 3040*f_12_2*f_8_2*J_0_2 + 3040*f_13_2*f_8_2*J_0_2 + 1050*f_14_2*f_8_2*J_0_2 + 1050*f_15_2*f_8_2*J_0_2 + 80*f_1_2*f_8_2*J_0_2 + 165*f_2_2*f_8_2*J_0_2 + 165*f_3_2*f_8_2*J_0_2 + 1240*f_4_2*f_8_2*J_0_2 + 1240*f_5_2*f_8_2*J_0_2 - 450*f_6_2*f_8_2*J_0_2 - 450*f_7_2*f_8_2*J_0_2 - 1200*f_8_4*J_0_2 + 80*f_0_2*f_9_2*J_0_2 - 720*f_10_2*f_9_2*J_0_2 - 720*f_11_2*f_9_2*J_0_2 + 3040*f_12_2*f_9_2*J_0_2 + 3040*f_13_2*f_9_2*J_0_2 + 1050*f_14_2*f_9_2*J_0_2 + 1050*f_15_2*f_9_2*J_0_2 + 80*f_1_2*f_9_2*J_0_2 + 165*f_2_2*f_9_2*J_0_2 + 165*f_3_2*f_9_2*J_0_2 + 1240*f_4_2*f_9_2*J_0_2 + 1240*f_5_2*f_9_2*J_0_2 - 450*f_6_2*f_9_2*J_0_2 - 450*f_7_2*f_9_2*J_0_2 - 2400*f_8_2*f_9_2*J_0_2 - 1200*f_9_4*J_0_2 + 60*(f_0_2 + f_10_2 + f_11_2 + f_12_2 + f_13_2 + f_14_2 + f_15_2 + f_1_2 + f_2_2 + f_3_2 + f_4_2 + f_5_2 + f_6_2 + f_7_2 + f_8_2 + f_9_2)*(10*f_10_2 + 10*f_11_2 + 3*(5*f_12_2 + 5*f_13_2 + 7*(f_14_2 + f_15_2)) + f_4_2 + f_5_2 + 3*f_6_2 + 3*f_7_2 + 6*f_8_2 + 6*f_9_2)*U_0_2 - 60*((f_0_2 + f_10_2 + f_11_2 + f_12_2 + f_13_2 + f_14_2 + f_15_2 + f_1_2 + f_2_2 + f_3_2 + f_4_2 + f_5_2 + f_6_2 + f_7_2 + f_8_2 + f_9_2)*(5*f_10_2 + 5*f_11_2 + 6*f_12_2 + 6*f_13_2 + 7*(f_14_2 + f_15_2) + f_2_2 + f_3_2 + 2*f_4_2 + 2*f_5_2 + 3*f_6_2 + 3*f_7_2 + 4*f_8_2 + 4*f_9_2)*mu + (6*f_10_2*f_12_2 + 6*f_11_2*f_12_2 + 6*f_10_2*f_13_2 + 6*f_11_2*f_13_2 + 7*f_12_2*f_14_2 + 7*f_13_2*f_14_2 + 7*(f_12_2 + f_13_2)*f_15_2 + f_0_2*(f_2_2 + f_3_2) + f_1_2*(f_2_2 + f_3_2) + 2*f_2_2*f_4_2 + 2*f_3_2*f_4_2 + 2*f_2_2*f_5_2 + 2*f_3_2*f_5_2 + 3*f_4_2*f_6_2 + 3*f_5_2*f_6_2 + 3*f_4_2*f_7_2 + 3*f_5_2*f_7_2 + 5*f_10_2*f_8_2 + 5*f_11_2*f_8_2 + 4*f_6_2*f_8_2 + 4*f_7_2*f_8_2 + 5*f_10_2*f_9_2 + 5*f_11_2*f_9_2 + 4*f_6_2*f_9_2 + 4*f_7_2*f_9_2)*J[0])*U[0])/(60.*Sqr(f_0_2 + f_10_2 + f_11_2 + f_12_2 + f_13_2 + f_14_2 + f_15_2 + f_1_2 + f_2_2 + f_3_2 + f_4_2 + f_5_2 + f_6_2 + f_7_2 + f_8_2 + f_9_2)*U[0]);
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

				if (n < nmax) {
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
				}
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
