#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//Author: Ketson R. M. dos Santos
//Institution: Columbia University
//e-mail: kmd2191@columbia.edu
//
//This is the GPU version of the program
//to solve the WPI problem in stochastic dynamics

double * rk4(double bc0, double bc1, double bc2, double bc3, double *resp, int pi, int pf, int Nx) {

	double m = 1;
	double c = 0.5;
	double k = 0.5;
	double e1 = 0.1;
	double e2 = 0.1;
	static double b[8];

	double xi = 0.0;
	double xf = 1.0;
	double h, t0, exp10, exp20;
	double y0, y1, y2, y3, y4, y5, y6, y7;
	double yaux0, yaux1, yaux2, yaux3, yaux4, yaux5, yaux6, yaux7;
	double k10, k11, k12, k13, k14, k15, k16, k17;
	double k20, k21, k22, k23, k24, k25, k26, k27;
	double k30, k31, k32, k33, k34, k35, k36, k37;
	double k40, k41, k42, k43, k44, k45, k46, k47;
	int N0 = Nx, cont;

	h = (xf - xi) / (N0 - 1);
	t0 = xi;
	y0 = 0; y1 = 0; y2 = bc0; y3 = bc1; y4 = 0; y5 = 0; y6 = bc2; y7 = bc3;
	
	resp[0 + pi] = y0;
	resp[1 + pi] = y1;
	resp[2 + pi] = y2;
	resp[3 + pi] = y3;
	resp[4 + pi] = y4;
	resp[5 + pi] = y5;
	resp[6 + pi] = y6;
	resp[7 + pi] = y7;

	cont = 8;
	for (int i = 1;i < N0;i++) {
		b[0] = y1;
		b[1] = y2;
		b[2] = y3;
		exp10 = pow(m, -2)*((-3)*e1*e1*k*k*pow(y0, 5) + 4 * c*e2*k*pow(y1, 3) + (-3)*c*e2*k*y1*y1*y5 + \
			5 * c*c*y2 + (-4)*k*m*y2 + 24 * c*c*e2*y1*y1*y2 + 15 * c*c*e2*e2*pow(y1, 4)*y2 + \
			(-6)*c*c*e2*y1*y5*y2 + 3 * e1*k*y0*y0*(k*y4 + 2 * c*e2*y1*y1*y1 + c * y5 + \
			(-2)*m*y2) + 2 * k*y4*(2 * k + (-3)*c*e2*y1*y2) + 2 * e1*k*y0*y0*y0*((-4)*k + 3 * c*e2*y1*y2)\
			+ k * y0*((-5)*k + (-6)*e1*m*y1*y1 + 12 * c*e2*y1*y2)\
			+ (-4)*c*c*y6 + 2 * k*m*y6 + (-3)*c*c*e2*y1*y1*y6);

		b[3] = exp10;
		b[4] = y5;
		b[5] = y6;
		b[6] = y7;
		exp20 = pow(m, -2)*(4 * k*k*y0 + e1 * k*k*y0*y0*y0 + (-5)*k*k*y4 + (-3)*c*e1*k*y0*y0*y1 + \
			c*e2*k*y1*y1*y1 + (-4)*c*c*y2 + 2 * k*m*y2 + (-3)*c*c*e2*y1*y1*y2 + 5 * c*c*y6 + (-4)*k*m*y6);
		b[7] = exp20;

		k10 = b[0] * h; k11 = b[1] * h; k12 = b[2] * h; k13 = b[3] * h;
		k14 = b[4] * h; k15 = b[5] * h; k16 = b[6] * h; k17 = b[7] * h;
		
		yaux0 = y0 + 0.5*k10;
		yaux1 = y1 + 0.5*k11;
		yaux2 = y2 + 0.5*k12;
		yaux3 = y3 + 0.5*k13;
		yaux4 = y4 + 0.5*k14;
		yaux5 = y5 + 0.5*k15;
		yaux6 = y6 + 0.5*k16;
		yaux7 = y7 + 0.5*k17;

		//b = OdeFun(t0+0.5*h,yaux0,yaux1,yaux2,yaux3); //k2
		b[0] = yaux1;
		b[1] = yaux2;
		b[2] = yaux3;
		exp10 = pow(m, -2)*((-3)*e1*e1*k*k*pow(yaux0, 5) + 4 * c*e2*k*pow(yaux1, 3) + (-3)*c*e2*k*yaux1*yaux1*yaux5 + \
			5 * c*c*yaux2 + (-4)*k*m*yaux2 + 24 * c*c*e2*yaux1*yaux1*yaux2 + 15 * c*c*e2*e2*pow(yaux1, 4)*yaux2 + \
			(-6)*c*c*e2*yaux1*yaux5*yaux2 + 3 * e1*k*yaux0*yaux0*(k*yaux4 + 2 * c*e2*yaux1*yaux1*yaux1 + c * yaux5 + \
			(-2)*m*yaux2) + 2 * k*yaux4*(2 * k + (-3)*c*e2*yaux1*yaux2) + 2 * e1*k*yaux0*yaux0*yaux0*((-4)*k + 3 * c*e2*yaux1*yaux2)\
			+ k * yaux0*((-5)*k + (-6)*e1*m*yaux1*yaux1 + 12 * c*e2*yaux1*yaux2)\
			+ (-4)*c*c*yaux6 + 2 * k*m*yaux6 + (-3)*c*c*e2*yaux1*yaux1*yaux6);

		b[3] = exp10;
		b[4] = yaux5;
		b[5] = yaux6;
		b[6] = yaux7;
		exp20 = pow(m, -2)*(4 * k*k*yaux0 + e1 * k*k*yaux0*yaux0*yaux0 + (-5)*k*k*yaux4 + (-3)*c*e1*k*yaux0*yaux0*yaux1 + \
			c*e2*k*yaux1*yaux1*yaux1 + (-4)*c*c*yaux2 + 2 * k*m*yaux2 + (-3)*c*c*e2*yaux1*yaux1*yaux2 + 5 * c*c*yaux6 + (-4)*k*m*yaux6);
		b[7] = exp20;

		k20 = b[0] * h; k21 = b[1] * h; k22 = b[2] * h; k23 = b[3] * h;
		k24 = b[4] * h; k25 = b[5] * h; k26 = b[6] * h; k27 = b[7] * h;

		yaux0 = y0 + 0.5*k20;
		yaux1 = y1 + 0.5*k21;
		yaux2 = y2 + 0.5*k22;
		yaux3 = y3 + 0.5*k23;
		yaux4 = y4 + 0.5*k24;
		yaux5 = y5 + 0.5*k25;
		yaux6 = y6 + 0.5*k26;
		yaux7 = y7 + 0.5*k27;

		//b = OdeFun(t0+0.5*h,yaux0,yaux1,yaux2,yaux3); //k3
		b[0] = yaux1;
		b[1] = yaux2;
		b[2] = yaux3;
		exp10 = pow(m, -2)*((-3)*e1*e1*k*k*pow(yaux0, 5) + 4 * c*e2*k*pow(yaux1, 3) + (-3)*c*e2*k*yaux1*yaux1*yaux5 + \
			5 * c*c*yaux2 + (-4)*k*m*yaux2 + 24 * c*c*e2*yaux1*yaux1*yaux2 + 15 * c*c*e2*e2*pow(yaux1, 4)*yaux2 + \
			(-6)*c*c*e2*yaux1*yaux5*yaux2 + 3 * e1*k*yaux0*yaux0*(k*yaux4 + 2 * c*e2*yaux1*yaux1*yaux1 + c * yaux5 + \
			(-2)*m*yaux2) + 2 * k*yaux4*(2 * k + (-3)*c*e2*yaux1*yaux2) + 2 * e1*k*yaux0*yaux0*yaux0*((-4)*k + 3 * c*e2*yaux1*yaux2)\
			+ k * yaux0*((-5)*k + (-6)*e1*m*yaux1*yaux1 + 12 * c*e2*yaux1*yaux2)\
			+ (-4)*c*c*yaux6 + 2 * k*m*yaux6 + (-3)*c*c*e2*yaux1*yaux1*yaux6);

		b[3] = exp10;
		b[4] = yaux5;
		b[5] = yaux6;
		b[6] = yaux7;
		exp20 = pow(m, -2)*(4 * k*k*yaux0 + e1 * k*k*yaux0*yaux0*yaux0 + (-5)*k*k*yaux4 + (-3)*c*e1*k*yaux0*yaux0*yaux1 + \
			c*e2*k*yaux1*yaux1*yaux1 + (-4)*c*c*yaux2 + 2 * k*m*yaux2 + (-3)*c*c*e2*yaux1*yaux1*yaux2 + 5 * c*c*yaux6 + (-4)*k*m*yaux6);
		b[7] = exp20;

		k30 = b[0] * h; k31 = b[1] * h; k32 = b[2] * h; k33 = b[3] * h;
		k34 = b[4] * h; k35 = b[5] * h; k36 = b[6] * h; k37 = b[7] * h;

		yaux0 = y0 + k30;
		yaux1 = y1 + k31;
		yaux2 = y2 + k32;
		yaux3 = y3 + k33;
		yaux4 = y4 + k34;
		yaux5 = y5 + k35;
		yaux6 = y6 + k36;
		yaux7 = y7 + k37;

		//b = OdeFun(t0+h,yaux0,yaux1,yaux2,yaux3); //k4
		b[0] = yaux1;
		b[1] = yaux2;
		b[2] = yaux3;
		exp10 = pow(m, -2)*((-3)*e1*e1*k*k*pow(yaux0, 5) + 4 * c*e2*k*pow(yaux1, 3) + (-3)*c*e2*k*yaux1*yaux1*yaux5 + \
			5 * c*c*yaux2 + (-4)*k*m*yaux2 + 24 * c*c*e2*yaux1*yaux1*yaux2 + 15 * c*c*e2*e2*pow(yaux1, 4)*yaux2 + \
			(-6)*c*c*e2*yaux1*yaux5*yaux2 + 3 * e1*k*yaux0*yaux0*(k*yaux4 + 2 * c*e2*yaux1*yaux1*yaux1 + c * yaux5 + \
			(-2)*m*yaux2) + 2 * k*yaux4*(2 * k + (-3)*c*e2*yaux1*yaux2) + 2 * e1*k*yaux0*yaux0*yaux0*((-4)*k + 3 * c*e2*yaux1*yaux2)\
			+ k * yaux0*((-5)*k + (-6)*e1*m*yaux1*yaux1 + 12 * c*e2*yaux1*yaux2)\
			+ (-4)*c*c*yaux6 + 2 * k*m*yaux6 + (-3)*c*c*e2*yaux1*yaux1*yaux6);

		b[3] = exp10;
		b[4] = yaux5;
		b[5] = yaux6;
		b[6] = yaux7;
		exp20 = pow(m, -2)*(4 * k*k*yaux0 + e1 * k*k*yaux0*yaux0*yaux0 + (-5)*k*k*yaux4 + (-3)*c*e1*k*yaux0*yaux0*yaux1 + \
			c*e2*k*yaux1*yaux1*yaux1 + (-4)*c*c*yaux2 + 2 * k*m*yaux2 + (-3)*c*c*e2*yaux1*yaux1*yaux2 + 5 * c*c*yaux6 + (-4)*k*m*yaux6);
		b[7] = exp20;

		k40 = b[0] * h; k41 = b[1] * h; k42 = b[2] * h; k43 = b[3] * h;
		k44 = b[4] * h; k45 = b[5] * h; k46 = b[6] * h; k47 = b[7] * h;

		y0 = y0 + (k10 + 2 * k20 + 2 * k30 + k40) / 6;
		y1 = y1 + (k11 + 2 * k21 + 2 * k31 + k41) / 6;
		y2 = y2 + (k12 + 2 * k22 + 2 * k32 + k42) / 6;
		y3 = y3 + (k13 + 2 * k23 + 2 * k33 + k43) / 6;
		y4 = y4 + (k14 + 2 * k24 + 2 * k34 + k44) / 6;
		y5 = y5 + (k15 + 2 * k25 + 2 * k35 + k45) / 6;
		y6 = y6 + (k16 + 2 * k26 + 2 * k36 + k46) / 6;
		y7 = y7 + (k17 + 2 * k27 + 2 * k37 + k47) / 6;

		t0 = t0 + h;

		resp[cont + 0 + pi] = y0;
		resp[cont + 1 + pi] = y1;
		resp[cont + 2 + pi] = y2;
		resp[cont + 3 + pi] = y3;
		resp[cont + 4 + pi] = y4;
		resp[cont + 5 + pi] = y5;
		resp[cont + 6 + pi] = y6;
		resp[cont + 7 + pi] = y7;
		cont = cont + 8;


	}
	
	return resp;

}

double * numjac(double *resp, double bc0, double bc1, double bc2, double bc3, double yb0, double yb1, double yb2, double yb3, int Nt0, int Nx, int pi, int pf) {

	double dy, v1, v2, bcaux;
	int num, i;
	static double  Jac[16];

	num = Nt0 - 1;
	dy = 0.0001;

	//=========================================================================
	i = 0;
	bcaux = bc0 - dy;
	resp = rk4(bcaux, bc1, bc2, bc3, resp, pi, pf, Nx);
	v1 = resp[num + i - 7 + pi] - yb0;

	bcaux = bc0 + dy;
	resp = rk4(bcaux, bc1, bc2, bc3, resp, pi, pf, Nx);
	v2 = resp[num + i - 7 + pi] - yb0;

	Jac[0] = (v2 - v1) / (2 * dy);
	//=========================================================================
	i = 0;
	bcaux = bc1 - dy;
	resp = rk4(bc0, bcaux, bc2, bc3, resp, pi, pf, Nx);
	v1 = resp[num + i - 7 + pi] - yb0;

	bcaux = bc1 + dy;
	resp = rk4(bc0, bcaux, bc2, bc3, resp, pi, pf, Nx);
	v2 = resp[num + i - 7 + pi] - yb0;

	Jac[1] = (v2 - v1) / (2 * dy);
	//=========================================================================
	i = 0;
	bcaux = bc2 - dy;
	resp = rk4(bc0, bc1, bcaux, bc3, resp, pi, pf, Nx);
	v1 = resp[num + i - 7 + pi] - yb0;

	bcaux = bc2 + dy;
	resp = rk4(bc0, bc1, bcaux, bc3, resp, pi, pf, Nx);
	v2 = resp[num + i - 7 + pi] - yb0;

	Jac[2] = (v2 - v1) / (2 * dy);
	//=========================================================================
	i = 0;
	bcaux = bc3 - dy;
	resp = rk4(bc0, bc1, bc2, bcaux, resp, pi, pf, Nx);
	v1 = resp[num + i - 7 + pi] - yb0;

	bcaux = bc3 + dy;
	resp = rk4(bc0, bc1, bc2, bcaux, resp, pi, pf, Nx);
	v2 = resp[num + i - 7 + pi] - yb0;

	Jac[3] = (v2 - v1) / (2 * dy);
	//=========================================================================
	i = 1;
	bcaux = bc0 - dy;
	resp = rk4(bcaux, bc1, bc2, bc3, resp, pi, pf, Nx);
	v1 = resp[num + i - 7 + pi] - yb1;

	bcaux = bc0 + dy;
	resp = rk4(bcaux, bc1, bc2, bc3, resp, pi, pf, Nx);
	v2 = resp[num + i - 7 + pi] - yb1;

	Jac[4] = (v2 - v1) / (2 * dy);
	//=========================================================================
	i = 1;
	bcaux = bc1 - dy;
	resp = rk4(bc0, bcaux, bc2, bc3, resp, pi, pf, Nx);
	v1 = resp[num + i - 7 + pi] - yb1;

	bcaux = bc1 + dy;
	resp = rk4(bc0, bcaux, bc2, bc3, resp, pi, pf, Nx);
	v2 = resp[num + i - 7 + pi] - yb1;

	Jac[5] = (v2 - v1) / (2 * dy);
	//=========================================================================
	i = 1;
	bcaux = bc2 - dy;
	resp = rk4(bc0, bc1, bcaux, bc3, resp, pi, pf, Nx);
	v1 = resp[num + i - 7 + pi] - yb1;

	bcaux = bc2 + dy;
	resp = rk4(bc0, bc1, bcaux, bc3, resp, pi, pf, Nx);
	v2 = resp[num + i - 7 + pi] - yb1;

	Jac[6] = (v2 - v1) / (2 * dy);
	//=========================================================================
	i = 1;
	bcaux = bc3 - dy;
	resp = rk4(bc0, bc1, bc2, bcaux, resp, pi, pf, Nx);
	v1 = resp[num + i - 7 + pi] - yb1;

	bcaux = bc3 + dy;
	resp = rk4(bc0, bc1, bc2, bcaux, resp, pi, pf, Nx);
	v2 = resp[num + i - 7 + pi] - yb1;

	Jac[7] = (v2 - v1) / (2 * dy);
	//=========================================================================
	i = 0;
	bcaux = bc0 - dy;
	resp = rk4(bcaux, bc1, bc2, bc3, resp, pi, pf, Nx);
	v1 = resp[num + i - 3 + pi] - yb2;

	bcaux = bc0 + dy;
	resp = rk4(bcaux, bc1, bc2, bc3, resp, pi, pf, Nx);
	v2 = resp[num + i - 3 + pi] - yb2;

	Jac[8] = (v2 - v1) / (2 * dy);
	//=========================================================================
	i = 0;
	bcaux = bc1 - dy;
	resp = rk4(bc0, bcaux, bc2, bc3, resp, pi, pf, Nx);
	v1 = resp[num + i - 3 + pi] - yb2;

	bcaux = bc1 + dy;
	resp = rk4(bc0, bcaux, bc2, bc3, resp, pi, pf, Nx);
	v2 = resp[num + i - 3 + pi] - yb2;

	Jac[9] = (v2 - v1) / (2 * dy);
	//=========================================================================
	i = 0;
	bcaux = bc2 - dy;
	resp = rk4(bc0, bc1, bcaux, bc3, resp, pi, pf, Nx);
	v1 = resp[num + i - 3 + pi] - yb2;

	bcaux = bc2 + dy;
	resp = rk4(bc0, bc1, bcaux, bc3, resp, pi, pf, Nx);
	v2 = resp[num + i - 3 + pi] - yb2;

	Jac[10] = (v2 - v1) / (2 * dy);
	//=========================================================================
	i = 0;
	bcaux = bc3 - dy;
	resp = rk4(bc0, bc1, bc2, bcaux, resp, pi, pf, Nx);
	v1 = resp[num + i - 3 + pi] - yb2;

	bcaux = bc3 + dy;
	resp = rk4(bc0, bc1, bc2, bcaux, resp, pi, pf, Nx);
	v2 = resp[num + i - 3 + pi] - yb2;

	Jac[11] = (v2 - v1) / (2 * dy);
	//=========================================================================
	i = 1;
	bcaux = bc0 - dy;
	resp = rk4(bcaux, bc1, bc2, bc3, resp, pi, pf, Nx);
	v1 = resp[num + i - 3 + pi] - yb3;

	bcaux = bc0 + dy;
	resp = rk4(bcaux, bc1, bc2, bc3, resp, pi, pf, Nx);
	v2 = resp[num + i - 3 + pi] - yb3;

	Jac[12] = (v2 - v1) / (2 * dy);

	//=========================================================================
	i = 1;
	bcaux = bc1 - dy;
	resp = rk4(bc0, bcaux, bc2, bc3, resp, pi, pf, Nx);
	v1 = resp[num + i - 3 + pi] - yb3;

	bcaux = bc1 + dy;
	resp = rk4(bc0, bcaux, bc2, bc3, resp, pi, pf, Nx);
	v2 = resp[num + i - 3 + pi] - yb3;

	Jac[13] = (v2 - v1) / (2 * dy);
	//=========================================================================
	i = 1;
	bcaux = bc2 - dy;
	resp = rk4(bc0, bc1, bcaux, bc3, resp, pi, pf, Nx);
	v1 = resp[num + i - 3 + pi] - yb3;

	bcaux = bc2 + dy;
	resp = rk4(bc0, bc1, bcaux, bc3, resp, pi, pf, Nx);
	v2 = resp[num + i - 3 + pi] - yb3;

	Jac[14] = (v2 - v1) / (2 * dy);
	//=========================================================================
	i = 1;
	bcaux = bc3 - dy;
	resp = rk4(bc0, bc1, bc2, bcaux, resp, pi, pf, Nx);
	v1 = resp[num + i - 3 + pi] - yb3;

	bcaux = bc3 + dy;
	resp = rk4(bc0, bc1, bc2, bcaux, resp, pi, pf, Nx);
	v2 = resp[num + i - 3 + pi] - yb3;

	Jac[15] = (v2 - v1) / (2 * dy);
	//=========================================================================

	return Jac;
}

double bvpsol(double *resp, double *J, double yb0, double yb1, double yb2, double yb3, int Nt0, int Nx, int pi, int pf) {

	static double  yr[8];
	double Deter;

	double J11, J12, J13, J14;
	double J21, J22, J23, J24;
	double J31, J32, J33, J34;
	double J41, J42, J43, J44;

	double J11i, J12i, J13i, J14i;
	double J21i, J22i, J23i, J24i;
	double J31i, J32i, J33i, J34i;
	double J41i, J42i, J43i, J44i;

	double D1, D2, D3, D4;
	double cof1, cof2, cof3, cof4;

	double r0, r1, r2, r3, bc00, bc10, bc20, bc30, bc0, bc1, bc2, bc3, erro1, erro2;
	double tol = 0.0001;
	int num = Nt0 - 1;
	int Nmax = 100, runs;

	//Modify!!!
	double m = 1;
	double c = 0.1;
	double k = 1;
	double e1 = 1;
	double e2 = 1;
	double myPI = 3.14159265358979323846;
	double soma, pdf, fu0, fu1, h;
	double v0, v1, v2, v3, v4, v5, v6, v7;
	double S0 = 1;
	double xi = 0.0;
	double xf = 1.0;
	int N0 = Nx, cont;

	bc0 = 0.0; bc1 = 0.0; bc2 = 0.0; bc3 = 0.0;
	J[0] = 0.0;J[1] = 0.0;J[2] = 0.0;J[3] = 0.0;
	J[4] = 0.0;J[5] = 0.0;J[6] = 0.0;J[7] = 0.0;
	J[8] = 0.0;J[9] = 0.0;J[10] = 0.0;J[11] = 0.0;
	J[12] = 0.0;J[13] = 0.0;J[14] = 0.0;J[15] = 0.0;

	runs = 1;
	for (int i = 0;i < Nmax;i++) {

		J = numjac(resp, bc0, bc1, bc2, bc3, yb0, yb1, yb2, yb3, Nt0, N0, pi, pf);

		J11 = J[0];  J12 = J[1];  J13 = J[2]; J14 = J[3];
		J21 = J[4];  J22 = J[5];  J23 = J[6]; J24 = J[7];
		J31 = J[8];  J32 = J[9];  J33 = J[10];J34 = J[11];
		J41 = J[12]; J42 = J[13]; J43 = J[14];J44 = J[15];

		//Calcule the determinant
		D1 = J22 * J33*J44 + J23 * J34*J42 + J24 * J32*J43 - J24 * J33*J42 - J23 * J32*J44 - J22 * J34*J43;
		D2 = J12 * J33*J44 + J13 * J34*J42 + J14 * J32*J43 - J14 * J33*J42 - J13 * J32*J44 - J12 * J34*J43;
		D3 = J12 * J23*J44 + J13 * J24*J42 + J14 * J22*J43 - J14 * J23*J42 - J13 * J22*J44 - J12 * J24*J43;
		D4 = J12 * J23*J34 + J13 * J24*J32 + J14 * J22*J33 - J14 * J23*J32 - J13 * J22*J34 - J12 * J24*J33;
		cof1 = J11; cof2 = J21; cof3 = J31; cof4 = J41;
	
		Deter = cof1 * D1 - cof2 * D2 + cof3 * D3 - cof4 * D4;

		J11i = (J22*J33*J44 + J23 * J34*J42 + J24 * J32*J43 - J24 * J33*J42 - J23 * J32*J44 - J22 * J34*J43) / Deter;
		J12i = (-J12 * J33*J44 - J13 * J34*J42 - J14 * J32*J43 + J14 * J33*J42 + J13 * J32*J44 + J12 * J34*J43) / Deter;
		J13i = (J12*J23*J44 + J13 * J24*J42 + J14 * J22*J43 - J14 * J23*J42 - J13 * J22*J44 - J12 * J24*J43) / Deter;
		J14i = (-J12 * J23*J34 - J13 * J24*J32 - J14 * J22*J33 + J14 * J23*J32 + J13 * J22*J34 + J14 * J24*J33) / Deter;

		J21i = (-J21 * J33*J44 - J23 * J34*J41 - J24 * J31*J43 + J24 * J33*J41 + J23 * J31*J44 + J21 * J34*J43) / Deter;
		J22i = (J11*J33*J44 + J13 * J34*J41 + J14 * J31*J43 - J14 * J33*J41 - J13 * J31*J44 - J11 * J34*J43) / Deter;
		J23i = (-J11 * J23*J44 - J13 * J24*J41 - J14 * J21*J43 + J14 * J23*J41 + J13 * J21*J44 + J11 * J24*J43) / Deter;
		J24i = (J11*J23*J34 + J13 * J24*J31 + J14 * J21*J33 - J14 * J23*J31 - J13 * J21*J34 - J11 * J24*J33) / Deter;

		J31i = (J21*J32*J44 + J22 * J34*J41 + J24 * J31*J42 - J24 * J32*J41 - J22 * J31*J44 - J21 * J34*J42) / Deter;
		J32i = (-J11 * J32*J44 - J12 * J34*J41 - J14 * J31*J42 + J14 * J32*J41 + J12 * J31*J44 + J11 * J34*J42) / Deter;
		J33i = (J11*J22*J44 + J12 * J24*J41 + J14 * J21*J42 - J14 * J22*J41 - J12 * J21*J44 - J11 * J24*J42) / Deter;
		J34i = (-J11 * J22*J34 - J12 * J24*J31 - J14 * J21*J32 + J14 * J22*J31 + J12 * J21*J34 + J11 * J24*J32) / Deter;

		J41i = (-J21 * J32*J43 - J22 * J33*J41 - J23 * J31*J42 + J23 * J32*J41 + J22 * J31*J43 + J21 * J33*J42) / Deter;
		J42i = (J11*J32*J43 + J12 * J33*J41 + J13 * J31*J42 - J13 * J32*J41 - J12 * J31*J43 - J11 * J33*J42) / Deter;
		J43i = (-J11 * J22*J43 - J12 * J33*J41 - J13 * J21*J42 + J13 * J22*J41 + J12 * J21*J43 + J11 * J23*J42) / Deter;
		J44i = (J11*J22*J33 + J12 * J23*J31 + J13 * J21*J32 - J13 * J22*J31 - J12 * J21*J33 - J11 * J23*J32) / Deter;

		resp = rk4(bc0, bc1, bc2, bc3, resp, pi, pf, N0);
		r0 = resp[num + 0 - 7 + pi] - yb0;
		r1 = resp[num + 1 - 7 + pi] - yb1;
		r2 = resp[num + 0 - 3 + pi] - yb2;
		r3 = resp[num + 1 - 3 + pi] - yb3;

		bc00 = bc0; bc10 = bc1; bc20 = bc2; bc30 = bc3;

		bc0 = bc0 - (J11i*r0 + J12i * r1 + J13i * r2 + J14i * r3);
		bc1 = bc1 - (J21i*r0 + J22i * r1 + J23i * r2 + J24i * r3);
		bc2 = bc2 - (J31i*r0 + J32i * r1 + J33i * r2 + J34i * r3);
		bc3 = bc3 - (J41i*r0 + J42i * r1 + J43i * r2 + J44i * r3);

		erro1 = sqrt((bc0 - bc00)*(bc0 - bc00) + (bc1 - bc10)*(bc1 - bc10) + (bc2 - bc20)*(bc2 - bc20) + (bc3 - bc30)*(bc3 - bc30));
		erro2 = sqrt(r0*r0 + r1*r1 + r2*r2 + r3*r3);

		if (erro1 <= tol) { break; }
		if (erro2 <= tol) { break; }
		runs = runs + 1;

	}
	//printf("%f, %f, %d\n", erro1, erro2, runs);
	h = (xf - xi) / (N0 - 1);
	cont = 0;
	v0 = resp[cont + pi];
	v1 = resp[cont + 1 + pi];
	v2 = resp[cont + 2 + pi];
	v3 = resp[cont + 3 + pi];
	v4 = resp[cont + 4 + pi];
	v5 = resp[cont + 5 + pi];
	v6 = resp[cont + 6 + pi];
	v7 = resp[cont + 7 + pi];
	fu0 = 0.25*(1.0 / (myPI*S0))*(pow(2*k*v0 + e1*k*v0*v0*v0 - k*v4 + 2*c*v1 + c*e2*v1*v1*v1 - c*v5 + m*v2, 2) +\
		pow(k*v0 - 2*k*v4 + c*v1 - 2*c*v5 - m*v6, 2));
	cont = 8;
	soma = 0.0;
	for (int i = 1;i < N0;i++) {
		v0 = resp[cont + pi];
		v1 = resp[cont + 1 + pi];
		v2 = resp[cont + 2 + pi];
		v3 = resp[cont + 3 + pi];
		v4 = resp[cont + 4 + pi];
		v5 = resp[cont + 5 + pi];
		v6 = resp[cont + 6 + pi];
		v7 = resp[cont + 7 + pi];

		fu1 = 0.25*(1.0 / (myPI*S0))*(pow(2 * k*v0 + e1 * k*v0*v0*v0 - k * v4 + 2 * c*v1 + c * e2*v1*v1*v1 - c * v5 + m * v2, 2) + \
			pow(k*v0 - 2 * k*v4 + c * v1 - 2 * c*v5 - m * v6, 2));

		soma = soma + 0.5*(fu0 + fu1)*h;

		fu0 = fu1;

		cont = cont + 8;
	}

	//pdf = exp(-soma);
	pdf = soma;

	return pdf;
}

//Kernel
void shoot(double *c, double *X1, double *X2, double *X3, double *X4, double *resp, double *J, int Nt, int Nx, int nrd, int N, int n1max, int n2max, int n3max, int n4max) {

	double x1, x2, x3, x4;
	int Nt0 = nrd*Nx;
	int pi, pf;
	int index;

	/*index = 0;
	x1 = X1[index];
	x2 = X2[index];
	x3 = X3[index];
	x4 = X4[index];
	pi = index * Nt0; pf = (index + 1)*Nt0 - 1;
	c[index] = bvpsol(resp, J, x1, x2, x3, x4, Nt0, Nx, pi, pf);//should be x1 and x2
	index = index + 1;*/

	index = 0;
	for (int i = 0;i < n4max;i++) {
		for (int j = 0;j < n3max;j++) {
			for (int k = 0;k < n2max;k++) {
				for (int s = 0;s < n1max;s++) {
					x1 = X1[index];
					x2 = X2[index];
					x3 = X3[index];
					x4 = X4[index];
					//printf("i = %f, j = %f, k = %f, s = %f\n", x1, x2, x3, x4);
					pi = index * Nt0; pf = (index + 1)*Nt0 - 1;
					c[index] = bvpsol(resp, J, x1, x2, x3, x4, Nt0, Nx, pi, pf);//should be x1 and x2
					index = index + 1;
				}
			}
		}
	}
}

//############################################################################################

int main(void) {

	double x1_max = 2, x2_max = 2, x3_max = 2, x4_max = 2;
	double dx1, dx2, dx3, dx4, x10, x20, x30, x40;
	int neval = 16, nevalt = 1, cont;
	int n1max = neval, n2max = neval, n3max = neval, n4max = neval;
	//========================================================================================
	int N = n1max*n2max*n3max*n4max; //Grid dimension
	int nrd = 8;         //System order
	int Nx = 300;          //Discretization BVP
	int Nt = Nx * nrd;     //Vector size

	int size_N = N * sizeof(double);
	int size_J = 16 * sizeof(double); //CHECK THIS IN THE SDOF
	int size_resp = N * Nt * sizeof(double);
	int size_int = sizeof(int);

	double *c, *X1, *X2, *X3, *X4, *resp, *J;

	FILE * fp;
	fp = fopen("CPU2DOF.txt", "w");

	//========================================================================================
	dx1 = 2 * x1_max / ((float)(n1max - 1));
	dx2 = 2 * x2_max / ((float)(n2max - 1));
	dx3 = 2 * x3_max / ((float)(n3max - 1));
	dx4 = 2 * x4_max / ((float)(n4max - 1));
	//========================================================================================

	J = (double *)malloc(size_J);
	c = (double *)malloc(size_N);
	X1 = (double *)malloc(size_N);
	X2 = (double *)malloc(size_N);
	X3 = (double *)malloc(size_N);
	X4 = (double *)malloc(size_N);
	resp = (double *)malloc(size_resp);

	for (int i = 0;i < 16;i++) {
		J[i] = 0.0;
	}

	cont = 0;
	x40 = -x4_max;
	for (int i = 0;i < n4max;i++) {
		x30 = -x3_max;
		for (int j = 0;j < n3max;j++) {
			x20 = -x2_max;
			for (int k = 0;k < n2max;k++) {
				x10 = -x1_max;
				for (int s = 0;s < n1max;s++) {

					X1[cont] = x10;
					X2[cont] = x20;
					X3[cont] = x30;
					X4[cont] = x40;

					x10 = x10 + dx1;
					cont = cont + 1;
				}
				x20 = x20 + dx2;
			}
			x30 = x30 + dx3;
		}
		x40 = x40 + dx4;
	}
	//========================================================================================

	clock_t tStart = clock();
	shoot(c, X1, X2, X3, X4, resp, J, Nt, Nx, nrd, N, n1max, n2max, n3max, n4max);
	printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

	for (int i = 0;i < N;i++) {
		printf("%f, %f, %f, %f, %f\n", X1[i], X2[i], X3[i], X4[i], c[i]);
		fprintf(fp, "%f, %f, %f, %f, %f\n", X1[i], X2[i], X3[i], X4[i], c[i]);

	}

	free(c); free(X1); free(X2); free(X3); free(X4); free(resp); free(J);

	fclose(fp);

	return 0;
}

