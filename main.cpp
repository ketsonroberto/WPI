#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


//Author: Ketson R. M. dos Santos
//Institution: Columbia University
//e-mail: kmd2191@columbia.edu
//
//This is the C/C++ version of the program
//to solve the WPI problem in stochastic dynamics (Single-degree-of-freedom)

double * rk4(double xx, double yy, double *resp, int pi, int pf) {

	double m = 1;
	double c = 0.1;
	double k = 1;
	double e = 0;
	static double  b[4];

	double xi = 0.0;
	double xf = 1.0;
	double h, t0;
	double y0, y1, y2, y3;
	double yaux0, yaux1, yaux2, yaux3;
	double k10, k11, k12, k13;
	double k20, k21, k22, k23;
	double k30, k31, k32, k33;
	double k40, k41, k42, k43;
	int N0 = 300, cont;

	h = (xf - xi) / (N0 - 1);
	t0 = xi;
	y0 = 0; y1 = 0; y2 = xx; y3 = yy;

	resp[0 + pi] = y0;
	resp[1 + pi] = y1;
	resp[2 + pi] = y2;
	resp[3 + pi] = y3;

	cont = 4;
	for (int i = 1;i < 300;i++) {
		//b = OdeFun(t0,y0,y1,y2,y3); //k1
		b[0] = y1;
		b[1] = y2;
		b[2] = y3;
		b[3] = pow(m, -2)*((-1)*k*k*y0 + (-4)*e*k*k*pow(y0, 3) + (-3)*pow(e, 2)*pow(k, 2)*pow(y0, 5) + (-6)*e*k*m*y0*pow(y1, 2) + pow(c, 2)*y2 + (-2)*k*m*y2 + (-6)*e*k*m*pow(y0, 2)*y2);

		k10 = b[0] * h; k11 = b[1] * h; k12 = b[2] * h; k13 = b[3] * h;

		yaux0 = y0 + 0.5*k10;
		yaux1 = y1 + 0.5*k11;
		yaux2 = y2 + 0.5*k12;
		yaux3 = y3 + 0.5*k13;

		//b = OdeFun(t0+0.5*h,yaux0,yaux1,yaux2,yaux3); //k2
		b[0] = yaux1;
		b[1] = yaux2;
		b[2] = yaux3;
		b[3] = pow(m, -2)*((-1)*k*k*yaux0 + (-4)*e*k*k*pow(yaux0, 3) + (-3)*pow(e, 2)*pow(k, 2)*pow(yaux0, 5) + (-6)*e*k*m*yaux0*pow(yaux1, 2) + pow(c, 2)*yaux2 + (-2)*k*m*yaux2 + (-6)*e*k*m*pow(yaux0, 2)*yaux2);

		k20 = b[0] * h; k21 = b[1] * h; k22 = b[2] * h; k23 = b[3] * h;

		yaux0 = y0 + 0.5*k20;
		yaux1 = y1 + 0.5*k21;
		yaux2 = y2 + 0.5*k22;
		yaux3 = y3 + 0.5*k23;

		//b = OdeFun(t0+0.5*h,yaux0,yaux1,yaux2,yaux3); //k3
		b[0] = yaux1;
		b[1] = yaux2;
		b[2] = yaux3;
		b[3] = pow(m, -2)*((-1)*k*k*yaux0 + (-4)*e*k*k*pow(yaux0, 3) + (-3)*pow(e, 2)*pow(k, 2)*pow(yaux0, 5) + (-6)*e*k*m*yaux0*pow(yaux1, 2) + pow(c, 2)*yaux2 + (-2)*k*m*yaux2 + (-6)*e*k*m*pow(yaux0, 2)*yaux2);

		k30 = b[0] * h; k31 = b[1] * h; k32 = b[2] * h; k33 = b[3] * h;

		yaux0 = y0 + k30;
		yaux1 = y1 + k31;
		yaux2 = y2 + k32;
		yaux3 = y3 + k33;

		//b = OdeFun(t0+h,yaux0,yaux1,yaux2,yaux3); //k4
		b[0] = yaux1;
		b[1] = yaux2;
		b[2] = yaux3;
		b[3] = pow(m, -2)*((-1)*k*k*yaux0 + (-4)*e*k*k*pow(yaux0, 3) + (-3)*pow(e, 2)*pow(k, 2)*pow(yaux0, 5) + (-6)*e*k*m*yaux0*pow(yaux1, 2) + pow(c, 2)*yaux2 + (-2)*k*m*yaux2 + (-6)*e*k*m*pow(yaux0, 2)*yaux2);

		k40 = b[0] * h; k41 = b[1] * h; k42 = b[2] * h; k43 = b[3] * h;

		y0 = y0 + (k10 + 2 * k20 + 2 * k30 + k40) / 6;
		y1 = y1 + (k11 + 2 * k21 + 2 * k31 + k41) / 6;
		y2 = y2 + (k12 + 2 * k22 + 2 * k32 + k42) / 6;
		y3 = y3 + (k13 + 2 * k23 + 2 * k33 + k43) / 6;

		t0 = t0 + h;

		resp[cont + 0 + pi] = y0;
		resp[cont + 1 + pi] = y1;
		resp[cont + 2 + pi] = y2;
		resp[cont + 3 + pi] = y3;
		cont = cont + 4;
	}

	return resp;

}

double * numjac(double *resp, double xx, double yy, double yb1, double yb2, int Nt0, int pi, int pf) {

	double dy, xx0, yy0, v1, v2;
	int num, i;
	static double  Jac[4];

	num = Nt0 - 1;
	dy = 0.0001;

	i = 0;
	xx0 = xx - dy;
	resp = rk4(xx0, yy, resp, pi, pf);
	v1 = resp[num + i - 3 + pi] - yb1;

	xx0 = xx + dy;
	resp = rk4(xx0, yy, resp, pi, pf);
	v2 = resp[num + i - 3 + pi] - yb1;

	Jac[0] = (v2 - v1) / (2 * dy);

	i = 0;
	yy0 = yy - dy;
	resp = rk4(xx, yy0, resp, pi, pf);
	v1 = resp[num + i - 3 + pi] - yb1;

	yy0 = yy + dy;
	resp = rk4(xx, yy0, resp, pi, pf);
	v2 = resp[num + i - 3 + pi] - yb1;

	Jac[1] = (v2 - v1) / (2 * dy);


	i = 1;
	xx0 = xx - dy;
	resp = rk4(xx0, yy, resp, pi, pf);
	v1 = resp[num + i - 3 + pi] - yb2;

	xx0 = xx + dy;
	resp = rk4(xx0, yy, resp, pi, pf);
	v2 = resp[num + i - 3 + pi] - yb2;

	Jac[2] = (v2 - v1) / (2 * dy);


	i = 1;
	yy0 = yy - dy;
	resp = rk4(xx, yy0, resp, pi, pf);
	v1 = resp[num + i - 3 + pi] - yb2;

	yy0 = yy + dy;
	resp = rk4(xx, yy0, resp, pi, pf);
	v2 = resp[num + i - 3 + pi] - yb2;

	Jac[3] = (v2 - v1) / (2 * dy);

	return Jac;
}

double bvpsol(double *resp, double *J, double yb1, double yb2, int Nt0, int pi, int pf) {

	static double  yr[4];
	double det;
	double J0, J1, J2, J3;
	double J0i, J1i, J2i, J3i;
	double r0, r1, x0, y0, xx, yy, erro1, erro2;
	double tol = 0.00001;
	int num = Nt0 - 1;
	int Nmax = 100, runs;

	//Modify!!!
	double m = 1;
	double c = 0.1;
	double k = 1;
	double e = 0;
	double myPI = 3.14159265358979323846;
	double soma, pdf, v0, v1, v2, v3, fu0, fu1, h;
	double S0 = 1;
	double xi = 0.0;
	double xf = 1.0;
	int N0 = 300, cont;

	xx = 0.0; yy = 0.0;
	J[0] = 0.0;J[1] = 0.0;J[2] = 0.0;J[3] = 0.0;

	runs = 1;
	for (int i = 0;i < Nmax;i++) {

		J = numjac(resp, xx, yy, yb1, yb2, Nt0, pi, pf);

		J0 = J[0]; J1 = J[1];J2 = J[2];J3 = J[3];

		det = J0 * J3 - J1 * J2;
		J0i = J3 / det;
		J3i = J0 / det;
		J1i = -J1 / det;
		J2i = -J2 / det;

		resp = rk4(xx, yy, resp, pi, pf);
		r0 = resp[num + 0 - 3 + pi] - yb1;
		r1 = resp[num + 1 - 3 + pi] - yb2;

		x0 = xx; y0 = yy;

		xx = xx - (J0i*r0 + J1i * r1);
		yy = yy - (J2i*r0 + J3i * r1);

		erro1 = sqrt((xx - x0)*(xx - x0) + (yy - y0)*(yy - y0));
		erro2 = sqrt(r0*r0 + r1 * r1);

		if (erro1 <= tol) { break; }
		if (erro2 <= tol) { break; }
		runs = runs + 1;

	}

	h = (xf - xi) / (N0 - 1);
	cont = 0;
	v0 = resp[cont + pi];
	v1 = resp[cont + 1 + pi];
	v2 = resp[cont + 2 + pi];
	v3 = resp[cont + 3 + pi];
	fu0 = 0.25*(1.0 / (myPI*S0))*pow(k*v0 + e * k*v0*v0*v0 + c * v1 + m * v2, 2);
	cont = 4;
	soma = 0.0;
	for (int i = 1;i < N0;i++) {
		v0 = resp[cont + pi];
		v1 = resp[cont + 1 + pi];
		v2 = resp[cont + 2 + pi];
		v3 = resp[cont + 3 + pi];

		fu1 = 0.25*(1.0 / (myPI*S0))*pow(k*v0 + e * k*v0*v0*v0 + c * v1 + m * v2, 2);

		soma = soma + 0.5*(fu0 + fu1)*h;

		fu0 = fu1;

		cont = cont + 4;
	}

	pdf = exp(-soma);

	return pdf;
}

//Kernel
void shoot(double *c, double *X1, double *X2, double *resp, double *J, int Nt, int n1max, int n2max) {

	double x1, x2;
	int Nt0 = 1200;
	int pi, pf;
	int index;

	index = 0;
	for (int i = 0;i < n2max;i++) {
		for (int j = 0;j < n1max;j++) {

			x1 = X1[index];
			x2 = X2[index];
			pi = index * Nt0; pf = (index + 1)*Nt0 - 1;
			c[index] = bvpsol(resp, J, x1, x2, Nt0, pi, pf);//should be x1 and x2
                        index = index + 1;

		}
	}

}

//############################################################################################

int main(void) {

	double x1_max = 7.0, x2_max = 10.0, dx1, dx2;
	double x10, x20;
	int n1max = 32, n2max = 32, cont;
	//========================================================================================
	int N = n1max * n2max; //Grid dimension
	int nrd = 4;         //System order
	int Nx = 300;          //Discretization BVP
	int Nt = Nx * nrd;     //Vector size

	int size_N = N * sizeof(double);
	int size_J = nrd * sizeof(double);
	int size_resp = N * Nt * sizeof(double);
	int size_int = sizeof(int);

	double *c, *X1, *X2, *resp, *J;

	FILE * fp;
	fp = fopen("CPU.txt", "w");

	//========================================================================================
	dx1 = 2 * x1_max / ((float)(n1max - 1));
	dx2 = 2 * x2_max / ((float)(n2max - 1));
	//========================================================================================

	J = (double *)malloc(size_J);
	c = (double *)malloc(size_N);
	X1 = (double *)malloc(size_N);
	X2 = (double *)malloc(size_N);
	resp = (double *)malloc(size_resp);

	J[0] = 0.0;
	J[1] = 0.0;
	J[2] = 0.0;
	J[3] = 0.0;

	cont = 0;
	x20 = -x2_max;
	for (int i = 0;i < n2max;i++) {
		x10 = -x1_max;
		for (int j = 0;j < n1max;j++) {

			X1[cont] = x10;
			X2[cont] = x20;

			x10 = x10 + dx1;
			cont = cont + 1;

		}
		x20 = x20 + dx2;
	}
	//========================================================================================

        clock_t tStart = clock();
	shoot(c, X1, X2, resp, J, Nt, n1max, n2max);
        printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

	for (int i = 0;i < N;i++) {
		printf("%f, %f, %f\n", X1[i], X2[i], c[i]);
		fprintf(fp, "%f, %f, %f\n", X1[i], X2[i], c[i]);

	}

	free(c); free(X1); free(X2); free(resp); free(J);

	fclose(fp);

	return 0;
}
