#ifndef DLM_BESSEL_H
#define DLM_BESSEL_H

double DLM_Bessel1(const double& order, const double& arg, const bool& Threadsafe=true);

//code taken from:
//Numerical Recipes
//The Art of Scientific Computing
//3rd Edition
//isbn: 9780521880688
//further modifications done by Dimitar Lubomirov Mihaylov
struct Bessel {
	static const int NUSE1=7, NUSE2=8;
	static const double c1[NUSE1],c2[NUSE2];
	double xo,nuo;
	double jo,yo,jpo,ypo;
	double io,ko,ipo,kpo;
	double aio,bio,aipo,bipo;
	double sphjo,sphyo,sphjpo,sphypo;
	int sphno;

	Bessel() : xo(9.99e99), nuo(9.99e99), sphno(-9999) {}

	void besseljy(const double nu, const double x);
	void besselik(const double nu, const double x);

	double jnu(const double nu, const double x);
	double ynu(const double nu, const double x);
	double inu(const double nu, const double x);
	double knu(const double nu, const double x);

	void airy(const double x);
	double airy_ai(const double x);
	double airy_bi(const double x);

	void sphbes(const int n, const double x);
	double sphbesj(const int n, const double x);
	double sphbesy(const int n, const double x);

	inline double chebev(const double *c, const int m, const double x);
};



#endif

