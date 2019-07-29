#pragma once

#include "math.h"

#define MATH_PI 3.14159265358979323846
#define MATH_PI_2 1.57079632679489661923
#define MATH_SQRT2 1.41421356237309504880
#define MATH_SQRT3 1.73205080756887729353
#define MACH_SQRT_DBL_EPSILON 1.4901161193847656e-08
#define MACH_LOG_DBL_MIN (-7.0839641853226408e+02)
#define MACH_LOG_FLT_MIN (-8.7336544750553108e+01)

static inline double hornerSchem(const double* a, const double x, const int n)
{
	int i;
	double f = a[n];
	for (i = n - 1; i >= 0; i--)
	{
		f = x * f + a[i];
	}
	return f;
}

static inline float hornerSchemF(const float* a, const float x, const int n)
{
	int i;
	float f = a[n];
	for (i = n - 1; i >= 0; i--)
	{
		f = x * f + a[i];
	}
	return f;
}

static double synchrotron11_coef[12] = {
	2.149528241534478485,
	0.4030365452877174690,
	0.6045548179314911165e-1,
	0.2361542257554632558e-2,
	0.4600406994651064777e-4,
	5.421908264461817201e-7,
	4.293365805443271182e-9,
	2.443775407411153507e-11,
	1.048462185580334617e-13,
	3.517673757866801084e-16,
	9.276203335336769246e-19,
	2.487518452145833156e-21
}; //6e-17
static double synchrotron12_coef[11] = {
	0.1423934037247282815,
	0.7628218056681864558e-2,
	0.1820370218072000925e-3,
	2.500508541288898197e-6,
	2.240804989687901660e-8,
	1.415245229181725731e-10,
	6.642707404690081889e-13,
	2.408835566075670457e-15,
	6.957217376253642924e-18,
	1.616072251127992819e-20,
	3.693431297197193414e-23
}; //3e-18   3e-15
static double synchrotron1a_coef[15] = {
	1.253314151318345289,
	0.95739073440816383,
	-1.22694931936470073,
	2.94131175393275699,
	-10.0160956662322282,
	42.2597994853969905,
	-195.974329891525252,
	893.719183427841709,
	-3662.99838165288478,
	12610.7593887738478,
	-34527.7964181278657,
	71190.6177419997158,
	-102965.357425859132,
	92643.3250058469047,
	-38888.8222805020856
}; //6e-15 4to20 >35 f=1e-16

static double synchrotron21_coef[12] = {
	1.07476412076723916,
	0.806073090575432292,
	0.755693522414381793e-1,
	0.269890543720422123e-2,
	0.506044769415638960e-4,
	5.83897812999689049e-7,
	4.56170118337595284e-9,
	2.57239500863771439e-11,
	1.09612063108607483e-13,
	3.65833441841716940e-16,
	9.60863798895216726e-19,
	2.56653370498025160e-21
}; //1e-16
static double synchrotron22_coef[12] = {
	1.26571914421980694,
	0.189857871632971079,
	0.889958773279540523e-2,
	0.202263357563665276e-3,
	2.70888425300440443e-6,
	2.39019199080360306e-8,
	1.49386992984264251e-10,
	6.95903096426353705e-13,
	2.50917148507849502e-15,
	7.21643605478982684e-18,
	1.66577631618849652e-20,
	3.85373079376373116e-23
}; //1.5e-18   4e-15
static double synchrotron2a_coef[13] = {
	1.25331413731553056,
	0.121849985551466331,
	-0.550017255754534718e-1,
	0.532189220291644753e-1,
	-0.785192304361570963e-1,
	0.155075556066478267,
	-0.377291378427079637,
	1.03064346511217399,
	-2.79297293304295070,
	6.54033276198134163,
	-11.5393857144838457,
	13.0543592851241271,
	-6.94238421837777902
}; //6e-15

double synchrotron_1(const double x)
{
	if (x < 0.0) {
		return 0.0;
	}
	else if (x < 0.001) {
		double px = pow(x, 1.0 / 3.0);
		return -1.8137993642342178506 * x +
			(2.1495282415344786368 + (0.40303654528771474439 - 0.14239340372472828097 * x * px) * x * x) * px;
	}
	else if (x <= 4.0) {
		const double c0 = MATH_PI / MATH_SQRT3;
		const double px = pow(x, 1.0 / 3.0);
		const double px11 = pow(px, 11);
		const double y = x * x;
		return px * hornerSchem(synchrotron11_coef, y, 11) - px11 * hornerSchem(synchrotron12_coef, y, 10) - c0 * x;
	}
	else if (x < 20.0) {
		const double y = 1.0 / x;
		return sqrt(x) * hornerSchem(synchrotron1a_coef, y, 14) * exp(-x);
	}
	else if (x < -8.0 * MACH_LOG_DBL_MIN / 7.0) {
		const double y = 1.0 / x;
		return sqrt(x) * hornerSchem(synchrotron1a_coef, y, 14) * exp(-x);
	}
	else {
		return 0.0;
	}
}
double synchrotron_2(const double x)
{
	if (x < 0.0) {
		return 0.0;
	}
	else if (x < 0.001) {
		double px = pow(x, 1.0 / 3.0);
		return (1.0747641207672393184 + (-1.2657191442198069419 + (0.80607309057542948877 - 0.18985787163297104129 * px * x)
			* px * px) * px * x) * px;
	}
	else if (x <= 4.0) {
		const double px = pow(x, 1.0 / 3.0);
		const double px5 = pow(px, 5);
		const double y = x * x;
		return px * hornerSchem(synchrotron21_coef, y, 11) - px5 * hornerSchem(synchrotron22_coef, y, 11);
	}
	else if (x < -8.0 * MACH_LOG_DBL_MIN / 7.0) {
		const double y = 1.0 / x;
		return sqrt(x) * exp(-x) * hornerSchem(synchrotron2a_coef, y, 12);
	}
	else {
		return 0.0;
	}
}