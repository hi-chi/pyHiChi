#pragma once
/*
This is an implicit solver of the equation of the Relativistic Electronic Spring (RES) model 
(for details see A. Gonoskov, Physics of Plasmas 25 (1), 013108 (2018)).
Note that the sign of the incidence angle theta here is the opposite to what is used in the 
article and the Ey_reflected is also having the opposite sign. The sign of Ey_reflected is 
chosen so that in the limit of high plasma density the reflected radiation tends to the 
incident one.

To run the solver from Python do the following steps:
	- initiate an instance of the solver (this allocates data for the computations
	and should be called once in case of series of computations that can vary
	in tersm of all orther settings (E_y(x), E_z(x), N(x), theta, etc.):
	X = res_solver(maxIndex_in, maxIndex_out, maxIndex_plasma, L_in, L_out, L_plasma)
	The first three parameters are the sizes of the grids used for sampling incoming 
	field, outgoing field and the density distribution at the plasma interface; the last 
	three parameters are the physical coordinate limits for this functions: x \in [-L_in, 0],
	x \in [-L_out, 0], x \in [0, L_plasma], respectively. Everything is assumed to be 
	specified in lab frame and in CGS units. Note that the density N_plasma must be 
	notably > 0 in all points (not equal or close to 0 at any point in the region X > 0).

	- set the grid values for the electric field components (Ey and Ez) and for the plasma
	density distribution by calling X.E_in_set(index, Ey, Ez) for all index values \in {0, ..., 
	maxIndex_in - 1} and X.N_plasma_set(index, N) for all index values \in {0, ..., 
	maxIndex_plasma - 1}. For x > L_plasma the solvers assumes N(x) = N[maxIndex_plasma - 1].

	- set the incidence angle (theta) by calling X.incidence_angle_set(value)

	- run the solver by calling X.compute()

	- read the values of the outgoing field by calling X.Ey_out_get(index), X.Ez_out_get(index),
	where index varies from 0(x = -L_out) to maxIndex_out-1 (x = 0; leading point)

The computation can be repeated using the same solver instance with new settings. 

Optional parameters to be configured after the initiation, the defaule value is given in 
square brackets:
	- time_step [0.001 * L_out / LightVelocity], time step
	- beta_bound [0.99], a value that is ought to prevent the singularity (close but < 1)
	- x_s_bound_factor [0.01], the factor in the estimate of the minimal accepted value
	of the sheet coordinate (acts also as a seed for the computation loop): 
	x_s_min = x_s_bound_factor (2 \pi |e|/\cos^2(\theta))*E_max/N_av, where E_max is the 
	estimated maximal value of field strength and N_av is the average plasma density.
	- max_time_factor [5], the factor that defines the limit for the time integration
	accoirding to time <= max_time_factor * L_in / LightVelocity

As an option one can initiate tracking of the sheet coordinate x_s by calling:
initiate_x_s_tracking(int maxIndex, double maxTime). This allocates memory and thus in case
of repeating computations with the same instance of solver this function should be called 
only once. To read the data one should use x_s_get_value(int index). 

Contacts: arkady.gonoskov@gmail.com
*/

#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//The following auxiliary function is a solver of cubic equation adopted from http://www.cplusplus.com/forum/beginner/234717/.
//Only real roots are of interest here.
int solvCubicEq(double a, double b, double c, double d, double &x1, double &x2, double &x3) // the function returns the number of real roots and place them sequentially in x1, x2, and x3 (some may be not used if the number of roots is smaller than 3).
{
	b /= a;
	c /= a;
	d /= a;
	double disc, q, r, dum1, s, t, term1, r13;
	q = (3.0*c - (b*b)) / 9.0;
	r = -(27.0*d) + b * (9.0*c - 2.0*(b*b));
	r /= 54.0;
	disc = q * q * q + r * r;
	term1 = (b / 3.0);
	if (disc > 0)
	{
		s = r + sqrt(disc);
		s = s < 0 ? -cbrt(-s) : cbrt(s);
		t = r - sqrt(disc);
		t = t < 0 ? -cbrt(-t) : cbrt(t);
		x1 = -term1 + s + t;
		return 1;
	}
	else if (disc == 0)
	{
		r13 = r < 0 ? -cbrt(-r) : cbrt(r);
		x1 = -term1 + 2.0*r13;
		x3 = x2 = -(r13 + term1);
		return 3;
	}
	else
	{
		q = -q;
		dum1 = q * q*q;
		dum1 = acos(r / sqrt(dum1));
		r13 = 2.0*sqrt(q);
		x1 = -term1 + r13 * cos(dum1 / 3.0);
		x2 = -term1 + r13 * cos((dum1 + 2.0*M_PI) / 3.0);
		x3 = -term1 + r13 * cos((dum1 + 4.0*M_PI) / 3.0);
		return 3;
	}
	return -1;
}

class res_solver
{
public:
	double *Ey_in, *Ez_in, *Ey_out, *Ez_out, *N_plasma; // lab frame, CGS units
	//The electric field compoents of the incoming ratioation are defined in [-L_in, 0] 
	//([-L_out, 0] for the outgoing radiation) and the radiation propagates in the 
	//positive X direction, i.e. X = 0 is the leading point.
	//The Y axis is in the incidence plane and is orientated so that at theta = Pi/4 coincides with the specular direction.
	//The density is defined at X > 0, note that N_plasma must be notably > 0 (not equal or close to 0 at any point in the region X > 0).
	double theta; // in radians
	int maxIndex_in, maxIndex_out, maxIndex_plasma; // the sizes of the grids used sampling incoming and outgoing fields, as well as the density distribution at the plasma interface
													// ! must be set on creation, must not be modified afterwards
	double L_in, L_out, L_plasma;
	double time_step;
	double beta_bound; // a value that is close but < 1 that is ought to prevent the singularity
	double max_time_factor;   // the time is advanced not longer than maxTimeFactor*L_in/LightVelocity 
							  //(if the outgoing pulse is not formed by this time, the plasma is not dense enough)
	double x_s_bound_factor; // a small value that is used as a factor in the estimate for the minimal allowed value of x_s 
							 // (prevents singularity)
	
	//auxilary diagnostics: tracking of x_s;
	bool x_s_tracking;
	double *x_s_value, *x_s_time; // data of x_s(time)
	int maxIndex_x_s;
	double maxTime_x_s;

private:
	double LightVelocity;
	double ElectronCharge;
	double *q_s_data;
	double coeff_q;
	double x_s_bound; // the minimal value for the coordinate of the electron sheet
	double E_max; // this is an estimate, not exact value! (lab frame)
public:
	res_solver(int maxIndex_in, int maxIndex_out, int maxIndex_plasma, double L_in, double L_out, double L_plasma) :
		maxIndex_in(maxIndex_in), maxIndex_out(maxIndex_out), maxIndex_plasma(maxIndex_plasma),
		L_in(L_in), L_out(L_out), L_plasma(L_plasma)
	{
		Ey_in = new double[maxIndex_in];
		Ez_in = new double[maxIndex_in];
		Ey_out = new double[maxIndex_out];
		Ez_out = new double[maxIndex_out];
		N_plasma = new double[maxIndex_plasma];
		q_s_data = new double[maxIndex_plasma];
		theta = 0; // default value
		LightVelocity = 2.99792458e+10;
		ElectronCharge = -4.80320427e-10;
		time_step = 0.001 * L_out / LightVelocity; // default value
		beta_bound = 0.99; // default value
		max_time_factor = 5; // default value
		x_s_bound_factor = 0.01; // default value
		x_s_tracking = false;
	}
	void initiate_x_s_tracking(int maxIndex, double maxTime)
	{
		maxIndex_x_s = maxIndex;
		maxTime_x_s = maxTime;
		x_s_value = new double[maxIndex_x_s];
		x_s_time = new double[maxIndex_x_s];
		for (int i = 0; i < maxIndex_x_s; i++) x_s_value[i] = 0;
		x_s_tracking = true;
	}
	double x_s_get(int index)
	{
		return x_s_value[index];
	}
	~res_solver()
	{
		delete[]Ey_in;
		delete[]Ez_in;
		delete[]Ey_out;
		delete[]Ez_out;
		delete[]N_plasma;
		delete[]q_s_data;
		if (x_s_tracking) delete[]x_s_value;
	}
	void E_in_set(int index_in, double Ey_value, double Ez_value)
	{
		Ey_in[index_in] = Ey_value;
		Ez_in[index_in] = Ez_value;
	}
	void N_plasma_set(int index_plasma, double N_plasma_value)
	{
		N_plasma[index_plasma] = N_plasma_value;
	}
	void incidence_angle_set(double value)
	{
		theta = value;
	}
	double Ey_out_get(int index_out)
	{
		return Ey_out[index_out];
	}
	double Ez_out_get(int index_out)
	{
		return Ez_out[index_out];
	}
	double X_in(int index_in)
	{
		return -L_in * (1 - index_in / double(maxIndex_in - 1));
	}
	double X_out(int index_out)
	{
		return -L_out * (1 - index_out / double(maxIndex_out - 1));
	}
	double X_plasma(int index_plasma)
	{
		return L_plasma * index_plasma / double(maxIndex_plasma);
	}
	int compute()
	{
		compute_q_s();
		compute_x_bound();

		double x_s, beta_x, Ey, Ez, Ey_out_, Ez_out_, X_out_;
		int out_index_previous = maxIndex_out - 1;
		double X_out_previous = 0, Ey_out_previous = 0, Ez_out_previous = 0;

		memset(Ey_out, 0, sizeof(double) * maxIndex_out);
		memset(Ez_out, 0, sizeof(double) * maxIndex_out);

		x_s = x_s_bound; // seeded value

		int x_s_tracking_last = 0; // last index value

		int counter = 0;
		for (double time = 0; ((time <= max_time_factor * L_in / LightVelocity) && (out_index_previous > 0));)
		{
			counter++;

			//perferct start
			beta_x = beta_x_(x_s, time); //explicite estimate
			double beta_next = beta_x_(x_s, time + time_step); //explicite estimate
			if ((beta_x < 0) && (beta_next > 0.001))
			{
				time += time_step * (-beta_x) / (beta_next - beta_x);
				beta_x = beta_x_(x_s, time);
			}
			
			beta_x = implicite_step(time, x_s, time_step);
			//cout << counter << ": x_s = " << x_s << ", beta_x = " << beta_x << ", time = " << time << endl;

			//compute the reflected field
			{
				double q_s = q_s_(x_s);
				incomingField(cos(theta)*(x_s - time * LightVelocity), Ey, Ez);
				double Ry = -sin(theta) - Ey / q_s;
				double Rz = -Ez / q_s;
				Ey_out_ = -q_s * (Ry*(1 - beta_x) / (1 + beta_x) + sin(theta));
				Ez_out_ = q_s * (Rz*(1 - beta_x) / (1 + beta_x));
				X_out_ = -cos(theta) * (x_s + LightVelocity * time);
				int out_index = int((maxIndex_out - 1)*(1 + X_out_ / L_out));
				if (out_index < 0) out_index = 0;

				if (out_index != out_index_previous)
				{
					for (int i = out_index_previous - 1; i >= out_index; i--)
					{
						double X = -L_out + i * L_out / double(maxIndex_out - 1);
						Ey_out[i] = (Ey_out_*(X_out_previous - X) + Ey_out_previous * (X - X_out_)) / (X_out_previous - X_out_);
						Ez_out[i] = (Ez_out_*(X_out_previous - X) + Ez_out_previous * (X - X_out_)) / (X_out_previous - X_out_);
					}
					Ey_out_previous = Ey_out_;
					Ez_out_previous = Ez_out_;
					X_out_previous = X_out_;
					out_index_previous = out_index;
				}
			}

			// tracking of x_s
			if (x_s_tracking)
			{
				int ind = int(maxIndex_x_s * time / maxTime_x_s);
				if (ind >= maxIndex_x_s) ind = maxIndex_x_s - 1;
				if (x_s_tracking_last < ind)
					for (int k = x_s_tracking_last + 1; k <= ind; k++)
						x_s_value[k] = (x_s_value[x_s_tracking_last] * (ind - k) + x_s * (k - x_s_tracking_last)) / double(ind - x_s_tracking_last);
				x_s_tracking_last = ind;
			}
		}
		return 0;
	}
private:
	double sqr(double x)
	{
		return x * x;
	}
	double implicite_step(double &time, double &x_s, double timeStep)
	{
		double beta_result;
		double beta_0 = beta_x_(x_s, time);
		double Ey_0, Ez_0;
		incomingField(cos(theta)*(x_s - (time + timeStep) * LightVelocity), Ey_0, Ez_0);
		double Ey_p, Ez_p;
		double dx_E = 0.001 * (L_in / double(maxIndex_in - 1)) / cos(theta);
		incomingField(cos(theta)*(x_s + dx_E - (time + timeStep) * LightVelocity), Ey_p, Ez_p);
		double Ey_1 = ((Ey_0 - Ey_p) / dx_E) * LightVelocity * timeStep;
		double Ez_1 = ((Ez_0 - Ez_p) / dx_E) * LightVelocity * timeStep;

		double q0 = q_s_(x_s);
		double dx_q = 0.001 * L_plasma / double(maxIndex_plasma - 1);
		double q1 = ((q_s_(x_s + dx_q) - q0) / dx_q) * LightVelocity * timeStep;

		double A = sqr(sin(theta) * q1) + sqr(Ey_1);
		double B = sqr(sin(theta)) * 2 * q0 * q1 + 2 * sin(theta) * q1 * Ey_0 + 2 * sin(theta) * q0 * Ey_1 + 2 * Ez_0 * Ez_1;
		double C = sqr(sin(theta) * q0) + 2 * sin(theta) * q0 * Ey_0 + sqr(Ey_0) + sqr(Ez_0);

		double a = A + sqr(q1);
		double b = B + 2 * q0 * q1 - A + sqr(q1) * beta_bound;
		double c = C + sqr(q0) - B + 2 * q0 * q1 * beta_bound;
		double d = -C + sqr(q0) * beta_bound;

		double x[3];

		int n = solvCubicEq(a, b, c, d, x[0], x[1], x[2]);
		if (n == 1) beta_result = x[0];

		if (n > 1)
		{
			// If there are more than one real root, choose the one that is closest to the explicite value
			double beta_explicite = beta_x_(x_s, time);
			double min_dif = 2;
			int i_min = 0;
			for (int i = 0; i < 3; i++)
				if (min_dif > abs(beta_explicite - x[i]))
				{
					min_dif = abs(beta_explicite - x[i]);
					i_min = i;
				}
			beta_result = x[i_min];
		}
		//verification and correction of the errors due to limitations of the numerical arithmetics
		{
			if (beta_result > beta_bound)beta_result = beta_bound;
			if (beta_result < -beta_bound)beta_result = -beta_bound;
		}

		x_s += timeStep * beta_result * LightVelocity;
		time += timeStep;
		if (x_s < x_s_bound)
		{
			x_s = x_s_bound;
			beta_result = 0;
		}
		return beta_result;
	}
	void compute_q_s()
	{
		coeff_q = 2 * M_PI * ElectronCharge / pow(cos(theta), 2); // (note that it is < 0)
		double q_s = 0;
		for (int i = 0; i <= maxIndex_plasma; i++)
		{
			q_s_data[i] = q_s;
			q_s += coeff_q * N_plasma[i] * L_plasma / double(maxIndex_plasma - 1);
		}
	}
	void compute_x_bound()
	{
		E_max = 0; // this is an estimate, not exact value! (lab frame)
		for (int i = 0; i < maxIndex_in; i++)
		{
			if (E_max < abs(Ey_in[i]))E_max = abs(Ey_in[i]);
			if (E_max < abs(Ez_in[i]))E_max = abs(Ez_in[i]);
		}
		double averageDensity = 0;
		for (int i = 0; i < maxIndex_plasma; i++)
			averageDensity += N_plasma[i] / double(maxIndex_plasma);

		coeff_q = 2 * M_PI * ElectronCharge / pow(cos(theta), 2);

		x_s_bound = - x_s_bound_factor *E_max / (coeff_q * averageDensity);
	}
	double q_s_(double x_s)
	{
		if (x_s < L_plasma)
			return LI(q_s_data, (maxIndex_plasma - 1) * x_s / L_plasma);
		return q_s_data[maxIndex_plasma - 1] + (x_s - L_plasma) * N_plasma[maxIndex_plasma - 1] * coeff_q;
	}
	double beta_x_(double x_s, double t)
	{
		double Ey, Ez;
		incomingField(cos(theta)*(x_s - t * LightVelocity), Ey, Ez);
		double Ry = -sin(theta) - Ey / q_s_(x_s);
		double Rz = -Ez / q_s_(x_s);
		return (Ry*Ry + Rz * Rz - beta_bound) / (Ry*Ry + Rz * Rz + 1);
	}
	double LI(double *array, double index) // linear interpolation
	{
		return (index - int(index))*array[int(index) + 1] + (1 - index + int(index))*array[int(index)];
	}
	double N(double X)// density, lab frame
	{
		if (X <= 0) return 0;
		if (X >= L_plasma) return N_plasma[maxIndex_plasma - 1];
		return LI(N_plasma, (maxIndex_plasma - 1)*X / L_plasma);
	}
	int incomingField(double X, double &Ey, double &Ez) // lab frame
	{
		Ey = 0;
		Ez = 0;
		if (X >= 0) return -1;
		if (X <= -L_in) return -1;
		Ey = LI(Ey_in, (maxIndex_in - 1) * (1 + X / L_in));
		Ez = LI(Ez_in, (maxIndex_in - 1) * (1 + X / L_in));
		return 0;
	}
};



