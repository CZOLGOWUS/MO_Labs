#include "lab10.hpp"


namespace lab10
{
	using namespace std;

#define func(t,yt,dyt) dyt + (10.0 * t * t + 20.0) / (t * t + 1.0) * (yt - 1.0)
#define analitic(t) 1.0 - exp(-10.0 * (t + atan(t)))

	// - might be wrong check later
	//metoda bezposrednia eulera
	double fMBE(double t, double yt)
	{
		return -((10.0 * t * t + 20.0) / (t * t + 1.0)) * (yt - 1.0);
	}


	double MBE(double dt,double tMax, double (*f)(double, double))
	{
		double yt = 0.0;
		for (double t = 0.0; t < tMax; t += dt)
			yt += dt*f(t, yt);

		return yt;
	}

	//metoda Posrednia Eulera

	double fMPE(double t, double h)
	{
		return (10.0 * (t + h) * (t + h) + 20.0) / ((t + h) * (t + h) + 1.0);
	}


	double MPE(double dt, double tMax ,double (*f)(double,double))
	{
		double yt = 0.0;
		for (double t = 0.0; t < tMax; t += dt)
			yt = (yt + dt * f(t, dt)) / (1 + dt * f(t, dt));
		return yt;
	}


	//metoda trapezów
	double f1MT(double t)
	{
		return (10.0 * t * t + 20.0) / (t * t + 1.0);
	}

	double f2MT(double t,double h)
	{
		return (10.0 * (t + h) * (t + h) + 20.0) / ((t + h) * (t + h) + 1.0);
	}

	double MT(double dt, double tMax)
	{
		double yt = 0.0;
		for (double t = 0.0; t < tMax; t += dt)
			yt = (-(dt * 0.5) * (f1MT(t) * (yt - 1.0) - f2MT(t, dt)) + yt) / (1.0 + (dt * 0.5) * f2MT(t, dt));
		return yt;
	}


	//calcualting errors
	double* calculateErrors(double dt, double numberOfIterations)
	{
		
		double* error = new double[3];

		for (int i = 0; i < 3; i++)
			error[i] = 0.0;


		//MBE
		double t = dt;
		double yt = 0.0;
		double err = 0.0;

		for (int i = 0; i < numberOfIterations; i++)
		{
			yt = yt + dt*fMBE(t, dt);
			err = abs(analitic(t) - yt);

			if (err > error[0])
				error[0] = err;

			t += dt;
		}


		//MPE
		t = dt;
		yt = 0.0;
		err = 0.0;

		for (int i = 0; i < numberOfIterations; i++)
		{
			yt = (yt + dt * fMPE(t, dt)) / (1 + dt * fMPE(t, dt));
			err = abs(analitic(t) - yt);

			if (err > error[1])
				error[1] = err;

			t += dt;
		}


		//MT
		t = dt;
		yt = 0.0;
		err = 0.0;

		for (int i = 0; i < numberOfIterations; i++)
		{
			yt = ( -(dt * 0.5) * (f1MT(t) * (yt - 1.0) - f2MT(t, dt)) + yt) / (1.0 + (dt * 0.5) * f2MT(t, dt));
			err = abs(analitic(t) - yt);

			if (err > error[2])
				error[2] = err;

			t += dt;
		}


		return error;

	}


#define tMax 3.0
#define iterations 500


	void test()
	{

		fstream errors, stable, notStable;

		errors.open("C:\\Users\\CZOLG\\Desktop\\lab10\\errors.txt",std::ios::out);
		stable.open("C:\\Users\\CZOLG\\Desktop\\lab10\\stable.txt",std::ios::out);
		notStable.open("C:\\Users\\CZOLG\\Desktop\\lab10\\notStable.txt",std::ios::out);
		
		if (errors.bad() || stable.bad() || notStable.bad()) throw exception("file bad");

		//errors
		double dt = 0.1;
		double* re = new double[3];
		for (double idt = 0.1; idt> 1e-16; idt= idt*  0.6)
		{
			re = calculateErrors(idt,iterations);
			errors << scientific << log10(idt) << " " << log10(re[0]) << " " << log10(re[1]) << " " << log10(re[2])<<endl;
		}
		errors.close();
		delete[] re;

		//stable dt is (0;1/10>
		dt = 0.01;
		for (double idt= 0; idt < 5.25; idt += 0.01)
		{
			stable << scientific << idt <<" "<< analitic(idt) << " " << MBE(dt,idt,fMBE) << " " << MPE(dt,idt,fMPE) << " " << MT(dt,idt)<<endl;
		}
		stable.close();

		//not stable (possible only for MBE and dt > 1/10)
		for (double idt = 0; idt < 5.0; idt += 0.1)
		{
			notStable << scientific << idt << " " << MBE(0.15, idt,fMBE)<<endl;
		}

		notStable.close();

	}

}