#ifndef RL_H
#define RL_H
#include <vector>
#include <cmath>
#include <string>

#include "gauss.h"
/* Regresion lineal simple */
using std::string;
using std::vector;

namespace regresion
{
	struct resultado_rl
	{
		double b0 = 0.0f;
		double b1 = 0.0f;
		double st = 0.0f;
		double sy = 0.0f;
		double sr = 0.0f;
		double syx = 0.0f;
		bool aceptable()
		{
			// si el error estadndar de aproximacion < dsv estandar
			return (syx < sy);
		}
		double r2;
		string ecuacion()
		{
			return "y = " + std::to_string(b1) + " x " + (string)((b0 > 0) ? "+ " : "- ") + std::to_string(fabs(b0));
		}
	};

	struct resultado_rl_potencia
	{
		double c;
		double a;
		double st;
		double sy;
		double sr;
		double syx;
		double r2;
		bool aceptable()
		{
			// si el error estadndar de aproximacion < dsv estandar
			return (syx < sy);
		}
		string ecuacion()
		{
			return "y= " + std::to_string(c) + "x^" + std::to_string(a);
		}
	};

	struct resultado_rl_exponencial
	{
		double c;
		double a;
		double st;
		double sy;
		double sr;
		double syx;
		double r2;
		bool aceptable()
		{
			// si el error estadndar de aproximacion < dsv estandar
			return (syx < sy);
		}
		string ecuacion()
		{
			return "y= " + std::to_string(c) + "e^(" + std::to_string(a) + "x)";
		}
	};

	class rl
	{
	public:
		static resultado_rl calcular(vector<double> x, vector<double> y)
		{
			resultado_rl r;
			double sum_xy = 0.0f;
			double sum_x = 0.0f;
			double sum_y = 0.0f;
			double sum_x2 = 0.0f;
			double xprom = 0.0f;
			double yprom = 0.0f;
			size_t n = x.size();
			if (n == 0)
			{
				return r;
			};
			for (size_t i = 0; i < n; i++)
			{
				sum_xy += x[i] * y[i];
				sum_x += x[i];
				sum_y += y[i];
				sum_x2 += x[i] * x[i];
			}
			yprom = sum_y / (double)n;
			xprom = sum_x / (double)n;
			r.b1 = (sum_xy - (yprom * sum_x)) / (sum_x2 - (xprom * sum_x));
			r.b0 = yprom - (r.b1 * xprom);
			for (size_t i = 0; i < n; i++)
			{
				r.st += pow(y[i] - yprom, 2);
				r.sr += pow(y[i] - (r.b1 * x[i] + r.b0), 2);
			}
			// Coeficientes de la recta
			// DEv. estandar
			r.sy = sqrt(r.st / (double(n - 1)));
			// Error de aproximacion
			r.syx = sqrt(r.sr / (double)(n - 2));
			// coeficiente de determinacion
			r.r2 = (r.st - r.sr) / r.st;

			return r;
		}
	};
	class rl_potencia
	{
	public:
		static resultado_rl_potencia calcular(vector<double> x, vector<double> y)
		{
			resultado_rl_potencia r;
			vector<double> X{x};
			vector<double> Y{y};
			// trnasformar los datos
			for (size_t i = 0; i < X.size(); i++)
			{
				X[i] = log10(X[i]);
				Y[i] = log10(Y[i]);
			}
			// realizar la regresion lineal sobre los datos tranformados
			resultado_rl r_rl = rl::calcular(X, Y);
			r.c = pow(10.0, r_rl.b0);
			r.a = r_rl.b1;
			r.st = r_rl.st;
			r.sy = r_rl.sy;
			r.syx = r_rl.syx;
			r.r2 = r_rl.r2;
			return r;
		}
	};
	class rl_exponencial
	{
	public:
		static resultado_rl_exponencial calcular(vector<double> x, vector<double> y)
		{
			resultado_rl_exponencial r;
			vector<double> Y{y};
			for (size_t i = 0; i < y.size(); i++)
			{
				Y[i] = log(Y[i]);
			}
			// realizar la regresion lineal sobre los datos transformados
			resultado_rl r_rl = rl::calcular(x, Y);
			r.c = exp(r_rl.b0);
			r.a = r_rl.b1;
			r.st = r_rl.st;
			r.sy = r_rl.sy;
			r.syx = r_rl.syx;
			r.r2 = r_rl.r2;
			return r;
		}
	};
	struct resultado_cuadratica
	{
		double a2;
		double a1;
		double a0;
		double st = 0.0f;
		double sy = 0.0f;
		double sr = 0.0f;
		double syx = 0.0f;
		bool aceptable()
		{
			// se considera aceptable si
			// el error estandar de aproximacion < desviacion estandar
			return (syx < sy);
		}
		double r2;
		string ecuacion()
		{
			return "y = " +
				   std::to_string(a2) + "x^2 + " + std::to_string(a1) + "x + " + std::to_string((a0));
		}
	};
	class cuadratica
	{
	public:
		static resultado_cuadratica cualcular(vector<double> x, vector<double> y)
		{
			resultado_cuadratica r;
			int n = x.size();
			double sx = 0.0f;
			double sx2 = 0.0f;
			double sx3 = 0.0f;
			double sx4 = 0.0f;
			double sy = 0.0f;
			double sxy = 0.0f;
			double sx2y = 0.0f;
			double y_prom = 0.0f;
			// Todo calcular los valores!!!!
			for (int i = 0; i < n; i++)
			{
				sx += x[i];
				sx2 += x[i] * x[i];
				sy += y[i];
				sx3 += x[i] * x[i] * x[i];
				sxy += x[i] * y[i];
				sx4 += x[i] * x[i] * x[i] * x[i];
				sx2y += (x[i] * x[i]) * y[i];
			}
			vector<vector<double>> matriz{
				{(double)n, sx, sx2, sy},
				{sx, sx2, sx3, sxy},
				{sx2, sx3, sx4, sx2y}};

			vector<double> coeficiente = gauss(matriz);
			r.a2 = coeficiente[2];
			r.a1 = coeficiente[1];
			r.a0 = coeficiente[0];

			y_prom = sy / (double)n;
			for (size_t i = 0; i < n; i++)
			{
				r.st += pow(y[i] - y_prom, 2);
				r.sr += pow(y[i] - r.a0 - r.a1 * x[i] - r.a2 * pow(x[i], 2), 2);
			}

			// Desviacion Estandar
			r.sy = sqrt(r.st / (double)(n - 1));

			// Error estandar de aproximacion
			r.syx = sqrt(r.sr / (double)(n - 3));

			// coeficiente de determinacion
			r.r2 = (r.st - r.sr) / r.st;
			return r;
		}
	};
};
#endif
