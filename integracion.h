#ifndef INTEGRACION_H
#define INTEGRACION_H
#define _USE_MATH_DEFINES
#include <cmath>
#include <functional>
#include <vector>
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::function;
using std::vector;

namespace integracion
{
	void imprimir_tabla(vector<double> x, vector<double> y);
	struct resultado_romber
	{
		double integral;
		double error;
	};
	class integral
	{
	public:
		/**
		*@brief Calcula la precision de k .
		*@param double resultado, primer parametro
		*@return la cantidad de cibras que tuvo la precision
		*/
		static int precisionIntegral(double resultado)
		{
			int cont = 0;
			while (resultado < 1.0f)
			{
				cont++;
				resultado *= 10.0f;
			}
			return cont - 1;
		}
		/**
		*@brief Calcula una integral definida por el metodo del trapecio.
		*@param funcion f, primer parametro
		*@param limite_inferior segundo parametro
		*@param limite_superior tercer parametro
		*@param segmentos cuarto parametro
		*@return devuelve el resultado de la integral
		*/
		static double trapecio(function<double(double)> f, double limite_inferior, double limite_superior, int segmentos)
		{
			double result = 0.0f;
			vector<double> x;
			vector<double> y;
			// Generar tabla
			generar_tabla(f, limite_inferior, limite_superior, segmentos, x, y);
			// imprimir los datos
			//imprimir_tabla(x, y);
			// calcular la integral de la tabla de datos
			result = trapecio(x, y);
			return result;
		}
		/**
		*@brief Calcula las operaciones del metodo del trapecio.
		*@param x, primer parametro
		*@param y segundo parametro
		*@return devuelve el resultado de la integral
		*/
		static double trapecio(vector<double> x, vector<double> y)
		{
			double result = 0.0f;
			int n = x.size() - 1;
			double fx0 = y[0];
			double fxn = y[n];
			double a = x[0];
			double b = x[n];

			double sum = 0.0f;

			for (int i = 1; i < n; i++)
			{
				sum += y[i];
			}
			result = ((b - a) * (fx0 + 2 * sum + fxn)) / (2.0 * (double)n);
			return result;
		}
		
		/**
		*@brief Calcula el error estimado de la integral definida por el metodo del trapecio.
		*@param funcion df2, primer parametro
		*@param a segundo parametro
		*@param b tercer parametro
		*@return devuelve el error estimado de la integral
		*/
		static double trapecioErrorEstimadoNoPolinomica(function<double(double)> df2, double a, double b)
		{
			double result = 0.0f;
			double h = b - a;
			double fmayor = mayorValorFuncion(df2, a, b, 100);
			result = (pow(h, 3) / 12) * fabs(fmayor);
			return fabs(result);
		}
		
		/**
		*@brief ajusta el resultado de la integral definida por el metodo del trapecio.
		*@param resultado integral, primer parametro
		*@param funcion f, segundo parametro parametro
		*@param a tercer parametro
		*@param b cuarto parametro
		*@return devuelve la integral ajustada
		*/
		static double trapecioAjusteIntegralNoPolinomica(double resultIntegral, function<double(double)> df2, double a, double b)
		{
			double result = 0.0f;
			double h = b - a;
			double fmayor = mayorValorFuncion(df2, a, b, 100);
			result = -(pow(h, 3) / 12) * fmayor;
			return resultIntegral + result;
		}
		/**
		*@brief Calcula una integral definida por el metodo de Simpson 1/3.
		*@param funcion f, primer parametro
		*@param limite_inferior segundo parametro
		*@param limite_superior tercer parametro
		*@param segmentos cuarto parametro
		*@return devuelve el resultado de la integral
		*/
		static double simpson13(function<double(double)> f, double limite_inferior, double limite_superior, int segmentos)
		{
			double result = 0.0f;
			vector<double> x;
			vector<double> y;

			// Generar la tabla de datos
			generar_tabla(f, limite_inferior, limite_superior, segmentos, x, y);
			// Integrar la tabla de datos
			//imprimir_tabla(x, y);
			result = simpson13(x, y, ((limite_superior - limite_inferior) / (double)segmentos));
			return result;
		}
		/**
		*@brief Calcula las operaciones del metodo de Simpson 1/3.
		*@param x, primer parametro
		*@param y segundo parametro
		*@param h tercer parametro
		*@return devuelve el resultado de la integral
		*/
		static double simpson13(vector<double> x, vector<double> y, double h)
		{
			double result = 0.0f;
			int n = x.size() - 1;
			double fx0 = y[0];
			double fxn = y[n];
			double sumXi = 0.0f;
			double sumXj = 0.0f;
			for (int i = 1; i < n; i += 2)
				sumXi += y[i];
			for (int j = 2; j < n - 1; j += 2)
				sumXj += y[j];
			result = (h / 3) * (fx0 + (4 * sumXi) + (2 * sumXj) + fxn);
			return result;
		}
		/**
		*@brief Calcula el error estimado de la integral definida no polinomica por el metodo Simpson 1/3.
		*@param funcion df4, primer parametro
		*@param a segundo parametro
		*@param b tercer parametro
		*@param n cuarto parametro
		*@return devuelve el error estimado de la integral
		*/
		static double simpson13ErrorEstimadoNoPolinomica(function<double(double)> df4, double a, double b, int n)
		{
			double result = 0.0f;
			double h = (b - a) / n;
			double fmayor = mayorValorFuncion(df4, a, b, 100);
			result = (pow(h, 5) / 90) * fabs(fmayor);
			return fabs(result);
		}
		/**
		*@brief Calcula la cantidad de segmentos que necesita el metodo Simpson 1/3.
		*@param funcion df4, primer parametro
		*@param a segundo parametro
		*@param b tercer parametro
		*@return devuelve la cantidad de segmentos que necesita la integral
		*/
		static int CantidadSegmentosSimpson13(function<double(double)> df4, double a, double b)
		{
			double fmayor = mayorValorFuncion(df4, a, b, 100);
			double operacion = ((pow(b - a, 5) / (90 * 1 * pow(10, -7))) * fabs(fmayor));
			double n = pow(operacion, 1.0 / 5.0);
			int numero=round(n);
			while (numero%2!=0)
				numero++;
			return numero;
		}
		/**
		*@brief ajusta el resultado de la integral definida por el metodo Simpson 1/3
		*@param resultado integral, primer parametro
		*@param funcion df4, segundo parametro parametro
		*@param a tercer parametro
		*@param b cuarto parametro
		*@param n quinto parametro
		*@return devuelve la integral ajustada
		*/
		static double simpson13AjusteIntegralNoPolinomica(double resultIntegral, function<double(double)> df4, double a, double b, int n)
		{
			double result = 0.0f;
			double h = (b - a) / n;
			double fmayor = mayorValorFuncion(df4, a, b, 100);
			result = -(pow(h, 5) / 90) * fmayor;
			return resultIntegral + result;
		}
		
		/**
		*@brief Calcula el error estimado de la integral definida polinomica por el metodo Simpson 1/3.
		*@param funcion df4, primer parametro
		*@param a segundo parametro
		*@param b tercer parametro
		*@param n cuarto parametro
		*@return devuelve el error estimado de la integral
		*/
		static double simpson13ErrorEstimadoPolinomica(function<double(double)> df4, double a, double b, int n)
		{
			double result = 0.0f;
			double h = (b - a) / n;
			result = -(pow(h, 4)) * ((b - a) / 180) * df4((a + b) / 2);
			return result;
		}
		
		/**
		*@brief Calcula una integral definida por el metodo de Simpson 3/8.
		*@param funcion f, primer parametro
		*@param limite_inferior segundo parametro
		*@param limite_superior tercer parametro
		*@param segmentos cuarto parametro
		*@return devuelve el resultado de la integral
		*/
		static double simpson38(function<double(double)> f, double limite_inferior, double limite_superior, int segmentos)
		{
			double result = 0.0f;
			vector<double> x;
			vector<double> y;
			// Generar la tabla de datos
			generar_tabla(f, limite_inferior, limite_superior, segmentos, x, y);
			// Integrar la tabla de datos
			//imprimir_tabla(x, y);
			result = simpson38(x, y, ((limite_superior - limite_inferior) / (double)segmentos));
			return result;
		}
		/**
		*@brief Calcula las operaciones del metodo de Simpson 3/8.
		*@param x, primer parametro
		*@param y segundo parametro
		*@param h tercer parametro
		*@return devuelve el resultado de la integral
		*/
		static double simpson38(vector<double> x, vector<double> y, double h)
		{
			double result = 0.0f;
			int n = x.size() - 1;
			double fx0 = y[0];
			double fxn = y[n];
			double sumXi = 0.0f;
			double sumXj = 0.0f;
			double sumXk = 0.0f;
			for (int i = 1; i < n - 1; i += 3)
				sumXi += y[i];
			for (int j = 2; j < n; j += 3)
				sumXj += y[j];
			for (int k = 3; k < n - 2; k += 3)
				sumXk += y[k];
			result = ((3 * h) / 8) * (fx0 + (3 * sumXi) + (3 * sumXj) + (2 * sumXk) + fxn);
			return result;
		}
		/**
		*@brief Calcula la cantidad de segmentos que necesita el metodo Simpson 3/8.
		*@param funcion df4, primer parametro
		*@param a segundo parametro
		*@param b tercer parametro
		*@return devuelve la cantidad de segmentos que necesita la integral
		*/
		static int CantidadSegmentosSimpson38(function<double(double)> df4, double a, double b)
		{
			double fmayor = mayorValorFuncion(df4, a, b, 100);
			double operacion = ((3*pow(b - a, 5) / (80 * 1 * pow(10, -7))) * fabs(fmayor));
			double n = pow(operacion, 1.0 / 5.0);
			int numero=round(n);
			while (numero%2==0 || numero%3!=0)
				numero++;
			return numero;
		}
		/**
		*@brief Calcula el error estimado de la integral definida no polinomica por el metodo Simpson 3/8.
		*@param funcion df4, primer parametro
		*@param a segundo parametro
		*@param b tercer parametro
		*@param n cuarto parametro
		*@return devuelve el error estimado de la integral
		*/
		static double simpson38ErrorEstimadoNoPolinomica(function<double(double)> df4, double a, double b, int n)
		{
			double result = 0.0f;
			double h = (b - a) / n;
			double fmayor = mayorValorFuncion(df4, a, b, 100);
			result = 3 * (pow(h, 5) / 80) * fabs(fmayor);
			return fabs(result);
		}
		/**
		*@brief ajusta el resultado de la integral definida por el metodo Simpson 3/8
		*@param resultado integral, primer parametro
		*@param funcion df4, segundo parametro parametro
		*@param a tercer parametro
		*@param b cuarto parametro
		*@param n quinto parametro
		*@return devuelve la integral ajustada
		*/
		static double simpson38AjusteIntegralNoPolinomica(double resultIntegral, function<double(double)> df4, double a, double b, int n)
		{
			double result = 0.0f;
			double h = (b - a) / n;

			double fmayor = mayorValorFuncion(df4, a, b, 100);
			result = -(3 * (pow(h, 5) / 80)) * fmayor;
			return resultIntegral + result;
		}

	private:
		/**
		*@brief Genera una tabla de la funcion 
		*@param funcion f, primer parametro
		*@param limite_inferior, segundo parametro parametro
		*@param limite_superior tercer parametro
		*@param segmentos cuarto parametro
		*@param x quinto parametro
		*@param y sexto parametro
		*/
		static void generar_tabla(function<double(double)> f, double limite_inferior, double limite_superior, int segmentos, vector<double> &x, vector<double> &y)
		{
			x.clear();
			y.clear();
			x.resize(segmentos + 1);
			y.resize(segmentos + 1);
			double h = (limite_superior - limite_inferior) / (double)segmentos;

			// generar la tabla de datos x,y
			double xi = limite_inferior;
			for (int i = 0; i <= segmentos; i++)
			{
				x[i] = xi;
				y[i] = f(xi);
				xi += h;
			}
		}
		/**
		*@brief Encuentra el mayor intervalo de una funcion
		*@param funcion f, primer parametro
		*@param limite_inferior, segundo parametro parametro
		*@param limite_superior tercer parametro
		*@param segmentos cuarto parametro
		*/
		static double mayorValorFuncion(function<double(double)> f, double limite_inferior, double limite_superior, int segmentos)
		{
			vector<double> x1;
			vector<double> y2;
			x1.clear();
			y2.clear();
			x1.resize(segmentos + 1);
			y2.resize(segmentos + 1);
			double h = (limite_superior - limite_inferior) / (double)segmentos;
			// generar la tabla de datos x,y
			double mayor = 0;
			double xi = limite_inferior;
			for (int i = 0; i <= segmentos; i++)
			{
				x1[i] = xi;
				y2[i] = f(xi);
				xi += h;
				if (fabs(y2[i]) > fabs(mayor))
					mayor = y2[i];
			}
			return mayor;
		}
	};
	/**
	*@brief Imprime los intervalos de una funcion
	*@param x, primer parametro
	*@param y, segundo parametro parametro
	*/
	void imprimir_tabla(vector<double> x, vector<double> y)
	{
		int n = x.size();
		cout << "-------------------" << endl;
		cout << "X" << std::setw(11) << "Y" << endl;
		cout << "-------------------" << endl;

		for (int i = 0; i < n; i++)
		{
			cout << x[i] << std::setw(11) << y[i] << endl;
		}
		cout << "-------------------" << endl;
	}
};
#endif
