#ifndef MULLER_H
#define MULLER_H
#include <functional>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "precision.h"
#include "raices.h"
using std::function;
using std::cout;
using std::isnan;
using std::endl;
using std::isnan;
using precision::error_relativo_porcentual;
using raices::resultado_raiz;
namespace raices{
	class muller{
	public:
		/**
		*@brief Calcula una solucion dada una funcion, utiliza una aproximaci�n cuadr�tica en lugar de una aproximaci�n lineal.
		*@param funcion f, primer parametro
		*@param x0 segundo parametro
		*@param x1 tercer parametro
		*@param x2 cuarto parametro
		*@param erp quinto parametro
		*@param maxIter sexto parametro
		*@return solucion raiz aproximada o resultado de tipo NAN
		*/
		static resultado_raiz calcular (function<double(double)> f, double x0, double x1, double x2, double erp, int maxIter){
			resultado_raiz resultado;
			
			double h1 = x1 - x0;
			double h2 = x2 - x1;
			double delta1 = (f(x1)-f(x0)) / h1;
			double delta2 = (f(x2)-f(x1)) / h2;
			double a = ( delta2 - delta1 ) / (h2 + h1);
			int i = 2;			
			while(i <= maxIter){
				double b = delta2  + (h2 * a);
				double D = sqrt( (b*b) - (4 * f(x2) * a) ); 
				if(isnan(D)){
					break;
				}
				double E;
				if(fabs(b - D) < fabs(b + D)){
					E = b + D;
				}else{
					E = b - D;
				}
				double h = (-2*f(x2)) / E;
				double p = x2+h;
				double er = error_relativo_porcentual(p,x2); 
				if(er < erp){
					resultado = resultado_raiz(p, er,i);
					break;
				}
				//Siguiente aproximacion
				x0 = x1;
				x1 = x2;
				x2 = p;				
				h1 = x1 - x0;
				h2 = x2 - x1;
				delta1 = (f(x1)-f(x0)) / h1;
				delta2 = (f(x2)-f(x1)) / h2;
				a = ( delta2 - delta1 ) / (h2 + h1);
				i++;
			}
			return resultado;
		};
	};
};
#endif
