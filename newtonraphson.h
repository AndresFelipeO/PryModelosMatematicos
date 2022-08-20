#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

#include <functional>
#include <cmath>
#include <iostream>
#include "precision.h"
#include "raices.h"

using std::function;
using std::isnan;
using std::cout;
using std::endl;
using precision::error_relativo_porcentual;
using raices::resultado_raiz;

namespace raices{
	class newton_raphson{
	public:
		/**
		*@brief Encuentra la funcion fx dada una aproximacion inical
		*@param funcion f, primer parametro
		*@param funcion df segundo parametro
		*@param po tercer parametro
		*@param erp cuarto parametro
		*@param maxIter quinto parametro
		*@return solucion aproximada p o resultado NAN
		*/
		static resultado_raiz calcular (function<double(double)> f, function<double(double)> df, double p0, double erp, int maxIter){
			resultado_raiz resultado;
			int i = 1;
			while(i <= maxIter){
				double p = p0 - ( f(p0) / df(p0) );
				//Validar la aproximacion
				if(isnan(p)){
					break;
				}
				cout << "p = " << p << endl;
				double er = error_relativo_porcentual(p, p0);
				if(er < erp){
					resultado = resultado_raiz(p, er,i);
					break;
				}else{
					i++;
					p0 = p;
				}
			}
			return resultado;
		};
		/**
		*@brief Encuentra la funcion fx dada una aproximacion inical de raices multiples
		*@param funcion f, primer parametro
		*@param funcion df segundo parametro
		*@param po tercer parametro
		*@param erp cuarto parametro
		*@param maxIter quinto parametro
		*@return solucion aproximada p o resultado NAN
		*/
		static resultado_raiz calcular_generalizado(function<double(double)> f, function<double(double)> df, function<double(double)> ddf, double p0, double erp, int maxIter){
			resultado_raiz resultado;
			int i = 1;
			while(i <= maxIter){
				double p = p0 - ( (f(p0) * df(p0)) / ((df(p0)*df(p0)) - (f(p0)*ddf(p0))));
				//Validar la aproximacion
				if(isnan(p)){
					break;
				}
				cout << "p = " << p << endl;
				double er = error_relativo_porcentual(p, p0);
				if(er < erp){
					resultado = resultado_raiz(p,er,i);
					break;
				}else{
					i++;
					p0 = p;
				}
			}
			return resultado;
		};
	};
};
#endif
