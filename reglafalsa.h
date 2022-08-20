#ifndef REGLAFALSA_H
#define REGLAFALSA_H

#include "tvi.h"
#include "precision.h"
#include <cmath>
#include <iomanip>
#include <functional>
#include "raices.h"
using raices::tvi;
using std::function;
using std::isnan;
//using precision::error_relativo;
using precision::es_cercano;
using precision::error_relativo_porcentual;
using raices::resultado_raiz;

namespace raices{
	class reglafalsa{
	public:
		/**
		*@brief metodo regla falsa que genera conjunto de aproximaciones a la solucion de una ecuacion
		*@param funcion f, primer parametro
		*@param xi segundo parametro
		*@param xs tercer paramentro 
		*@param erp cuarto parametro
		*@param maxIter quinto parametro
		*@return solucion con una aproximacion de la raiz en la tolerancia requeria o NAN
		*/
		static resultado_raiz calcular(function<double(double)> f, double xi, double xs, double erp, int maxIter){
			resultado_raiz resultado;
			//			double tolerancia = erp / 100.f;
			//			double cero = 0.000001f;
			//Xi = Xa
			//Xs = Xb
			int i = 1;
			double xant = xs - ((f(xs)*(xi - xs))/(f(xi)-f(xs))); //Xr 
			if(!tvi::existeRaiz(f, xi, xant)){
				xi = xant;
			}else{
				xs = xant;
			}
			while(i<=maxIter){
				double xnueva = xs - ((f(xs)*(xi - xs))/(f(xi)-f(xs)));
				if(isnan(xnueva)){
					break;
				}
				double er = error_relativo_porcentual(xnueva, xant);
				
				if(es_cercano(f(xnueva), (double)0.0f) || er < erp){
					resultado = resultado_raiz(xnueva,er,i);
					return resultado;
				}
				i++;
				if(!tvi::existeRaiz(f, xi, xnueva)){
					xi = xnueva;
				}else{
					xs = xnueva;
				}
				xant = xnueva;
				
			}
			return resultado;
		};
			
	};
};
#endif
