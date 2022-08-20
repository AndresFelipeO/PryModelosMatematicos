#ifndef BISECCION_H
#define BISECCION_H
#pragma once

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
	class biseccion{
	public:
		/**
		*@brief Calcula una solucion dada una funcion para la aproximacion de una raiz
		*@param funcion f, primer parametro
		*@param xa segundo parametro
		*@param xb tercer parametro
		*@param erp cuarto parametro
		*@param maxIter quinto parametro
		*@return solucion raiz aproximada o resultado de tipo NAN
		*/
		static resultado_raiz calcular(function<double(double)> f, double xa, double xb, double erp, int maxIter){
			resultado_raiz resultado;
			//			double tolerancia = erp / 100.f;
			//			double cero = 0.000001f;			
			int i = 1;
			double xant = (xa + xb) / 2.0f;
			if(!tvi::existeRaiz(f, xa, xant)){
				xa = xant;
			}else{
				xb = xant;
			}
			while(i<=maxIter){
				//Aproximacion
				double xnueva = (xa + xb) / 2.0f;
				if(isnan(xnueva)){
					break;
				}
				double er = error_relativo_porcentual(xnueva, xant);
				
				if(es_cercano(f(xnueva), (double)0.0f) || er < erp){
					resultado = resultado_raiz(xnueva,er,i);
					return resultado;
				}
				i++;
				if(!tvi::existeRaiz(f, xa, xnueva)){
					//No existe raiz, tomar xa = xnueva
					xa = xnueva;
				}else{
					//Si existe raiz, tomar xb = xnueva
					xb = xnueva;
				}
				xant = xnueva;
			}
			//No se encontro la raiz
			return resultado;
		}
			
	};
};
#endif
