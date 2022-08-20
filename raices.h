#ifndef RAICES_H
#define RAICES_H
#include <cmath>

namespace raices{
	struct resultado_raiz{
		/**
		*@brief Constructor por defecto de la structura resultado_raiz
		*/
		resultado_raiz() = default;
		/**
		*@brief Constructor parametrizado de la structura resultado_raiz
		*@param p_raiz primer parametro
		*@param p_erp segundo parametro
		*@param p_iter tercer parametro
		*/
		resultado_raiz(double p_raiz, double p_erp, int p_iter) : 
			raiz(p_raiz), 
			erp(p_erp), 
			iteraciones(p_iter){
			
		};
		double raiz = NAN;
		double erp = 100.0f;
		int iteraciones = 0;
		/**
		*@brief Metodo para determinar si se encontro el resultado de una ecuacion
		*@return si se encontro el resultado
		*/
		bool encontrada(){
			return !std::isnan(raiz);
		}
	};
};
#endif
