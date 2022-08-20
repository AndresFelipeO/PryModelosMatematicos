#ifndef SECANTE_H
#define SECANTE_H	

#define SECANTE_H
#include <functional>
#include <cmath>
#include "precision.h"
#include "raices.h"

using std::function;
using std::cout;
using std::isnan;
using std::endl;
using precision::error_relativo_porcentual;
using raices::resultado_raiz;

namespace raices{
	class secante{
	public:
		/**
		*@brief Proceso iteratico del metodo secante
		*@param funcion f, primer parametro
		*@param xo segundo parametro
		*@param x2 tercer parametro
		*@param erp cuarto parametro
		*@param maxIter quinto parametro
		*@return solucion aproximada p o resultado NAN
		*/
		static resultado_raiz calcular (function<double(double)> f, double x0, double x1, double erp, int maxIter){
			resultado_raiz resultado;
			int i = 1;
			while(i <= maxIter){
				double x2 = x1 - (( f(x1) * (x0 - x1)) / (f(x0) - f(x1)) );
				//validar si la aproximacion es valida
				if(isnan(x2)){
					break;
				}
				double er = error_relativo_porcentual(x2,x1); 
				if(er < erp){
					resultado = resultado_raiz(x2,er,i);
					break;
				}else{
					i++;
					x0 = x1;
					x1 = x2;
				}
			}
			return resultado;
		};
	};
};
#endif
