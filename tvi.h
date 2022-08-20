#ifndef TVI_H
#define TVI_H
#include <functional>
using std:: function;
namespace raices{
	class tvi{
	public:
		/**
		*@brief Calcula si existe almenos una raiz en dos puntos
		*@param f funcion, primer parametro
		*@param a segundo parametro
		*@param b tercer parametro
		*@return verdadero si existe almenos una raiz o falso si no existe
		*/
		static bool existeRaiz(function<double(double)> f, double a, double b){
			//aplicar la definicion del tvi pdf
			if(f(a) * f(b)<0.0f){
				return true;		
			}else{
				return false;
			}		
		}
		/**
		*@brief Calcula si existe una unica raiz
		*@param funcion f, primer parametro
		*@param a segundo parametro
		*@param b tercer parametro
		*@param n cuarto parametro
		*@return verdadero si tiene una unica raiz o falso si no.
		*/
		bool existeUnicaRaiz(function<double(double)> df,double a, double b,int n){
			double paso=(b-a)/(double)n;
			double r=df(a);
			double x=a+paso;
			for(int i=1;i<=n;i++){
				double r_p=df(x);
				if(r_p*r<0.0f)
					return false;
				x+=paso;
			}
			return true;
		}
	};
};
#endif
