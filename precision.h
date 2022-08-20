#ifndef PRECISION_H
#define PRECISION_H

#include <cmath>
#include <cfloat>
#include <limits>


namespace precision{
	
	
	template<typename T>
	/**
	* @brief calcula el error absoluto entre dos valores
	* @param  x primer valor
	* @param  y segundo valor
	*return valor abdosluto de la diferencia
	*/
	T error_absoluto(T x, T y){
		return fabs(x-y);
	}

	template<typename T>
	/**
	*@brief Calcula el error relativo entre dos variables
	*@param xNuevo primer valor
	*@param xAnterior segundo valor
	*@return error relativo 
	*/
	T error_relativo(T xnuevo, T xant){
		return ( fabs(xnuevo-xant)/fabs(xnuevo));
	}
		
	template<typename T>
	/**
	*@brief Calcula el error relativo porcentual entre dos variables
	*@param xNuevo primer valor
	*@param xAnterior segundo valor
	*@return error relativo en porcentaje 
	*/
	T error_relativo_porcentual(T xnuevo, T xant){
		return ( fabs(xnuevo-xant)/fabs(xnuevo))*100.0f;
	}	
	
	template<typename T>
	/**
	*@brief Calcula si un numero es cercado a otro con epsilon como entrada por el usuario
	*@param x primer parametro
	*@param y segundo parametro
	*@param epsilon tercer parametro
	*@return verdadero si el numero es cercano al epsilon p falso si no 
	*/
	bool es_cercano(T x, T y, T epsilon){
		constexpr T eps = std::numeric_limits<T>::epsilon();
		//tome ese valor y vuelvalo constante constexpr
		return(error_absoluto(x,y) < std::max(epsilon, eps));
	}
		
	template<typename T>
	/**
	*@brief cacula si un numero es cercano a otro
	*@param x primer parametro
	*@param y segundo parametro
	*@return verdadero si el numero es cercano al epsilon p falso si no
	*/
	bool es_cercano(T x, T y){
		return es_cercano(x,y, std::numeric_limits<T>::epsilon());
	}
	
};
#endif
