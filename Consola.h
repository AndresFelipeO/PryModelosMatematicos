#ifndef CONSOLA_H
#define CONSOLA_H

#include<iostream>
#include<iomanip>

using std::cout;
using std::cerr;
using std::cin;
using std::endl;
using std::string;
using std::to_string;
using std::setprecision;

namespace toolconsol{
	struct Consola{
		/**
		*@brief Escribe en consola 
		*@param mensaje, primer parametro
		*/
		void Write(string mensaje){
			cout<<mensaje;
		}
		/**
		*@brief Escribe en consola con salto de linea 
		*@param mensaje, primer parametro
		*/
		void WriteLn(string mensaje){
			cout<<mensaje<<endl;
		}
		/**
		*@brief Escribe en consola con salto de linea y un numero flotante 
		*@param mensaje, primer parametro
		*@param numero, segundo parametro
		*/
		void WriteLn(string mensaje,float num){
			cout<<mensaje<<num<<endl;
		}
		/**
		*@brief Escribe en consola un numero flotante 
		*@param mensaje, primer parametro
		*@param numero, segundo parametro
		*/
		void Write(string mensaje,float num){
			cout<<mensaje<<num;
		}
		/**
		*@brief Escribe en consola un mensaje de error 
		*@param mensaje, primer parametro
		*/
		void WriteError(string mensaje){
			cerr<<mensaje<<endl;
		}
		/**
		*@brief lee una variable de tipo entero 
		*@param mensaje, primer parametro
		*@return entero que lee en consola 
		*/
		int ScanInt(string mensaje){
			int entero;
			Write(mensaje);
			cin>>entero;
			return entero;
		}
		/**
		*@brief lee una variable de tipo flotante 
		*@param mensaje, primer parametro
		*@return flotante que lee en consola 
		*/
		float ScanFloat(string mensaje){
			float flotante;
			Write(mensaje);
			cin>>flotante;
			return flotante;
		}
		/**
		*@brief Pausa la consola
		*/
		void Pausar(){
			system("pause");
		}
		/**
		*@brief limpia la consola
		*/
		void Limpiar(){
			system("cls");
		}
	};
};
#endif
