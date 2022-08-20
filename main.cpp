/*
Elaborado por:
-Andres Felipe Ocampo Chaguendo
*/

#include "Consola.h"
#include "biseccion.h"
#include "precision.h"
#include "newtonraphson.h"
#include "reglafalsa.h"
#include "secante.h"
#include "muller.h"
#include "interpolacion.h"
#include "rl.h"
#include "integracion.h"

using raices::biseccion;
using raices::reglafalsa;
using raices::newton_raphson;
using raices::secante;
using raices::muller;
using toolconsol::Consola;

using interpolacion::inversa;
using interpolacion::lagrange;
using regresion::rl;
using regresion::rl_exponencial;
using regresion::rl_potencia;
using regresion::cuadratica;
using regresion::resultado_rl;
using regresion::resultado_rl_exponencial;
using regresion::resultado_rl_potencia;
using regresion::resultado_cuadratica;

using integracion::integral;

void Menu(int);
void ProcesarOpcion(int);

void Funciones(double,int,double,double,double);
void MetodoBiseccion(function<double(double)>,double,double,double,int);
void MetodoReglaFalsa(function<double(double)>,double,double,double,int);
void MetodoNewton(function<double(double)>,function<double(double)>,double,double,double,int);
void MetodoSecante(function<double(double)>,double,double,double,int);
void MetodoMuller(function<double(double)>,double,double,double,double,int);
void EvaluarRaiz(resultado_raiz);

void InterInversa();
int intervalo(vector<double>, double);
void Regresion();
double RegresionLinealSimple(vector<double>,vector<double>,double);
double RegresionLinealPotencia(vector<double>,vector<double>,double);
double CalcularEstimacion(char,double,double,double,double);
void MostrarInformacionRegresion(string,double,double,double,double,double,double,double,char,double,double);
double RegresionLinealExponencial(vector<double>,vector<double>,double);
double RegresionCuadratica(vector<double>,vector<double>,double);
int Mayor(vector<double>);
void InterpolacionLagrange();
void ImprimirDatosInter(vector<double>,vector<double>,double,double,double);

void FuncionNoPolinomica(function<double(double)>,function<double(double)>,function<double(double)>,double,double,double,double,double,string);
void trapecioNoPoliomica(function<double(double)> ,function<double(double)> ,double,double,double,string);
void Simpson13NoPoliomica(function<double(double)>,function<double(double)>,double,double,double,string);
void Simpson38NoPoliomica(function<double(double)>,function<double(double)>,double,double,double,string);
void ImprimirResultadosNoPolinimoca(string,string,double,double,double,double,double,double,double);
void Simpson13Polinomica(function<double(double)>,function<double(double)>,double,double,double,string);
void ImrpimirInformacionPolinomicas(string,string,double,double,double,double,double);

Consola consola;

int main (int argc, char *argv[]) {
	Menu(0);
	return 0;
}

void Menu(int opc)
{
	if(opc!=4){
		system("cls");
		consola.WriteLn("   Menu Modelos Matematicos  ");
		consola.WriteLn("--------------------------------------");
		consola.WriteLn("1.Raices De Ecuaciones");
		consola.WriteLn("2.Ajuste De Curvas");
		consola.WriteLn("3.Integracion");
		consola.WriteLn("4.Salir");
		consola.WriteLn("--------------------------------------");
		opc=consola.ScanInt("Digite una opcion: ");
		ProcesarOpcion(opc);
		Menu(opc);
	}
}

void ProcesarOpcion(int opc){
	consola.Limpiar();
	switch (opc)
	{
	case 1:
		consola.WriteLn("---Parte I. Raices De Ecuaciones---");
		consola.WriteLn("Funcion de ejemplo: f(x)= 5-((21750 sqrt(x))/(138 sqrt(x)+275)) (x*0.001)^(((1)/(2))) ");
		Funciones(1*pow(10,-7),consola.ScanInt("Digite el numero de iteracciones: "),consola.ScanFloat("Digite x2: "),consola.ScanFloat("Digite x1: "),consola.ScanFloat("Digite x0: "));
		break;
	case 2:
		consola.WriteLn("---Parte II. Ajuste De Curvas---");
		consola.WriteLn("II.I Regresion");
		Regresion();
		consola.Pausar();
		consola.Limpiar();
		consola.WriteLn("II.II Interpolacion");
		InterpolacionLagrange();
		break;
	case 3:
		consola.WriteLn("---Parte III. Integracion---");
		consola.WriteLn("III.I No polinomica");
		FuncionNoPolinomica([](double x)->double {return 5*x*sin(x);},[](double x)->double {return 10*cos(x)-5*x*sin(x);},[](double x)->double {return 5*x*sin(x)-20*cos(x);},1,3,31,32,33,"5xsen(x)");
		consola.Pausar();
		consola.Limpiar();
		consola.WriteLn("III.II Polinomica");
		Simpson13Polinomica([](double x)->double {return 1-x-4*pow(x,3)+3*pow(x,5);} ,[](double x)->double {return 360*x;} ,0,5,8,"1-x-4x^3+3x^5");
		break;
	case 4:
		consola.WriteLn("Saliendo del programa");
		break;
	default:
		consola.WriteLn("Opcion Incorrecta......");
		break;
	}
	consola.Pausar();
}

//------------------------region raices de una ecuacion------------------------

void Funciones(double erp,int maxIter,double x2,double x1,double x0){
	function<double(double)> f1 = [](double x)->double{ return  sin(0.5*x)-pow((pow(x,3)-1),2);};
	function<double(double)> df1 = [](double x)->double{ return -(2175*(69*sqrt(x)+275))/(sqrt(10)*pow(138*sqrt(x)+275,2));};
	MetodoBiseccion(f1,x0,x1,erp,maxIter);
	MetodoReglaFalsa(f1,x0,x1,erp,maxIter);
	MetodoNewton(f1,df1,x0,x1,erp,maxIter);
	MetodoSecante(f1,x0,x1,erp,maxIter);
	MetodoMuller(f1,x0,x1,x2,erp,maxIter);
}
//Metodo que calcula la funcion por Biseccion
void MetodoBiseccion(function<double(double)> f1,double x0,double x1,double erp,int maxIter){
	consola.WriteLn("\n ------- CASO BISECCION -------");
	resultado_raiz resultado= biseccion::calcular(f1,x0,x1,erp,maxIter);;
	EvaluarRaiz(resultado);
}
//Metodo que calcula la funcion por Regla falsa
void MetodoReglaFalsa(function<double(double)> f1,double x0,double x1,double erp,int maxIter){
	consola.WriteLn("\n ------- CASO REGLA FALSA -------");
	resultado_raiz resultado= reglafalsa::calcular(f1,x0,x1,erp,maxIter);
	EvaluarRaiz(resultado);
}
//Metodo que calcula la funcion por Newton
void MetodoNewton(function<double(double)> f1,function<double(double)> df1,double x0,double x1,double erp,int maxIter){
	consola.WriteLn("\n ------- CASO NEWTON RAPHSON -------");
	consola.WriteLn(" -- PUNTO X0 --");
	resultado_raiz resultadonewton= newton_raphson::calcular(f1,df1,x0,erp,maxIter);
	EvaluarRaiz(resultadonewton);
}
//Metodo que calcula la funcion por Secante
void MetodoSecante(function<double(double)> f1,double x0,double x1,double erp,int maxIter){
	consola.WriteLn("\n ------- CASO SECANTE -------");
	resultado_raiz resultado= secante::calcular(f1,x0,x1,erp,maxIter);
	EvaluarRaiz(resultado);
}
//Metodo que calcula la funcion por Muller
void MetodoMuller(function<double(double)> f1,double x0,double x1,double x2,double erp,int maxIter){
	consola.WriteLn("\n ------- CASO MULLER -------");
	resultado_raiz resultado = muller::calcular(f1,x0,x1,x2,erp,maxIter);
	EvaluarRaiz(resultado);
}
//Metodo que evalua la raiz de los metodos de raices de ecuaciones
void EvaluarRaiz(resultado_raiz resultado){
	if(!resultado.encontrada()){
		consola.WriteError("No se encontro la raiz");
		return;
	}
	consola.WriteLn("La raiz es: "+to_string(resultado.raiz));
	consola.WriteLn("Se encontro en "+to_string(resultado.iteraciones)+" Iteraciones");
	consola.WriteLn("Con un erp de ",resultado.erp);
}

//------------------------region Ajuste De Curvas------------------------

//metodo de regresion
void Regresion(){
	vector<double>x={40.0f,42.0f,43.0f,44.0f,46.0f,48.0f,49.0f,52.0f,53.0f,54.0f,57.0f,58.0f};
	vector<double>y={825.0f,830.0f,960.0f,840.0f,895.0f,910.0f,890.0f,1010.0f,990.0f,1012.0f,1030.0f,1050.0f};
	double xint=50;
	integracion::imprimir_tabla(x,y);
	consola.WriteLn("Para x=50");
	double resultadoRegreLinealS=RegresionLinealSimple(x,y,xint);
	double resultadoRegreLinealP= RegresionLinealPotencia(x,y,xint);
	double resultadoRegreLinealE= RegresionLinealExponencial(x,y,xint);
	double resultadoRegreCuadratica=RegresionCuadratica(x,y,xint);
	vector<double>resultados={resultadoRegreLinealS,resultadoRegreLinealP,resultadoRegreLinealE,resultadoRegreCuadratica};
	vector<string>metodos={"metodo regresion simple","metodo regresion potencial","metodo regresion exponencial","metodo regresion cuadratica"};
	int posMetodo=Mayor(resultados);
	consola.WriteLn("\nEn x=50 el metodo mas eficiente fue el \n'"+(metodos[posMetodo])+"' con una incertidumbre de "+to_string(resultados[posMetodo]*100)+"%");
}
//metodo que retorna la posion de un numero mayor en un vector
int Mayor(vector<double> vec)
{
    int may=vec[0];
    int pos=0;
    for(size_t i=1;i<vec.size();i++)
        if (vec[i]>may)
        {
            may=vec[i];
            pos=i;
        }
	return pos;
}
//metodo de regresion lineal simple
double RegresionLinealSimple(vector<double> x, vector<double> y,double xint)
{
	consola.WriteLn("\nRegresion lineal simple. ");
	resultado_rl r = rl::calcular(x, y);
	MostrarInformacionRegresion(r.ecuacion(),r.sy,r.syx,r.r2,r.b0,r.b1,0,xint,'S',r.st,r.sr);
	return r.r2;
}
//metodo de regresion Lineal Potencial
double RegresionLinealPotencia(vector<double> x, vector<double> y,double xint)
{
	consola.WriteLn("\nRegresion lineal funcion potencia. ");
	resultado_rl_potencia r = rl_potencia::calcular(x, y);
	MostrarInformacionRegresion(r.ecuacion(),r.sy,r.syx,r.r2,r.a,r.c,0,xint,'P',r.st,r.sr);
	return r.r2;
}
//metodo de regresion Lineal exponencial
double RegresionLinealExponencial(vector<double> x, vector<double> y,double xint)
{
	consola.WriteLn("\nRegresion Lineal Exponencial" );
	resultado_rl_exponencial r = rl_exponencial::calcular(x, y);
	MostrarInformacionRegresion(r.ecuacion(),r.sy,r.syx,r.r2,r.a,r.c,0,xint,'E',r.st,r.sr);
	return r.r2;
}
//metodo de regresion cuadratica
double RegresionCuadratica(vector<double> x, vector<double>y,double xint)
{
	consola.WriteLn("\nRegresion Cuadratica");
	resultado_cuadratica r = cuadratica::cualcular(x,y);
	MostrarInformacionRegresion(r.ecuacion(),r.sy,r.syx,r.r2,r.a0,r.a1,r.a2,xint,'C',r.st,r.sr);
	return r.r2;
}
//Metodo que imprimie los resultados de los metodos de regresion
void MostrarInformacionRegresion(string ecuacion,double sy,double syx,double r2,double aux0,double aux1,double aux2,double xint,char tipo,double st,double sr){
	consola.WriteLn("Ecuacion de regresion: "+ecuacion);
	consola.WriteLn("Sumatoria: " ,st);
	consola.WriteLn("Discrepancia: " ,sr );
	consola.WriteLn("Desviacion estandar: ", sy);
	consola.WriteLn("Error estandar de aproximacion: " ,syx);
	if (syx < sy)
		consola.WriteLn("La regresion se considera aceptable");
	else
		consola.WriteLn("La regresion no se considera aceptable");
	consola.WriteLn("En una estimacion de x="+to_string(xint)+" es igual a "+to_string(CalcularEstimacion(tipo,aux0,aux1,aux2,xint)));
	consola.WriteLn("Cofieciente de determinacion: ",r2 );
}
//Metodo que calcula la estimacion dependiendo del metodo
double CalcularEstimacion(char tipo,double aux0,double aux1,double aux2,double xint){
	if(tipo=='S')
		return aux1*xint+aux0;
	else if (tipo=='P')
		return aux1*pow(xint,aux0);
	else if(tipo=='E')
		return aux1*exp(aux0*xint);
	else 
		return aux2*xint*xint + aux1*xint + aux0;
}
//metodo de interpolacion de Lagrange
void InterpolacionLagrange()
{
	consola.WriteLn("Tabla con los intervalos");
	vector<double> x{1.0f,1.1f,1.2f,1.3f,1.4f,1.5f};
	vector<double> y{-0.0778461f,-0.16174304f,-0.2486851f,-0.33795622f,-0.42861239f,-0.51944702f};
	integracion::imprimir_tabla(x,y);
	consola.WriteLn("Estimacion de interpolacion para x=1.24");
	//primeros intervalos
	vector<double> x1{1.1f,1.2f,1.3f};
	vector<double> y1{-0.16174304f,-0.2486851f,-0.33795622f};
	//segundo intervalo
	vector<double> x2{1.2f,1.3f,1.4f};
	vector<double> y2{-0.2486851f,-0.33795622f,-0.42861239f};
	//--------Primeros Puntos-----------------------
	lagrange l(x1, y1);
	double xint = 1.24f;
	double yint = l.interpolar(xint);
	double erroEstimado=l.ErrorInterpolacion(x,y,xint);
	ImprimirDatosInter(x1,y1,xint,yint,erroEstimado);
	//--------Segundos Puntos-----------------------
	lagrange l2(x2, y2);
	double yint2 = l2.interpolar(xint);
	double erroEstimado2=l2.ErrorInterpolacion(x,y,xint);
	ImprimirDatosInter(x2,y2,xint,yint2,erroEstimado2);
	consola.WriteLn( "-----------------------------------------------------------------------------");
	cout<<"Dado los resultados obtenidos el "<<(fabs(erroEstimado)<fabs(erroEstimado2)?"primer":"segundo")<<" intervalo es el mas adecuado"<<endl;
	consola.WriteLn( "-----------------------------------------------------------------------------" );
}

void ImprimirDatosInter(vector<double> x,vector<double> y,double xint,double yint,double errorEstimado){
	consola.WriteLn("\nInterpolacion de Lagrnage");
	integracion::imprimir_tabla(x,y);
	consola.Write("p1(", xint);
	consola.WriteLn(")=",yint);
	consola.WriteLn("Error estimado: ",errorEstimado);
}

//------------------------region Integracion------------------------

//dependiendo de la funcion no polinomica la calcula por todos los metodos de integrales
void FuncionNoPolinomica(function<double(double)> f,function<double(double)> df2,function<double(double)> df4,double limite_inf,double limite_sup,double segmentoT,double segmento13,double segmento38,string funcion){
	trapecioNoPoliomica(f,df2,limite_inf,limite_sup,segmentoT,funcion);
	Simpson13NoPoliomica(f,df4,limite_inf,limite_sup,segmento13,funcion);
	Simpson38NoPoliomica(f,df4,limite_inf,limite_sup,segmento38,funcion);
}
//metodo trapecio para funciones no polinomicas
void trapecioNoPoliomica(function<double(double)> f,function<double(double)> df2,double limite_inf,double limite_sup,double segmentos,string funcion){
	double resultado = integral::trapecio(f , limite_inf, limite_sup, segmentos);
	double error=integral::trapecioErrorEstimadoNoPolinomica(df2 , limite_inf, limite_sup);
	int k=integral::precisionIntegral(error);
	double ajuste=integral::trapecioAjusteIntegralNoPolinomica(resultado,df2,limite_inf,limite_sup);
	ImprimirResultadosNoPolinimoca("Metodo Trapecio",funcion,limite_inf,limite_sup,segmentos,resultado,error,k,ajuste);
}
//metodo simpson 1/3 para funciones no polinomicas
void Simpson13NoPoliomica(function<double(double)> f,function<double(double)> df4,double limite_inf,double limite_sup,double segmentos,string funcion){
	double resultado = integral::simpson13(f, limite_inf, limite_sup, segmentos);
	double error=integral::simpson13ErrorEstimadoNoPolinomica(df4 , limite_inf, limite_sup,segmentos);	
	int k=integral::precisionIntegral(error);
	double ajuste=integral::simpson13AjusteIntegralNoPolinomica(resultado,df4,limite_inf,limite_sup,segmentos);
	ImprimirResultadosNoPolinimoca("Metodo Simpson 1/3",funcion,limite_inf,limite_sup,segmentos,resultado,error,k,ajuste);
	consola.WriteLn("Los segmentos que requerimos para obtener un error (Et)<1*10^-7");
	segmentos=integral::CantidadSegmentosSimpson13([](double x)->double {return ((5*x)*sin(x))-(20*cos(x));},limite_inf,limite_sup);
	error=integral::simpson13ErrorEstimadoNoPolinomica([](double x)->double {return ((5*x)*sin(x))-(20*cos(x));} , limite_inf, limite_sup,segmentos);	
	k=integral::precisionIntegral(error)-1;
	consola.WriteLn("Usando el metodo simpsion 1/3 cuantos segmentos se requiere para obtener un error (ET)<1*10^-7?");
	consola.WriteLn("Teniendo en cuenta lo que exige el metodo, la cantidad requeridas de segmentos debe ser");
	consola.WriteLn("Segmentos: ",segmentos);
	consola.WriteLn("entonces, k: ",k);
	consola.WriteLn("ErrorEstimado: ",error);
}
//metodo simpson 1/3 para funciones no polinomicas
void Simpson38NoPoliomica(function<double(double)> f,function<double(double)> df4,double limite_inf,double limite_sup,double segmentos,string funcion){
	double resultado = integral::simpson38(f , limite_inf, limite_sup, segmentos);
	double error=integral::simpson38ErrorEstimadoNoPolinomica(df4, limite_inf, limite_sup,segmentos);	
	int k=integral::precisionIntegral(error);
	double ajuste=integral::simpson13AjusteIntegralNoPolinomica(resultado,df4,limite_inf,limite_sup,segmentos);
	ImprimirResultadosNoPolinimoca("Metodo Simpson 3/8",funcion,limite_inf,limite_sup,segmentos,resultado,error,k,ajuste);
	segmentos=integral::CantidadSegmentosSimpson38([](double x)->double {return ((5*x)*sin(x))-(20*cos(x));},limite_inf,limite_sup);
	error=integral::simpson38ErrorEstimadoNoPolinomica([](double x)->double {return ((5*x)*sin(x))-(20*cos(x));} , limite_inf, limite_sup,segmentos);	
	k=integral::precisionIntegral(error)-1;
	consola.WriteLn("Usando el metodo simpsion 3/8 cuantos segmentos se requiere para obtener un error (ET)<1*10^-7?");
	consola.WriteLn("Teniendo en cuenta lo que exige el metodo, la cantidad requeridas de segmentos debe ser");
	consola.WriteLn("Segmentos: ",segmentos);
	consola.WriteLn("Entonces, k: ",k);
	consola.WriteLn("ErrorEstimado: ",error);
}
//Imprime los datos obtenidos de cualquiera de los metodos de funciones no polinimicas
void ImprimirResultadosNoPolinimoca(string metodo,string funcion,double limite_inf,double limite_sup,double segmentos,double resultado,double error,double k,double ajuste){
	consola.WriteLn("\n--------"+metodo+"--------");
	consola.WriteLn("Funcion: "+funcion);
	consola.Write("Limite inferior (a): ",limite_inf);
	consola.Write(", Limite superior (b): ",limite_sup);
	consola.WriteLn(", Numero de segmentos ",segmentos);
	consola.WriteLn("Integral: ", resultado);
	consola.WriteLn("Error estimado: ",error);
	consola.WriteLn("precision k= ",k);
	if(k>0)
		consola.WriteLn("lo que asegura una precision de por lo menos ("+to_string(k)+") cifras decimales exactas en la aproximacion obtenida"
		+"Final mente ajustamos la integral\n"
		+"Integral ajustada: ",ajuste);
	else if(k!=-2)
		consola.WriteLn("como k=0 no garantiza que el valor obtenido aproxime al valor exacto con alguna cibra decimal exacta");	
}
//metodo simpson 1/3 para funciones polinomicas
void Simpson13Polinomica(function<double(double)> f,function<double(double)> df4,double limite_inf,double limite_sup,double segmentos,string funcion){
	double resultado = integral::simpson13(f, limite_inf, limite_sup, segmentos);
	double error=integral::simpson13ErrorEstimadoPolinomica(df4 , limite_inf, limite_sup,segmentos);	
	ImrpimirInformacionPolinomicas("Metodo Simpson 1/3",funcion,limite_inf,limite_sup,segmentos,resultado,error);	
}
//Imprime los datos obtenidos de cualquiera de los metodos de funciones polinimicas
void ImrpimirInformacionPolinomicas(string metodo,string funcion,double limite_inf,double limite_sup,double segmentos,double resultado,double error){
	consola.WriteLn("\n--------"+metodo+"--------");
	consola.WriteLn("Funcion: "+funcion);
	consola.Write("Limite inferior (a): ",limite_inf);
	consola.Write(", Limite superior (b): ",limite_sup);
	consola.WriteLn(", Numero de segmentos ",segmentos);
	consola.WriteLn("Integral: ", resultado);
	consola.WriteLn("Error estimado: ",error);
	if(error!=0)
		consola.WriteLn("Integral ajustada: ",resultado+error);
}
