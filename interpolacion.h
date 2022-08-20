#ifndef INTERPOLACION_H
#define INTERPOLACION_H

#include <vector>
#include <iostream>

#include "biseccion.h"
#include "raices.h"

using raices::biseccion;
using raices::resultado_raiz;
using std::vector;

using namespace std;

namespace interpolacion
{
    class newton
    {
    public:
        newton() = delete; // no se acepta crear insancias sin datos
        newton(vector<double> p_x, vector<double> p_y) : x(p_x), y(p_y)
        {
            calcular_coeficientes();
        }
        double interpolar(double xint)
        {
            double yint;
            yint = b[0];
            // double p = 1.0f;
            for (size_t i = 1; i < b.size(); i++)
            {
                double p = 1.0f;
                for (size_t j = 0; j < 1; j++)
                {
                    p *= (xint - x[j]);
                }
                yint += b[i] * p;
            }
            return yint;
        }
        vector<double> coeficientes()
        {
            return b;
        }

    private:
        vector<double> x;
        vector<double> y;
        vector<double> b;
        void calcular_coeficientes()
        {
            size_t n = x.size();
            vector<vector<double>> f;
            // dimensionar la matriz
            f.resize(n); // redimensiona la cantidad de filas
            for (size_t i = 0; i < n; i++)
            {
                f[i].resize(n);
            }
            // llenar la primera columna
            for (size_t i = 0; i < n; i++)
            {
                f[i][0] = y[i];
            }

            // calcular las demas columnas
            for (size_t j = 1; j < n; j++)
            {
                for (size_t i = 0; i < n - j; i++)
                {
                    f[i][j] = (f[i + 1][j - 1] - f[i][j - 1]) / (x[i + j] - x[i]);
                }
            }
            // obtener las coeficientes de la primera columna
            b = f[0];
        }
    };

    class lagrange
    {
    public:
        lagrange() = delete;
        lagrange(vector<double> p_x, vector<double> p_y) : x(p_x), y(p_y)
        {
        }
        double interpolar(double xint)
        {
            double yint = 0.0f;
            size_t n = x.size();
            for (size_t j = 0; j < n; j++)
            {
                double lj = 1.0f;
                for (size_t k = 0; k < n; k++)
                {
                    if (k != j)
                    {
                        lj *= (xint - x[k]) / (x[j] - x[k]);
                    }
                }
                yint += y[j] * lj;
            }
            return yint;
        }
        double ErrorInterpolacion(vector<double> x1, vector<double> y1,double xint)
        {
            double resultado;
            //---------------newton----------------
            newton n(x1, y1);
            // obtener los coeficientes
            vector<double> b = n.coeficientes();
            resultado=b[b.size()-1];
            for(size_t i=0;i<x.size();i++)
                resultado*=(xint-x[i]);

            return resultado;
        }

    private:
        vector<double> x;
        vector<double> y;
    };
    class inversa
    {
    public:
        inversa() = delete;
        inversa(vector<double> p_x, vector<double> p_y) : x(p_x), y(p_y)
        {
        }
        double calcular_biccesion(double yint, int intervalo) // valor de la y, la x es la que se quiere hallar
        {
            
            //encuentra la raiz del primer interval
            double primerRaiz = NAN;
            vector<double> x2{x[intervalo-1],x[intervalo],x[intervalo+1]};
            vector<double> y2{y[intervalo-1],y[intervalo],y[intervalo+1]};
            double primerRango=fabs(y2[intervalo+1]-y2[intervalo-1]);
            primerRaiz=EncontrarRaizBiccesion(x2,y2,yint);

            //encuentra la raiz del segundo intervalo
            double segundaRaiz = NAN;
            vector<double> x3{x[intervalo],x[intervalo+1],x[intervalo+2]};
            vector<double> y3{y[intervalo],y[intervalo+1],y[intervalo+2]};
            double segundoRango=fabs(y3[intervalo+1]-y3[intervalo-1]);
            
            segundaRaiz=EncontrarRaizBiccesion(x3,y3,yint);
            //Comprueba cual es la raiz mas "Util"
            if(primerRango<segundoRango)
                return primerRaiz;

            return segundaRaiz;
        }

    private:
        vector<double> x;
        vector<double> y;

        double EncontrarRaizBiccesion(vector<double> xAux, vector<double> yAux, double yint)
        {
            resultado_raiz resultadoRaiz;
            double resultado = NAN;

            double x0 = xAux[0];
            double x1 = xAux[1];
            double x2 = xAux[2];
            // determinar el polinimio
            newton n(xAux, yAux);
            // usar newton para calculcar el valor de coeficientes
            vector<double> coeficientes = n.coeficientes();

            if (coeficientes.size() != 3)
            {
                return resultado;
            }
            // obtenemos los coeficientes del polinimio
            double b0 = coeficientes[0];
            double b1 = coeficientes[1];
            double b2 = coeficientes[2];
            double erp = 0.01f;
            double maxIter = 100;

            resultadoRaiz = biseccion::calcular([&](double x) -> double
                                                { return (b0 + b1 * (x - x0) + b2 * ((x - x0) * (x - x1)) - yint); },
                                                x0, x2, erp, maxIter);
            resultado = resultadoRaiz.raiz;
            return resultado;
        }
    };
};

#endif
