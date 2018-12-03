#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define tiempo (float)(4E4)//Tiempo final aunque aquí en realidad es adimensional
#define Temperatura (float)(1)//Esto en realidad es energía pues hago T*k_b
#define dT (float)(0.01)//Paso de T*k_b
#define N 1

/*
->ini_ran(int SEMILLA) inicia el array "Wheel" usado para obtener los números aleatorios
->float RandomC(float Max,float Min) implementa el método Parisi-rapuano  descorrelacionando los números guardados en Wheel
->box_muller(float m, float s) devuelve un número asociado con distribución gaussiana
*/
void ini_ran(int SEMILLA);
float RandomC(float Max,float Min);
float Gauss(float m, float s);
float F(float x);
void Evoluciona(float t,int i, int j);

/*Definición de los parámetros usados en float ini_ran(int SEMILLA) y RandomC(float Max,float Min)  para descorrelacionar los números generados con rand()*/
unsigned char ind_ran,ig1,ig2,ig3;
unsigned int Wheel[256],ir1;
/*
Ficheros
->En D1 se guardan los valores de temperatura, energía cinética y energía potencial
Este fichero esta pensado para comprobar a partir de el la relación de Einstein y el teorema de equipartición de la energía
->En D2 se guarda
Esta pensado para comprobar la relación de fluctuación disipación
*/
FILE *D1;
FILE *D2;

float Eta,M,K,chi,dt;
float r[3][1],v[3][1],D[2];
int main()
{

    scanf("%f",&K);
    scanf("%f",&Eta);
    scanf("%f",&dt);
    M=1;
    chi=Eta/2/sqrt(K*M);
    ini_ran(123456789);
    D1=fopen("Eta0.01_Equi.csv","w");
    D2=fopen("Eta0.01.csv","w");                                                              /*Posición, velocidad, array de promedios para posición y velocidad, array de numeros aleatorios para el box_muller (ahorra tiempo guardarlo en memoria)*/
    fprintf(D2,"Tiempo,Posicion,Velocidad\n");
    for(int i=0;i<N;i++)for(int j=0;j<1;j++){r[i][j]=0;  v[i][j]=1;}

    for(int t=0;t<(int)(tiempo/dt);t++)
    {
      for(int i=0;i<N;i++)
      {
        for(int j=0;j<1;j++)
        {
          Evoluciona(t,i,j);
        }
      }
    }


    fprintf(D1,"Temperatura,Cinetica,Potencial\n");
    for(float T=0;T<=Temperatura;T+=dT)/*Para varias temperaturas*/
    {
      //Guardo la temperatura, el promedio de energía cinética y energía potencial
      fprintf(D1,"%.3f, %.3f, %.3f\n",T,0.5*T*D[0]/(int)(tiempo/dt),0.5*T*D[1]/(int)(tiempo/dt));
    }

    fclose(D1);
    fclose(D2);
    return 0;
}

float F(float x){return -x;}

void Evoluciona(float t,int i,int j)
{                                                                                     /*Integración mediante algoritmo rk-estocástico 2 orden*/

        register float g11,g12,g21,g22,z; z=sqrt(4*chi*dt)*Gauss(0,1);
        g11=v[i][j]+z;                                                                         /*Calculo las funciones g11,g12,... para el rk estoc�stico*/
        g12=-2*chi*g11+F(r[i][j]);
        g21=v[i][j]+g12*dt;
        g22=-2*chi*(v[i][0]+g12*dt)+F(r[i][j]+g11*dt);
        r[i][j]=r[i][j]+0.5*dt*(g11+g21);                                                            /*Calculo posiciones y v[0][0]elocidades en cada momento*/
        v[i][j]=v[i][j]+0.5*dt*(g12+g22)+z;
        //fprintf(D2,"%.3f, %.3f, %.3f, %.3f\n",t*dt,r[0][0],v[0][0],z);
        //fprintf(D2,"%.3f, %.3f, %.3f\n",t*dt,r[0][0],v[0][0]);
        D[0]+=v[i][j]*v[i][j]; D[1]+=r[i][j]*r[i][j];
}

void ini_ran(int SEMILLA)
{
    srand(SEMILLA);
    for(int k=0;k<256;k++)Wheel[k]=(rand()<<16)+rand();
    ind_ran=ig1=ig2=ig3=0;
}

float RandomC(float Max,float Min)/*Modifico rand() con Parisi-Rapuano*/
{
    float r;
    ig1=ind_ran-11;ig2=ind_ran-29;ig3=ind_ran-17;
    Wheel[ind_ran]=Wheel[ig1]+Wheel[ig2];
    ir1=Wheel[ind_ran]^Wheel[ig3];
    ind_ran++;
    r=Min+(Max-Min)*ir1*(2.3283063671E-10F);
    return r;
}

float Gauss(float m, float s)/*Con media M y varianza s*/
{
    register float x1, x2, w, y1;
    static float y2;
    static int use_last = 0;

    if (use_last){y1 = y2; use_last = 0;}
    else
    {
        do {
            x1 = RandomC(1,-1);
            x2 = RandomC(1,-1);
            w = x1 * x1 + x2 * x2;
        } while(w>=1.0);

        w = sqrt( (-2.0*log(w))/w);
        y1 = x1*w;
        y2 = x2*w;
        use_last=1;
    }
    return(m + y1*s);
}
