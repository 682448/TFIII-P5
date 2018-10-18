#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define tiempo (float)(100)//Tiempo final aunque aquí en realidad es adimensional
#define dt (float)(0.01)//Paso en tiempo
#define Temperatura (float)(1)//Esto en realidad es energía pues hago T*k_b
#define dT (float)(0.001)//Paso de T*k_b
#define trayec (int)(10000)

#define K   (float)(1)//Coeficiente del "muelle" para el oscilador
#define M   (float)(1)//Masa
#define Eta (float)(1)//Coeficiente viscosidad del medio
#define chi (float)(Eta/2/sqrt(K*M))//El coeficiente
/*
->ini_ran(int SEMILLA) inicia el array "Wheel" usado para obtener los números aleatorios
->float RandomC(float Max,float Min) implementa el método Parisi-rapuano  descorrelacionando los números guardados en Wheel
->box_muller(float m, float s) devuelve un número asociado con distribución gaussiana
*/
void ini_ran(int SEMILLA);
float RandomC(float Max,float Min);
float box_muller(float m, float s);
float F(float x);

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
int main()
{
    ini_ran(123456789);
    D1=fopen("Rk_Equiparticion.txt","w");
    D2=fopen("Rk_txv.txt","w");
    /*Posición, velocidad, array de promedios para posición y velocidad, array de numeros aleatorios para el box_muller (ahorra tiempo guardarlo en memoria)*/
    float x,v,D[2],**zz;

    /*Para cada temperatura uso los mismos números aletorios en el box_muller(0,1) por ello guardo en memoria con malloc()*/
  	zz = (float **)malloc(trayec*sizeof(float*));
  	for (int i=0;i<trayec;i++)zz[i]=(float*)malloc((int)(tiempo/dt)*sizeof(float));
    for(int k=0;k<trayec;k++)for(int t=0;t<(int)(tiempo/dt);t++)zz[k][t]=sqrt(4*chi*dt)*box_muller(0,1);


    D[0]=D[1]=0;
    for(int k=0;k<trayec;k++)
    {
        x=0;  v=1;

        for(int t=0;t<(int)(tiempo/dt);t++)/*Integración mediante algoritmo rk-estocástico 2 orden*/
        {
            register float g11,g12,g21,g22;
            g11=v+zz[k][t];  /*Calculo las funciones g11,g12,... para el rk estoc�stico*/
            g12=-2*chi*g11-F(x);
            g21=v+g12*dt;
            g22=-2*chi*(v+g12*dt)-F(x+g11*dt);
            x=x+0.5*dt*(g11+g21);    /*Calculo posiciones y velocidades en cada momento*/
            v=v+0.5*dt*(g12+g22)+zz[k][t];

            if(k==0){fprintf(D2,"%.3f\t %.3f\t %.3f\n",t*dt,x,v);}
        }
        D[0]+=v*v; D[1]+=x*x;
    }

    for(float T=0;T<=Temperatura;T+=dT)/*Para varias temperaturas*/
    {
      //Guardo la temperatura, el promedio de energía cinética y energía potencial
      fprintf(D1,"%.3f\t %.3f\t %.3f\n",T,0.5*T*D[0]/(float)(trayec),0.5*T*D[1]/(float)(trayec));
    }

    fclose(D1);
    fclose(D2);
    return 0;
}

float F(float x){return x;}

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

float box_muller(float m, float s)/*Con media M y varianza s*/
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
