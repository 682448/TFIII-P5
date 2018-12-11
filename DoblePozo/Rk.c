#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define tiempo (float)(4E3)
#define dt (float)(1E-2)//Paso en tiempo
#define dT (float)(0.001)//Paso de T*k_b


/*
->ini_ran(int SEMILLA) inicia el array "Wheel" usado para obtener los números aleatorios
->float RandomC(float Max,float Min) implementa el método Parisi-rapuano  descorrelacionando los números guardados en Wheel
->box_muller(float m, float s) devuelve un número asociado con distribución gaussiana
*/
void ini_ran(int SEMILLA);
float RandomC(float Max,float Min);
float Gauss(float m, float s);
float F(float x, float KOnM);

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
    float B,Eta,M,T,K,chi;
    scanf("%f",&T);
    scanf("%f",&B);
    scanf("%f",&Eta);
    M=1;
    K=2*B;
    chi=2*Eta*T/M/M;
    ini_ran(123456789);
    D1=fopen("Equiparticion.csv","a+");
    D2=fopen("Data_Basico.csv","w");

    fprintf(D2,"Tiempo,Posicion,Velocidad\n");
    /*Posición, velocidad, array de promedios para posición y velocidad, array de numeros aleatorios para el box_muller (ahorra tiempo guardarlo en memoria)*/
    float x,v,D[3],EtaOnM,KOnM,PositiveX,NegativeX;
    EtaOnM=Eta/M; KOnM=K/M;
    D[0]=D[1]=D[2]=0;
    x=0;  v=sqrt(T/M);

    for(int t=0;t<(int)(tiempo/dt);t++ )                                                 /*Integración mediante algoritmo rk-estocástico 2 orden*/
    {
        register float g11,g12,g21,g22,z; z=sqrt(chi*dt)*Gauss(0,1);
        g11=v+z;                                                                         /*Calculo las funciones g11,g12,... para el rk estoc�stico*/
        g12=-EtaOnM*g11/M+F(x,KOnM);
        g21=v+g12*dt;
        g22=-EtaOnM*(v+g12*dt)+F(x+g11*dt,KOnM);
        x=x+0.5*dt*(g11+g21);                                                            /*Calculo posiciones y velocidades en cada momento*/
        v=v+0.5*dt*(g12+g22)+z;

        fprintf(D2,"%.3f, %.3f, %.3f\n",t*dt,x,v);
        D[0]+=v*v; D[1]+=x*x; D[2]=x*x*x*x;
    }
    fprintf(D1,"%.3f, %.3f, %.3f\n",T,0.5*M*D[0]/(int)(tiempo/dt),0.5*B*(2*D[1]+1-D[2])/(int)(tiempo/dt));


    fclose(D1);
    fclose(D2);
    return 0;
}

float F(float x,float KOnM){return -KOnM*x*(x*x-1);}

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
