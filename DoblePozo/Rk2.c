#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define tiempo (float)(4E4)
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
float F(float x, float KOnM,float cte);

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
FILE *D3;
FILE *D4;
int main()
{
    float B,Eta,M,T,K,chi,cte;
    scanf("%f",&T);
    scanf("%f",&B);
    scanf("%f",&Eta);
    scanf("%f",&cte);
    M=1;
    K=2*B;
    chi=2*Eta*T/M/M;
    ini_ran(123456789);
    D1=fopen("Data/Equiparticion.csv","a+");
    D2=fopen("Data/Probabilidades.csv","a+");
    D3=fopen("Distribucion_t_pos.csv","w");
    D4=fopen("Distribucion_t_neg.csv","w");

    /*Posición, velocidad, array de promedios para posición y velocidad, array de numeros aleatorios para el box_muller (ahorra tiempo guardarlo en memoria)*/
    float x,v,D[3],EtaOnM,KOnM;
    float t_salto_pos,t_salto_neg,x_ant;
    double prob;
    prob=t_salto_pos=t_salto_neg=0;
    EtaOnM=Eta/M; KOnM=K/M;
    D[0]=D[1]=D[2]=0;
    x=0;  v=sqrt(T/M);
    // Pongo etiquetas para la distribucion de tiempos
    fprintf(D3, "t\n");
    fprintf(D4, "t\n");
    // Hago evolucionar el sistema en el tiempo
    for(int t=0;t<(int)(tiempo/dt);t++)                                                 /*Integración mediante algoritmo rk-estocástico 2 orden*/
    {
        // Guardo la posición para comparar entre un paso y el anterior en la distribucion de tiempos
        x_ant=x;
        // Integro con el Runge Kutta
        register float g11,g12,g21,g22,z; z=sqrt(chi*dt)*Gauss(0,1);
        g11=v+z;                                                                         /*Calculo las funciones g11,g12,... para el rk estoc�stico*/
        g12=-EtaOnM*g11/M+F(x,KOnM,cte);
        g21=v+g12*dt;
        g22=-EtaOnM*(v+g12*dt)+F(x+g11*dt,KOnM,cte);
        x=x+0.5*dt*(g11+g21);                                                            /*Calculo posiciones y velocidades en cada momento*/
        v=v+0.5*dt*(g12+g22)+z;
        D[0]+=v*v; D[1]+=x*x; D[2]=x*x*x*x;
        // La ocupacion en el lado positivo del eje x aumenta en uno
        if(x>=0)prob+=1;
        // Calculo la distribución de tiempos
        if(x_ant>=0&&x>=0)t_salto_pos+=dt;
        else if(x_ant>=0&&x<0)
        {
          //printf(D3,"%f\n",t_salto_pos);
          t_salto_pos=0;
        }
        if(x_ant<=0&&x<=0)t_salto_neg+=dt;
        else if(x_ant>=0&&x<0)
        {
          //fprintf(D4,"%f\n",t_salto_neg);
          t_salto_neg=0;
        }
    }
    // Escribo en fichero la temperatura, B, Eta, probabilidades de ocupacion en -x y en +x
    fprintf(D2,"%.3f, %.3f, %.3f,  %.3f, %.3f, %.3f\n",T,cte,B,Eta,1-prob/(tiempo/dt), prob/(tiempo/dt));
    // Escribo en fichero la temperatura, B, Eta, energia cinética y potencial
    fprintf(D1,"%.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n",T,cte,B,Eta,0.5*M*D[0]/(tiempo/dt),0.5*B*(2*D[1]+1-D[2])/(tiempo/dt));


    fclose(D1);
    fclose(D2);
    fclose(D3);
    fclose(D4);
    return 0;
}

float F(float x,float KOnM,float cte){return -KOnM*x*(x*x-1)+cte;}

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
