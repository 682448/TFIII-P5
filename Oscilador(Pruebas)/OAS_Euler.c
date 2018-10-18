#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define tiempo (float)(1000)//Tiempo final
#define dt (float)(0.1)//Paso en tiempo
#define Temperatura (float)(10)//Temperatura final
#define dT (float)(0.5)//Paso en Temperatura
#define trayec (int)(20000)
#define M   (float)(1.0)//Masa
#define Eta (float)(1.0)//Coeficiente viscosidad del medio
#define K   (float)(0.0)//Coeficiente del "muelle" para el oscilador
//Se tomará que la constante de Bolztmann tiene valor unitario

/*
->ini_ran(int SEMILLA) inicia el array "Wheel" usado para obtener los números aleatorios
->float RandomC(float Max,float Min) implementa el método Parisi-rapuano  descorrelacionando los números guardados en Wheel
->box_muller(float m, float s) devuelve un número asociado con distribución gaussiana
*/
void ini_ran(int SEMILLA);
float RandomC(float Max,float Min);
float box_muller(float m, float s);

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
//FILE *D2;
int main()
{
    ini_ran(123456789);
    D1=fopen("OAS_Euler_T-E_cin-D.txt","w");
    //D2=fopen("2.txt","w");
    /*Posición, velocidad, array de promedios para posición y velocidad, array de numeros aleatorios para el box_muller (ahorra tiempo guardarlo en memoria)*/
    float x,v,D[2],epsilon,**zz;
    /*Colección de constantes físicas, ahorra tiempo guardarlas en memoria y usarlas cada que calcularlas (producto a producto, suma a suma)  miles de veces*/
    register float temp,temp1,temp2;temp=K/M;temp1=Eta/M;

    /*Para cada temperatura uso los mismos números aletorios en el box_muller(0,1) por ello guardo en memoria con malloc()*/
	zz = (float **)malloc(trayec*sizeof(float*));
	for (int i=0;i<trayec;i++)zz[i]=(float*)malloc((int)(tiempo/dt)*sizeof(float));
    for(int k=0;k<trayec;k++)for(int t=0;t<(int)(tiempo/dt);t++)zz[k][t]=box_muller(0,1);

    for(float T=0;T<=Temperatura;T+=dT)/*Para varias temperaturas*/
    {
        D[0]=D[1]=0;epsilon=pow((dt*2*T*temp1/M),0.5);
        for(int k=0;k<trayec;k++)
        {
            x=v=0;
            for(int t=0;t<(int)(tiempo/dt);t++)/*Integración mediante algoritmo Euler*/
            {
                x=x+dt*v;
                v=v-dt*temp1*v-dt*temp*x+epsilon*zz[k][t];
            }
            D[0]+=v*v; D[1]+=x*x;
        }
        //Guardo la temperatura, el promedio de energía cinética y energía potencial
        fprintf(D1,"%f\t %f\t %f\n",T,0.5*M*D[0]/(float)(trayec),0.5*M*D[1]/(float)(trayec)/2/tiempo);
    }

    fclose(D1);
    //fclose(D2);
    return 0;
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
