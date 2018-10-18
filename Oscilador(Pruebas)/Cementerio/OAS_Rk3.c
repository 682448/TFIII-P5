#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Num_t (float)(100000)//Tiempo final
#define dt (float)(0.001)//Paso en tiempo
#define Num_T (float)(10)//Temperatura final
#define dT (float)(0.5)//Paso en Temperatura
#define trayec (float)(1000)

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
float RandomC(void);
float box_muller(void);

float F(float *xv);

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
float **zz;
float* p_g;
float* p_xvee;

int main()
{
    ini_ran(123456789);
    D1=fopen("OAS_Rk_T-E_cin-D.txt","w");
    //D2=fopen("2.txt","w");
    /*Posición, velocidad, array de promedios para posición y velocidad, array de numeros aleatorios para el box_muller (ahorra tiempo guardarlo en memoria)*/
    float x,v,D[2],epsilon;
    /*Colección de constantes físicas, ahorra tiempo guardarlas en memoria y usarlas cada que calcularlas (producto a producto, suma a suma)  miles de veces*/
    register float temp,temp1,temp2;temp=K/M;temp1=Eta/M,temp2=0.5*dt;

    /*Para cada temperatura uso los mismos números aletorios en el box_muller(0,1) por ello guardo en memoria con malloc()*/
	zz = (float **)malloc(trayec*sizeof(float*));
	for (int i=0;i<trayec;i++)zz[i]=(float*)malloc(Num_t*sizeof(float));
    for(int k=0;k<trayec;k++)for(int t=0;t<Num_t;t++)zz[k][t]=box_muller();

    float g[4];//g[0]=g_11,g[1]=g_21,g[2]=g_12,g[3]=g_22
    float xvee[4];
    p_xvee = &xvee[0];
    p_g = &g[0];

    for(float T=0;T<=dT*Num_T;T+=dT)/*Para varias temperaturas*/
    {
        *(p_xvee + 2) = *(p_xvee + 3) = 0;  epsilon=pow((dt*2*T*temp1/M),0.5);
        for(int k=0;k<trayec;k++)
        {
            *p_xvee = *(p_xvee+1) = 0;
            for(int t=0;t<Num_t;t++)/*Integración mediante algoritmo rk-estocástico 2 orden*/
            {
                *(p_xvee + 1) += epsilon*zz[k][t];    
                *p_g = *(p_xvee + 1);//g_11
                *(p_g + 2) = F(p_xvee);//g_12


                *p_xvee += *p_g*dt;
                *(p_xvee + 1) += *(p_g+2)*dt;
                *(p_g + 1) = *(p_xvee + 1);//g_21
                *(p_g + 3) = F(p_xvee);//g_22

                *p_xvee += 0.5*dt*(*p_g + *(p_g + 1));
                *(p_xvee + 1) += 0.5*dt*(*(p_g + 2) + *(p_g + 3));
                
                




                /*register float g11,g12,g21,g22;
                g11=v+epsilon*zz[k][t];  //Calculo las funciones g11,g12,... para el rk estoc�stico
                g12=-temp1*g11-temp*x;
                g21=v+g12*dt;
                g22=-temp1*g21-temp*(x+g11*dt);
                x=x+temp2*(g11+g21);    //Calculo posiciones y velocidades en cada momento
                v=v+temp2*(g12+g22)+epsilon*zz[k][t];*/
            }

            //D[0]+=v*v; D[1]+=x*x;
            //printf("%f\t %f\t %f\n",T,*(p_xvee),*(p_xvee+1));
            *(p_xvee+2) = *(p_xvee)*(*(p_xvee))*0.5*M/trayec;
            *(p_xvee+3) = *(p_xvee+1)*(*(p_xvee+1))*0.5*M/trayec/2/(dt*Num_t);
            //printf("%f\t %f\t %f\n",T,*(p_xvee+2),*(p_xvee+3));
            //getchar();           
        }
        //Guardo la temperatura, el promedio de energía cinética y energía potencial
        printf("%f\t %f\t %f\n",T,*(p_xvee+2),*(p_xvee+3));
        fprintf(D1,"%f\t %f\t %f\n",T,*(p_xvee+2),*(p_xvee+3));
    }

    fclose(D1);
    //fclose(D2);
    return 0;
}

float F(float *xv){return *xv*K + *(xv + 1)*Eta;}
void ini_ran(int SEMILLA)
{
    srand(SEMILLA);
    for(int k=0;k<256;k++)Wheel[k]=(rand()<<16)+rand();
    ind_ran=ig1=ig2=ig3=0;
}

float RandomC(void)/*Modifico rand() con Parisi-Rapuano*/
{
    float r;
    ig1=ind_ran-11;ig2=ind_ran-29;ig3=ind_ran-17;
    Wheel[ind_ran]=Wheel[ig1]+Wheel[ig2];
    ir1=Wheel[ind_ran]^Wheel[ig3];
    ind_ran++;
    //r=Min+(Max-Min)*ir1*(2.3283063671E-10F);
    r=2.0*ir1*(2.3283063671E-10F)-1;
    return r;
}

float box_muller(void)/*Con media M y varianza s*/
{
    register float x1, x2, w, y1;
    static float y2;
    static int use_last = 0;

    if (use_last){y1 = y2; use_last = 0;}
    else
    {
        do {
            x1 = RandomC();
            x2 = RandomC();
            w = x1 * x1 + x2 * x2;
        } while(w>=1.0);

        w = sqrt( (-2.0*log(w))/w);
        y1 = x1*w;
        y2 = x2*w;
        use_last=1;
    }
    return(y1);//numero:media+y1*varianza, en nuestro caso 0+y1*1
}
