#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define tiempo (float)(4E4)//Tiempo final aunque aquí en realidad es adimensional
#define Temperatura (float)(1)//Esto en realidad es energía pues hago T*k_b
#define dT (float)(0.01)//Paso de T*k_b
#define N 5
#define D 3

/*
->ini_ran(int SEMILLA) inicia el array "Wheel" usado para obtener los números aleatorios
->float RandomC(float Max,float Min) implementa el método Parisi-rapuano  descorrelacionando los números guardados en Wheel
->box_muller(float m, float s) devuelve un número asociado con distribución gaussiana
*/
void ini_ran(int SEMILLA);
float RandomC(float Max,float Min);
float Gauss(float m, float s);
float F(float r, float r_ant, float r_pos);
void Evoluciona(float t,int i, int j);
void crea_copia();
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

float Eta,M,K,B,chi,dt;
float r[N][D],v[N][D],Acumulador[N][D][2];
float r_copy[N][D],Long_ant,Long_pos;
int main()
{

    //scanf("%f",&K);
    //scanf("%f",&B);
    //scanf("%f",&Eta);
    //scanf("%f",&dt);
    K=1;
    B=1;
    Eta=1;
    dt=0.01;
    M=1;
    chi=Eta/2/sqrt(K*M);
    ini_ran(123456789);
    D1=fopen("EquiparticionOnParticle.csv","w");
    //D2=fopen("Eta0.01.csv","w");
    //fprintf(D2,"Tiempo,Posicion,Velocidad\n");

    //Posiciones y velocidades iniciales de cada patícula del polimero
    for(int i=0;i<N;i++)for(int j=0;j<D;j++)
    {
      r[i][j]=i*B; v[i][j]=1;
      Acumulador[i][j][0]=0;
      Acumulador[i][j][1]=0;
    }

    for(int t=0;t<(int)(tiempo/dt);t++)
    {
      crea_copia();
      for(int i=0;i<N;i++)
      {
        //printf("-------Particula:%d-----\n",i);
        //Long_pos=sqrt(pow(r_copy[i][0]-r_copy[i+1][0],2)+pow(r_copy[i][1]-r_copy[i+1][1],2)+pow(r_copy[i][2]-r_copy[i+1][2],2));
        //Long_ant=sqrt(pow(r_copy[i][0]-r_copy[i-1][0],2)+pow(r_copy[i][1]-r_copy[i-1][1],2)+pow(r_copy[i][2]-r_copy[i-1][2],2));
        for(int j=0;j<D;j++)
        {
          //printf("%f\n", r_copy[i][j]);
          Evoluciona(t,i,j);
        }
        //printf("Long_pos: %f",Long_pos);
        //getchar();
      }
    }
    fprintf(D1,"Temperatura,Cinetica,Potencial\n");
    float r2_media,v2_media;
    r2_media=v2_media=0;

    for(int i=0;i<N;i++)
    {
      for(int j=0;j<D;j++)
        {
          r2_media+=Acumulador[i][j][0];
          v2_media+=Acumulador[i][j][1];
        }
    }
    r2_media/=(int)(N*tiempo/dt);
    v2_media/=(int)(N*tiempo/dt);
    for(float T=0;T<=Temperatura;T+=dT)
    {
      //Guardo la temperatura, el promedio de energía cinética y energía potencial
      fprintf(D1,"%.3f, %.3f, %.3f\n",T,0.5*T*r2_media,0.5*T*v2_media);
    }

    fclose(D1);
    fclose(D2);
    return 0;
}


float F(float r, float r_ant, float r_pos)
{
  float long_ant,long_pos;
  long_ant=r_ant-r;
  long_pos=r-r_pos;
  return -K*r;
  //return -K*(2*r-r_ant-r_pos)-K*B*(long_ant/Long_ant-long_pos/Long_pos);
}

void Evoluciona(float t,int i,int j)
{
  register float g11,g12,g21,g22,z; z=sqrt(4*chi*dt)*Gauss(0,1);

  g11=v[i][j]+z;                                                                         /*Calculo las funciones g11,g12,... para el rk estoc�stico*/
  g12=-2*chi*g11+F(r_copy[i][j],r_copy[i-1][j],r_copy[i+1][j]);
  g21=v[i][j]+g12*dt;
  g22=-2*chi*(v[i][j]+g12*dt)+F(r_copy[i][j]+g11*dt,r_copy[i-1][j],r_copy[i+1][j]);

  r[i][j]=r[i][j]+0.5*dt*(g11+g21);                                                            /*Calculo posiciones y v[0][0]elocidades en cada momento*/
  v[i][j]=v[i][j]+0.5*dt*(g12+g22)+z;
  //fprintf(D2,"%.3f, %.3f, %.3f, %.3f\n",t*dt,r[0][0],v[0][0],z);
  //fprintf(D2,"%.3f, %.3f, %.3f\n",t*dt,r[0][0],v[0][0]);
  Acumulador[i][j][0]+=v[i][j]*v[i][j]; Acumulador[i][j][1]+=r[i][j]*r[i][j];
  //printf("%d%d) VAlor:%f\n",i,j, Acumulador[i][j][0]);
}

void crea_copia()
{
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<D;j++)
    {
      r_copy[i][j]=r[i][j];
    }
  }
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
