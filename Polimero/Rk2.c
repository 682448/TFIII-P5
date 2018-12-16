#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define tiempo (float)(1E4)//Tiempo final aunque aquí en realidad es adimensional
#define Temperatura (float)(1)//Esto en realidad es energía pues hago T*k_b
#define dT (float)(0.01)//Paso de T*k_b
#define N 4
#define D 3

/*
->ini_ran(int SEMILLA) inicia el array "Wheel" usado para obtener los números aleatorios
->float RandomC(float Max,float Min) implementa el método Parisi-rapuano  descorrelacionando los números guardados en Wheel
->box_muller(float m, float s) devuelve un número asociado con distribución gaussiana
*/
void ini_ran(int SEMILLA);
float RandomC(float Max,float Min);
float Gauss(float m, float s);
float F(float r, float r_ant, float r_pos, int i,int flag);
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

float EtaOnM,KOnM,B,T,chi,dt;
float r[N][D][10],v[N][D][10],z[N][D];
float Long_ant[N][2],Long_pos[N][2],Long_N;
int main()
{

    //scanf("%f",&K);
    //scanf("%f",&B);
    //scanf("%f",&Eta);
    //scanf("%f",&dt);
    B=1;
    KOnM=1;
    EtaOnM=1;
    dt=1E-2;
    T=0.2;


    chi=2*EtaOnM*T;
    ini_ran(123456789);
    D1=fopen("EquiparticionOnParticle.csv","w");

    // Doy las condiciones iniciales para la posición y la velocidad
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<D;j++)
        {
          r[i][j][0]=i*B/sqrt(3); v[i][j][0]=sqrt(T);
          v[i][j][1]=0;
        }
    }
    Long_N=0;

    for(int t=0;t<(int)(tiempo/dt);t++)
    {
      for(int i=0;i<N;i++){// Creo los números de distribución gaussiana multiplicados por el factor asociada con el rk
        for(int j=0;j<D;j++)z[i][j]=sqrt(chi*dt)*Gauss(0,1);
      }
      // Este bule inicializa las variable Long_ que luego utilizaré en la fuerza
      for(int i=0;i<N-1;i++)Long_pos[i][0]=Long_pos[i][1]=Long_ant[i][0]=Long_ant[i][1]=0;

      for(int i=0;i<N;i++){// Con este bucle calculo las longitudes relevantes para la función fuerza
        for(int k;k<D;k++){// Aqui hago la suma a cada una de las dimensiones
          if(i<=N-1){// Lon_pos no esta bien definido para N
            Long_pos[i][0]+=pow(r[i][k][1]-r[i+1][k][1],2);
            Long_pos[i][1]+=pow(r[i][k][1]+dt*(v[i][k][0]+z[i][k])-r[i+1][k][1]-dt*(v[i+1][k][0]+z[i+1][k]),2);
          }
          if(i!=0){// Long_ant no esta bien definido para la primera posición
            Long_ant[i][0]+=pow(r[i][k][1]-r[i-1][k][1],2);
            Long_ant[i][1]+=pow(r[i][k][1]+dt*(v[i][k][0]+z[i][k])-r[i-1][k][1]-dt*(v[i-1][k][0]+z[i-1][k]),2);
          }
        // Deberíamos observar que Long_pos[N]=Long_ant[0]=0
        }
        // Despues de hacer la raiz tiene dimensión de longitud
        Long_pos[i][0]=sqrt(Long_pos[i][0]);
        Long_ant[i][0]=sqrt(Long_ant[i][0]);
        Long_pos[i][1]=sqrt(Long_pos[i][1]);
        Long_ant[i][1]=sqrt(Long_ant[i][1]);
      }

      crea_copia();

      for(int i=0;i<N;i++){// Hago que recorra todas las partículas
        //printf("-------Particula:%d-----\n",i);
        for(int j=0;j<D;j++){// Hago que recorra para cada partícula, cada componente
          Evoluciona(t,i,j);


        }
        if(i<N-1){// Aqui sumo las longitudes de los muelles, al final la longitud será la suma de longitudes del polimero a lo largo del tiempo
          Long_N+=pow((Long_pos[i][0]-B),2);
        }
      }
    }

    fprintf(D1,"Temperatura,Cinetica,Potencial\n");
    float r2_media,v2_media;
    r2_media=v2_media=0;

    for(int i=0;i<N;i++)
    {
      for(int j=0;j<D;j++)
        {
          v2_media+=v[i][j][1];
        }
    }
    Long_N/=((N-1)*tiempo/dt);
    v2_media/=(int)(N*tiempo/dt);
    fprintf(D1,"%.3f, %.3f, %.3f\n",T,0.5*v2_media,0.5*Long_N);


    fclose(D1);
    return 0;
}


float F(float r, float r_ant, float r_pos,int i,int flag)
{
  float long_ant,long_pos;
  long_ant=r_ant-r;
  long_pos=r-r_pos;
  //return -K*r;
  if(i==0)
  {
    return -KOnM*(r-r_pos)+KOnM*B*long_pos/Long_pos[i][flag];
  }
  else if(i==N-1)
  {
    return -KOnM*(r-r_ant)-KOnM*B*long_ant/Long_ant[i][flag];
  }
  else
  {
    return -KOnM*(2*r-r_ant-r_pos)-KOnM*B*(long_ant/Long_ant[i][flag]-long_pos/Long_pos[i][flag]);
  }
}

void Evoluciona(float t,int i,int j)
{
  register float g11,g12,g21,g22;

  g11=v[i][j][0]+z[i][j];
  g12=-EtaOnM*g11+F(r[i][j][1],r[i-1][j][1],r[i+1][j][1],i,0);
  g21=v[i][j][0]+g12*dt;
  g22=-EtaOnM*(v[i][j][0]+g12*dt)+F(r[i][j][1]+g11*dt,r[i-1][j][1]+(v[i-1][j][0]+z[i-1][j])*dt,r[i+1][j][1]+(v[i+1][j][0]+z[i+1][j])*dt,i,1);

  r[i][j][0]=r[i][j][0]+0.5*dt*(g11+g21);
  v[i][j][0]=v[i][j][0]+0.5*dt*(g12+g22)+z[i][j];
  //printf("%f\n", v[i][j][1]);getchar();
  v[i][j][1]+=v[i][j][0]*v[i][j][0];
}

void crea_copia()
{
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<D;j++)
    {
      r[i][j][1]=r[i][j][0];
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

float Gauss(float m, float s)// Con media M y varianza s
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
