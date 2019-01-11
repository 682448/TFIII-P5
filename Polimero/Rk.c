#include <stdio.h>
#include <stdlib.h>
#include <math.h>




#define tiempo (float)(3E3)//Tiempo final aunque aquí en realidad es adimensional
#define N 65
#define D 3
#define cte (float)(2.2)
float RandomC(float Max,float Min);
float Gauss(float m, float s);
float F(float r, float r_ant, float r_pos, int i,int j,int k);
void Evoluciona(float t,float dt,int i, int j);
void ini_ran(int SEMILLA);
void crea_copia();
// Definición de los parámetros usados en float ini_ran(int SEMILLA) y RandomC(float Max,float Min)  para descorrelacionar los números generados con rand()
unsigned char ind_ran,ig1,ig2,ig3;
unsigned int Wheel[256],ir1;
// Fichero
FILE *D1;

float EtaOnM,KOnM,B,T,chi,dt;
float r[N][D][2],v[N][D][3],z[N][D];
float Long_ant[2],Long_pos[2],Long_N,R_G,Ree;
double err,err1;
float Fconst[3];
int main()
{
    Fconst[1]=Fconst[2]=0;
    Fconst[0]=cte;
    B=1;
    EtaOnM=1;
    dt=1E-3;
    scanf("%f",&T);
    scanf("%f",&KOnM);
    chi=2*EtaOnM*T;
    ini_ran(123456789);
    char noum[50];
    scanf("%s\n",noum);

    D1=fopen(noum,"a+");
    // Doy las condiciones iniciales para la posición y la velocidad
    for(int i=0;i<N;i++){
        for(int j=0;j<D;j++){
          r[i][j][0]=v[i][j][2]=0;
          v[i][j][0]=0;
          //printf("%f\n",pow(-1,(int)(RandomC(20,1)))); getchar();
        }
        r[i][0][0]=i*B;
    }
    v[0][0][0]=+10;
    v[N-1][0][0]=-10;
    Long_N=R_G=Ree=err=err1=0;
    for(int t=0;t<(int)(tiempo/dt);t++){
      // Ahora paso a calcuar el centro de masas
      float R_CM[3];
      R_CM[0]=R_CM[1]=R_CM[2]=0;
      crea_copia();
      for(int i=0;i<N;i++)for(int j=0;j<D;j++)z[i][j]=sqrt(chi*dt)*Gauss(0,1);
      for(int j=0;j<D;j++){
        for(int i=0;i<N;i++){
          R_CM[j]+=r[i][j][1];
        }
      }
      R_CM[0]/=N;R_CM[1]/=N;R_CM[2]/=N;
      //if(t==0)printf("%f\t%f\t%f\t\n",R_CM[0],R_CM[1],R_CM[2]);getchar();
      // Calculo las longitudes que luego uso para calcular la fuerza
      for(int i=0;i<N;i++){
        Long_pos[0]=Long_pos[1]=Long_ant[0]=Long_ant[1]=0;
        for(int k=0;k<D;k++){
            Long_pos[0]+=pow(r[i][k][1]-r[i+1][k][1],2);
            Long_ant[0]+=pow(r[i][k][1]-r[i-1][k][1],2);
            Long_pos[1]+=pow(r[i][k][1]+(v[i][k][1]+z[i][k])*dt-r[i+1][k][1]-(v[i+1][k][1]+z[i+1][k])*dt,2);
            Long_ant[1]+=pow(r[i][k][1]+(v[i][k][1]+z[i][k])*dt-r[i-1][k][1]-(v[i-1][k][1]+z[i-1][k])*dt,2);
        }
        Long_pos[0]=sqrt(Long_pos[0]);
        Long_ant[0]=sqrt(Long_ant[0]);
        Long_pos[1]=sqrt(Long_pos[1]);
        Long_ant[1]=sqrt(Long_ant[1]);
        for(int j=0;j<D;j++){// Hago que el sistema evolucione
          Evoluciona(t,dt,i,j);
        }
        // Sumo las longitudes del polimero para todos los tiempos
        if(i<N-1)Long_N+=pow((Long_pos[0]-B),2);
      }
      for(int i=0;i<N;i++){// Esto aculumla los términos cuadráticos del radio de giro
        for(int j=0;j<D;j++){
          R_G+=pow(r[i][j][1]-R_CM[j],2);
        }
      }
      // Distacia extremo a extremo
    Ree+=r[N-1][0][1]-r[0][0][1];
    err1+=(r[N-1][0][1]-r[0][0][1])*(r[N-1][0][1]-r[0][0][1]);
    //printf("%f\n%f\n\n",Ree,Ree2);getchar();
    }
    R_G/=(N*tiempo/dt);// Este es el radio de giro promedio en el tiempo
    R_G=sqrt(R_G);
    Ree/=(tiempo/dt);// Esta es la distacia extremo a extremo promedio a lo largo del tiempo
    for(int t=0;t<(int)(tiempo/dt);t++){
      err+=pow(r[N-1][0][1]-r[0][0][1]-Ree,2);
    }
    err/=(tiempo/dt);
    err1/=(tiempo/dt);
    err=sqrt(err);
    float r2_media,v2_media;
    r2_media=v2_media=0;

    for(int i=0;i<N;i++){
      for(int j=0;j<D;j++){
          v2_media+=v[i][j][2];
        }
    }
    Long_N/=((N-1)*tiempo/dt);// Esta es la longitud promedio en el tiempo de cada muelle
    v2_media/=(int)(N*tiempo/dt);// Esta es la velocidad promedio en el tiempo de cada partícula del polímero
    //fprintf(D1,"\n%d, %.2f, %.3f, %.2f, %.2f, %.2f, %.1e, %.1e, %.2f, %.2f, %.2f, %.2f\n",N,cte,T,B,EtaOnM,KOnM,tiempo,dt,R_G,Ree,0.5*v2_media/T,0.5*KOnM*Long_N/T);
    fprintf(D1,"\n%d, %.2f, %.3f, %.2f, %.2f, %.2f, %.1e, %.1e, %.2f, %.2f, %.2f, %.2f\n",N,cte,T,B,EtaOnM,KOnM,tiempo,dt,Ree,Ree-(N-1)*cte/KOnM,err,sqrt(err1-Ree*Ree));

    fclose(D1);
    return 0;
}


float F(float r, float r_ant, float r_pos,int i,int j,int k)
{
  float long_ant,long_pos;
  long_ant=r_ant-r;
  long_pos=r-r_pos;
  //return -K*r;
  if(i==0){
    return -KOnM*(r-r_pos)+KOnM*B*long_pos/Long_pos[k]-Fconst[j];
  }
  else if(i==N-1){
    return -KOnM*(r-r_ant)-KOnM*B*long_ant/Long_ant[k]+Fconst[j];
  }
  else{
    return -KOnM*(2*r-r_ant-r_pos)-KOnM*B*(long_ant/Long_ant[k]-long_pos/Long_pos[k]);
  }
}

void Evoluciona(float t,float dt,int i,int j){
  register float g11,g12,g21,g22;

  g11=v[i][j][0]+z[i][j];
  g12=-EtaOnM*g11+F(r[i][j][1],r[i-1][j][1],r[i+1][j][1],i,j,0);
  g21=v[i][j][0]+g12*dt;
  g22=-EtaOnM*(v[i][j][0]+g12*dt)+F(r[i][j][1]+(v[i][j][0]+z[i][j])*dt,r[i-1][j][1]+(v[i-1][j][0]+z[i-1][j])*dt,r[i+1][j][1]+(v[i+1][j][0]+z[i+1][j])*dt,i,j,1);

  r[i][j][0]=r[i][j][0]+0.5*dt*(g11+g21);
  v[i][j][0]=v[i][j][0]+0.5*dt*(g12+g22)+z[i][j];
  v[i][j][2]+=v[i][j][0]*v[i][j][0];
}

void crea_copia(){
  for(int i=0;i<N;i++){
    for(int j=0;j<D;j++){
      r[i][j][1]=r[i][j][0];
      v[i][j][1]=v[i][j][0];
    }
  }
}

void ini_ran(int SEMILLA){
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
