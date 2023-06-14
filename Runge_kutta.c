#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define p 10.0
#define r 28.0
#define b 8.0/3.0

#define dt 0.02
#define tmax 100

double f(double t,double x,double y,double z){
    return -p*x+p*y;
}

double g(double t,double x,double y,double z){
    return -x*z+r*x-y;
}

double h(double t,double x,double y,double z){
    return x*y-b*z;
}

int main(void){
    double x=0.0,y=1,z=2;
    FILE *fp;
    fp=fopen("Lolentz.txt","w");

    if(fp==NULL){printf("File error!\n");exit(1);}

    for(double t=0;t<tmax;t+=dt){
        double K1=dt*f(t,x,y,z);
        double L1=dt*g(t,x,y,z);
        double M1=dt*h(t,x,y,z);

        double K2=dt*f(t+dt/2.0,x+K1/2.0,y+L1/2.0,z+M1/2.0);
        double L2=dt*g(t+dt/2.0,x+K1/2.0,y+L1/2.0,z+M1/2.0);
        double M2=dt*h(t+dt/2.0,x+K1/2.0,y+L1/2.0,z+M1/2.0);

        double K3=dt*f(t+dt/2.0,x+K2/2.0,y+L2/2.0,z+M2/2.0);
        double L3=dt*g(t+dt/2.0,x+K2/2.0,y+L2/2.0,z+M2/2.0);
        double M3=dt*h(t+dt/2.0,x+K2/2.0,y+L2/2.0,z+M2/2.0);

        double K4=dt*f(t+dt,x+K3,y+L3,z+M3);
        double L4=dt*g(t+dt,x+K3,y+L3,z+M3);
        double M4=dt*h(t+dt,x+K3,y+L3,z+M3);

        x+=(K1+2.0*K2+2.0*K3+K4)/6.0;
        y+=(L1+2.0*L2+2.0*L3+L4)/6.0;
        z+=(M1+2.0*M2+2.0*M3+M4)/6.0;

        printf("x=%f,y=%f,z=%f\n",x,y,z);
        // fprintf(fp,"%f %f %f\n",x,y,z);
        fprintf(fp,"%f\n",x);
    }
    fclose(fp);

    return 0;
}

// int main(void) {
// double x = 0.0;
// double y = 1.0;
// double z = 2.0;
// double k1, k2, k3, k4;
// double l1, l2, l3, l4;
// double m1, m2, m3, m4;

//     for (double t = 0.0; t <= tmax; t += dt) {
//         k1 = dt * f(t, x, y, z);
//         l1 = dt * g(t ,x, y, z);
//         m1 = dt * h(t ,x, y, z);   

//         k2 = dt * f(t + dt / 2.0, x + k1 / 2.0, y + l1 / 2.0, z + m1 / 2.0);
//         l2 = dt * g(t + dt / 2.0, x + k1 / 2.0, y + l1 / 2.0, z + m1 / 2.0);
//         m2 = dt * h(t + dt / 2.0, x + k1 / 2.0, y + l1 / 2.0, z + m1 / 2.0);
//         k3 = dt * f(t + dt / 2.0, x + k2 / 2.0, y + l2 / 2.0, z + m2 / 2.0);
//         l3 = dt * g(t + dt / 2.0, x + k2 / 2.0, y + l2 / 2.0, z + m2 / 2.0);
//         m3 = dt * h(t + dt / 2.0, x + k2 / 2.0, y + l2 / 2.0, z + m2 / 2.0);

//         k4 = dt * f(t + dt, x + k3, y + l3, z + m3);
//         l4 = dt * g(t + dt, x + k3, y + l3, z + m3);
//         m4 = dt * h(t + dt, x + k3, y + l3, z + m3);

//         x += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
//         y += (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0;
//         z += (m1 + 2.0 * m2 + 2.0 * m3 + m4) / 6.0;

//         printf("x=%f,y=%f,z=%f\n",x,y,z);
//     }

// return(0);
// }
