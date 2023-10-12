#include "Include/functions.h"

// #define nu 1.0/100.0//拡散係数
// #define Re 1/(nu)//Reynolds数
#define dt 0.01//時間の刻み幅
#define N 500//空間方向の分割数
#define T 2100//時間方向の分割数
#define a -1.0
#define b 1.0//区間の端
#define dx (b-a)/(double)N//空間の刻み幅
#define al 1./(dx)//なんか逆数定義しないと上手くいかん

//初期条件
double f(double x){
    return 1.0;//exp(-x*x)+1.0;
}

//次のステップの数値解
void *u_new(double *up,double rand){
    double nu=0.01;
    /*==============係数行列の作成==============*/
    //漸化式au_{i+1}+bui+cu_{i-1}から，行列の作成
    double **At=dmatrix(0,N,0,N);//時刻tにおける係数行列
    double alpha=dt*al,beta=(nu*al)*(dt*al);//刻み幅に応じたパラメータ
    
    At[0][0]=1.0,At[N][N]=1.0;//境界条件用の係数考慮
    //漸化式によって，離散化した各方向の解を求める
    for(int i=1;i<N;i++){
        double coef1=-beta;//u_{i+1}の係数
        double coef0=1+alpha*up[i]+2*beta;//u_{i}の係数
        double coef2=-alpha*up[i]-beta;//u_{i-1}の係数

        At[i][i+1]=coef1,At[i][i]=coef0,At[i][i-1]=coef2;
    }
    /*========================================*/

    /*==============連立方程式によって，次のステップの解を計算==============*/
    double **L=dmatrix(0,N,0,N),**U=dmatrix(0,N,0,N);
    LU(At,N+1,L,U);
    up[0]=rand,up[N]=up[N-1];

    double *u=LU_Decomp(L,U,up,N+1);

    //解の更新
    for(int i=0;i<=N;i++){
        // double xi=a+dx*i;
        // if(i==0){
        //     up[i]=f(a);//ディリクレ境界条件の考慮
        // }else if(i==N){
        //     up[i]=u[i-1];//ノイマン条件の考慮
        // }else{
        //     up[i]=u[i];
        // }
        up[i]=u[i];
    }
    /*===============================================================*/
    free_dmatrix(At,0,N,0,N);
    free_dmatrix(L,0,N,0,N),free_dmatrix(U,0,N,0,N);
    free_dvector(u,0,N+1);
}

double *input(int input){
    double *ret=dvector(0,T-1);

    FILE *fp;
    char str[200];
    sprintf(str,"inputdata/randomset%d.txt",input);
    fp=fopen(str,"r");
    if (fp == NULL) {
        printf("入力読み込み失敗\n");
        exit(1);
    }

    double data;
    //解の表示
    int i=0;
    while(i<T && fscanf(fp,"%lf",&data)!=EOF){
        fscanf(fp,"%lf\n",&data);
        ret[i]=data;
        i++;
    }
    fclose(fp);

    return ret;
}

//近似関数の値を得る
double u_app(double x,double *up){
    double ret=0.0;
    for(int i=0;i<N;i++){
        double r=a+i*dx;
        if(r<=x && x<=r+dx){//xが入る要素の確認
            double coef=(up[i+1]-up[i])/dx;
            ret=coef*(x-r)+up[i];
        }
    }
    return ret;
}

int main(int argc,char *argv[]){
    int Nr=232;//格子点の総数
    int input_num=0;//入力の番号
    double tau=5;//データをとる感覚

    double Re=atof(argv[1]);
    double nu=1./Re;

    if(argc<2){
        printf("please input the reynold number.\n");
        exit(1);
    }
    double *Input=input(input_num);
    // double *Input=dvector(0,T-1);
    // for(int i=0;i<T;i++)Input[i]=0.5*cos(i)+0.5;

    //時刻tnにおける解を計算する．
    double *up=dvector(0,N);//１ステップ前の数値解
    for(int i=0;i<=N;i++){
        double xi=a+dx*i;
        up[i]=f(xi);
    }
    
    FILE *res_mat;
    char str_res[40];
    sprintf(str_res,"Reservoir_data/Re=%f.txt",Re);
    res_mat=fopen(str_res,"w");

    int k=0,n=0;//リザバーの時間
    double t=0;

    while(k<T){
        double rand=Input[k];//入力データ
        // /*==============解のファイル保存==============*/
        // FILE *fp;//,*res_mat;
        // char str[50];
        // sprintf(str,"solution/sol%d.txt",n);
        // fp=fopen(str,"w");

        printf("Re=%f:t=%f,k=%d\n",Re,t,k);fflush(stdout);
        // //解の表示
        // for(int i=0;i<=N;i++){
        //     double xi=a+dx*i;
        //     fprintf(fp,"%f %f\n",xi,up[i]);
        // }
        // fclose(fp);
        // n+=1;
        // /*========================================*/
        
        //リザバーに用いるデータの取得
        if(fabs((k+1)*tau-t)<0.01){
            //リザバーの行列作成
            for(int i=0;i<=Nr;i++){
                //up[i]にはu(a+i*dx)の数値解が入る．Nr>Nになったらセグフォ
                double a1=a+dx,b1=b-dx;//格子点を取る区間の両端
                double xi=a1+i*(b1-a1)/Nr;//格子点の座標
                double data;
                if(i==0){
                    data=1.0;
                }else{
                    data=u_app(xi,up);
                }
                fprintf(res_mat,"%f ",data);
            }
            fprintf(res_mat,"\n");
            k+=1;
        }

        // u_new(up,rand);//解の更新
        
        /*==============係数行列の作成==============*/
        //漸化式au_{i+1}+bui+cu_{i-1}から，行列の作成
        double **At=dmatrix(0,N,0,N);//時刻tにおける係数行列
        double alpha=dt*al,beta=(nu*al)*(dt*al);//刻み幅に応じたパラメータ
        
        At[0][0]=1.0,At[N][N]=1.0;//境界条件用の係数考慮
        //漸化式によって，離散化した各方向の解を求める
        for(int i=1;i<N;i++){
            double coef1=-beta;//u_{i+1}の係数
            double coef0=1+alpha*up[i]+2*beta;//u_{i}の係数
            double coef2=-alpha*up[i]-beta;//u_{i-1}の係数

            At[i][i+1]=coef1,At[i][i]=coef0,At[i][i-1]=coef2;
        }
        /*========================================*/

        /*==============連立方程式によって，次のステップの解を計算==============*/
        double **L=dmatrix(0,N,0,N),**U=dmatrix(0,N,0,N);
        LU(At,N+1,L,U);
        up[0]=rand,up[N]=up[N-1];

        double *u=LU_Decomp(L,U,up,N+1);

        //解の更新
        for(int i=0;i<=N;i++){
            // double xi=a+dx*i;
            // if(i==0){
            //     up[i]=f(a);//ディリクレ境界条件の考慮
            // }else if(i==N){
            //     up[i]=u[i-1];//ノイマン条件の考慮
            // }else{
            //     up[i]=u[i];
            // }
            up[i]=u[i];
        }
        /*===============================================================*/
        free_dmatrix(At,0,N,0,N);
        free_dmatrix(L,0,N,0,N),free_dmatrix(U,0,N,0,N);
        free_dvector(u,0,N+1);

        t+=dt;
    }

    free_dvector(up,0,N);
    free_dvector(Input,0,T-1);
    fclose(res_mat);

    return 0;
}