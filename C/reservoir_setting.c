#include "functions.h"

#define N 10//ネットワークサイズ
#define T0 100//過渡時間
#define Tr (5000+T0)//データの数
#define Tl 800//評価用の時間
#define ALL_Length (Tr+Tl)//全データの長さ

double **SPmatrix(int col,int row,double p);
double *SPvector(int col,double p);
double max_eigenvalue(double **A,int M);
void normalize(double *y,int n,int m);
double **W(double a);

double *SPvector(int col,double p){//スパースなベクトル
    double *Sp=dvector(1,col);

    int Nonzero_ele=(int)(p*col);
    while(Nonzero_ele>=0){
        //全要素を確認
        int m_col=rand_I(1,col);//代入する要素位置を適当に決める
        if(Sp[m_col]==0.0){
            Sp[m_col]=rand_z(-1,1);
            Nonzero_ele-=1;//残り入れる個数を確認
        }
    }
    return Sp;
}

double **SPmatrix(int col,int row,double p){//スパースな行列
    double **Sp=dmatrix(1,col,1,row);

    int Nonzero_ele=(int)(p*col*row);
    while(Nonzero_ele>0){
        //全要素を確認
        int m_col=rand_I(1,col),m_row=rand_I(1,row);//代入する要素位置を適当に決める
        if(Sp[m_col][m_row]==0.0){
            Sp[m_col][m_row]=rand_z(-1,1);
            Nonzero_ele-=1;//残り入れる個数を確認
        }
    }
    return Sp;
}

//最大固有値を返す関数
double max_eigenvalue(double **A,int M){
    printf("Search eigen\n");
    /*ながさ1のベクトル*/
    double *y=dvector(1,M);for(int i=1;i<=M;i++)y[i]=1.0/pow(M,0.5);//初期ベクトル

    double lambda=0.0;
    for(int t=0;t<=pow(10,3);t++){
        double *Ay=matrix_vector_product_CRS(A,y,M,M);

        double lambda1=inner_product(1,M,Ay,y);

        /*========収束判定=========*/
        // double sum=0.0;
        // for(int i=1;i<=M;i++){
        //     sum+=pow(Ay[i]-lambda*y[i],2);
        // }

        // if(sum<0.0001){
        //     break;
        // }
        // printf("%d,sum=%f\n",t,sum);
        double th=fabs(lambda-lambda1)/fabs(lambda);
        if(th<0.001){
            break;
        }

        if(t>=pow(10,3)){
            printf("can't eigenval\n");
            exit(1);
        }

        lambda=lambda1;
        for(int i=1;i<=M;i++)y[i]=Ay[i];
        normalize(y,1,M);

        free_dvector(Ay,1,M);
    }

    free_dvector(y,1,M);

    return fabs(lambda);
}

double **W(double a){
    printf("W\n");
    // double **ret=dmatrix(1,N,1,N);

    double **Sp=SPmatrix(N,N,0.2);
    double l_max=max_eigenvalue(Sp,N);
    // printf("l=%f\n",l_max);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            // printf("%f,",Sp[i][j]);
            Sp[i][j]=(a/l_max)*(Sp[i][j]);
        }
        // printf("\n");
    }
    // free_dmatrix(Sp,1,N,1,N);

    // printf("%f\n",max_eigenvalue(Sp,N));

    double l=max_eigenvalue(Sp,N);

    if(fabs(l-a)>0.0001){
        printf("Reservoir is not succsess\n");
        printf("λ=%f\n",l);
        exit(1);
    }else{
        printf("ρ=%f\n",l);
    }

    return Sp;
}


double *Wout(double **M,int N1,int M1,double *T,double lambda){
    printf("Wout\n");
    printf("Mt,");
    double **Mt=trans(M,M1,N1);//MtにはM1,N1型の行列が入る
    printf("OK\n");

    printf("MTM,");
    double **MTM=AB(Mt,M1,N1,M,N1,M1);//係数行列,正定値対称行列
    for(int i=1;i<=M1;i++){
        MTM[i][i]+=lambda;
    }
    printf("OK\n");

    double *MtT=matrix_vector_product_CRS(Mt,T,M1,N1);//右辺

    // for(int i=1;i<=M1;i++){
    //     for(int j=1;j<=M1;j++){
    //         printf("%0.2f,",MTM[i][j]-MTM[j][i]);
    //     }
    //     printf(":%0.2f\n",MtT[i]);
    // }
    double *WOUT=CG_CRS(MTM,MtT,M1);

    // double **L=dmatrix(1,M1,1,M1),**U=dmatrix(1,M1,1,M1);
    // LU(MTM,M1,L,U);

    // double *WOUT=LU_Decomp(L,U,MtT,M1);
    // free_dmatrix(L,1,M1,1,M1);
    // free_dmatrix(U,1,M1,1,M1);

    double *ans=matrix_vector_product_CRS(MTM,WOUT,M1,M1);
    double sum=0.0;
    for(int i=1;i<=M1;i++)sum+=pow(ans[i]-MtT[i],2)/M1;
    if(sum>=0.01){
        printf("solution can't get\n");
        exit(1);
    }else{
        printf("solved\n");
        printf("sum=%f\n",sum);   
    }

    free_dvector(ans,1,M1);
    
    free_dmatrix(Mt,1,M1,1,N1);
    free_dmatrix(MTM,1,M1,1,M1);
    free_dvector(MtT,1,M1);

    return WOUT;
}
