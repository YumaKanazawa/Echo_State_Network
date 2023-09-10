#include "reservoir_setting.c"

typedef double (* ODE_Sol)(double t,double *x,int i,int n);//被積分関数の定義
typedef double** (* Jacobi)(double t,double *x,int n);//時刻tにおけるJacobi行列

#define dt 0.01//時間ステップ
#define Sym_L 100000.0//ALL_Length

#define sigma 10.0
#define gamma 28.0
#define beta (8./3.)

//1ステップ先の微分方程式の解
double *x_dt(ODE_Sol df,double t,double *x,int n,double dt_1){
  //時刻tにおける解ベクトルがxで与えられるとき，t+Δtでのベクトルを返す

  double *x_dt=dvector(1,n);//返す用の配列

  double *K1=dvector(1,n);
  double *K2=dvector(1,n);
  double *K3=dvector(1,n);
  double *K4=dvector(1,n);

  //時刻tからdt_1ステップ先の解を計算s

  /*K1=Δt*F(t,x)*/
  for(int i=1;i<=n;i++)K1[i]=dt_1*df(t,x,i,n);

  /*K2=Δt*F(t+Δt/2,x+K1/2)*/
  double *x_K=dvector(1,n);
  for(int i=1;i<=n;i++)x_K[i]=x[i]+K1[i]/2.;//x+K1/2
  for(int i=1;i<=n;i++)K2[i]=dt_1*df(t+dt_1/2.,x_K,i,n);

  /*K3=Δt*F(t+Δt/2,x+K2/2)*/
  for(int i=1;i<=n;i++)x_K[i]=x[i]+K2[i]/2.;//x+K2/2
  for(int i=1;i<=n;i++)K3[i]=dt_1*df(t+dt_1/2.,x_K,i,n);

  /*K4=Δt*F(t+Δt,x+K3)*/
  for(int i=1;i<=n;i++)x_K[i]=x[i]+K3[i];//x+K3
  for(int i=1;i<=n;i++)K4[i]=dt_1*df(t+dt_1,x_K,i,n);

  for(int j=1;j<=n;j++){
    x_dt[j]=x[j]+(K1[j]+2*K2[j]+2*K3[j]+K4[j])/6.;
  }

  free_dvector(K1,1,n);
  free_dvector(K2,1,n);
  free_dvector(K3,1,n);
  free_dvector(K4,1,n);


  return x_dt;
}


double *Lyapunov_exponent(ODE_Sol df,Jacobi J_l,int n){//n次元力学系のリアプノフ指数
  double T=100.0;//Sym_L-1-num;

  double *riap=dvector(1,n);//返すべきリアプノフ指数の値の格納
 
  //行列の微分方程式を解く
  double **X=I(1,n,1,n);//微分方程式の解となる行列(初期行列)

  double *x_sol=dvector(1,n);//元の常微分方程式の初期値
  for(int i=1;i<=n;i++)x_sol[i]=1.0;//rand_z(-1.0,1.0)

  FILE *fp;
  fp=fopen("lorenz.txt","w");
  if(fp==NULL){
    printf("err");
    exit(1);
  }

  double t=0;
  while(1){
    printf("%f:",t);fflush(stdout);
    /*======リアプノフ指数の計算用のグラムシュミットの直交化===========*/
    //行列Xを正規直交化する．

    // //グラムシュミットの直交化による各ベクトルの直交化
    // double **W_orth=dmatrix(1,n,1,n);

    // //L列目までを考える
    // //ベクトル行列積による計算
    // for(int L=1;L<=n;L++){
    //   double *al=dvector(1,n);//直交化をするベクトルの切り出し
    //   for(int i=1;i<=n;i++)al[i]=X[i][L];

    //   double **Q_s_i=dmatrix(1,n,1,n);
    //   for(int i=1;i<=n;i++){
    //     for(int j=1;j<L;j++){
    //       Q_s_i[i][j]=W_orth[i][j];
    //     }
    //   }
    //   double **Q_T=trans(Q_s_i,n,n);
      
    //   double *w=matrix_vector_product_CRS(Q_T,al,n,n);

    //   double *q_hat=matrix_vector_product_CRS(Q_s_i,w,n,n);
    //   for(int i=1;i<=n;i++)al[i]-=q_hat[i];
    //   free_dvector(q_hat,1,n);

    //   double norm_l=sqrt(inner_product(1,n,al,al));//直交化のみをした際のベクトルの大きさ

    //   if(t>=T){
    //     riap[L]+=log(norm_l);
    //   }

    //   normalize(al,1,n);//ベクトルの正規化

    //   for(int i=1;i<=n;i++)W_orth[i][L]=al[i];//正規直交基底の切り出し

    //   free_dvector(w,1,n);
    //   free_dmatrix(Q_T,1,n,1,n);
    //   free_dmatrix(Q_s_i,1,n,1,n); 
    //   free_dvector(al,1,n);
    // }
    int sta=1,M=n;

    double **W_orth=dmatrix(sta,sta+M-1,sta,sta+M-1);
    /*------------------Step1 sta+(M/2)-1までを直交化する---------------*/
    int P=sta+(M/2)-1;

    // int P=sta+M-1;
    //P列目までを考える
    //ベクトル行列積による計算
    double **Q_s_i=dmatrix(sta,sta+M-1,sta,sta+M-1);//直交化済みのベクトルが入った行列

    for(int L=sta;L<=P;L++){//sta+(M/2)-1までを直交化する
      // printf("%d\n",L);
      double *al=dvector(sta,sta+M-1);//直交化をするベクトルの切り出し
      for(int i=sta;i<=sta+M-1;i++)al[i]=X[i][L];//L列目を切り出し

      for(int i=sta;i<=sta+M-1;i++){//直交化したベクトルの格納，L-1までのベクトルを格納
        for(int j=sta;j<=L-1;j++){
          Q_s_i[i][j]=W_orth[i][j];
        }
      }

      double **Q_T=trans(Q_s_i,M,M);//直交化したベクトルの転置行列
      double *w=matrix_vector_product_CRS(Q_T,al,M,M);//内積計算

      double *q_hat=matrix_vector_product_CRS(Q_s_i,w,M,M);

      for(int i=sta;i<=sta+M-1;i++)al[i]-=q_hat[i];
      free_dvector(q_hat,sta,sta+M-1);

      double norm_l=vector_norm1(al,sta,sta+M-1,2.0);
      if(t>=T){
        riap[L]+=log(norm_l);
      }

      normalize(al,sta,sta+M-1);//ベクトルの正規化

      for(int i=sta;i<=sta+M-1;i++)W_orth[i][L]=al[i];//正規直交基底の切り出し

      free_dvector(w,sta,sta+M-1);
      free_dmatrix(Q_T,sta,sta+M-1,sta,sta+M-1);

      free_dvector(al,sta,sta+M-1);

      //Q_S_iに直交化済みのベクトル全てを格納する
      if(L==P){
        for(int i=sta;i<=sta+M-1;i++)Q_s_i[i][L]=W_orth[i][L];
      }
    }
    /*----------------------------------------------------------------*/
    //Q_S_iにはsta~P-1までの直交化済みのベクトルが入った行列

    /*-------------------------step2 行列積を用いたグラムシュミットの直交化法-----------------------------*/
    double **A_res=dmatrix(sta,sta+M-1,sta,sta+M-1);//直交化が施されていないベクトルをまとめた行列(P+1~Last)
    for(int i=sta;i<=sta+M-1;i++){
      for(int j=P+1;j<=sta+M-1;j++){
        A_res[i][j]=X[i][j];
      }
    }

    double **Q_T=trans(Q_s_i,M,M);//Q_s_iはP列目までの直交化が済んだもの
    double **Product_matrix=AB(Q_T,M,M,A_res,M,M);//直交化が済んでいないベクトルと直交ベクトルとの内積
    double **Matrix_minus=AB(Q_s_i,M,M,Product_matrix,M,M);
    for(int i=sta;i<=sta+M-1;i++){
      for(int j=sta;j<=sta+M-1;j++){
        A_res[i][j]-=Matrix_minus[i][j];
      }
    }


    //A_resの各方向を正規化する
    for(int L=P+1;L<=sta+M-1;L++){
      double *al=dvector(sta,sta+M-1);//直交化をするベクトルの切り出し
      for(int i=sta;i<=sta+M-1;i++)al[i]=A_res[i][L];//L列目を切り出し
    
      double **Q_s_i_1=dmatrix(sta,sta+M-1,sta,sta+M-1);//L列目以降の直交化済みのベクトルを格納
      for(int i=sta;i<=sta+M-1;i++){
        for(int j=P+1;j<=L-1;j++){
          Q_s_i_1[i][j]=W_orth[i][j];//Q_i_sに直交化済みのベクトルを格納
        }
      }

      double **Q_1_T=trans(Q_s_i_1,M,M);

      double *w=matrix_vector_product_CRS(Q_1_T,al,M,M);
      double *w1=matrix_vector_product_CRS(Q_s_i_1,w,M,M);
      for(int i=sta;i<=sta+M-1;i++)al[i]-=w1[i];

      free_dvector(w1,sta,sta+M-1);
      free_dvector(w,sta,sta+M-1);
      free_dmatrix(Q_1_T,sta,sta+M-1,sta,sta+M-1);
      free_dmatrix(Q_s_i_1,sta,sta+M-1,sta,sta+M-1);

      double norm_l=vector_norm1(al,sta,sta+M-1,2.0);
      if(t>=T){
        riap[L]+=log(norm_l);
      }

      normalize(al,sta,sta+M-1);

      for(int i=sta;i<=sta+M-1;i++)W_orth[i][L]=al[i];//正規直交基底の切り出し
      free_dvector(al,sta,sta+M-1);
    }  

    free_dmatrix(Matrix_minus,sta,sta+M-1,sta,sta+M-1);
    free_dmatrix(Product_matrix,sta,sta+M-1,sta,sta+M-1);
    free_dmatrix(Q_T,sta,sta+M-1,sta,sta+M-1);
    free_dmatrix(A_res,sta,sta+M-1,sta,sta+M-1);
    free_dmatrix(Q_s_i,sta,P,sta,sta+M-1); 

  
    /*==============行列の数値解X_{t+Δt}================*/
    /*======行列微分方程式をここで解く===========*/
    //時間依存の行列微分方程式をルンゲ・クッタ法で解く
    //必要なものの定義
    double *x_next_2=x_dt(df,t,x_sol,n,dt/2.);//t+Δt/2秒後の解
    double *x_next_1=x_dt(df,t,x_sol,n,dt);//t+Δt秒後の解
    
    double **Jt=J_l(t,x_sol,n);//時刻tにおけるヤコビ行列
    double **Jt_2=J_l(t+dt/2.,x_next_2,n);//t+Δt/2のヤコビ行列
    double **Jt_1=J_l(t+dt,x_next_1,n);//t+Δtのヤコビ行列

    /*======行列K1の計算======*/
    double **K1=AB(Jt,n,n,W_orth,n,n);//JtX
    /*======================*/

    /*======行列K2の計算======*/
    double **K1_X=dmatrix(1,n,1,n);//X+ΔtK1/2
    for(int i=1;i<=n;i++){
      for(int j=1;j<=n;j++){
        K1_X[i][j]=W_orth[i][j]+dt*K1[i][j]/2.;
      }
    }
    double **K2=AB(Jt_2,n,n,K1_X,n,n);//J(t+Δt/2)(X+ΔtK1/2)
    free_dmatrix(K1_X,1,n,1,n);
    /*======================*/

    /*======行列K3の計算======*/
    double **K2_X=dmatrix(1,n,1,n);//X+ΔtK2/2
    for(int i=1;i<=n;i++){
      for(int j=1;j<=n;j++){
        K2_X[i][j]=W_orth[i][j]+dt*K2[i][j]/2.;
      }
    }
    double **K3=AB(Jt_2,n,n,K2_X,n,n);//J(t+Δt/2)(X+ΔtK2/2)
    free_dmatrix(K2_X,1,n,1,n);
    /*======================*/

    /*======行列K4の計算======*/
    double **K3_X=dmatrix(1,n,1,n);//X+ΔtK3
    for(int i=1;i<=n;i++){
      for(int j=1;j<=n;j++){
        K3_X[i][j]=W_orth[i][j]+dt*K3[i][j];
      }
    }
    double **K4=AB(Jt_1,n,n,K3_X,n,n);//J(t+Δt)(X+ΔtK3)
    free_dmatrix(K3_X,1,n,1,n);
    /*======================*/

    //X=(dt/6)*(K1+2K2+2K3+K4);
    for(int i=1;i<=n;i++){
      for(int j=1;j<=n;j++){
        W_orth[i][j]+=(dt/6.)*(K1[i][j]+2.*K2[i][j]+2.*K3[i][j]+K4[i][j]);
      }
    }

    for(int i=1;i<=n;i++){
      for(int j=1;j<=n;j++){
        X[i][j]=W_orth[i][j];
      }
    }

    free_dmatrix(K1,1,n,1,n);
    free_dmatrix(K2,1,n,1,n);
    free_dmatrix(K3,1,n,1,n);
    free_dmatrix(K4,1,n,1,n);

    free_dvector(x_next_1,1,n);
    free_dvector(x_next_2,1,n);
    free_dmatrix(Jt_1,1,n,1,n);
    free_dmatrix(Jt_2,1,n,1,n);
    free_dmatrix(Jt,1,n,1,n);
    /*=======================================*/
  
    /*================収束確認=======================*/
    if(t>=T){
      double trA=0.0;//リアプノフ指数の収束の確認
      double **Jt=J_l(t,x_sol,n);//時刻tにおけるヤコビ行列
      for(int i=1;i<=n;i++)trA+=Jt[i][i];
      free_dmatrix(Jt,1,n,1,n);

      double sum_riap=0.0;
      for(int i=1;i<=n;i++)sum_riap+=riap[i]/(t-T+1);

      // printf("%f:",t);
      for(int i=1;i<=n;i++){
        printf("%f,",riap[i]/(t-T+1));
      }
      printf("%f\n",fabs(trA-sum_riap));
      

      // 収束判定
      if(t>=Sym_L){
        printf("finished,%f\n",t);
        for(int i=1;i<=n;i++)riap[i]=riap[i]/(t-T+1);
        break;
      }
    }
    

    /*============解の更新=============*/
    double *x_new=x_dt(df,t,x_sol,n,dt);//xは時刻tでの解
    for(int i=1;i<=n;i++)x_sol[i]=x_new[i];//時刻の更新
    free_dvector(x_new,1,n);
    /*===============================*/

    t+=dt;
    free_dmatrix(W_orth,1,n,1,n);
    printf("\n");
  }
  fclose(fp);

  free_dvector(x_sol,1,n);
  free_dmatrix(X,1,n,1,n);

  return riap;
}
  

/*=============扱う方程式についての記述===========================*/
double **J(double t,double *x,int n){//xは時刻tにおける解,ローレンツ方程式のヤコビ行列
  double **Ret=dmatrix(1,n,1,n);

  Ret[1][1]=-sigma,Ret[1][2]=sigma,Ret[1][3]=0;
  Ret[2][1]=gamma-x[3],Ret[2][2]=-1.,Ret[2][3]=-x[1];
  Ret[3][1]=x[2],Ret[3][2]=x[1],Ret[3][3]=-beta;

  // Ret[1][1]=x[1],Ret[1][2]=0.0,Ret[1][3]=0.0;
  // Ret[2][1]=0.0,Ret[2][2]=x[2],Ret[2][3]=0.0;
  // Ret[3][1]=0.0,Ret[3][2]=0.0,Ret[3][3]=x[3];

  return Ret;
}

double Lorenz(double t,double *x,int i,int n){//ローレンツ方程式の右辺
  double a=0.0;
  if(i==1){
    a=-sigma*x[1]+sigma*x[2];
  }else if(i==2){
    a=-x[1]*x[3]+gamma*x[1]-x[2];
  }else if(i==3){
    a=x[1]*x[2]-beta*x[3];
  }else{
    a=0.0;
  }
  return a;
}

int main(void){
  int n=3;
  double *lambda=Lyapunov_exponent(&Lorenz,J,n);
  for(int i=1;i<=n;i++){
    printf("%f\n",lambda[i]);
  }
  free_dvector(lambda,1,n);
  return 0;
}