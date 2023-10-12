#include "Include/functions.h"
#include<glsc3d_3.h>//glsc3Dを用いるときに使用

//標準座標系のパラメータ
#define back_x 1500
#define back_y 1000

//仮想座標系の範囲指定
// #define x0 10
// #define x1 100//グラフの横幅の長さ，最大値と最小値の差
#define base_x 10//グラフの開始位置x
#define base_y 10//グラフの開始位置y

//仮想座標系のサイズ
#define size_x 0.9*back_x
#define size_y 0.9*back_y

void plot(double t);

void plot_set(double b_x,double b_y){
    /*--------ここは脳死でコピペ----------------*/
    // g_enable_highdpi();
    // g_set_antialiasing(4);//ubuntuの場合は必要ない
    g_init("Burgers eq",b_x,b_y);//標準座標系の作成。上部の文字日本語対応
    g_scr_color(1,1,1);//標準座標系の背景カラーrgb
    g_cls();//標準座標系の初期化
    g_capture_set("Plot_fig");//画像保存
}

//パラメータ param，上下が[y_min,y_max]であるようなグラフの描画
void glaph_outline(double x_min,double x_max,double y_min,double y_max,char *title,char *xlabel,char *ylabel){
    
    int j=0;//グラフの配置の仕方．基本0でOK

    /*=================枠の描画=================================*/
    g_line_color(0,0,0,1);//グラフの枠の色
    g_line_width(2);//文字の太さ
    g_def_scale_2D(j,x_min,x_max,y_min,y_max,base_x,base_y,size_x,size_y);//標準座標系の定義
    g_sel_scale(j);
    g_box_2D(x_min,x_max,y_min,y_max,G_YES,G_NO);//箱の描画
    /*=========================================================*/

    /*========================縦軸，横軸の描画=================================*/
    g_text_size(30);
    //縦軸の範囲設定
    char y_M[50];
    sprintf(y_M,"%0.1f",y_max);
    g_text_standard(base_x+size_x,base_y+20,y_M);//縦軸の最大値の描画
    char y_m[50];
    sprintf(y_m,"%0.1f",y_min);
    g_text_standard(base_x+size_x,base_y+size_y,y_m);//縦軸の最小値の描画

    //横軸の範囲設定
    // char x_M[50];
    sprintf(y_M,"%0.1f",x_max);
    g_text_standard(base_x+0.98*size_x,base_y+1.05*size_y,y_M);//縦軸の最大値の描画
    // char x_m[50];
    sprintf(y_m,"%0.1f",x_min);
    g_text_standard(base_y,base_y+1.05*size_y,y_m);//縦軸の最小値の描画
    
    g_text_size(50);
    g_line_width(10);//線の太さ
    g_text_standard(base_y+30,base_x+50,title);//タイトル

    g_text_standard(base_x+size_x,base_y+size_y/2,ylabel);//縦軸のラベル
    g_text_standard(base_x+size_x/2,base_y+size_y+70,xlabel);//横軸のラベル
    /*======================================================================*/

    // for(int k=sta;k<length;k+=1){
    //     double Re=param[k];
    //     char str[5];
    //     snprintf(str,5,"%lf",Re);   
    //     g_text_standard(base_x+(param[k]-x_min)*size_x/(x_max-x_min),0+size_y+40,str);//横軸の描画
    // }

    // //注釈の表示
    // int base=50;
    // int delta=50;
    // g_text_color(1,0,0,1);
    // g_text_standard(back_x-300,base,"RedLine：Circle");//赤線の描画
    // // g_text_standard(back_x-300,base+delta,"Error Value");//赤線の描画
    // g_text_color(0,1,0,1);
    // g_text_standard(back_x-300,base+delta,"GreenLine：Curl");//緑線の描画

    // g_text_color(0,0,1,1);
    // g_text_standard(back_x-300,base+2*delta,"BlueLine：narrow");//青線の描画
    // g_text_color(0,0,0,1);
    // // g_text_standard(back_x-300,base+2*delta,"黒色：");//黒線の描画
    // // g_text_standard(back_x-300,base+3*delta,"Training Data");//黒線の描画
    // g_text_standard(back_x-300,base+3*delta,"Re：Reynolds Number");//kの描画
}

void legend(char *legend,int *color){
    g_text_color(color[0],color[1],color[2],1);
    double legend_x=base_x+0.8*size_x;
    double legend_y=base_y+0.2*size_y;
    g_text_standard(legend_x,legend_y,legend);//赤線の描画
}

//形状に応じたグラフの色の変化
void color_change(char *Answer_data){
    int base=50;
    int delta=50;
    if(strcmp(Answer_data,"Dataset")==0){
        g_line_color(1,0,0,1);//赤が円柱
        g_marker_color(1,0,0,1);//赤が円柱
        g_text_color(1,0,0,1);
        g_text_standard(back_x-300,base,"RedLine：Circle");//赤線の描画
    }else if(strcmp(Answer_data,"narrow")==0){
        g_line_color(0,0,1,1);//赤が円柱
        g_marker_color(0,0,1,1);

        g_text_color(0,0,1,1);
        g_text_standard(back_x-300,base+2*delta,"BlueLine：narrow");//青線の描画
       
    }else if(strcmp(Answer_data,"stressfree")==0){
        g_line_color(0,1,0,1);
        g_marker_color(0,1,0,1);

        g_text_color(0,1,0,1);
        g_text_standard(back_x-300,base+delta,"GreenLine：Curl");//緑線の描画
    }else{
        g_line_color(0,0,0,1);
    }
}

//param,Answerのグラフ
void data_plot(double *Answer,double *param,int start,int length){
    double min_err=1000;
    int min_Re=1000;

    for(int i=start;i<length-1;i++){
        g_move_2D(param[i],Answer[i]);
        g_plot_2D(param[i+1],Answer[i+1]);
    
        //点のプロット
        g_marker_size(20);//点の大きさ指定

        g_marker_2D(param[i],Answer[i]);//点のプロットはg_marker_2D
        if(Answer[i]<min_err){
            min_err=Answer[i];
            min_Re=param[i];
        }
    }

    g_marker_color(0,0,0,1);//点の色指定RGB
    g_marker_2D(min_Re,min_err);//最小値の色変換
}

void fileplot_2D(char *filename){//ファイル名
    
    //ファイル読み込み
    FILE *fp;
    fp=fopen(filename,"r");
    if(fp==NULL){
        printf("To open plotting file is failed\n");
        exit(1);
    }
    /*==========仮想座標系の縦横の幅=============*/
    double x_min=-1.0,x_max=1.0;
    double y_min=0.95,y_max=1.05;

    // int j=0;
    // g_line_color(0,0,0,1);//グラフの枠の色
    // g_line_width(2);//文字の太さ
    // g_def_scale_2D(j,x_min,x_max,y_min,y_max,base_x,base_y,size_x,size_y);//標準座標系の定義
    // g_sel_scale(j);
    // g_box_2D(x_min,x_max,y_min,y_max,G_YES,G_NO);//箱の描画

    g_text_size(50);
    glaph_outline(x_min,x_max,y_min,y_max,"Burgers","x","ux");

    g_text_size(20);
    g_text_standard(base_x+0.02*size_x,base_y+0.2*size_y,filename);//縦軸のラベル
    /*========================================*/

    double x,ux;//格納するための変数
    fscanf(fp,"%lf %lf\n",&x,&ux);//初期条件の読み込み

    double x_new,ux_new;//次のステップの解
    while(fscanf(fp,"%lf %lf\n",&x_new,&ux_new)!=EOF){
        g_move_2D(x,ux);
        g_plot_2D(x_new,ux_new);
        // printf("%f %f\n",x,ux);
        x=x_new,ux=ux_new;
    }

    fclose(fp);
}

void plot(double t){
    g_finish();//ここまでのg_関数を実行
    g_sleep(t);//描画の時間，-1の時は無限に
}

int main(void){

    int T=0;
    plot_set(back_x,back_y);
    
    int k=0;
    while(k<2100){
        g_cls();
        char str[50];
        
        sprintf(str,"solution/sol%d.txt",T);
       
        fileplot_2D(str);

        plot(0.01);
        g_capture();
        
        T++;
        if(T%5==0){
            k+=1;
        }
    }
    
    return 0;
}