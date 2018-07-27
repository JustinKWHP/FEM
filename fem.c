# -*- coding:UTF-8 -*-
/*             平面有限元法计算程序            */
/*  nn:节点总数；ne:单元总数；nd:总自由度数；E:弹性模量；
    nfix:节点位移约束的总自由度数；anu:泊松比;
    T:结构的厚度；gm:结构材料的重度;
    NTYPE:NTYPE=1表示平面应力，其他为平面应变。
    loc(I,1),loc(I,2),loc(I,3)：第I个单元逆时针排列的节点号
    cx(J),cy(J)：第J个节点的坐标
    ifix(K)：第K个给定约束的自由度号
    f(I)：第I个自由度上的约束力
    stres(I,J)：各单元的应力
*/              

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define  NN  6
#define  NE  4
#define  ND  12
#define  IFIX  6

/*   高斯消去法解线性方程组   */  
float GAUSS(float gk[ND][ND],float f[ND]){
	
	    int i,i1,j,m;
	    for(i=0;i<ND;i++){
	    	    i1=i+1;
	    	    for(j=i1;j<ND;j++){
	    	    	    gk[i][j]=gk[i][j]/gk[i][i]; 
	    	    }
	    	    f[i]=f[i]/gk[i][i];
	    	    gk[i][i]=1.0;
	    	    for(j=i1;j<ND;j++){
	    	    	    for(m=i1;m<ND;m++){
	    	    	    	    gk[j][m]=gk[j][m]-gk[j][i]*gk[i][m];
	    	    	    }
	    	    	    f[j]=f[j]-gk[j][i]*f[i];
	    	    }
	    }
	    
	    for(i=ND-2;i>=0;i--){
	    	    for(j=i+1;j<ND;j++){
	    	    	    f[i]=f[i]-gk[i][j]*f[j];
	    	    }
	    }	   
	    
}

/*   计算刚度矩阵，节点位移和单元应力   */
float CST(int loc[NE][3],float cx[NN],float cy[NN],int ifix[IFIX],float f[ND],float gk[ND][ND],float stres[NE][3],float bak[NE][3][6],int nn,int ne,int nd,int nfix,int NTYPE,float E,float anu,float t,float gm,char resultfile[10]){
	
	    FILE *result;
	    char *str;
         int m,n,l,unit,i1,i2,i3; 	
         float s,s2,b1;
	    //d[][] 弹性矩阵；bb[][] 单元应变矩阵；ba[3][6] 单元应力矩阵；ek[][]单元刚度矩阵；xx[]单元节点位移  ；be[] b(i,j,m)；ce[] c(i,j,m)；
	    float d[3][3],bb[3][6],ba[3][6],ek[6][6],be[3],ce[3],xx[6];
	    
	    //现将整体刚度矩阵置0
	    for(m=0;m<ND;m++){
	    	    for(n=0;n<ND;n++){
	    	    	    gk[m][n]=0.0;
	    	    }
	    }
	    
	    //计算弹性矩阵d[],NTYPE=1为平面应力，否则为平面应变
	    for(m=0;m<3;m++){
	    	    for(n=0;n<3;n++){
	    	    	    d[m][n]=0.0;
	    	    }
	    }
	    if(NTYPE!=1){               //平面应变问题
	    	  E=E/(1.0-anu*anu);
	    	  anu=anu/(1.0-anu);
	    }
	    s=E/(1.0-anu*anu);
	    d[0][0]=s;
	    d[0][1]=s*anu;
	    d[1][1]=s;
	    d[1][0]=d[0][1];
	    d[2][2]=0.5*s*(1.0-anu);

         //计算整体刚度矩阵  
         for(unit=0;unit<NE;unit++){
          	 
          	//首先置bb[]为0
              for(m=0;m<3;m++){
	    	         for(n=0;n<6;n++){
	    	    	         bb[m][n]=0.0;
	    	         }
	         }
	          
	         // 单元i,j,m
	         i1=loc[unit][0]-1;
	         i2=loc[unit][1]-1;
	         i3=loc[unit][2]-1;	          
	          
	         // 计算 b(i,j,m） c(i,j,m)
	         be[0]=cy[i2]-cy[i3];
	         be[1]=cy[i3]-cy[i1];
	         be[2]=cy[i1]-cy[i2];
          	ce[0]=cx[i3]-cx[i2];
          	ce[1]=cx[i1]-cx[i3];
              ce[2]=cx[i2]-cx[i1];
          	 
          	// 2倍单元面积
          	s2=cx[i1]*be[0]+cx[i2]*be[1]+cx[i3]*be[2];
          	 
          	// 计算单元应变矩阵bb[][]
          	for(m=0;m<3;m++){
               	n=m*2+1;
          	 	l=n-1;
          	 	bb[0][l]=be[m]/s2;
          	 	bb[1][n]=ce[m]/s2;
          	 	bb[2][l]=bb[1][n];
          	 	bb[2][n]=bb[0][l];
          	}
          	  
          	// 计算[ba]=[d][bb],并保存在bak[][][]中
          	for(m=0;m<3;m++){
          		for(n=0;n<6;n++){
          			ba[m][n]=0;
          			for(l=0;l<3;l++){
          				ba[m][n]=ba[m][n]+d[m][l]*bb[l][n];
          			}
          			bak[unit][m][n]=ba[m][n];
          		}
          	} 
          	
          	//计算重力引起的节点力并叠加到f[]
          	if(gm!=0){
          		for(m=0;m<3;m++){
          			n=loc[unit][m];
          			l=n*2-1;
          			f[l]=f[l]-t*gm*(0.5*s2)/3;
          		}
          	}
          	
          	// 计算单元刚度矩阵ek[][]
          	for(m=0;m<6;m++){
          		for(n=0;n<6;n++){
          			b1=0;
          			for(l=0;l<3;l++){
          				b1=b1+bb[l][m]*ba[l][n];
          			}
          			ek[m][n]=0.5*s2*b1*t;
          		}
          	}        	
          	
          	// 组装整体刚度矩阵
          	int inode,nodei,idofn,nrows,nrowe,jnode,jdofn,nodej,ncols,ncole;
          	for(inode=0;inode<3;inode++){
          		nodei=loc[unit][inode];
          		for(idofn=0;idofn<2;idofn++){
          			nrows=(nodei-1)*2+idofn;
          			nrowe=inode*2+idofn;
          			for(jnode=0;jnode<3;jnode++){
          				nodej=loc[unit][jnode];
          				for(jdofn=0;jdofn<2;jdofn++){
          					ncols=(nodej-1)*2+jdofn;
          					ncole=jnode*2+jdofn;
          					gk[nrows][ncols]=gk[nrows][ncols]+ek[nrowe][ncole];
          				}
          			}
          		}
          	}
          	 
          }
	    
	    //输出整体载荷列阵和整体刚度矩阵
	    if((result=fopen(resultfile,"a"))==NULL){
	  	   printf("cannot open %s\n",resultfile);
	  	   exit(0);
	    }
	    // 输出整体载荷列阵
	    str="  整体载荷列阵\n";
	    fputs(str,result);
	    str="  NODE                 X—LOAD                         Y—LOAD\n";
	    fputs(str,result);
	    n=0;
	    for(m=0;m<NN;m++){	  	                  //写入整体载荷列阵
	  	    n=m+1; 
	         fprintf(result,"      %d",n);
	         fprintf(result,"                 %9.6E",f[m*2]);
	         fprintf(result,"                %9.6E\n",f[m*2+1]);
	    }
	    str="\n";
	    fputs(str,result);
	    // 输出整体刚度矩阵
	    str="  整体刚度矩阵\n";
	    fputs(str,result);
	    for(m=0;m<ND;m++){
	    	    for(n=0;n<ND;n++){
	    	    	    fprintf(result,"  %9.6E        ",gk[m][n]);
	    	    }
	    	    str="\n";
	         fputs(str,result);
	    }
	    str="\n";
	    fputs(str,result);
	    fclose(result);
	    
	    // 将约束对应gk[][]的主对角线元素乘大数
	    for(m=0;m<nfix;m++){
	    	    n=ifix[m]-1;
	    	    gk[n][n]=gk[n][n]*1.0e15;
	    }
	    
	    // 调用解方程组程序
	    GAUSS(gk,f);
	    
	    // 将单元应力矩阵置0
	    for(m=0;m<NE;m++){
	    	    for(n=0;n<3;n++){
	    	    	    stres[m][n]=0;
	    	    }
	    }	    
	    
	    //引入单元节点位移
	    for(m=0;m<NE;m++){
	    	    for(n=0;n<3;n++){
	    	    	    xx[0]=f[2*loc[m][0]-2];
	    	    	    xx[1]=f[2*loc[m][0]-1];
	    	    	    xx[2]=f[2*loc[m][1]-2];
	    	    	    xx[3]=f[2*loc[m][1]-1];
	    	    	    xx[4]=f[2*loc[m][2]-2];
	    	    	    xx[5]=f[2*loc[m][2]-1];
	    	    	    // 计算单元应力
	    	    	    for(l=0;l<6;l++){
	    	    	    	    stres[m][n]=stres[m][n]+bak[m][n][l]*xx[l];
	    	    	    }
	    	    }
	    }
	    
}

/*   主程序   */
void main(){
	  
	  printf("\n");
	  printf("\n");
	  printf("****************本程序用于三角形单元的有限元计算******************");
	  printf("\n");
	  printf("\n");
	  printf("************************编写者:       ****************************");
	  printf("\n");
	  printf("\n");
	
	  //定义数组及变量
	  int nn,ne,nd,nfix,NTYPE;
	  float E,anu,t,gm;
	  int loc[NE][3],ifix[IFIX];
	  float cx[NN],cy[NN],f[ND],gk[ND][ND],stres[NE][3],bak[NE][3][6];
	  int i,j,k;
	  
	  FILE *data,*result;                         //定义文件指针
	  char datafile[10],resultfile[10];
	  printf("请输入数据文件的文件名:\n");
	  printf("\n");
	  scanf("%s",datafile);
	  printf("\n");
	  printf("\n");
	  printf("请输入结果文件的文件名:\n");
	  printf("\n");
	  scanf("%s",resultfile);
	  printf("\n");
	  printf("\n");
	  
	  //从数据文件中读入基本参数、单元编号、节点坐标、约束自由度及节点载荷
	  if((data=fopen(datafile,"r"))==NULL){
	  	   printf("cannot open %s\n",datafile);
	  	   exit(0);
	  }
       while(!feof(data)){   	   	  
	           fscanf(data,"%d",&nn);          //读入基本参数
	           fscanf(data,"%d",&ne);
	           fscanf(data,"%d",&nd);
	           fscanf(data,"%d",&nfix);
	           fscanf(data,"%f",&E);
	           fscanf(data,"%f",&anu);
	           fscanf(data,"%f",&t);
	           fscanf(data,"%f",&gm);
	           fscanf(data,"%d",&NTYPE);                           
                for(i=0;i<NE;i++){                              
	           	  for(j=0;j<3;j++){ 
	           	  	  fscanf(data,"%d",&loc[i][j]);                  //读入单元节点编号
	           	  }
	           }          
		      for(i=0;i<NN;i++){
	           	  fscanf(data,"%f",&cx[i]);                           //读入节点坐标
	           	  fscanf(data,"%f",&cy[i]);
	           }    
	           for(j=0;j<IFIX;j++){                                       //读入约束自由度
	          	 fscanf(data,"%d",&ifix[j]);  
	           }  
	           for(k=0;k<ND;k++){                                       //读入节点载荷
	           	 fscanf(data,"%f",&f[k]);
	          }
       }	  
	  fclose(data);
	  
	  //将基本数据写入到输出文件中
	  if((result=fopen(resultfile,"a"))==NULL){
	  	   printf("cannot open %s\n",resultfile);
	  	   exit(0);
	  }
	  char  *str;
	  str="  基本参数\n";        //写入基本参数
	  fputs(str,result);
	  str="  NN   NE     ND       NFIX            E                   ANU                     T                      GM             NTYPE\n";
	  fputs(str,result);
	  fprintf(result,"  %3d    %3d       %3d         %3d         %f         %f          %f          %f              %d\n\n",nn,ne,nd,nfix,E,anu,t,gm,NTYPE);
	  str="  单元节点编号\n";                  //写入单元节点编号  
	  fputs(str,result);
	  str="  UNIT         I             J            M\n";
	  fputs(str,result);
	  k=0;
	  for(i=0;i<NE;i++){  
	  	  k=i+1;
	  	  fprintf(result,"      %d",k);                             
	       for(j=0;j<3;j++){ 
	            fprintf(result,"            %d",loc[i][j]);    	            	            	                  
	       }
	       str="\n";
	       fputs(str,result); 
	  }
	  str="\n";
	  fputs(str,result);
	  str="  节点坐标\n";                         //写入节点坐标
	  fputs(str,result);
	  str="  NODE                 X                        Y\n";
	  fputs(str,result); 
	  k=0;
	  for(i=0;i<NN;i++){
	  	  k=k+1;
	  	  fprintf(result,"      %d",k);
	       fprintf(result,"                 %7.2f                %7.2f\n",cx[i],cy[i]);                          	          	 
	  }
	  str="\n";
	  fputs(str,result);    	   	 
	  str="  约束自由度\n";                      //写入约束自由度
	  fputs(str,result); 
	  for(j=0;j<IFIX;j++){
	       fprintf(result,"  %d      ",ifix[j]);
	  }
	  str="\n";
	  fputs(str,result);
	  str="\n";
	  fputs(str,result);	  
	  fclose(result);
	  
	  //调CST子程序，计算刚度矩阵、节点位移和单元应力
	  CST(loc,cx,cy,ifix,f,gk,stres,bak,nn,ne,nd,nfix,NTYPE,E,anu,t,gm,resultfile);
	  
	  //将计算结果写入到输出文件中
	  if((result=fopen(resultfile,"a"))==NULL){
	  	   printf("cannot open %s\n",resultfile);
	  	   exit(0);
	  }
	  //写入节点位移
	  str="  节点位移\n";
	  fputs(str,result);
	  str="  NODE                  X—DISP                            Y—DISP\n";
	  fputs(str,result);
	  j=0;
	  for(i=0;i<NN;i++){	  	                 
	  	  j=i+1; 
	       fprintf(result,"      %d",j);
	       fprintf(result,"                 %9.6E",f[i*2]);
	       fprintf(result,"                %9.6E\n",f[i*2+1]);
	    }
	    str="\n";
	    fputs(str,result);                         
	  
	  //写入单元应力
	  str="  单元应力\n";
	  fputs(str,result);
	  str="  ELEMENT              X—STR                              Y—STR                            XY—STR\n";
	  fputs(str,result);
	  for(i=0;i<NE;i++){
	  	  k=i+1;
	  	  fprintf(result,"          %d",k);
	  	  for(j=0;j<3;j++){
	  	  	  fprintf(result,"                %9.6E",stres[i][j]);
	  	  }
	  	  str="\n";
	  	  fputs(str,result);
	  }
	  fclose(result);  
	
	  //程序结束
	  printf("按任意键结束程序!");                     	 
	  getchar(); 
	  getchar();  
}
