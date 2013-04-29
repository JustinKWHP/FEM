/*             ƽ������Ԫ���������            */
/*  nn:�ڵ�������ne:��Ԫ������nd:�����ɶ�����E:����ģ����
    nfix:�ڵ�λ��Լ���������ɶ�����anu:���ɱ�;
    T:�ṹ�ĺ�ȣ�gm:�ṹ���ϵ��ض�;
    NTYPE:NTYPE=1��ʾƽ��Ӧ��������Ϊƽ��Ӧ�䡣
    loc(I,1),loc(I,2),loc(I,3)����I����Ԫ��ʱ�����еĽڵ��
    cx(J),cy(J)����J���ڵ������
    ifix(K)����K������Լ�������ɶȺ�
    f(I)����I�����ɶ��ϵ�Լ����
    stres(I,J)������Ԫ��Ӧ��
*/              

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define  NN  6
#define  NE  4
#define  ND  12
#define  IFIX  6

/*   ��˹��ȥ�������Է�����   */  
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

/*   ����նȾ��󣬽ڵ�λ�ƺ͵�ԪӦ��   */
float CST(int loc[NE][3],float cx[NN],float cy[NN],int ifix[IFIX],float f[ND],float gk[ND][ND],float stres[NE][3],float bak[NE][3][6],int nn,int ne,int nd,int nfix,int NTYPE,float E,float anu,float t,float gm,char resultfile[10]){
	
	    FILE *result;
	    char *str;
         int m,n,l,unit,i1,i2,i3; 	
         float s,s2,b1;
	    //d[][] ���Ծ���bb[][] ��ԪӦ�����ba[3][6] ��ԪӦ������ek[][]��Ԫ�նȾ���xx[]��Ԫ�ڵ�λ��  ��be[] b(i,j,m)��ce[] c(i,j,m)��
	    float d[3][3],bb[3][6],ba[3][6],ek[6][6],be[3],ce[3],xx[6];
	    
	    //�ֽ�����նȾ�����0
	    for(m=0;m<ND;m++){
	    	    for(n=0;n<ND;n++){
	    	    	    gk[m][n]=0.0;
	    	    }
	    }
	    
	    //���㵯�Ծ���d[],NTYPE=1Ϊƽ��Ӧ��������Ϊƽ��Ӧ��
	    for(m=0;m<3;m++){
	    	    for(n=0;n<3;n++){
	    	    	    d[m][n]=0.0;
	    	    }
	    }
	    if(NTYPE!=1){               //ƽ��Ӧ������
	    	  E=E/(1.0-anu*anu);
	    	  anu=anu/(1.0-anu);
	    }
	    s=E/(1.0-anu*anu);
	    d[0][0]=s;
	    d[0][1]=s*anu;
	    d[1][1]=s;
	    d[1][0]=d[0][1];
	    d[2][2]=0.5*s*(1.0-anu);

         //��������նȾ���  
         for(unit=0;unit<NE;unit++){
          	 
          	//������bb[]Ϊ0
              for(m=0;m<3;m++){
	    	         for(n=0;n<6;n++){
	    	    	         bb[m][n]=0.0;
	    	         }
	         }
	          
	         // ��Ԫi,j,m
	         i1=loc[unit][0]-1;
	         i2=loc[unit][1]-1;
	         i3=loc[unit][2]-1;	          
	          
	         // ���� b(i,j,m�� c(i,j,m)
	         be[0]=cy[i2]-cy[i3];
	         be[1]=cy[i3]-cy[i1];
	         be[2]=cy[i1]-cy[i2];
          	ce[0]=cx[i3]-cx[i2];
          	ce[1]=cx[i1]-cx[i3];
              ce[2]=cx[i2]-cx[i1];
          	 
          	// 2����Ԫ���
          	s2=cx[i1]*be[0]+cx[i2]*be[1]+cx[i3]*be[2];
          	 
          	// ���㵥ԪӦ�����bb[][]
          	for(m=0;m<3;m++){
               	n=m*2+1;
          	 	l=n-1;
          	 	bb[0][l]=be[m]/s2;
          	 	bb[1][n]=ce[m]/s2;
          	 	bb[2][l]=bb[1][n];
          	 	bb[2][n]=bb[0][l];
          	}
          	  
          	// ����[ba]=[d][bb],��������bak[][][]��
          	for(m=0;m<3;m++){
          		for(n=0;n<6;n++){
          			ba[m][n]=0;
          			for(l=0;l<3;l++){
          				ba[m][n]=ba[m][n]+d[m][l]*bb[l][n];
          			}
          			bak[unit][m][n]=ba[m][n];
          		}
          	} 
          	
          	//������������Ľڵ��������ӵ�f[]
          	if(gm!=0){
          		for(m=0;m<3;m++){
          			n=loc[unit][m];
          			l=n*2-1;
          			f[l]=f[l]-t*gm*(0.5*s2)/3;
          		}
          	}
          	
          	// ���㵥Ԫ�նȾ���ek[][]
          	for(m=0;m<6;m++){
          		for(n=0;n<6;n++){
          			b1=0;
          			for(l=0;l<3;l++){
          				b1=b1+bb[l][m]*ba[l][n];
          			}
          			ek[m][n]=0.5*s2*b1*t;
          		}
          	}        	
          	
          	// ��װ����նȾ���
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
	    
	    //��������غ����������նȾ���
	    if((result=fopen(resultfile,"a"))==NULL){
	  	   printf("cannot open %s\n",resultfile);
	  	   exit(0);
	    }
	    // ��������غ�����
	    str="  �����غ�����\n";
	    fputs(str,result);
	    str="  NODE                 X��LOAD                         Y��LOAD\n";
	    fputs(str,result);
	    n=0;
	    for(m=0;m<NN;m++){	  	                  //д�������غ�����
	  	    n=m+1; 
	         fprintf(result,"      %d",n);
	         fprintf(result,"                 %9.6E",f[m*2]);
	         fprintf(result,"                %9.6E\n",f[m*2+1]);
	    }
	    str="\n";
	    fputs(str,result);
	    // �������նȾ���
	    str="  ����նȾ���\n";
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
	    
	    // ��Լ����Ӧgk[][]�����Խ���Ԫ�س˴���
	    for(m=0;m<nfix;m++){
	    	    n=ifix[m]-1;
	    	    gk[n][n]=gk[n][n]*1.0e15;
	    }
	    
	    // ���ýⷽ�������
	    GAUSS(gk,f);
	    
	    // ����ԪӦ��������0
	    for(m=0;m<NE;m++){
	    	    for(n=0;n<3;n++){
	    	    	    stres[m][n]=0;
	    	    }
	    }	    
	    
	    //���뵥Ԫ�ڵ�λ��
	    for(m=0;m<NE;m++){
	    	    for(n=0;n<3;n++){
	    	    	    xx[0]=f[2*loc[m][0]-2];
	    	    	    xx[1]=f[2*loc[m][0]-1];
	    	    	    xx[2]=f[2*loc[m][1]-2];
	    	    	    xx[3]=f[2*loc[m][1]-1];
	    	    	    xx[4]=f[2*loc[m][2]-2];
	    	    	    xx[5]=f[2*loc[m][2]-1];
	    	    	    // ���㵥ԪӦ��
	    	    	    for(l=0;l<6;l++){
	    	    	    	    stres[m][n]=stres[m][n]+bak[m][n][l]*xx[l];
	    	    	    }
	    	    }
	    }
	    
}

/*   ������   */
void main(){
	  
	  printf("\n");
	  printf("\n");
	  printf("****************���������������ε�Ԫ������Ԫ����******************");
	  printf("\n");
	  printf("\n");
	  printf("************************��д��:       ****************************");
	  printf("\n");
	  printf("\n");
	
	  //�������鼰����
	  int nn,ne,nd,nfix,NTYPE;
	  float E,anu,t,gm;
	  int loc[NE][3],ifix[IFIX];
	  float cx[NN],cy[NN],f[ND],gk[ND][ND],stres[NE][3],bak[NE][3][6];
	  int i,j,k;
	  
	  FILE *data,*result;                         //�����ļ�ָ��
	  char datafile[10],resultfile[10];
	  printf("�����������ļ����ļ���:\n");
	  printf("\n");
	  scanf("%s",datafile);
	  printf("\n");
	  printf("\n");
	  printf("���������ļ����ļ���:\n");
	  printf("\n");
	  scanf("%s",resultfile);
	  printf("\n");
	  printf("\n");
	  
	  //�������ļ��ж��������������Ԫ��š��ڵ����ꡢԼ�����ɶȼ��ڵ��غ�
	  if((data=fopen(datafile,"r"))==NULL){
	  	   printf("cannot open %s\n",datafile);
	  	   exit(0);
	  }
       while(!feof(data)){   	   	  
	           fscanf(data,"%d",&nn);          //�����������
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
	           	  	  fscanf(data,"%d",&loc[i][j]);                  //���뵥Ԫ�ڵ���
	           	  }
	           }          
		      for(i=0;i<NN;i++){
	           	  fscanf(data,"%f",&cx[i]);                           //����ڵ�����
	           	  fscanf(data,"%f",&cy[i]);
	           }    
	           for(j=0;j<IFIX;j++){                                       //����Լ�����ɶ�
	          	 fscanf(data,"%d",&ifix[j]);  
	           }  
	           for(k=0;k<ND;k++){                                       //����ڵ��غ�
	           	 fscanf(data,"%f",&f[k]);
	          }
       }	  
	  fclose(data);
	  
	  //����������д�뵽����ļ���
	  if((result=fopen(resultfile,"a"))==NULL){
	  	   printf("cannot open %s\n",resultfile);
	  	   exit(0);
	  }
	  char  *str;
	  str="  ��������\n";        //д���������
	  fputs(str,result);
	  str="  NN   NE     ND       NFIX            E                   ANU                     T                      GM             NTYPE\n";
	  fputs(str,result);
	  fprintf(result,"  %3d    %3d       %3d         %3d         %f         %f          %f          %f              %d\n\n",nn,ne,nd,nfix,E,anu,t,gm,NTYPE);
	  str="  ��Ԫ�ڵ���\n";                  //д�뵥Ԫ�ڵ���  
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
	  str="  �ڵ�����\n";                         //д��ڵ�����
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
	  str="  Լ�����ɶ�\n";                      //д��Լ�����ɶ�
	  fputs(str,result); 
	  for(j=0;j<IFIX;j++){
	       fprintf(result,"  %d      ",ifix[j]);
	  }
	  str="\n";
	  fputs(str,result);
	  str="\n";
	  fputs(str,result);	  
	  fclose(result);
	  
	  //��CST�ӳ��򣬼���նȾ��󡢽ڵ�λ�ƺ͵�ԪӦ��
	  CST(loc,cx,cy,ifix,f,gk,stres,bak,nn,ne,nd,nfix,NTYPE,E,anu,t,gm,resultfile);
	  
	  //��������д�뵽����ļ���
	  if((result=fopen(resultfile,"a"))==NULL){
	  	   printf("cannot open %s\n",resultfile);
	  	   exit(0);
	  }
	  //д��ڵ�λ��
	  str="  �ڵ�λ��\n";
	  fputs(str,result);
	  str="  NODE                  X��DISP                            Y��DISP\n";
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
	  
	  //д�뵥ԪӦ��
	  str="  ��ԪӦ��\n";
	  fputs(str,result);
	  str="  ELEMENT              X��STR                              Y��STR                            XY��STR\n";
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
	
	  //�������
	  printf("���������������!");                     	 
	  getchar(); 
	  getchar();  
}
