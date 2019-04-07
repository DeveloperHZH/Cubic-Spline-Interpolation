#include<stdio.h>
#include<math.h>
#define N 20
double *Gauss(int n,double Matrix[N][N],double Value[N],double *solve)
{
	double tmp;
	for(int i=0;i<n;i++) { //In this case, it is no need to confirm that Matrix[i][i]!=0 for Matrix[i][i] always equals 2 initially.
		for(int j=i+1;j<=n;j++) {
			tmp=Matrix[j][i]/Matrix[i][i];
			Value[j]-=Value[i]*tmp;
			for(int k=i;k<=n;k++) Matrix[j][k]-=Matrix[i][k]*tmp;
		}
	}
	for(int i=n;i>=0;i--) {
		solve[i]=Value[i];
		for(int j=n;j>i;j--) solve[i]-=solve[j]*Matrix[i][j];
		solve[i]/=Matrix[i][i];
	}
	//for(int i=0;i<=n;i++) printf("%7.2lf",solve[i]);printf("\n");
	return(solve);
}
int POutput(int n,double M[N],double h[N],double x[N],double y[N])
{
	printf("\n");
	double A,B,C,D;
	A=B=C=D=0;
	for(int i=0;i<n;i++) {
		printf("In condition that %7.2lf<x<%.2lf  ,S(x)=",x[i],x[i+1]);
		A=(M[i+1]-M[i])/(6*h[i]);
		B=(x[i+1]*M[i]-x[i]*M[i+1])/(2*h[i]);
		C=(x[i]*x[i]*M[i+1]-x[i+1]*x[i+1]*M[i])/(2*h[i])-y[i]/h[i]+h[i]*M[i]/6+y[i+1]/h[i]-h[i]*M[i+1]/6;
		D=(y[i]/h[i]-h[i]*M[i]/6)*x[i+1]-x[i]*(y[i+1]/h[i]-h[i]*M[i+1]/6)+x[i+1]*x[i+1]*x[i+1]*M[i]/(6*h[i])-x[i]*x[i]*x[i]*M[i+1]/(6*h[i]);
		printf("%10.4lfx^3  ",A);
		if(fabs(B)<1e-8) {} else {
			if(B>0) {
				printf("+%10.4lfx^2  ",B);
			} else printf("-%10.4lfx^2  ",-B);
		}
		if(fabs(C)<1e-8) {} else {
			if(C>0) {
				printf("+%10.4lfx  ",C);
			} else printf("-%10.4lfx  ",-C);
		}
		if(fabs(D)<1e-8) {} else {
			if(D>0) {
				printf("+%10.4lf  ",D);
			} else printf("-%10.4lf  ",-D);
			printf("\n");
		}
	}
	return(0);
}
int main()
{
	double x[N],y[N],h[N],miu[N],lambda[N],d[N],solve[N],M[N],value[N];
	double matrix[N][N]={0};
	double m0,mn;
	/*for(int i=0;i<=N;i++) {
		for(int j=0;j<=N;j++) printf("%5.1f",matrix[i][j]);
		printf("\n");
	}*/	 
	int n=0,flag=1;
	printf("Please input the data as \'x y¨L\', and end with \'666 666\'.\n");
	while(flag) {
		scanf("%lf %lf",&x[n],&y[n]);
		if((fabs(x[n]-666)<1e-8)&&(fabs(y[n]-666)<1e-8)){
			flag=0;
		} else {n++;}
	} 
	n--;
	if(n==-1){
		printf("Invalid data!\n");
	} else {
		for(int i=0;i<n;i++) { //Initiate
			h[i]=x[i+1]-x[i];
			if(i!=0) {
				lambda[i]=h[i]/(h[i]+h[i-1]);
				miu[i]=1-lambda[i];
				d[i]=6*((y[i+1]-y[i])/h[i]-(y[i]-y[i-1])/h[i-1])/(h[i]+h[i-1]);
				//printf("%.2lf %.2lf %.2lf\n",h[i],lambda[i],d[i]);
			}	
		}
		printf("Please input the type of boundary conditions(with number 1,2 or 3).\n"); 
		printf("Type 1: Giving M0 and Mn.\n");
		printf("Type 2: Giving m0 and mn.\n");
		printf("Type 3: Unfinished.\n");
		int type;scanf("%d",&type);
		switch(type) {
			case 1: {
				printf("Please input M0 and Mn:\n");
				scanf("%lf%lf",&M[0],&M[n]);
				matrix[0][0]=2;matrix[0][1]=lambda[1];
				matrix[n-2][n-3]=miu[n-1];matrix[n-2][n-2]=2;
				for(int i=1;i<n-2;i++) {
					matrix[i][i-1]=miu[i+1];matrix[i][i+1]=lambda[1];matrix[i][i]=2;
				}
				value[0]=d[1]-miu[1]*M[0];value[n-2]=d[n-1]-lambda[n-1]*M[n];
				for(int i=1;i<n-2;i++) value[i]=d[i+1];
	 		 	/*for(int i=0;i<N;i++) {
					for(int j=0;j<N;j++) printf("%6.1lf",matrix[i][j]);
					printf("\n");
				}
				for(int i=0;i<=n-2;i++) printf("%7.2lf",value[i]=d[i+1]);printf("\n"); */
				Gauss(n-2,matrix,value,M+sizeof(double)/8);
				//for(int i=0;i<=n;i++) printf("%7.2lf",M[i]);
				POutput(n,M,h,x,y);
				break;
			}
			case 2: {
				printf("Please input m0 and mn:\n");
				scanf("%lf%lf",&m0,&mn);
				matrix[0][0]=2;matrix[0][1]=1;
				matrix[n][n]=2;matrix[n][n-1]=1;
				for(int i=1;i<n;i++) {
					matrix[i][i-1]=miu[i];
					matrix[i][i]=2;
					matrix[i][i+1]=lambda[i];
				}
				d[0]=((y[0]-y[1])/(x[0]-x[1])-m0)*6/h[0];d[n]=(mn-(y[n-1]-y[n])/(x[n-1]-x[n]))*6/h[n-1];
				Gauss(n,matrix,d,M);
				POutput(n,M,h,x,y);
				break;
			}
			case 3: {
				printf("Unfinished function.\n");
				break;
			}
			default:{
				printf("Invalid input!\n");
				break;
			}
		}
	}
	return 0;
}
