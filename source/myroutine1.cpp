#include "myroutine1.h"

void Fsc(double Ksc,double ds,double sx2,double sz2,double x,double z,double &Fscx, double &Fscz)
{
	int flag; //flag=0:sx>sz, 1:sz>sx, 2:sx~sz, 3:sx~sz, w<<1
	double sx=sqrt(sx2);
	double sz=sqrt(sz2);
	double eps=1-sz/sx;
	if(fabs(eps) < 0.01){
		double p=x*x/(sx*(sx+sz));
		double q=z*z/(sx*(sx+sz));
//swap or not!!
		double w=p+q;
		double f0=0, f1=0, f2=0, f3=0, f4=0, f5=0, f6=0, f7=0, f8=0, f9=0, f10=0;
		if(w>0.001){
			flag=2; //A15
			/*
			   double f0 =      1./pow(w,1) *(1-exp(-w)*1);
			   double f1 =      1./pow(w,2) *(1-exp(-w)*(1+w));
			   double f2 =      2./pow(w,3) *(1-exp(-w)*(1+w+pow(w,2)/2));
			   double f3 =      6./pow(w,4) *(1-exp(-w)*(1+w+pow(w,2)/2+pow(w,3)/6));
			   double f4 =     24./pow(w,5) *(1-exp(-w)*(1+w+pow(w,2)/2+pow(w,3)/6+pow(w,4)/24));
			   double f5 =    120./pow(w,6) *(1-exp(-w)*(1+w+pow(w,2)/2+pow(w,3)/6+pow(w,4)/24+pow(w,5)/120));
			   double f6 =    720./pow(w,7) *(1-exp(-w)*.... lazy)
			   */
			double sum=0,wn=1,wn1=w;
			int n=0,nfactorial=1;
			f0=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f1=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f2=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f3=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f4=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f5=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f6=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f7=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f8=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f9=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f10=nfactorial/wn1*(1-exp(-w)*sum);
		} else {
			flag=3; //A19
			for(int k=0;k<=10;k++) f0+=pow(-w,k)/(0+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f1+=pow(-w,k)/(1+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f2+=pow(-w,k)/(2+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f3+=pow(-w,k)/(3+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f4+=pow(-w,k)/(4+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f5+=pow(-w,k)/(5+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f6+=pow(-w,k)/(6+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f7+=pow(-w,k)/(7+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f8+=pow(-w,k)/(8+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f9+=pow(-w,k)/(9+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f10+=pow(-w,k)/(10+k+1)/factorial(k);
		}
		//cout<<"f0="<<f0<<endl<<"f1="<<f1<<endl<<"f2="<<f2<<endl<<"f3="<<f3<<endl<<"f4="<<f4<<endl<<"f5="<<f5<<endl<<"f6="<<f6<<endl;
		double fx=sx/(sx+sz)*(f0+eps*(p-3*q)/2*f2
				+pow(eps,2)*(-2*q*f3+pow((p-3*q),2)/8*f4)
				+pow(eps,3)*(-5*q/2*f4-(p-3*q)*q*f5+pow((p-3*q),3)/48*f6)
				+pow(eps,4)*(-3*q*f5-(5*p*q-23*q*q)/4*f6-pow((p-3*q),3)/4*q*f7+pow((p-3*q),4)/384*f8)
				+pow(eps,5)*(-7/2*q*f6-(3*p*q-19*q*q)/2*f7-(5*p*p*q-46*p*q*q+93*q*q*q)/16*f8-pow((p-3*q),3)*q/24*f9+pow((p-3*q),5)/3840*f10));
		double fz=sz2/(sx*(sx+sz))*(f0+eps*(2*f1+(p-3*q)/2*f2)
				+pow(eps,2)*(3*f2+(p-5*q)*f3+pow((p-3*q),2)/8*f4)
				+pow(eps,3)*(4*f3+(3*p-22*q)/2*f4+(p*p-10*p*q+21*q*q)/2*f5+pow((p-3*q),3)/48*f6)
				+pow(eps,4)*(5*f4+(2*p-20*q)*f5+(3*p*p-44*p*q+121*q*q)/8*f6+(pow((p-3*q),3)-6*pow((p-3*q),2))/24*f7+pow((p-3*q),4)/384*f8)
				+pow(eps,5)*(6*f5+(5*p-65*q)/2*f6+(p*p-20*p*q-69*q*q)/2*f7+(p*p*p-22*p*p*q+121*p*q*q-192*q*q*q)/16*f8+(pow((p-3*q),4)-8*pow((p-3*q),3)*q)/192*f9+pow((p-3*q),5)/3840*f10));
		Fscx=Ksc*ds*x/sx2*fx;
		Fscz=Ksc*ds*z/sz2*fz;
	}else{
		flag=0; //sx2>sz2
		if(sx2<sz2){
			flag=1;
			swap(sx2,sz2);
		}

		complex<double> J(0,-1);
		double r=sqrt(sz2/sx2);
		double a=x/sqrt(2*(sx2-sz2));
		double b=z/sqrt(2*(sx2-sz2));
		if(flag==1){
			swap(a,b);
			b*=-1;
		}

		//complex<double> Fsc=J*ds*Ksc*sqrt(M_PI/2/(sx2-sz2))*(Faddeeva::w(a+J*b)-exp(-pow(a+J*b,2)+pow(a*r+J*b/r,2))*Faddeeva::w(a*r+J*b/r));
		complex<double> Fsc=J*ds*Ksc*sqrt(M_PI/2/(sx2-sz2))*exp(-pow(a+J*b,2))*(Faddeeva::erfc(J*(a+J*b))-Faddeeva::erfc(J*(a*r+J*b/r)));
		Fscx=real(Fsc);
		Fscz=imag(Fsc);

		if(flag==1){
			swap(Fscx,Fscz);
			Fscx*=-1;
		}
		if(isnan(Fscx)) Fscx=0;
		if(isnan(Fscz)) Fscz=0;
		if(isinf(Fscx)) Fscx=0;
		if(isinf(Fscz)) Fscz=0;
	}
//cout<<"flag="<<flag<<endl;
}

//use w function
void Fsc2(double Ksc,double ds,double sx2,double sz2,double x,double z,double &Fscx, double &Fscz)
{
	int flag; //flag=0:sx>sz, 1:sz>sx, 2:sx~sz, 3:sx~sz, w<<1
	double sx=sqrt(sx2);
	double sz=sqrt(sz2);
	double eps=1-sz/sx;
	if(fabs(eps) < 0.01){
		double p=x*x/(sx*(sx+sz));
		double q=z*z/(sx*(sx+sz));
//swap or not!!
		double w=p+q;
		double f0=0, f1=0, f2=0, f3=0, f4=0, f5=0, f6=0, f7=0, f8=0, f9=0, f10=0;
		if(w>0.001){
			flag=2; //A15
			/*
			   double f0 =      1./pow(w,1) *(1-exp(-w)*1);
			   double f1 =      1./pow(w,2) *(1-exp(-w)*(1+w));
			   double f2 =      2./pow(w,3) *(1-exp(-w)*(1+w+pow(w,2)/2));
			   double f3 =      6./pow(w,4) *(1-exp(-w)*(1+w+pow(w,2)/2+pow(w,3)/6));
			   double f4 =     24./pow(w,5) *(1-exp(-w)*(1+w+pow(w,2)/2+pow(w,3)/6+pow(w,4)/24));
			   double f5 =    120./pow(w,6) *(1-exp(-w)*(1+w+pow(w,2)/2+pow(w,3)/6+pow(w,4)/24+pow(w,5)/120));
			   double f6 =    720./pow(w,7) *(1-exp(-w)*.... lazy)
			   */
			double sum=0,wn=1,wn1=w;
			int n=0,nfactorial=1;
			f0=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f1=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f2=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f3=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f4=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f5=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f6=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f7=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f8=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f9=nfactorial/wn1*(1-exp(-w)*sum);
			n+=1;nfactorial*=n;wn=wn1;wn1*=w;sum+=wn/nfactorial;
			f10=nfactorial/wn1*(1-exp(-w)*sum);
		} else {
			flag=3; //A19
			for(int k=0;k<=10;k++) f0+=pow(-w,k)/(0+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f1+=pow(-w,k)/(1+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f2+=pow(-w,k)/(2+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f3+=pow(-w,k)/(3+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f4+=pow(-w,k)/(4+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f5+=pow(-w,k)/(5+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f6+=pow(-w,k)/(6+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f7+=pow(-w,k)/(7+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f8+=pow(-w,k)/(8+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f9+=pow(-w,k)/(9+k+1)/factorial(k);
			for(int k=0;k<=10;k++) f10+=pow(-w,k)/(10+k+1)/factorial(k);
		}
		//cout<<"f0="<<f0<<endl<<"f1="<<f1<<endl<<"f2="<<f2<<endl<<"f3="<<f3<<endl<<"f4="<<f4<<endl<<"f5="<<f5<<endl<<"f6="<<f6<<endl;
		double fx=sx/(sx+sz)*(f0+eps*(p-3*q)/2*f2
				+pow(eps,2)*(-2*q*f3+pow((p-3*q),2)/8*f4)
				+pow(eps,3)*(-5*q/2*f4-(p-3*q)*q*f5+pow((p-3*q),3)/48*f6)
				+pow(eps,4)*(-3*q*f5-(5*p*q-23*q*q)/4*f6-pow((p-3*q),3)/4*q*f7+pow((p-3*q),4)/384*f8)
				+pow(eps,5)*(-7/2*q*f6-(3*p*q-19*q*q)/2*f7-(5*p*p*q-46*p*q*q+93*q*q*q)/16*f8-pow((p-3*q),3)*q/24*f9+pow((p-3*q),5)/3840*f10));
		double fz=sz2/(sx*(sx+sz))*(f0+eps*(2*f1+(p-3*q)/2*f2)
				+pow(eps,2)*(3*f2+(p-5*q)*f3+pow((p-3*q),2)/8*f4)
				+pow(eps,3)*(4*f3+(3*p-22*q)/2*f4+(p*p-10*p*q+21*q*q)/2*f5+pow((p-3*q),3)/48*f6)
				+pow(eps,4)*(5*f4+(2*p-20*q)*f5+(3*p*p-44*p*q+121*q*q)/8*f6+(pow((p-3*q),3)-6*pow((p-3*q),2))/24*f7+pow((p-3*q),4)/384*f8)
				+pow(eps,5)*(6*f5+(5*p-65*q)/2*f6+(p*p-20*p*q-69*q*q)/2*f7+(p*p*p-22*p*p*q+121*p*q*q-192*q*q*q)/16*f8+(pow((p-3*q),4)-8*pow((p-3*q),3)*q)/192*f9+pow((p-3*q),5)/3840*f10));
		Fscx=Ksc*ds*x/sx2*fx;
		Fscz=Ksc*ds*z/sz2*fz;
	}else{
		flag=0; //sx2>sz2
		if(sx2<sz2){
			flag=1;
			swap(sx2,sz2);
		}

		complex<double> J(0,-1);
		double r=sqrt(sz2/sx2);
		double a=x/sqrt(2*(sx2-sz2));
		double b=z/sqrt(2*(sx2-sz2));
		if(flag==1){
			swap(a,b);
			b*=-1;
		}


		//test with divergence
		complex<double> probe1;
		int flag1=0;//0: converge, 1:diverge at x->+inf, 2:at z->+inf, 3:x->-inf, 4:z->-inf
		double a1,b1;
		a1=1/sqrt(2*(sx2-sz2));
		b1=0;
		probe1=(Faddeeva::w(a1+J*b1)-exp(-pow(a1+J*b1,2)+pow(a1*r+J*b1/r,2))*Faddeeva::w(a1*r+J*b1/r));
		if(isnan(real(probe1)) or isinf(real(probe1))) flag1=1;
		a1=0;
		b1=1/sqrt(2*(sx2-sz2));
		probe1=(Faddeeva::w(a1+J*b1)-exp(-pow(a1+J*b1,2)+pow(a1*r+J*b1/r,2))*Faddeeva::w(a1*r+J*b1/r));
		if(isnan(real(probe1)) or isinf(real(probe1))) flag1=2;
		a1=-1/sqrt(2*(sx2-sz2));
		b1=0;
		probe1=(Faddeeva::w(a1+J*b1)-exp(-pow(a1+J*b1,2)+pow(a1*r+J*b1/r,2))*Faddeeva::w(a1*r+J*b1/r));
		if(isnan(real(probe1)) or isinf(real(probe1))) flag1=3;
		a1=0;
		b1=-1/sqrt(2*(sx2-sz2));
		probe1=(Faddeeva::w(a1+J*b1)-exp(-pow(a1+J*b1,2)+pow(a1*r+J*b1/r,2))*Faddeeva::w(a1*r+J*b1/r));
		if(isnan(real(probe1)) or isinf(real(probe1))) flag1=4;

		complex<double> Fsc;
		if(flag1==1){
			if(a>0){
				a=-a;
				Fsc=J*ds*Ksc*sqrt(M_PI/2/(sx2-sz2))*(Faddeeva::w(a+J*b)-exp(-pow(a+J*b,2)+pow(a*r+J*b/r,2))*Faddeeva::w(a*r+J*b/r));
				Fscx=-real(Fsc);
				Fscz=imag(Fsc);
			} else{
				Fsc=J*ds*Ksc*sqrt(M_PI/2/(sx2-sz2))*(Faddeeva::w(a+J*b)-exp(-pow(a+J*b,2)+pow(a*r+J*b/r,2))*Faddeeva::w(a*r+J*b/r));
				Fscx=real(Fsc);
				Fscz=imag(Fsc);
			}
		}else if(flag1==2){
			if(b>0){
				b=-b;
				Fsc=J*ds*Ksc*sqrt(M_PI/2/(sx2-sz2))*(Faddeeva::w(a+J*b)-exp(-pow(a+J*b,2)+pow(a*r+J*b/r,2))*Faddeeva::w(a*r+J*b/r));
				Fscx=real(Fsc);
				Fscz=-imag(Fsc);
			} else{
				Fsc=J*ds*Ksc*sqrt(M_PI/2/(sx2-sz2))*(Faddeeva::w(a+J*b)-exp(-pow(a+J*b,2)+pow(a*r+J*b/r,2))*Faddeeva::w(a*r+J*b/r));
				Fscx=real(Fsc);
				Fscz=imag(Fsc);
			}
		}else if(flag1==3){
			if(a<0){
				a=-a;
				Fsc=J*ds*Ksc*sqrt(M_PI/2/(sx2-sz2))*(Faddeeva::w(a+J*b)-exp(-pow(a+J*b,2)+pow(a*r+J*b/r,2))*Faddeeva::w(a*r+J*b/r));
				Fscx=-real(Fsc);
				Fscz=imag(Fsc);
			} else{
				Fsc=J*ds*Ksc*sqrt(M_PI/2/(sx2-sz2))*(Faddeeva::w(a+J*b)-exp(-pow(a+J*b,2)+pow(a*r+J*b/r,2))*Faddeeva::w(a*r+J*b/r));
				Fscx=real(Fsc);
				Fscz=imag(Fsc);
			}
		}else if(flag1==4){
			if(b<0){
				b=-b;
				Fsc=J*ds*Ksc*sqrt(M_PI/2/(sx2-sz2))*(Faddeeva::w(a+J*b)-exp(-pow(a+J*b,2)+pow(a*r+J*b/r,2))*Faddeeva::w(a*r+J*b/r));
				Fscx=real(Fsc);
				Fscz=-imag(Fsc);
			} else{
				Fsc=J*ds*Ksc*sqrt(M_PI/2/(sx2-sz2))*(Faddeeva::w(a+J*b)-exp(-pow(a+J*b,2)+pow(a*r+J*b/r,2))*Faddeeva::w(a*r+J*b/r));
				Fscx=real(Fsc);
				Fscz=imag(Fsc);
			}
		}else{
			//cout<<"flag1=0! ??? divergence test is converged"<<endl;
		}

		if(flag==1){
			swap(Fscx,Fscz);
			Fscx*=-1;
		}

	}
//cout<<"flag="<<flag<<endl;
}

