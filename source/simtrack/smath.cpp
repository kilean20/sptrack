#include "smath.h"
using namespace std;
//===========================================
//
//          math functions
//
//===========================================
long int  fac(int n)
{
  if (n<=0) {
    cout<<"n should be non-negative."<<endl;
    exit(1);
  }
  if(n==1) return 1;
  if(n==2) return 2;
  if(n==3) return 6;
  if(n==4) return 24;
  if(n==5) return 120;
  if(n==6) return 720;
  if(n==7) return 5040;
  if(n==8) return 40320;
  if(n==9) return 362880;
  if(n==10) return 3628800;  
  if(n >10) {
    cout<<" fac(n), n>10. exit."<<endl;
    exit(1);
  }   
  return 0;
}

//---random number generator
double rnd(const double & r)
{
  int m;
  double s,u,v,p,r1;
  s=65536.0; u=2053.0; v=13849.0;
  m=int(r/s);  r1=r-m*s;  r1=u*r1+v; 
  m=(int)(r1/s); r1=r1-m*s; p=r1/s;
  return(p);
}

double guass(double u,double g, double & r)
{ int i,m;
 double s,w,v,t;
 s=65536.0; w=2053.0; v=13849.0;
 t=0.0;
 for (i=1; i<=12; i++)
   { r=r*w+v; m=(int)(r/s);
   r=r-m*s; t=t+(r)/s;
   }
 t=u+g*(t-6.0);
 return(t);
}

//---polynominal fitting
void pfit(double x[],double y[], int n, double a[], int m)
{ 
  int m1,i,j,l,ii,k,im,ix[21];
  double h[21],ha,hh,y1,y2,h1,h2,d,hm;
  for (i=0; i<=m; i++) a[i]=0.0;
  if (m>=n) m=n-1;
  if (m>=20) m=19;
  m1=m+1;
  ha=0.0;
  ix[0]=0; ix[m]=n-1;
  l=(n-1)/m; j=l;
  for (i=1; i<=m-1; i++)
    { ix[i]=j; j=j+l;}
  while (1==1)
    { hh=1.0;
    for (i=0; i<=m; i++)
      { a[i]=y[ix[i]]; h[i]=-hh; hh=-hh;}
    for (j=1; j<=m; j++)
      { ii=m1; y2=a[ii-1]; h2=h[ii-1];
      for (i=j; i<=m; i++)
	{ d=x[ix[ii-1]]-x[ix[m1-i-1]];
	y1=a[m-i+j-1];
	h1=h[m-i+j-1];
	a[ii-1]=(y2-y1)/d;
	h[ii-1]=(h2-h1)/d;
	ii=m-i+j; y2=y1; h2=h1;
	}
      }
    hh=-a[m]/h[m];
    for (i=0; i<=m; i++)
      a[i]=a[i]+h[i]*hh;
    for (j=1; j<=m-1; j++)
      { ii=m-j; d=x[ix[ii-1]];
      y2=a[ii-1];
      for (k=m1-j; k<=m; k++)
	{ y1=a[k-1]; a[ii-1]=y2-d*y1;
	y2=y1; ii=k;
	}
      }
    hm=fabs(hh);
    if (hm<=ha) { a[m]=-hm; return;}
    a[m]=hm; ha=hm; im=ix[0]; h1=hh;
    j=0;
    for (i=0; i<=n-1; i++)
      { if (i==ix[j])
	{ if (j<m) j=j+1;}
      else
	{ h2=a[m-1];
	for (k=m-2; k>=0; k--)
	  h2=h2*x[i]+a[k];
	h2=h2-y[i];
	if (fabs(h2)>hm)
	  { hm=fabs(h2); h1=h2; im=i;}
	}
      }
    if (im==ix[0]) return;
    i=0;l=1;
    while (l==1)
      { l=0;
      if (im>=ix[i])
	{ i=i+1;
	if (i<=m) l=1;
	}
      }
    if (i>m) i=m;
    if (i==(i/2)*2) h2=-hh;
    else h2=hh;
    if (h1*h2>=0.0) ix[i]=im;
    else
      { if (im<ix[0])
	{ for (j=m-1; j>=0; j--)
	  ix[j+1]=ix[j];
	ix[0]=im;
	}
      else
	{ if (im>ix[m])
	  { for (j=1; j<=m; j++)
	    ix[j-1]=ix[j];
	  ix[m]=im;
	  }
	else ix[i-1]=im;
	}
      }
    }
}

//---matrix kits
void mat_mult(double A[],  double B[], double C[], int m, int n, int k)
{
  int i, j, l;

  for(i=0; i<m;i++)
    for(j=0;j<k;j++) C[i*m+j]=0.;

  for(i=0; i<m;i++)
    for(j=0;j<k;j++)
      for (l=0;l<n;l++) C[i*m+j]=C[i*m+j] + A[i*m+l]*B[l*n+j];
}

double mat_det(double a[],int n)
{
  int i,j,k,is=0,js=0,l,u,v;
  double f,det,q,d;
  f=1.0; det=1.0;
  for (k=0; k<=n-2; k++)
    { q=0.0;
    for (i=k; i<=n-1; i++)
      for (j=k; j<=n-1; j++)
	{ l=i*n+j; d=fabs(a[l]);
	if (d>q) { q=d; is=i; js=j;}
	}
    if (q+1.0==1.0)
      { det=0.0; return(det);}
    if (is!=k)
      { f=-f;
      for (j=k; j<=n-1; j++)
	{ u=k*n+j; v=is*n+j;
	d=a[u]; a[u]=a[v]; a[v]=d;
	}
      }
    if (js!=k)
      { f=-f;
      for (i=k; i<=n-1; i++)
	{ u=i*n+js; v=i*n+k;
	d=a[u]; a[u]=a[v]; a[v]=d;
	}
      }
    l=k*n+k;
    det=det*a[l];
    for (i=k+1; i<=n-1; i++)
      { d=a[i*n+k]/a[l];
      for (j=k+1; j<=n-1; j++)
	{ u=i*n+j;
	a[u]=a[u]-d*a[k*n+j];
	}
      }
    }
  det=f*det*a[n*n-1];
  return(det);
}

int mat_inv(double a[],int n)
{
  int *is,*js,i,j,k,l,u,v;
  double d,p;
  is=new int[n];
  js=new int[n];
  for (k=0; k<=n-1; k++)
    { d=0.0;
    for (i=k; i<=n-1; i++)
      for (j=k; j<=n-1; j++)
	{ l=i*n+j; p=abs(a[l]);
	if (p>d) { d=p; is[k]=i; js[k]=j;}
	}
    if (d+1.0==1.0)
      { free(is); free(js); printf("err**not inv\n"); exit(1);
      return(0);
      }
    if (is[k]!=k)
      for (j=0; j<=n-1; j++)
	{ u=k*n+j; v=is[k]*n+j;
	p=a[u]; a[u]=a[v]; a[v]=p;
	}
    if (js[k]!=k)
      for (i=0; i<=n-1; i++)
	{ u=i*n+k; v=i*n+js[k];
	p=a[u]; a[u]=a[v]; a[v]=p;
	}
    l=k*n+k;
    a[l]=1.0/a[l];
    for (j=0; j<=n-1; j++)
      if (j!=k)
	{ u=k*n+j; a[u]=a[u]*a[l];}
    for (i=0; i<=n-1; i++)
      if (i!=k)
	for (j=0; j<=n-1; j++)
	  if (j!=k)
	    { u=i*n+j;
	    a[u]=a[u]-a[i*n+k]*a[k*n+j];
	    }
    for (i=0; i<=n-1; i++)
      if (i!=k)
	{ u=i*n+k; a[u]=-a[u]*a[l];}
    }
  for (k=n-1; k>=0; k--)
    { if (js[k]!=k)
      for (j=0; j<=n-1; j++)
	{ u=k*n+j; v=js[k]*n+j;
	p=a[u]; a[u]=a[v]; a[v]=p;
	}
    if (is[k]!=k)
      for (i=0; i<=n-1; i++)
	{ u=i*n+k; v=i*n+is[k];
	p=a[u]; a[u]=a[v]; a[v]=p;
	}
    }
  free(is); free(js);
  return(1);
}

//----eigen tune solver
void mat_change_hessenberg(double a[], int n)
{ 
  int i,j,k,u,v;
  double d,t;
  for (k=1; k<=n-2; k++)
    { d=0.0;
    for (j=k; j<=n-1; j++)
      { u=j*n+k-1; t=a[u];
      if (fabs(t)>fabs(d))
	{ d=t; i=j;}
      }
    if (fabs(d)+1.0!=1.0)
      { if (i!=k)
	{ for (j=k-1; j<=n-1; j++)
	  { u=i*n+j; v=k*n+j;
	  t=a[u]; a[u]=a[v]; a[v]=t;
	  }
	for (j=0; j<=n-1; j++)
	  { u=j*n+i; v=j*n+k;
	  t=a[u]; a[u]=a[v]; a[v]=t;
	  }
	}
      for (i=k+1; i<=n-1; i++)
	{ u=i*n+k-1; t=a[u]/d; a[u]=0.0;
	for (j=k; j<=n-1; j++)
	  { v=i*n+j;
	  a[v]=a[v]-t*a[k*n+j];
	  }
	for (j=0; j<=n-1; j++)
	  { v=j*n+k;
	  a[v]=a[v]+t*a[j*n+i];
	  }
	}
      }
    }
  return;
}

int mat_root_hessenberg(double a[],int n,double u[],double v[],double eps,int jt)
{ 
  int m,it,i,j,k,l,ii,jj,kk,ll;
  double b,c,w,g,xy,p,q,r,x,s,e,f,z,y;
  it=0; m=n;
  while (m!=0)
    { l=m-1;
    while ((l>0)&&(fabs(a[l*n+l-1])>eps*
		   (fabs(a[(l-1)*n+l-1])+fabs(a[l*n+l])))) l=l-1;
    ii=(m-1)*n+m-1; jj=(m-1)*n+m-2;
    kk=(m-2)*n+m-1; ll=(m-2)*n+m-2;
    if (l==m-1)
      { u[m-1]=a[(m-1)*n+m-1]; v[m-1]=0.0;
      m=m-1; it=0;
      }
    else if (l==m-2)
      { b=-(a[ii]+a[ll]);
      c=a[ii]*a[ll]-a[jj]*a[kk];
      w=b*b-4.0*c;
      y=sqrt(fabs(w));
      if (w>0.0)
	{ xy=1.0;
	if (b<0.0) xy=-1.0;
	u[m-1]=(-b-xy*y)/2.0;
	u[m-2]=c/u[m-1];
	v[m-1]=0.0; v[m-2]=0.0;
	}
      else
	{ u[m-1]=-b/2.0; u[m-2]=u[m-1];
	v[m-1]=y/2.0; v[m-2]=-v[m-1];
	}
      m=m-2; it=0;
      }
    else
      { if (it>=jt)
	{ printf("fail\n");
	return(-1);
	}
      it=it+1;
      for (j=l+2; j<=m-1; j++)
	a[j*n+j-2]=0.0;
      for (j=l+3; j<=m-1; j++)
	a[j*n+j-3]=0.0;
      for (k=l; k<=m-2; k++)
	{ if (k!=l)
	  { p=a[k*n+k-1]; q=a[(k+1)*n+k-1];
	  r=0.0;
	  if (k!=m-2) r=a[(k+2)*n+k-1];
	  }
	else
	  { x=a[ii]+a[ll];
	  y=a[ll]*a[ii]-a[kk]*a[jj];
	  ii=l*n+l; jj=l*n+l+1;
	  kk=(l+1)*n+l; ll=(l+1)*n+l+1;
	  p=a[ii]*(a[ii]-x)+a[jj]*a[kk]+y;
	  q=a[kk]*(a[ii]+a[ll]-x);
	  r=a[kk]*a[(l+2)*n+l+1];
	  }
	if ((fabs(p)+fabs(q)+fabs(r))!=0.0)
	  { xy=1.0;
	  if (p<0.0) xy=-1.0;
	  s=xy*sqrt(p*p+q*q+r*r);
	  if (k!=l) a[k*n+k-1]=-s;
	  e=-q/s; f=-r/s; x=-p/s;
	  y=-x-f*r/(p+s);
	  g=e*r/(p+s);
	  z=-x-e*q/(p+s);
	  for (j=k; j<=m-1; j++)
	    { ii=k*n+j; jj=(k+1)*n+j;
	    p=x*a[ii]+e*a[jj];
	    q=e*a[ii]+y*a[jj];
	    r=f*a[ii]+g*a[jj];
	    if (k!=m-2)
	      { kk=(k+2)*n+j;
	      p=p+f*a[kk];
	      q=q+g*a[kk];
	      r=r+z*a[kk]; a[kk]=r;
	      }
	    a[jj]=q; a[ii]=p;
	    }
	  j=k+3;
	  if (j>=m-1) j=m-1;
	  for (i=l; i<=j; i++)
	    { ii=i*n+k; jj=i*n+k+1;
	    p=x*a[ii]+e*a[jj];
	    q=e*a[ii]+y*a[jj];
	    r=f*a[ii]+g*a[jj];
	    if (k!=m-2)
	      { kk=i*n+k+2;
	      p=p+f*a[kk];
	      q=q+g*a[kk];
	      r=r+z*a[kk]; a[kk]=r;
	      }
	    a[jj]=q; a[ii]=p;
	    }
	  }
	}
      }
    }
  return(1);
}

void LinearEquations(double a1, double b1, double c1, double a2, double b2, double c2, double & x, double & y  )
// --the equations to be solved
//     a1*x +b1*y= c1
//     a2*x +b2*y= c2
{
  if(a1*b2-a2*b1 != 0. ){
    x=(c1*b2-c2*b1) / ( a1*b2-a2*b1);
    y=(c1*a2-c2*a1) / ( a2*b1-a1*b2);
  }
}

//---eigen solver from GSL
void  EigenSolver(double Matrix[4][4] , double wr[4], double  wi[4], double vr[4][4], double vi[4][4]) 
{
  int i, j;
  double data[16];

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      data[i*4+j]= Matrix[i][j];
  
  gsl_matrix_view m  = gsl_matrix_view_array (data, 4, 4);
  gsl_vector_complex *eval = gsl_vector_complex_alloc (4);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (4, 4);
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (4);
  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
  gsl_eigen_nonsymmv_free (w);
  gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
  
  for (i = 0; i < 4; i++)
    {
      gsl_complex eval_i = gsl_vector_complex_get (eval, i);
      gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);
      wr[i]=  	GSL_REAL(eval_i);
      wi[i]=  	GSL_IMAG(eval_i);
      
      for (j = 0; j < 4; j++)
	{
	  gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);	
	  vr[i][j]=GSL_REAL(z);
	  vi[i][j]=GSL_IMAG(z);
	}
    }
  
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
}

void  EigenSolver_6D(double Matrix[6][6] , double wr[6], double  wi[6], double vr[6][6], double vi[6][6]) 
{
  int i, j;
  double data[36];

  for(i=0;i<6;i++)
    for(j=0;j<6;j++)
      data[i*6+j]= Matrix[i][j];
  
  gsl_matrix_view m  = gsl_matrix_view_array (data, 6, 6);
  gsl_vector_complex *eval = gsl_vector_complex_alloc (6);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (6, 6);
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (6);
  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
  gsl_eigen_nonsymmv_free (w);
  gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
  
  for (i = 0; i < 6; i++)
    {
      gsl_complex eval_i = gsl_vector_complex_get (eval, i);
      gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);
      wr[i]=  	GSL_REAL(eval_i);
      wi[i]=  	GSL_IMAG(eval_i);
      
      for (j = 0; j < 6; j++)
	{
	  gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);	
	  vr[i][j]=GSL_REAL(z);
	  vi[i][j]=GSL_IMAG(z);
	}
    }
  
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
}
//---SVD stuffs
void brmul(double a[], double b[] , int m, int n, int k, double c[] )
{ 
  int i,j,l,u;
  for (i=0; i<=m-1; i++)
    for (j=0; j<=k-1; j++)
      { u=i*k+j; c[u]=0.0;
        for (l=0; l<=n-1; l++)
          c[u]=c[u]+a[i*n+l]*b[l*k+j];
      }
  return;
}

static void ppp(double *a,double *e,double *s,double *v,int m,int n)
{ 
  int i,j,p,q;
  double d;
  if (m>=n) i=n;
 else i=m;
 for (j=1; j<=i-1; j++)
   { a[(j-1)*n+j-1]=s[j-1];
   a[(j-1)*n+j]=e[j-1];
      }
    a[(i-1)*n+i-1]=s[i-1];
    if (m<n) a[(i-1)*n+i]=e[i-1];
    for (i=1; i<=n-1; i++)
    for (j=i+1; j<=n; j++)
      { p=(i-1)*n+j-1; q=(j-1)*n+i-1;
        d=v[p]; v[p]=v[q]; v[q]=d;
      }
    return;
  }

static void sss(double fg[2],double cs[2])
{ double r,d;
 if ((fabs(fg[0])+fabs(fg[1]))==0.0)
      { cs[0]=1.0; cs[1]=0.0; d=0.0;}
    else 
      { d=sqrt(fg[0]*fg[0]+fg[1]*fg[1]);
        if (fabs(fg[0])>fabs(fg[1]))
          { d=fabs(d);
            if (fg[0]<0.0) d=-d;
          }
        if (fabs(fg[1])>=fabs(fg[0]))
          { d=fabs(d);
            if (fg[1]<0.0) d=-d;
          }
        cs[0]=fg[0]/d; cs[1]=fg[1]/d;
      }
    r=1.0;
    if (fabs(fg[0])>fabs(fg[1])) r=cs[1];
    else
      if (cs[0]!=0.0) r=1.0/cs[0];
    fg[0]=d; fg[1]=r;
    return;
  }

int bmuav(double a[],int m,int n,double *u,double *v,double eps,int ka)
{ int i,j,k,l,it,ll,kk,ix,iy,mm,nn,iz,m1,ks;
 double d,dd,t,sm,sm1,em1,sk,ek,b,c,shh,fg[2],cs[2];
 double *s,*e,*w;
 
 s=new  double[ka];
 e=new  double[ka];
 w=new  double[ka];
 
 it=60; k=n;
 if (m-1<n) k=m-1;
 l=m;
 if (n-2<m) l=n-2;
 if (l<0) l=0;
 ll=k;
 if (l>k) ll=l;
 if (ll>=1)
   { for (kk=1; kk<=ll; kk++)
     { if (kk<=k)
       { d=0.0;
       for (i=kk; i<=m; i++)
	 { ix=(i-1)*n+kk-1; d=d+a[ix]*a[ix];}
       s[kk-1]=sqrt(d);
       if (s[kk-1]!=0.0)
	 { ix=(kk-1)*n+kk-1;
	 if (a[ix]!=0.0)
	   { s[kk-1]=fabs(s[kk-1]);
	   if (a[ix]<0.0) s[kk-1]=-s[kk-1];
	   }
	 for (i=kk; i<=m; i++)
	   { iy=(i-1)*n+kk-1;
	   a[iy]=a[iy]/s[kk-1];
	   }
	 a[ix]=1.0+a[ix];
	 }
       s[kk-1]=-s[kk-1];
       }
     if (n>=kk+1)
       { for (j=kk+1; j<=n; j++)
	 { if ((kk<=k)&&(s[kk-1]!=0.0))
	   { d=0.0;
	   for (i=kk; i<=m; i++)
	     { ix=(i-1)*n+kk-1;
	     iy=(i-1)*n+j-1;
	     d=d+a[ix]*a[iy];
	     }
	   d=-d/a[(kk-1)*n+kk-1];
	   for (i=kk; i<=m; i++)
	     { ix=(i-1)*n+j-1;
	     iy=(i-1)*n+kk-1;
	     a[ix]=a[ix]+d*a[iy];
	     }
	   }
                    e[j-1]=a[(kk-1)*n+j-1];
	 }
              }
            if (kk<=k)
              { for (i=kk; i<=m; i++)
                  { ix=(i-1)*m+kk-1; iy=(i-1)*n+kk-1;
                    u[ix]=a[iy];
                  }
              }
            if (kk<=l)
              { d=0.0;
                for (i=kk+1; i<=n; i++)
                  d=d+e[i-1]*e[i-1];
                e[kk-1]=sqrt(d);
                if (e[kk-1]!=0.0)
                  { if (e[kk]!=0.0)
                      { e[kk-1]=fabs(e[kk-1]);
                        if (e[kk]<0.0) e[kk-1]=-e[kk-1];
                      }
                    for (i=kk+1; i<=n; i++)
                      e[i-1]=e[i-1]/e[kk-1];
                    e[kk]=1.0+e[kk];
                  }
                e[kk-1]=-e[kk-1];
                if ((kk+1<=m)&&(e[kk-1]!=0.0))
                  { for (i=kk+1; i<=m; i++) w[i-1]=0.0;
                    for (j=kk+1; j<=n; j++)
                      for (i=kk+1; i<=m; i++)
                        w[i-1]=w[i-1]+e[j-1]*a[(i-1)*n+j-1];
                    for (j=kk+1; j<=n; j++)
                      for (i=kk+1; i<=m; i++)
                        { ix=(i-1)*n+j-1;
                          a[ix]=a[ix]-w[i-1]*e[j-1]/e[kk];
                        }
                  }
                for (i=kk+1; i<=n; i++)
                  v[(i-1)*n+kk-1]=e[i-1];
              }
          }
      }
    mm=n;
    if (m+1<n) mm=m+1;
    if (k<n) s[k]=a[k*n+k];
    if (m<mm) s[mm-1]=0.0;
    if (l+1<mm) e[l]=a[l*n+mm-1];
    e[mm-1]=0.0;
    nn=m;
    if (m>n) nn=n;
    if (nn>=k+1)
      { for (j=k+1; j<=nn; j++)
          { for (i=1; i<=m; i++)
              u[(i-1)*m+j-1]=0.0;
            u[(j-1)*m+j-1]=1.0;
          }
      }
    if (k>=1)
      { for (ll=1; ll<=k; ll++)
          { kk=k-ll+1; iz=(kk-1)*m+kk-1;
            if (s[kk-1]!=0.0)
              { if (nn>=kk+1)
                  for (j=kk+1; j<=nn; j++)
                    { d=0.0;
                      for (i=kk; i<=m; i++)
                        { ix=(i-1)*m+kk-1;
                          iy=(i-1)*m+j-1;
                          d=d+u[ix]*u[iy]/u[iz];
                        }
                      d=-d;
                      for (i=kk; i<=m; i++)
                        { ix=(i-1)*m+j-1;
                          iy=(i-1)*m+kk-1;
                          u[ix]=u[ix]+d*u[iy];
                        }
                    }
                  for (i=kk; i<=m; i++)
                    { ix=(i-1)*m+kk-1; u[ix]=-u[ix];}
                  u[iz]=1.0+u[iz];
                  if (kk-1>=1)
                    for (i=1; i<=kk-1; i++)
                      u[(i-1)*m+kk-1]=0.0;
              }
            else
              { for (i=1; i<=m; i++)
                  u[(i-1)*m+kk-1]=0.0;
                u[(kk-1)*m+kk-1]=1.0;
              }
          }
      }
    for (ll=1; ll<=n; ll++)
      { kk=n-ll+1; iz=kk*n+kk-1;
        if ((kk<=l)&&(e[kk-1]!=0.0))
          { for (j=kk+1; j<=n; j++)
              { d=0.0;
                for (i=kk+1; i<=n; i++)
                  { ix=(i-1)*n+kk-1; iy=(i-1)*n+j-1;
                    d=d+v[ix]*v[iy]/v[iz];
                  }
                d=-d;
                for (i=kk+1; i<=n; i++)
                  { ix=(i-1)*n+j-1; iy=(i-1)*n+kk-1;
                    v[ix]=v[ix]+d*v[iy];
                  }
              }
          }
        for (i=1; i<=n; i++)
          v[(i-1)*n+kk-1]=0.0;
        v[iz-n]=1.0;
      }
    for (i=1; i<=m; i++)
    for (j=1; j<=n; j++)
      a[(i-1)*n+j-1]=0.0;
    m1=mm; it=60;
    while (1==1)
      { if (mm==0)
          { ppp(a,e,s,v,m,n);
            delete s; delete e; delete w; return(1);
          }
        if (it==0)
          { ppp(a,e,s,v,m,n);
            delete s; delete e;  delete w; return(-1);
          }
        kk=mm-1;
	while ((kk!=0)&&(fabs(e[kk-1])!=0.0))
          { d=fabs(s[kk-1])+fabs(s[kk]);
            dd=fabs(e[kk-1]);
            if (dd>eps*d) kk=kk-1;
            else e[kk-1]=0.0;
          }
        if (kk==mm-1)
          { kk=kk+1;
            if (s[kk-1]<0.0)
              { s[kk-1]=-s[kk-1];
                for (i=1; i<=n; i++)
                  { ix=(i-1)*n+kk-1; v[ix]=-v[ix];}
              }
            while ((kk!=m1)&&(s[kk-1]<s[kk]))
              { d=s[kk-1]; s[kk-1]=s[kk]; s[kk]=d;
                if (kk<n)
                  for (i=1; i<=n; i++)
                    { ix=(i-1)*n+kk-1; iy=(i-1)*n+kk;
                      d=v[ix]; v[ix]=v[iy]; v[iy]=d;
                    }
                if (kk<m)
                  for (i=1; i<=m; i++)
                    { ix=(i-1)*m+kk-1; iy=(i-1)*m+kk;
                      d=u[ix]; u[ix]=u[iy]; u[iy]=d;
                    }
                kk=kk+1;
              }
            it=60;
            mm=mm-1;
          }
        else
          { ks=mm;
            while ((ks>kk)&&(fabs(s[ks-1])!=0.0))
              { d=0.0;
                if (ks!=mm) d=d+fabs(e[ks-1]);
                if (ks!=kk+1) d=d+fabs(e[ks-2]);
                dd=fabs(s[ks-1]);
                if (dd>eps*d) ks=ks-1;
                else s[ks-1]=0.0;
              }
            if (ks==kk)
              { kk=kk+1;
                d=fabs(s[mm-1]);
                t=fabs(s[mm-2]);
                if (t>d) d=t;
                t=fabs(e[mm-2]);
                if (t>d) d=t;
                t=fabs(s[kk-1]);
                if (t>d) d=t;
                t=fabs(e[kk-1]);
                if (t>d) d=t;
                sm=s[mm-1]/d; sm1=s[mm-2]/d;
                em1=e[mm-2]/d;
                sk=s[kk-1]/d; ek=e[kk-1]/d;
                b=((sm1+sm)*(sm1-sm)+em1*em1)/2.0;
                c=sm*em1; c=c*c; shh=0.0;
                if ((b!=0.0)||(c!=0.0))
                  { shh=sqrt(b*b+c);
                    if (b<0.0) shh=-shh;
                    shh=c/(b+shh);
                  }
                fg[0]=(sk+sm)*(sk-sm)-shh;
                fg[1]=sk*ek;
                for (i=kk; i<=mm-1; i++)
                  { sss(fg,cs);
                    if (i!=kk) e[i-2]=fg[0];
                    fg[0]=cs[0]*s[i-1]+cs[1]*e[i-1];
                    e[i-1]=cs[0]*e[i-1]-cs[1]*s[i-1];
                    fg[1]=cs[1]*s[i];
                    s[i]=cs[0]*s[i];
                    if ((cs[0]!=1.0)||(cs[1]!=0.0))
                      for (j=1; j<=n; j++)
                        { ix=(j-1)*n+i-1;
                          iy=(j-1)*n+i;
                          d=cs[0]*v[ix]+cs[1]*v[iy];
                          v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
                          v[ix]=d;
                        }
                    sss(fg,cs);
                    s[i-1]=fg[0];
                    fg[0]=cs[0]*e[i-1]+cs[1]*s[i];
                    s[i]=-cs[1]*e[i-1]+cs[0]*s[i];
                    fg[1]=cs[1]*e[i];
                    e[i]=cs[0]*e[i];
                    if (i<m)
                      if ((cs[0]!=1.0)||(cs[1]!=0.0))
                        for (j=1; j<=m; j++)
                          { ix=(j-1)*m+i-1;
                            iy=(j-1)*m+i;
                            d=cs[0]*u[ix]+cs[1]*u[iy];
                            u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
                            u[ix]=d;
                          }
                  }
                e[mm-2]=fg[0];
                it=it-1;
              }
            else
              { if (ks==mm)
                  { kk=kk+1;
                    fg[1]=e[mm-2]; e[mm-2]=0.0;
                    for (ll=kk; ll<=mm-1; ll++)
                      { i=mm+kk-ll-1;
                        fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        if (i!=kk)
                          { fg[1]=-cs[1]*e[i-2];
                            e[i-2]=cs[0]*e[i-2];
                          }
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
                          for (j=1; j<=n; j++)
                            { ix=(j-1)*n+i-1;
                              iy=(j-1)*n+mm-1;
                              d=cs[0]*v[ix]+cs[1]*v[iy];
                              v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
                              v[ix]=d;
                            }
                      }
                  }
                else
                  { kk=ks+1;
                    fg[1]=e[kk-2];
                    e[kk-2]=0.0;
                    for (i=kk; i<=mm; i++)
                      { fg[0]=s[i-1];
                        sss(fg,cs);
                        s[i-1]=fg[0];
                        fg[1]=-cs[1]*e[i-1];
                        e[i-1]=cs[0]*e[i-1];
                        if ((cs[0]!=1.0)||(cs[1]!=0.0))
                          for (j=1; j<=m; j++)
                            { ix=(j-1)*m+i-1;
                              iy=(j-1)*m+kk-2;
                              d=cs[0]*u[ix]+cs[1]*u[iy];
                              u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
                              u[ix]=d;
                            }
                      }
                  }
              }
          }
      }
  }

//------FFT, tune finders
void fft(int m, double*x, double*y)
{
   short int dir=1;
   int n,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   n = 1;
   for (i=0;i<m;i++) 
      n *= 2;

   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   c1 = -1.0; 
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0; 
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1; 
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1) 
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
   }
}

void  Sin2FFT(int n, double*xr, double*xi, int & nfft, double*famp, double*fphi)
{
  int i, m;
  double temp1,temp2;

  nfft= 1;
  m=0;
  do {
    m=m+1;
    nfft *= 2;} while (nfft <= n);  
  m=m-1;
  nfft=nfft/2;
  
  for (i=0;i<=nfft-1;i++){
    famp[i]=xr[i]*2.0*pow( sin(M_PI*(i+1)/nfft),2.); 
    fphi[i]=xi[i]*2.0*pow( sin(M_PI*(i+1)/nfft),2.);
  }
  
  fft(m,famp,fphi);
  
  for (i=0;i<nfft-1;i++){
    temp1=famp[i];
    temp2=fphi[i];
    famp[i]=sqrt(temp1*temp1+temp2*temp2);
    fphi[i]=atan2(temp2,temp1);
  }
}

void  FindMaxPeak(int n, double*xr, double*xi, int & nfft, double*famp, double*fphi, double& peakf, double& peakamp, double& peakphi, int flag)
{
  int i, k;
  double temp1,temp2;
  double peak1;
  double a,b,c;

  peak1=0.0;
  k=0;
  for (i=600;i<=750;i++){
    if(famp[i] > peak1){
      k=i;
      peak1=famp[i];
    }  
  }
  
  a=famp[k];
  if( famp[k+1] >  famp[k-1] ) {
    b=famp[k+1];
  }
  else{
    b=famp[k-1];
  }
  c=cos(2.0*M_PI/nfft);
  
  temp1=pow( c*(a+b),2.0)-2*a*b*(2*c*c-c-1);
  temp1=( -(a+b*c)*(a-b)+b*sqrt(temp1) )/( a*a+b*b+2.0*a*b*c );
  peakf=1.0*k/nfft+(1.0/2.0/M_PI)*asin(temp1 *sin(2.0*M_PI/nfft));
  //if(flag==1) peakf=1.-peakf;
  
  peakamp=0.;
  peakphi=0.;
  for (i=0;i<n;i++){
    peakamp=peakamp+xr[i]*cos(-2.0*M_PI*i*peakf) -xi[i]*sin(-2.0*M_PI*i*peakf);
    peakphi=peakphi+xr[i]*sin(-2.0*M_PI*i*peakf) + xi[i]*cos(-2.0*M_PI*i*peakf);
  }
  peakamp=peakamp/n;
  peakphi=peakphi/n;
  temp1=sqrt(peakamp*peakamp + peakphi*peakphi)*2.0;
  temp2=atan2(peakphi, peakamp);
  peakamp=temp1;
  peakphi=temp2;
}

void FineTuneFinder(int nfft, double *x, double & fpeak)
{
  int i;
  double xr[nfft], xi[nfft], famp[nfft], fphi[nfft];
  double peakfr0,peakamp, peakphi;
  
  for(i=0;i<nfft;i++){
    xr[i]=x[i];
    xi[i]=0.;
  }
  Sin2FFT(nfft, xr, xi, nfft, famp, fphi);
  FindMaxPeak(1024, xr, xi, nfft, famp, fphi, peakfr0, peakamp, peakphi,1);
  fpeak=peakfr0;  
}
