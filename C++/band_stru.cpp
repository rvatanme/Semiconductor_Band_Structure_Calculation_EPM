#include <iostream>
#include <cmath>
#include <fstream>
#include"lapacke.h"

/*Compile the code using : g++ band_stru.cpp -o band_ex_cpp -llapack
  OR
  g++ band_stru.cpp -o band_ex_cpp -L/usr/local/lib64 -llapack -lblas -lgfortran -lm */

using namespace std;

typedef struct ig_vec{
    int igx_vec, igy_vec, igz_vec;
} Ig_vec;

typedef struct pseudo {
    double Vs3, Vs4, Vs8, Vs11;
    double Va3, Va4, Va8, Va11;
} Ps;

typedef struct dk {
    double dkx, dky, dkz;
} Dk;

typedef struct crys {
    double ao, t;
} Crys;

typedef struct comple {
    double re;
    double im;
} Comple;

void eigenpro(Comple **A, double *ek);
void eigenpro(Comple **A, double *ek) {

    int size = 137;
    Comple b[size], WORK[2*size]; //RWORK[2*size];

    Comple w[size], vl[1][size], vr[1][size];

  double AT[2*size*size], tmp, RWORK[2*size];                 /* for transformed matrix */
  int i, j, ok;
  char jobvl, jobvr;
  int n, lda, ldvl, ldvr, lwork;


  n=size;
  jobvl='N';
  jobvr='N';
  lda=size;
  ldvl=1;
  ldvr=1;
  lwork=2*size;


  for (i=0; i<size; i++)          /* to call a Fortran routine from C we */
  {                               /* have to transform the matrix */
    for(j=0; j<size; j++)
    {
       AT[2*(j+size*i)]=A[j][i].re;
       AT[2*(j+size*i)+1]=A[j][i].im;
    }
  }

  /* find solution using LAPACK routine ZGEEV, all the non-array arguments
  have to be pointers 
  http://www.netlib.org/lapack/explore-html/db/d55/group__complex16_g_eeigen_ga0eb4e3d75621a1ce1685064db1ac58f0.html 
  http://www.netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_gaabef68a9c7b10df7aef8f4fec89fddbe.html
  https://kb.iu.edu/d/aqpn*/
  zgeev_(&jobvl, &jobvr,&n, reinterpret_cast <__complex__ double*> (AT), &lda, reinterpret_cast <__complex__ double*> (w),reinterpret_cast <__complex__ double*> (vl), &ldvl,reinterpret_cast <__complex__ double*> (vr), &ldvr,reinterpret_cast <__complex__ double*> (WORK), &lwork,
  RWORK, &ok);

  if (ok==0) {                             /* output of eigenvalues */
  for (i = 0; i < size; i++) {
      for (j = i + 1; j < size; j++) {
          if (w[i].re > w[j].re) {
             tmp = w[i].re;
             w[i].re = w[j].re;
             w[j].re = tmp;
          }
      }
  }
  for (i = 0; i < 8; i++) ek[i] = w[i].re;
  //for (i = 0; i < 8; i++) printf("%dth = %f\t%f\n", i, w[i].re, w[i].im);
     /*for (i=0; i<size; i++)
     {
        printf("%dth = %f\t%f\n", i, w[i].re, w[i].im);
     }*/
  }
  else printf("An error occurred");

  return;
}

void generate_Hamil(Crys *Inp, Dk *dk1, Ig_vec *Ig_Vec0, Ps *pseu, double *ek, Comple **H_EPM);
void generate_Hamil(Crys *Inp, Dk *dk1, Ig_vec *Ig_Vec0, Ps *pseu, double *ek, Comple **H_EPM) {
    Ig_vec *Ig_Vec, *Ig_Vec1;
    double pi = 4.0*atan(1.0), tau_fac, ct, st, Hreal, Himag;
    double k0 = 2.0*pi/Inp->ao, V_s, V_a;
    double gx_i, gy_i, gz_i, gx_j, gy_j, gz_j, gkx2, gky2, gkz2, dgx, dgy, dgz;
    int N_pw = 137;
    int i, j, id2, idgx, idgy, idgz;

    for (i =  0, Ig_Vec = Ig_Vec0; i < N_pw; i++, Ig_Vec++) {
        gx_i = (double) Ig_Vec->igx_vec;
        gy_i = (double) Ig_Vec->igy_vec;
        gz_i = (double) Ig_Vec->igz_vec;
        for (j = 0, Ig_Vec1 = Ig_Vec0; j < N_pw; j++, Ig_Vec1++) {
            gx_j = (double) Ig_Vec1->igx_vec;
            gy_j = (double) Ig_Vec1->igy_vec;
            gz_j = (double) Ig_Vec1->igz_vec;
            if (i == j) {
                gkx2 = (gx_i+dk1->dkx)*(gx_i+dk1->dkx);
                gky2=(gy_i+dk1->dky)*(gy_i+dk1->dky);
                gkz2=(gz_i+dk1->dkz)*(gz_i+dk1->dkz);
                Hreal= Inp->t*k0*k0*(gkx2+gky2+gkz2);
                H_EPM[i][j].re=Hreal;
                H_EPM[i][j].im=0.0;
            }
            else {
                dgx=gx_i-gx_j;
                dgy=gy_i-gy_j;
                dgz=gz_i-gz_j;
                idgx=round(dgx);
                idgy=round(dgy);
                idgz=round(dgz);
                id2=idgx*idgx+idgy*idgy+idgz*idgz;
                tau_fac=2.00*pi/8.00;
                if (id2==3) {
                    V_s=pseu->Vs3;
                    V_a=pseu->Va3;
                }
                else if (id2==4) {
                    V_s=pseu->Vs4;
                    V_a=pseu->Va4;
                }
                else if (id2==8) {
                    V_s=pseu->Vs8;
                    V_a=pseu->Va8;
                }
                else if (id2==11) {
                    V_s=pseu->Vs11;
                    V_a=pseu->Va11;
                }
                else {
                    V_s=0.00;
                    V_a=0.00;
                }
            ct=cos((dgx+dgy+dgz)*tau_fac);
            st=sin((dgx+dgy+dgz)*tau_fac);
            Hreal= V_s*ct;
            Himag= V_a*st;
            H_EPM[i][j].re= Hreal;
            H_EPM[i][j].im= Himag;
            }
        }
    }
    eigenpro(H_EPM,ek);
    return;
}

int main() {

    Ps *pseu;
    Ig_vec *Ig_Vec, *Ig_Vec0;
    Crys *Inp;
    Dk *dk1;
    double *ek, **ekm;
    double  Ry2eV=13.6056981;
    double hb=1.054e-34;
    double fh=6.625e-34;
    double q=1.602e-19;
    double mo=9.11e-31;
    double t=hb*hb/2.0/mo/q;
    double a0[14], str_fac[14][6];
    double a,X0,Y0,Z0,sumx,dkx,dky,dkz;
    int igx, igx2, jgy, jgy2, kgz, kgz2, igxy2, igxyz2;
    int namelen;
    char *filename;
    string txt = "_cpp.txt";
    const char *mat[14]={"Si", "Ge", "Sn", "GaP", "GaAs", "AlSb", "InP", "GaSb", "InAs", "InSb", "ZnS", "ZnSe", "ZnTe", "CdTe"};

    int i = 0, j = 0, icount = 0, ig2, ia, ic, ik, ik_curr, ijk, lh, kh, ii;
    int Num_semi=14;
    int NN_neighbours=10;
    int N_pw=137;
    int N_bands=12;
    int N_ek=40;
    int ig2_mag[10]={0, 3, 4, 8, 11, 12, 16, 19, 20, 24};   // Magnitude of 10 nearest G plane |(G-G0)2|
    int Lwork=2*N_pw-1;
    Comple **H_EPM;

    ifstream ofp;
    ofstream myfile;

    ofp.open("a0.txt");
    for (i = 0; i < 14; i++) {
        ofp >> a0[i];
    }
    ofp.close();
    ofp.open("str_fac.txt");
    for (i = 0; i < 14; i++){
        for (j = 0; j < 6; j++) {
            ofp >> str_fac[i][j];
        }
    }
    ofp.close();

    pseu = new Ps[1];
    Ig_Vec0 = new Ig_vec[137];
    Ig_Vec = Ig_Vec0;
    i = 0;

    icount = 0;
    for (j = 0; j < NN_neighbours; j++) {  // Total number of G planes with a given magnitude of |G2|
        ig2 = ig2_mag[j];
        a = (double) sqrt(ig2);
        ia =  floor(a);
        ic = 0;
        for (igx = -ia; igx <= ia; igx++) {
            igx2=igx*igx;
            for (jgy = -ia; jgy <= ia; jgy++) {
                jgy2=jgy*jgy;
                igxy2=igx2+jgy2;
                for (kgz = -ia; kgz <= ia; kgz++) {
                    kgz2=kgz*kgz;
                    igxyz2=igxy2+kgz2;
                    if (igxyz2==ig2) {
                        ic++;
                        Ig_Vec->igx_vec=igx;
                        Ig_Vec->igy_vec=jgy;
                        Ig_Vec->igz_vec=kgz;
                        Ig_Vec++;
                    }
                }
            }
        }
        icount=icount+ic;
            //printf("%d %d %d %d\n",ig2,ia,ic,icount);
    }

    dk1 = new Dk[1]; 
    Inp = new Crys[1]; 
    Inp->t = t;
    H_EPM = new Comple*[137];
    for (ii = 0; ii < 137; ii++) H_EPM[ii] = new Comple[137];
    if (H_EPM == NULL) printf("H_EPM : Memory allocation failed!\n");
    ekm = new double*[41]; 
    for (ii = 0; ii < 41; ii++) ekm[ii] = new double[8];
    ek = new double[8];    
    i = 0;
    while (i < 1) {
        namelen = txt.length() + 1;
        txt = mat[i];
        namelen += txt.length();
        filename =  new char[namelen];
        sprintf(filename, "%s%s", mat[i], "_cpp.txt");
	myfile.open(filename);
	myfile << "Band Structure of " << mat[i] << "\n";
	myfile.close();

        pseu->Vs3=str_fac[i][0]*Ry2eV;
        pseu->Vs4=0.0;
        pseu->Vs8=str_fac[i][1]*Ry2eV;
        pseu->Vs11=str_fac[i][2]*Ry2eV;

        pseu->Va3=str_fac[i][3]*Ry2eV;
        pseu->Va4=str_fac[i][4]*Ry2eV;
        pseu->Va8=0.0;
        pseu->Va11=str_fac[i][5]*Ry2eV;

        Inp->ao = a0[i]*1.0e-10;	

        /*Symmetry path L->Gamma*/
        ik_curr=0;
        sumx=0.0;
        X0=0.5;
        Y0=0.5;
        Z0=0.5;
        for (ik = 1; ik <= N_ek + 1; ik++) {
            dk1->dkx=0.500*(1.00-(double) (ik-1)/N_ek);
            dk1->dky=0.500*(1.00-(double) (ik-1)/N_ek);
            dk1->dkz=0.500*(1.00-(double) (ik-1)/N_ek);
            generate_Hamil(Inp, dk1, Ig_Vec0, pseu, ek, H_EPM);
            for (ijk = 0; ijk < 8; ijk++) ekm[ik-1][ijk] = ek[ijk];
        }

        /*Dumping data for Symmetry path L->Gamma*/
        myfile.open(filename, std::ios_base::app);
        for (lh = 0; lh < N_ek + 1; lh++){
            myfile << ik_curr + lh << "\t";
	    for (kh = 0; kh < 8; kh++) myfile << ekm[lh][kh] << "\t";
            myfile << "\n";
        }
        myfile.close(); 

	/*Symmetry path Gamma->X*/
	ik_curr += ik - 1;
	X0=0.0;
	Y0=0.0;
	Z0=0.0;
        for (ik = 1; ik <= N_ek; ik++) {
	    dk1->dkx=(double) ik/N_ek;
	    dk1->dky=0.0;
	    dk1->dkz=0.0;
            generate_Hamil(Inp, dk1, Ig_Vec0, pseu, ek, H_EPM);
            for (ijk = 0; ijk < 8; ijk++) ekm[ik-1][ijk] = ek[ijk];
        }

        /*Dumping data for Symmetry path Gamma->X*/
        myfile.open(filename, std::ios_base::app);
        for (lh = 0; lh < N_ek; lh++){
            myfile << ik_curr + lh << "\t";
	    for (kh = 0; kh < 8; kh++) myfile << ekm[lh][kh] << "\t";
            myfile << "\n";
        }
        myfile.close();


        /*Symmetry path X->U*/
        ik_curr += ik -1;
        X0=1.0;
        Y0=0.0;
        Z0=0.0;
        for (ik = 1; ik <= N_ek; ik++) {
            dk1->dkx=1.0;
            dk1->dky=(double) (ik)/(N_ek)*1.0/4.0;
            dk1->dkz=(double) (ik)/(N_ek)*1.0/4.0;
            generate_Hamil(Inp, dk1, Ig_Vec0, pseu, ek, H_EPM);
            for (ijk = 0; ijk < 8; ijk++) ekm[ik-1][ijk] = ek[ijk];
        }

        /*Dumping data for Symmetry path X->U*/
        myfile.open(filename, std::ios_base::app);
        for (lh = 0; lh < N_ek; lh++){
            myfile << (double) (ik_curr + lh/2.0) << "\t";
	    for (kh = 0; kh < 8; kh++) myfile <<  ekm[lh][kh] << "\t";
            myfile << "\n";
        }
        myfile.close();

        /*Symmetry path K->Gamma*/
        ik_curr += ik/2;
        X0=0.75;
        Y0=0.75;
        Z0=0.0;
        for (ik = 1; ik <= N_ek; ik++) {
            dk1->dkx=0.75*(1.0-(double) ik/N_ek);
            dk1->dky=0.75*(1.0-(double) ik/N_ek);
            dk1->dkz=0.0;
            generate_Hamil(Inp, dk1, Ig_Vec0, pseu, ek, H_EPM);
            for (ijk = 0; ijk < 8; ijk++) ekm[ik-1][ijk] = ek[ijk];
        }

        /*Dumping data for Symmetry path K->Gamma*/
        myfile.open(filename, std::ios_base::app);
        for (lh = 0; lh < N_ek; lh++){
            myfile << ik_curr + lh << "\t";
	    for (kh = 0; kh < 8; kh++) myfile << ekm[lh][kh] << "\t";
            myfile << "\n";
        }
        myfile.close();

        printf("%s\n", filename);

        i++;
    }
    return 0;
}
