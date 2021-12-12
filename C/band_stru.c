#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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

typedef struct complex {
    double re;
    double im;
} Complex;

void eigenpro(Complex **A, double *ek);
void eigenpro(Complex **A, double *ek) {

    int size = 137;
    Complex b[size], WORK[2*size], RWORK[2*size];

    Complex w[size], vl[1][size], vr[1][size];

  double AT[2*size*size], tmp;                 /* for transformed matrix */
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
  zgeev_(&jobvl, &jobvr,&n, AT, &lda, w, vl, &ldvl, vr, &ldvr, WORK, &lwork,
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
     /*for (i=0; i<size; i++)
     {
        printf("%dth = %f\t%f\n", i, w[i].re, w[i].im);
     }*/
  for (i = 0; i < 8; i++) ek[i] = w[i].re;
  }
  else printf("An error occurred");

  return;
}

void generate_Hamil(Crys *Inp, Dk *dk1, Ig_vec *Ig_Vec0, Ps *pseu, double *ek, Complex **H_EPM);
void generate_Hamil(Crys *Inp, Dk *dk1, Ig_vec *Ig_Vec0, Ps *pseu, double *ek, Complex **H_EPM) {
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
    FILE *myfile, *ofp;
    const char *mat[14]={"Si", "Ge", "Sn", "GaP", "GaAs", "AlSb", "InP", "GaSb", "InAs", "InSb", "ZnS", "ZnSe", "ZnTe", "CdTe"};

    int i = 0, j = 0, icount = 0, ig2, ia, ic, ik, ik_curr, ijk, lh, kh, ii;
    int Num_semi=14;
    int NN_neighbours=10;
    int N_pw=137;
    int N_bands=12;
    int N_ek=40;
    int ig2_mag[10]={0, 3, 4, 8, 11, 12, 16, 19, 20, 24};   // Magnitude of 10 nearest G plane |(G-G0)2|
    int Lwork=2*N_pw-1;
    Complex **H_EPM;

    ofp = fopen("a0.txt","r");
    for (i = 0; i < 14; i++) {
        fscanf(ofp, "%lf", &a0[i]);
    }
    fclose(ofp);
    ofp = fopen("str_fac.txt","r");
    for (i = 0; i < 14; i++){
        for (j = 0; j < 6; j++) {
            fscanf(ofp,"%lf",&str_fac[i][j]);
        }
    }
    fclose(ofp);

    pseu = (Ps *)calloc(1, sizeof(Ps));
    Ig_Vec0 = (Ig_vec *)calloc(137, sizeof(Ig_vec));
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

    dk1 = (Dk *)calloc(1, sizeof(Dk));
    Inp = (Crys *)calloc(1, sizeof(Crys));
    Inp->t = t;
    H_EPM = (Complex **)calloc(137,sizeof(Complex));
    for (ii = 0; ii < 137; ii++) H_EPM[ii] = (Complex *)calloc(137,sizeof(Complex));
    if (H_EPM == NULL) printf("H_EPM : Memory allocation failed!\n");
    ekm = (double **)calloc(41, sizeof(double *));
    for (ii = 0; ii < 41; ii++) ekm[ii] = (double *)calloc(8, sizeof(double));
    ek = (double *)calloc(8, sizeof(double ));    
    i = 0;
    while (i < 1) {
	
        namelen = strlen(mat[i]) + strlen(".txt") + 1;
        filename =  (char *)calloc(namelen, sizeof(char));
        sprintf(filename, "%s%s", mat[i], ".txt");
	myfile = fopen(filename, "w");
	fprintf(myfile, "Band Structure of %s\n", mat[i]);
	fclose(myfile);

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
        for (ik = 1; ik <= N_ek+1; ik++) {
            dk1->dkx=0.500*(1.00-(double) (ik-1)/N_ek);
            dk1->dky=0.500*(1.00-(double) (ik-1)/N_ek);
            dk1->dkz=0.500*(1.00-(double) (ik-1)/N_ek);
            generate_Hamil(Inp, dk1, Ig_Vec0, pseu, ek, H_EPM);
            for (ijk = 0; ijk < 8; ijk++) ekm[ik-1][ijk] = ek[ijk];
        }

        /*Dumping data for Symmetry path L->Gamma*/
        myfile = fopen(filename, "a");
        for (lh = 0; lh < N_ek + 1; lh++){
            fprintf(myfile, "%d\t", ik_curr + lh);
	    for (kh = 0; kh < 8; kh++) fprintf(myfile, "%f\t", ekm[lh][kh]);
            fprintf(myfile,"\n");
        }
        fclose(myfile);

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
        myfile = fopen(filename, "a");
        for (lh = 0; lh < N_ek; lh++){
            fprintf(myfile, "%d\t", ik_curr + lh);
	    for (kh = 0; kh < 8; kh++) fprintf(myfile, "%f\t", ekm[lh][kh]);
            fprintf(myfile,"\n");
        }
        fclose(myfile);

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
        myfile = fopen(filename, "a");
        for (lh = 0; lh < N_ek; lh++){
            fprintf(myfile, "%0.2f\t", (double) (ik_curr + lh/2.0));
	    for (kh = 0; kh < 8; kh++) fprintf(myfile, "%f\t", ekm[lh][kh]);
            fprintf(myfile,"\n");
        }
        fclose(myfile);

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
        myfile = fopen(filename, "a");
        for (lh = 0; lh < N_ek; lh++){
            fprintf(myfile, "%d\t", ik_curr + lh);
	    for (kh = 0; kh < 8; kh++) fprintf(myfile, "%f\t", ekm[lh][kh]);
            fprintf(myfile,"\n");
        }
        fclose(myfile);

        i++;
    }

    free(pseu);
    free(Ig_Vec0);
    free(dk1);
    free(Inp);
    free(H_EPM);
    free(ek);
    free(ekm);
    return 0;
}
