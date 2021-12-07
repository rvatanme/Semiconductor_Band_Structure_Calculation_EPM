#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>

typedef struct ig_vec{
    double igx_vec, igy_vec, igz_vec;
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

void generate_Hamil(Crys *Inp, Dk *dk1, Ig_vec *Ig_Vec0, Ps *pseu, double **ek);
void generate_Hamil(Crys *Inp, Dk *dk1, Ig_vec *Ig_Vec0, Ps *pseu, double **ek) {
    Ig_vec *Ig_Vec, *Ig_Vec1;
    double pi = 4.0*atan(1.0), tau_fac, ct, st, Hreal, Himag;
    double k0 = 2.0*pi/Inp->ao, V_s, V_a;
    double gx_i, gy_i, gz_i, gx_j, gy_j, gz_j, gkx2, gky2, gkz2, dgx, dgy, dgz;
    int N_pw = 137;
    int i, j, id2, idgx, idgy, idgz;
    double complex H_EPM[137][137];

    for (i =  1, Ig_Vec = Ig_Vec0; i <= N_pw; i++, Ig_Vec++) {
        gx_i = (double) Ig_Vec->igx_vec;
        gy_i = (double) Ig_Vec->igy_vec;
        gz_i = (double) Ig_Vec->igz_vec;
        for (j = 1, Ig_Vec1 = Ig_Vec0; j <= N_pw; j++, Ig_Vec1++) {
            gx_j = (double) Ig_Vec1->igx_vec;
            gy_j = (double) Ig_Vec1->igy_vec;
            gz_j = (double) Ig_Vec1->igz_vec;
            if (i == j) {
                gkx2 = (gx_i+dk1->dkx)*(gx_i+dk1->dkx);
                gky2=(gy_i+dk1->dky)*(gy_i+dk1->dky);
                gkz2=(gz_i+dk1->dkz)*(gz_i+dk1->dkz);
                Hreal= Inp->t*k0*k0*(gkx2+gky2+gkz2);
                H_EPM[i][j]=Hreal + 0.0*I;
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
            H_EPM[i][j]= Hreal + Himag*I;
            }
        }
    }
    //E=sort(eig(H_EPM));
    //ek(ik_curr+ik,1:8)=E(1:8)';
    return ek;
}

int main() {

    Ps *pseu;
    Ig_vec *Ig_Vec, *Ig_Vec0;
    Crys *Inp;
    Dk *dk1;
    double **ek;
    double Vs3, Vs8, Vs11, Vs12, Va3, Va4, Va11, Va12, ao, a;
    double  Ry2eV=13.6056981;
    double hb=1.054e-34;
    double fh=6.625e-34;
    double q=1.602e-19;
    double mo=9.11e-31;
    double t=hb*hb/2.0/mo/q;
    double a0[14], str_fac[14][6];
    double X0,Y0,Z0,sumx,dkx,dky,dkz, k_vec[200];
    int igx, igx2, jgy, jgy2, kgz, kgz2, igxy2, igxyz2;
    int igx_vec[200], igy_vec[200], igz_vec[200];
    int namelen;
    char *filename;
    const char *mat[14]={"Si", "Ge", "Sn", "GaP", "GaAs", "AlSb", "InP", "GaSb", "InAs", "InSb", "ZnS", "ZnSe", "ZnTe", "CdTe"};

    int i = 0, j = 0, icount = 0, ig2, ia, ic, ik, ik_curr;
    int Num_semi=14;
    int NN_neighbours=10;
    int N_pw=137;
    int N_bands=12;
    int N_ek=40;
    int ig2_mag[10]={0, 3, 4, 8, 11, 12, 16, 19, 20, 24};   // Magnitude of 10 nearest G plane |G2| to (0,0,0)
    int Lwork=2*N_pw-1;

    FILE *ofp = NULL;

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
                        Ig_Vec++;
                        ic++;
                        Ig_Vec->igx_vec=igx;
                        Ig_Vec->igy_vec=jgy;
                        Ig_Vec->igz_vec=kgz;
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
    ek = (double **)calloc(40, sizeof(double *));
    for (i = 0; i < 8; i++) *ek = (double *)calloc(8, sizeof(double));
    i = 0;
    while (i < 14) {

        ao=a0[i];
        pseu->Vs3=str_fac[i][0]*Ry2eV;
        pseu->Vs4=0.0;
        pseu->Vs8=str_fac[i][1]*Ry2eV;
        pseu->Vs11=str_fac[i][2]*Ry2eV;

        pseu->Va3=str_fac[i][3]*Ry2eV;
        pseu->Va4=str_fac[i][4]*Ry2eV;
        pseu->Va8=0.0;
        pseu->Va11=str_fac[i][5]*Ry2eV;

        Inp->ao = ao = ao*1.0e-10;

        /*Symmetry path L->Gamma*/
        ik_curr=0;
        sumx=0.0;
        X0=0.5;
        Y0=0.5;
        Z0=0.5;
        for (ik = 1; ik <= N_ek+1; ik++) {
            ik+ik_curr;
            dk1->dkx=0.500*(1.00-(double) (ik-1)/N_ek);
            dk1->dky=0.500*(1.00-(double) (ik-1)/N_ek);
            dk1->dkz=0.500*(1.00-(double) (ik-1)/N_ek);
            //generate_Hamil;
            sumx=sqrt((X0-dk1->dkx)*(X0-dk1->dkx)+(Y0-dk1->dky)*(Y0-dk1->dky)+(Z0-dk1->dkz)*(Z0-dk1->dkz));
            k_vec[ik+ik_curr]=sumx;
        }

        namelen = strlen(mat[i]) + strlen(".txt") + 1;
        filename =  (char *)calloc(namelen, sizeof(char));
        sprintf(filename, "%s%s", mat[i], ".txt");
        //printf("%s\n", filename);
        i++;
    }
    printf("%f\n", 4*atan(1.0));

    return 0;
}
