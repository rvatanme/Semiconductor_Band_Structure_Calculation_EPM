#include <stdio.h>
#define MAX_BUF 200

int main() {

    FILE *ofp = NULL, *ofp1 = NULL;
    char path[MAX_BUF];
    double Ry2eV = 13.6056981e0, hb = 1.054e-34, h = 6.625e-34;
    double q=1.602e-19, mo = 9.11e-31, t = hb*hb/2.0e0/mo/q;
    int Num_semi=14, NN_neighbor = 10, N_pw = 137, N_bands = 12, N_ek = 40, Lwork = 2*N_pw - 1, ichoose;
    int ig2_mag[10] = {0, 3, 4, 8, 11, 12, 16, 19, 20, 24};
    float **str_fac = NULL, num;
    float **dptr, *sptr;
    int k = 0, l = 0, m = 0, n = 0;

    getcwd(path, MAX_BUF);
    printf("Initial directory was \"%s\"\n", path);
    chdir("/home/reza/courses/EEE 533 Semiconductor Process & Device Simulation");
    getcwd(path, MAX_BUF);
    printf("The current directory is \"%s\"\n", path);
    ofp = fopen("str_fac.txt","r");
    ofp1 = fopen("a0.txt","r");
    for (m = 0; m < 14 ; m++) {
        for (n = 0; n < 6; n++) {
            fscanf(ofp, "%f", &num);
            printf("%f  ", num);
            if ((n + 1) % 3 == 0) printf("\n");
        }
    }

    str_fac = (float **)calloc(14, sizeof(float *));
    for (dptr = str_fac; dptr < str_fac + 14; dptr++) *dptr = (float *)calloc(8, sizeof(float));
    for (dptr = str_fac, m = 0; dptr < str_fac + 14; dptr++, m++)
        for (sptr = *dptr, n = 0; sptr < *dptr + 8; sptr++, n++) {
            *sptr = (k++)*0.1;  //*(*(arr+i)+j)
            //printf("%f ", str_fac[m][n]);
        }
    printf("This program calculates the bandstructure of semiconductors based on the Empirical");
    printf("Pseudopotential Method (EPM)\n");
    printf("The following semiconductors are included\n");
    printf("1 -> Si 2 -> Ge 3 -> Sn 4 -> GaP 5 -> GaAs 6 -> AlSb 7 -> InP");
    printf("8 -> GaSb 9 -> InAs 10 -> InSb 11 -> ZnS 12 -> ZnSe 13 -> ZnTe 14 -> CdTe\n\n");
    printf("Enter the semiconductor of your choice: ");
    scanf("%d", &ichoose);
    free(str_fac);
    fclose(ofp);
    fclose(ofp1);

    return 0;
}
