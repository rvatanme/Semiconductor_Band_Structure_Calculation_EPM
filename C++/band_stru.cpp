#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

int main() {

    double  Ry2eV=13.6056981;
    double hb=1.054D-34;
    double fh=6.625D-34;
    double q=1.602d-19;
    double mo=9.11D-31;
    double t=hb*hb/2.0/mo/q;
    double a0[14], str_fac[14][6], b0;
    int namelen;
    char *filename;
    string txt = ".txt";
    const char *mat[14]={"Si", "Ge", "Sn", "GaP", "GaAs", "AlSb", "InP", "GaSb", "InAs", "InSb", "ZnS", "ZnSe", "ZnTe", "CdTe"};

    int i = 0, j = 0;
    int Num_semi=14;
    int NN_neighbours=10;
    int N_pw=137;
    int N_bands=12;
    int N_ek=40;
    int ig2_mag[10]={0, 3, 4, 8, 11, 12, 16, 19, 20, 24};
    int Lwork=2*N_pw-1;

    ifstream myfile;

    myfile.open("a0.txt");
    for (i = 0; i < 14; i++) {
        myfile >> a0[i];
    }
    myfile.close();
    myfile.open("str_fac.txt");
    for (i = 0; i < 14; i++){
        for (j = 0; j < 6; j++) {
            myfile >> str_fac[i][j];
        }
    }
    myfile.close();
    i = 0;
    while (i < 14) {
        namelen = txt.length() + 1;
        txt = mat[i];
        namelen += txt.length();
        filename =  new char[namelen];
        sprintf(filename, "%s%s", mat[i], ".txt");
        printf("%s\n", filename);
        i++;
    }

    return 0;
}
