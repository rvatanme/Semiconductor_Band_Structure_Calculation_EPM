# Electronic Band Structure Calculation
Electronic band structure of a solid provides a set of possible energies that an electron might have with a wavevector of k in that solid. The variety of methods including Emprical Psuedopotential Methods (EPM), k.p perturbation theory, tight binding model, density functional theory, Green's function, etc have been proposed for band structure calculation. Among the many methods available, the EPM method, due to its short calculation runtime, has provided invaluable and practical tool in the design and modelling of modern electronic devices and semiconductor nanostrucre. Here three scripts written in C,C++ and python based on EPM model are provided that calculte the electronic band structure of a solid in diamond or zincblend crystal structure. 

## EPM Method
Consider a collection of noninteracting electrons in a box with L<sub>x</sub>,L<sub>y</sub>,L<sub>z</sub>, then the solution to the time independent Schroedinger’s equation in a priodic framework (ψ(x+L<sub>x</sub>, y+L<sub>y</sub>, z+L<sub>z</sub>) = ψ(x,y,z)) is:

![](https://latex.codecogs.com/svg.latex?%5Cdpi%7B300%7D%20%5CLARGE%20-%5Cfrac%7B%5Chbar%5E2%7D%7B2m%7D%5Ctriangledown%5E2%5Cpsi%20_n%3D%20%5Cepsilon%20_n%5Cpsi%20_n%20%5C%3B%5C%3B%20%5Crightarrow%20%5C%3B%5C%3B%20%5Cpsi%20_n%20%3D%20%5Cfrac%7B1%7D%7B%5Csqrt%7BV%7D%7De%5E%7Bik_xx%7De%5E%7Bik_yy%7De%5E%7Bik_zz%7D)

where V is the volume of the box and n is any set of whole numbers (n<sub>x</sub>, n<sub>y</sub>, n<sub>z</sub>). The exp(ik<sub>x</sub>x) describes an electron traveling along x axis with the average momentum of ћk<sub>x</sub> where k<sub>x</sub> is a real number 2n<sub>x</sub>π/L<sub>x</sub>. Defining r=xi+yj+zk and k=k<sub>x</sub>i+k<sub>y</sub>j+k<sub>z</sub>k, then ψ<sub>k</sub>(r) describes a free electron moving along r direction:

![](https://latex.codecogs.com/svg.latex?%5Cdpi%7B300%7D%20%5CLARGE%20%5Cpsi%20_%7B%5Cvec%7Bk%7D%7D%28%5Cvec%7Br%7D%29%20%3D%20%5Cfrac%7B1%7D%7B%5Csqrt%7BV%7D%7De%5E%7Bi%5Cvec%7Bk%7D.%5Cvec%7Br%7D%7D%20%5C%3B%5C%3B%5Crightarrow%20%5C%3B%5C%3BE_k%3D%5Cfrac%7B%5Chbar%20%5E2k%5E2%7D%7B2m%7D)

Noting the symmetry of ψ<sub>k</sub>(r), the probability density of an electron with any k wavevector is uniform in the space. Following the same procedure, the solution to the time dependent Schroedinger’s equation of free electron is:

![](https://latex.codecogs.com/svg.latex?%5Cdpi%7B300%7D%20%5CLARGE%20-%5Cfrac%7B%5Chbar%5E2%7D%7B2m%7D%5Ctriangledown%5E2%5Cpsi%20%28%5Cvec%7Br%7D%2Ct%29%3D%20i%5Chbar%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20t%7D%5Cpsi%20%28%5Cvec%7Br%7D%2Ct%29%20%5C%5C%20%5C%5C%20%5Cpsi%20_k%28%5Cvec%7Br%7D%2Ct%29%3DAe%5E%7Bi%28%5Cvec%7Bk%7D.%5Cvec%7Br%7D-%5Comega%20t%29%7D%20%3D%20Ae%5E%7Bi%28%5Cvec%7Bp%7D.%5Cvec%7Br%7D-E%20t%29/%5Chbar%7D%20%5C%5C%20%5C%5C%20%5Cvec%7Bp%7D%3D%5Chbar%20%5Cvec%7Bk%7D%20%5C%3B%5C%3B%5C%3B%5C%3B%20E%20%3D%20%5Chbar%20%5Comega)

Now let's imagine an array of atoms or cores (Coloumb potetial centers) are added to the mentioned box with the periodicity of the Bravais lattice. A Bravais lattice is an infinite array of discrete points generated by a set of discrete translation operations described in 3D space by **R** = n<sub>1</sub>**a<sub>1</sub>** + n<sub>2</sub>**a<sub>2</sub>** + n<sub>2</sub>**a<sub>2</sub>**, where n<sub>1</sub>, n<sub>2</sub>, and n<sub>3</sub> are integer numbers and **a<sub>1</sub>**, **a<sub>2</sub>**, **a<sub>3</sub>** are basis vectors.

The electrons in a crystal (basis + lattice) exprience a potential with the same priodicity as the crystal priodicity which means that U(**r+R**)=U(**r**). According to Bloch's theorem, the one electron Schroedinger’s equation for a periodic potential has the form of:

![](https://latex.codecogs.com/svg.latex?%5Cdpi%7B300%7D%20%5CLARGE%20%5Cpsi%20_%7B%5Cvec%7Bk%7D%7D%28%5Cvec%7Br%7D%29%3De%5E%7Bi%5Cvec%7Bk%7D.%5Cvec%7Br%7D%7Du_%7B%5Cvec%7Bk%7D%7D%28%5Cvec%7Br%7D%29%2C%20%5C%3B%5C%3B%5C%3B%5C%3B%20u_%7B%5Cvec%7Bk%7D%7D%28%5Cvec%7Br%7D&plus;%5Cvec%7BR%7D%29%3Du_%7B%5Cvec%7Bk%7D%7D%28%5Cvec%7Br%7D%29)

Since u<sub>k</sub>(r) and U(r) function have the priodicity of the crystal, they can be both expanded as a Fourier sum of M crystal plane wave: 




