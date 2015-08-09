# FFT and NTT comparison
University of Campinas
Institute of Computing
Laboratory of Security and Applied Cryptography (LASCA) 

Author: Pedro Alves

The intent of this repository is to study and compare a very simple implementation of a **FFT** (**Fast Fourier Transform**) algorithm and a very simple implementation of a **NTT** (**Number Theoretic Transform**) algorithm. For this, we based our work in a [existent](http://rosettacode.org/wiki/Fast_Fourier_transform#Python:_Recursive) FFT implementation in python of a recursive version of the algorithm.

As can be seen in [this reference](http://www.apfloat.org/ntt.html), the main (or the only) difference of FFT and NTT is that, the first operates over reals, while the second operates over a finite field. This way, FFT uses wN = exp(2j*pi/N) as primitive root, while NTT is able to use a generator of the finite field.

We apply the same code for FFT and NTT, but in the last one we replace all operators by modular operators, and the primitive root by a runtime computed value. For now, the algorithm used to generate this primitive root is a little dump.
