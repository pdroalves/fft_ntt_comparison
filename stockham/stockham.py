#!/usr/bin/python
#coding: utf-8

import numpy as np
from math import pi,cos,sin

is_power2 = lambda N: (N & (N - 1)) == 0 
MAX_ERROR = 0.0001
#####################
# Complex
def complex_add(a,b):
	return [a[0]+b[0],
			a[1]+b[1]]
def complex_sub(a,b):
	return [a[0]-b[0],
			a[1]-b[1]]
def complex_mul(a,b):
	return [a[0]*b[0]-a[1]*b[1],
			a[0]*b[1] + a[1]*b[0]]

######################
# Stockham
def CPU_FFT(data):
	N = len(data)
	assert N > 0 and is_power2(N)
	R = 2
	Ns = 1
	a = [[(float)(x),0.0] for x in data] + [[0.0,0.0]]*(N-len(data))
	b = [[0.0,0.0]]*N
	while Ns < N:
		for j in range(N/R):
			FftIteration(j,N,R,Ns,a,b)
		a,b = b,a
		Ns = Ns*R
	return a

def CPU_IFFT(data):
	N = len(data)
	assert N > 0 and is_power2(N)
	R = 2
	Ns = 1
	a = list(data)
	b = [[0.0,0.0]]*N
	while Ns < N:
		for j in range(N/R):
			IfftIteration(j,N,R,Ns,a,b)
		a,b = b,a
		Ns = Ns*R
	return a


def FftIteration(j,N,R,Ns,data0,data1):
	v = [[0,0]]*R
	idxS = j;
	angle = -2*pi*(j%Ns)/(Ns*R)
	for r in range(R):
		v[r] = complex_mul( data0[idxS+r*N/R],
							[cos(r*angle),sin(r*angle)]
							)
	v = FFT(R,v)
	idxD = expand(j,Ns,R)
	for r in range(R):
		data1[idxD+r*Ns] = v[r]
	return

def IfftIteration(j,N,R,Ns,data0,data1):
	v = [[0,0]]*R
	idxS = j;
	angle = 2*pi*(j%Ns)/(Ns*R)
	for r in range(R):
		v[r] = complex_mul( data0[idxS+r*N/R],
							[cos(r*angle),sin(r*angle)]
							)
	v = FFT(R,v)
	idxD = expand(j,Ns,R)
	for r in range(R):
		data1[idxD+r*Ns] = v[r]
	return

def FFT(R,v):
	if R == 2:
		v0 = v[0]
		return [complex_add(v0,v[1]),
				complex_sub(v0,v[1])]
	else:
		raise "Error!"

def expand(idxL,N1,N2):
	return (idxL/N1)*N1*N2 + (idxL%N1)

################################
# Test
from random import randint
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def test_fft_numpy():
	print "Test fft_numpy"
	powers2 = [pow(2,i) for i in range(1,19)]
	equal = lambda a,b: isclose(a[0],b.real,abs_tol=0.0001) and isclose(a[1],b.imag,abs_tol=0.0001)

	for N in powers2:
		print N
		a = [randint(0,2**16) for _ in xrange(N)]

		my_fft_a = CPU_FFT(a)
		np_fft_a = list(np.fft.fft(a))
		assert len(my_fft_a) == len(np_fft_a)
		l = len(my_fft_a)
		for i in range(l):
			assert( equal(
							my_fft_a[i],
							np_fft_a[i]
						)
					)
	print "We are good"

def test_fft_ifft():
	print "Test fft_ifft"
	powers2 = [pow(2,i) for i in range(1,19)]
	equal = lambda a,b: isclose(a[0],b.real,abs_tol=MAX_ERROR) and isclose(a[1],b.imag,abs_tol=0.0001)

	for N in powers2:
		print N
		a = [randint(0,2**16) for _ in xrange(N)]

		my_fft_a = CPU_FFT(a)
		fft_ifft_a = CPU_IFFT(my_fft_a)

		assert len(fft_ifft_a) == len(a)
		for i,v in enumerate(a):
			assert( isclose(v,fft_ifft_a[i][0]/N,abs_tol=MAX_ERROR))
	print "We are good"	

def test_fft_mul():
	print "test_fft_mul"
	powers2 = [pow(2,i) for i in range(1,19)]
	equal = lambda a,b: isclose(a[0],b.real,abs_tol=MAX_ERROR) and isclose(a[1],b.imag,abs_tol=0.0001)

	for N in powers2:
		print N

		a = [randint(0,2**16) for _ in xrange(N)]+[0]*N
		b = [randint(0,2**16) for _ in xrange(N)]+[0]*N

		polynomial_a = Polynomial()	
		polynomial_a.coef = a
		polynomial_b = Polynomial()
		polynomial_b.coef = b

		polynomial_c = polynomial_a*polynomial_b

		my_fft_a = CPU_FFT(a)
		my_fft_b = CPU_FFT(b)
		my_fft_c = [complex_mul(x[0],x[1]) for x in zip(my_fft_a,my_fft_b)]

		fft_ifft_c = CPU_IFFT(my_fft_c)

		# assert len(fft_ifft_c) == len(polynomial_c.coef)
		for i,v in enumerate(polynomial_c.coef):
			v2 = fft_ifft_c[i][0]/(2*N)
			# print v,v2
			assert isclose(v,v2)
	print "We are good"	

def main():
	test_fft_ifft()
	test_fft_numpy()
	test_fft_mul()

if __name__ == "__main__":
	main()
