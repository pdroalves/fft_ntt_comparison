#!/usr/bin/python
#coding: utf-8

import numpy as np
from math import pi,cos,sin
from polynomial import Polynomial

P = 2**64-2**32+1
w = []
wInv = []
is_power2 = lambda N: (N & (N - 1)) == 0 
MAX_ERROR = 0.0001

######################
# Stockham
def CPU_NTT(data):
	N = len(data)
	assert N > 0 and is_power2(N)
	R = 2
	Ns = 1
	a = list(data)+[0]*(N-len(data))
	b = [0]*N
	while Ns < N:
		for j in range(N/R):
			NTTIteration(j,N,R,Ns,a,b)
		a,b = b,a
		Ns = Ns*R
	return a

def CPU_INTT(data):
	N = len(data)
	assert N > 0 and is_power2(N)
	R = 2
	Ns = 1
	a = list(data)+[0]*(N-len(data))
	b = [0]*N
	while Ns < N:
		for j in range(N/R):
			INTTIteration(j,N,R,Ns,a,b)
		a,b = b,a
		Ns = Ns*R
	return a


def NTTIteration(j,N,R,Ns,data0,data1):
	v = [0]*R
	idxS = j;
	w_index = ((j%Ns)*N)/(Ns*R)
	assert( ((j%Ns)*N)%(Ns*R) == 0)

	for r in range(R):
		#print "Acessando data0[%d] e w[%d]" %(idxS+r*(N/R),r*w_index)
		v[r] = data0[idxS+r*(N/R)]*w[r*w_index] % P

	v = NTT(R,v)
	idxD = expand(j,Ns,R)
	for r in range(R):
		data1[idxD+r*Ns] = v[r]
	return

def INTTIteration(j,N,R,Ns,data0,data1):
	v = [0]*R
	idxS = j;
	w_index = ((j%Ns)*N)/(Ns*R)
	assert( ((j%Ns)*N)%(Ns*R) == 0)

	for r in range(R):
		v[r] = data0[idxS+r*(N/R)]*wInv[r*w_index] % P

	v = NTT(R,v)
	idxD = expand(j,Ns,R)
	for r in range(R):
		data1[idxD+r*Ns] = v[r]
	return

def NTT(R,v):
	if R == 2:
		return [ (v[0] + v[1]) % P,
				 (v[0] - v[1]) % P]
	else:
		raise "Error!"

def expand(idxL,N1,N2):
	return (idxL/N1)*N1*N2 + (idxL%N1)

################################
# Test
from random import randint
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def test_ntt_intt():
	global w,wInv
	print "Test ntt_intt"

	powers2 = [pow(2,i) for i in range(1,19)]

	for N in powers2:
		print N
		
		k = (P-1)/N
		assert (P-1)%N == 0
		wN = pow(7,k,P)
		assert pow(wN,N,P) == 1
		
		w = []
		for j in range(N):
			w.append(pow(wN,j,P))
		wInv = []
		for j in range(N):
			wInv.append(pow(w[j],P-2,P))

		a = [randint(0,2**32) for _ in xrange(N)]

		my_ntt_a = CPU_NTT(a)
		ntt_intt_a = CPU_INTT(my_ntt_a)

		assert len(ntt_intt_a) == len(a)
		for i,v in enumerate(a):
			v2 = ntt_intt_a[i]/N
			# print v,v2
			assert v == v2
	print "We are good"	

def test_ntt_mul():
	global w,wInv
	print "test_ntt_mul"

	powers2 = [pow(2,i) for i in range(1,19)]

	for N in powers2:
		print N
		
		k = (P-1)/(2*N)
		assert (P-1)%(2*N) == 0
		wN = pow(7,k,P)
		assert pow(wN,(2*N),P) == 1
		
		w = []
		for j in range(2*N):
			w.append(pow(wN,j,P))
		wInv = []
		for j in range(2*N):
			wInv.append(pow(w[j],P-2,P))

		a = [randint(0,2**16) for _ in xrange(N)]+[0]*N
		b = [randint(0,2**16) for _ in xrange(N)]+[0]*N

		polynomial_a = Polynomial()	
		polynomial_a.coef = a
		polynomial_b = Polynomial()
		polynomial_b.coef = b

		polynomial_c = polynomial_a*polynomial_b

		my_ntt_a = CPU_NTT(a)
		my_ntt_b = CPU_NTT(b)
		my_ntt_c = [x[0]*x[1] % P for x in zip(my_ntt_a,my_ntt_b)]

		ntt_intt_c = CPU_INTT(my_ntt_c)

		# assert len(ntt_intt_c) == len(polynomial_c.coef)
		for i,v in enumerate(polynomial_c.coef):
			v2 = ntt_intt_c[i]/(2*N)
			assert v == v2
	print "We are good"	

def test_big_integers():
	global w,wInv
	print "Test ntt_intt"

	powers2 = [pow(2,i) for i in range(1,19)]
	N = 1024

	for size in range(10,64):
		k = (P-1)/N
		assert (P-1)%N == 0
		wN = pow(7,k,P)
		assert pow(wN,N,P) == 1
		
		w = []
		for j in range(N):
			w.append(pow(wN,j,P))
		wInv = []
		for j in range(N):
			wInv.append(pow(w[j],P-2,P))

		a = [randint(2**(size-1),2**size) for _ in xrange(N)]

		my_ntt_a = CPU_NTT(a)
		ntt_intt_a = CPU_INTT(my_ntt_a)

		assert len(ntt_intt_a) == len(a)
		if [x/N for x in ntt_intt_a] == a:
			print "Works for %d bits" % size
		else:
			print "Doesn't works for %d bits" % size

	print "We are good"	

def main():
	test_ntt_intt()
	test_ntt_mul()
	test_big_integers()

if __name__ == "__main__":
	main()
