#!/usr/bin/python

from ntt_iterative import ntt,intt
from math import log
import sys
sys.path.append("../")
import generate_prime as Prime

print "Comparison between recursive NTT and FFT algorithm "

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def get_a_proper_prime(N):
	# We need a special prime p of the form p = k*N + 1, for some integer k
	k = 1
	while not Prime.is_prime(k * N + 1):
		k = k + 1
	return k * N + 1

equal = lambda a,b: isclose(a[0],b[0],abs_tol=0.0001) and isclose(a[1],b[1],abs_tol=0.0001)

# Test the transform
for i in range(2,9):
	N = 2**i
	A = range(N)
	q = get_a_proper_prime(N)

	print "Testing for N == %d) NTT: %s" % (
		N,
		A == intt(ntt(A,q),q))

try:
	from numpy import convolve
	from numpy import poly1d
except Exception as e:
	print "An error occured."
	print "We need numpy to validate the polynomial multiplication. Please check if it is installed."
	print e

# Test the polynomial multiplication:
for i in range(2,9):
	N = 2**i
	# q = 2**64-2**32+1
	q = get_a_proper_prime(N)

	A = range(N/2) + [0]*(N/2)
	B = [1]*(N/2) + [0]*(N/2)

	A_ntt = ntt(A, q)
	B_ntt = ntt(B, q)
	C_ntt = [x*y % q for x,y in zip(A_ntt, B_ntt)]
	# print poly1d(intt(C_ntt, q)[::-1])
	# print poly1d(convolve(A, B)[::-1] % q) 
	# print ""

	# We use convolve() to execute the polynomial multiplication through NumPy
	# and poly1d() to compare both objects
	print "Testing multiplication for N == %d) %s " % (
		N,
		poly1d(intt(C_ntt, q)[::-1]) == poly1d(convolve(A, B)[::-1] % q))
