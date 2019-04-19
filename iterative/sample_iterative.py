#!/usr/bin/python
import sys
sys.path.append("../")
import recursive.ntt_recursive as nttR
import ntt_iterative as nttI
from math import log
import generate_prime as Prime

print "Comparison between recursive and iterative NTT algorithm"

def get_a_proper_prime(N):
	# We need a special prime p of the form p = k*N + 1, for some integer k
	k = 1
	while not Prime.is_prime(k * N + 1):
		k = k + 1
	return k * N + 1

for i in range(2,10):
	N = 2**i
	A = range(N)
	q = get_a_proper_prime(N)

	print "Testing for N == %d) FFT: %s - NTT: %s" % (N, A == nttI.intt(nttI.ntt(A,q),q), A == nttR.intt(nttR.ntt(A,q),q))