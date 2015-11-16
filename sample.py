#!/usr/bin/python

from fft_recursive import fft,ifft
from ntt_recursive import ntt,intt
from math import log
import generate_prime as Prime

print "Comparison between recursive NTT and FFT algorithm "

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

equal = lambda a,b: isclose(a[0],b[0],abs_tol=0.0001) and isclose(a[1],b[1],abs_tol=0.0001)

for i in range(2,10):
    N = 2**i
    A = range(N)
    q = Prime.generate_large_prime(log(N,2)+1)
    print "Testing for N == %d) FFT: %s - NTT: %s" % (N,equal(A,ifft(fft(A))),A == intt(ntt(A,q),q))
    # if (A == ifft(fft(A))) is False:
    #     print "Incorrect FFT: %s == %s " %(A,ifft(fft(A)))
    #     print "Incorrect NTT: %s == %s " %(A,intt(ntt(A,q),q))
