#!/usr/bin/python
import sys
sys.path.append("../recursive")
import ntt_recursive as nttR
import ntt_iterative as nttI
from math import log
sys.path.append("../"
import generate_prime as Prime

print "Comparison between recursive and iterative NTT algorithm"

for i in range(2,10):
    N = 2**i
    A = range(N)
    q = Prime.generate_large_prime(log(N,2)+1)
    print "Testing for N == %d) FFT: %s - NTT: %s" % (N,A == nttI.intt(nttI.ntt(A,q),q),A == nttR.intt(nttR.ntt(A,q),q))
    # if (A == ifft(fft(A))) is False:
    #     print "Incorrect FFT: %s == %s " %(A,ifft(fft(A)))
    #     print "Incorrect NTT: %s == %s " %(A,intt(ntt(A,q),q))
