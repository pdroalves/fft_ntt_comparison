#!/usr/bin/python

from fft import fft,ifft
from ntt import ntt,intt
from math import log
import generate_prime as Prime


for i in range(2,10):
    N = 2**i
    A = range(N)
    q = Prime.generate_large_prime(log(N,2)+1)
    print "Testing for N == %d) FFT: %s - NTT: %s" % (N,A == ifft(fft(A)),A == intt(ntt(A,q),q))
    # if (A == ifft(fft(A))) is False:
    #     print "Incorrect FFT: %s == %s " %(A,ifft(fft(A)))
    #     print "Incorrect NTT: %s == %s " %(A,intt(ntt(A,q),q))
