#!/usr/bin/python

from fft import fft,ifft

A = range(8)
print "%s == %s " %(A,ifft(fft(A)))
