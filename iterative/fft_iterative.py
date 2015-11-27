#!/usr/bin/python
from cmath import exp, pi

def bitrev_shuffle(x):
    N = len(x)
    j = 0
    for i in xrange(1, N):
        b = N >> 1
        while j >= b:
            j -= b
            b >>= 1
        j += b
        if j > i:
            x[i], x[j] = x[j], x[i]


def fft_in_place(x):
    N = len(x)
    bitrev_shuffle(x)
    trans_size = 2
    while trans_size <= N:
        wb = 1+0j
        wb_step = exp(2j * pi / trans_size)
        for t in xrange(trans_size >> 1):
            for trans in xrange(N / trans_size):
                i = trans * trans_size + t
                j = i + (trans_size >> 1)
                a = x[i]
                b = x[j] * wb
                x[i] = a + b
                x[j] = a - b
            wb *= wb_step
        trans_size <<= 1


def fft(x):
    x = list(x)
    fft_in_place(x)
    return x

def ifft_in_place(x):
    N = len(x)
    bitrev_shuffle(x)
    trans_size = 2
    while trans_size <= N:
        wb = 1+0j
        wb_step = exp(-2j * pi / trans_size)
        for t in xrange(trans_size >> 1):
            for trans in xrange(N / trans_size):
                i = trans * trans_size + t
                j = i + (trans_size >> 1)
                a = x[i]
                b = x[j] * wb
                x[i] = a + b
                x[j] = a - b
            wb *= wb_step
        trans_size <<= 1
def ifft(x):
    N = len(x)
    x = list(x)
    ifft_in_place(x)
    return [X.real/N for X in x]
