#!/usr/bin/python
from math import log,floor
import sys
sys.path.append("../")
import generate_prime as Prime

invMod = lambda y,q:pow(y,q-2,q)
#################################################################

def generate_primitive_root(n,q):
    not_approved = []

    for i in xrange(q):
        if Prime.is_prime(i):
            s = set()
            for j in range(q):
                s.add(i**j % q)
            if s == set(range(q))-{0}:
                return i
    return None
#################################################################

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

#################################################################
def ntt_in_place(x,wN,q):
    N = len(x)
    bitrev_shuffle(x)
    trans_size = 2
    for trans_size in [2**i for i in range(1,int(log(N,2))+1)]:
        wb = 1
        # wb_step = exp(2j * pi / trans_size)
        wb_step = wN**(N/trans_size)

        for t in xrange(trans_size >> 1):
            for trans in xrange(N / trans_size):
                i = trans * trans_size + t
                j = i + (trans_size >> 1)
                a = x[i] % q
                b = x[j] * wb
                x[i] = a + b % q
                x[j] = a - b % q
            wb *= wb_step


def ntt(x,q):
    x = list(x)
    N = len(x)
    assert (q-1) % N == 0 # If this is not true, we won't find a proper k

    k = (q-1) / N
    r = generate_primitive_root(N,q)
    wN = r**k
    assert pow(wN, len(x), q) == 1

    ntt_in_place(x,wN,q)
    return  [y % q for y in x]

def intt_in_place(x,wN,q):
    N = len(x)
    bitrev_shuffle(x)
    trans_size = 2
    for trans_size in [2**i for i in range(1,int(log(N,2))+1)]:
        wb = 1
        # wb_step = exp(2j * pi / trans_size)
        wb_step = invMod(wN**(N/trans_size) % q,q)

        for t in xrange(trans_size >> 1):
            for trans in xrange(N / trans_size):
                i = trans * trans_size + t
                j = i + (trans_size >> 1)
                a = x[i]
                b = x[j] * wb
                x[i] = a + b % q
                x[j] = a - b % q
            wb *= wb_step

def intt(x,q):
    N = len(x)
    x = list(x)

    k = (q-1) / N
    r = generate_primitive_root(N,q)
    wN = r**k
    assert pow(wN, len(x), q) == 1
    intt_in_place(x,wN,q)
    return  [y*invMod(len(x),q)%q for y in x]
