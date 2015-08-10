#!/usr/bin/python
from math import log,floor
import generate_prime as Prime

invMod = lambda y,q:pow(y,q-2,q)
#################################################################

def generate_k(n,not_approved=[]):
    k = 1
    while (not Prime.is_prime(k*n+1)) or (k in not_approved):
        k+=1
    return k

def generate_primitive_root(n,q):
    not_approved = []
    r = None

    while r is None:
        k = generate_k(n,not_approved)

        for i in range(q):
            if Prime.is_prime(i):
                s = set()
                for j in range(q):
                    s.add(i**j % q)
                if s == set(range(q))-{0}:
                    r = i
        if r is None:
            # print "%d not approved" % k
            not_approved.append(k)
    return r,k
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

    r,k = generate_primitive_root(len(x),q)
    wN = r**k
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

    r,k = generate_primitive_root(len(x),q)
    wN = r**k
    intt_in_place(x,wN,q)
    return  [y*invMod(len(x),q)%q for y in x]
