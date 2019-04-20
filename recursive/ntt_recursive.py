from cmath import exp, pi
import sys
sys.path.append("../")
import generate_prime as Prime

invMod = lambda y,q:pow(y,q-2,q)
def generate_primitive_root(n,q):
    if q == 2**64-2**32+1:
        return 7
    not_approved = []

    for i in xrange(q):
        if Prime.is_prime(i):
            s = set()
            for j in range(q):
                s.add(i**j % q)
            if s == set(range(q))-{0}:
                return i
    return None

def ntt(x, q, flag_inverse = False):
    N = len(x)
    assert (q-1) % N == 0 # If this is not true, we won't find a proper k

    # Compute wN at each step
    k = (q-1) / N
    r = generate_primitive_root(N,q)
    wN = pow(r, k, q)
    assert pow(wN, N, q) == 1

    if flag_inverse: # INTT
        wN = invMod(wN, q)

    # End recursion
    if N <= 1: return x

    # Start recursion
    even = ntt(x[0::2], q, flag_inverse)
    odd =  ntt(x[1::2], q, flag_inverse)

    # Apply the transform at each step
    T= [pow(wN, k, q) * odd[k] % q for k in xrange(N/2)]
    return [(even[k] + T[k]) % q for k in xrange(N/2)] + [(even[k] - T[k]) % q for k in xrange(N/2)]

def intt(x,q):
    N = len(x)
    assert (q-1) % N == 0 # If this is not true, we won't find a proper k

    return [y*invMod(N,q)%q for y in ntt(x, q, -1)]

def mul(a, b, q):
    return [x * y % q for x,y in zip(a,b)]