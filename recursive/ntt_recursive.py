from cmath import exp, pi
import sys
sys.path.append("../")
import generate_prime as Prime

invMod = lambda y,q:pow(y,q-2,q)
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

def ntt_aux(x,wN,q):
    N = len(x)
    if N <= 1: return x
    even = ntt_aux(x[0::2],wN,q)
    odd =  ntt_aux(x[1::2],wN,q)
    # import pdb;pdb.set_trace()
    T= [((wN**(k/N))*odd[k])%q for k in xrange(N/2)]
    return [(even[k] + T[k]) for k in xrange(N/2)] + [(even[k] - T[k]) for k in xrange(N/2)]

def ntt(x,q):
    N = len(x)
    assert (q-1) % N == 0 # If this is not true, we won't find a proper k

    k = (q-1) / N
    r = generate_primitive_root(N,q)
    wN = r**k
    assert pow(wN, len(x), q) == 1

    return [y % q for y in ntt_aux(x,wN,q)]

def intt_aux(x,wN,q):
    N = len(x)
    if N <= 1: return x
    even = intt_aux(x[0::2],wN,q)
    odd =  intt_aux(x[1::2],wN,q)
    # T= [invMod(wN,q)*odd[k] % q for k in xrange(N/2)]
    T= [(invMod(wN**(k/N) %q,q)*odd[k])% q for k in xrange(N/2)]
    return [(even[k] + T[k]) for k in xrange(N/2)] + [(even[k] - T[k]) for k in xrange(N/2)]

def intt(x,q):
    N = len(x)
    assert (q-1) % N == 0 # If this is not true, we won't find a proper k

    k = (q-1) / N
    r = generate_primitive_root(len(x),q)
    wN = r**k
    return [y*invMod(len(x),q)%q for y in intt_aux(x,wN,q)]
