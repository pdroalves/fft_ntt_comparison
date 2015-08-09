from cmath import exp, pi
import generate_prime as Prime

invMod = lambda y,q:pow(y,q-2,q)

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

def ntt_aux(x,wN,q):
    N = len(x)
    if N <= 1: return x
    even = ntt_aux(x[0::2],wN,q)
    odd =  ntt_aux(x[1::2],wN,q)
    # import pdb;pdb.set_trace()
    T= [((wN**(k/N))*odd[k])%q for k in xrange(N/2)]
    return [(even[k] + T[k]) for k in xrange(N/2)] + [(even[k] - T[k]) for k in xrange(N/2)]

def ntt(x,q):
    r,k = generate_primitive_root(len(x),q)
    wN = r**k
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
    r,k = generate_primitive_root(len(x),q)
    wN = r**k
    return [y*invMod(len(x),q)%q for y in intt_aux(x,wN,q)]
