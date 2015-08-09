from cmath import exp, pi

def fft(x):
    N = len(x)
    if N <= 1: return x
    even = fft(x[0::2])
    odd =  fft(x[1::2])
    T= [exp(2j*pi*k/N)*odd[k] for k in xrange(N/2)]
    return [even[k] + T[k] for k in xrange(N/2)] + [even[k] - T[k] for k in xrange(N/2)]

def ifft_aux(x):
    N = len(x)
    if N <= 1: return x
    even = ifft_aux(x[0::2])
    odd =  ifft_aux(x[1::2])
    T= [exp(-2j*pi*k/N)*odd[k] for k in xrange(N/2)]
    return [even[k] + T[k] for k in xrange(N/2)] + [even[k] - T[k] for k in xrange(N/2)]

def ifft(x):
    return [int(y.real)/len(x) for y in ifft_aux(x)]
