#!/usr/bin/python
# coding: utf-8
#from numpy.polynomial.polynomial import polyadd
#from numpy.polynomial.polynomial import polysub
#from numpy.polynomial.polynomial import polymul
#from numpy.polynomial.polynomial import polydiv
#from numpy.polynomial.polynomial import polypow
#from numpy.polynomial import Polynomial as NPPolynomial
import itertools
from thereisnoinverseexception import ThereIsNoInverseException
import sys; sys.path.append("../")
import generate_prime as Prime
from math import log
from math import pi as PI
from fractions import gcd
import cmath
import operator
import functools
import json

# see http://jeremykun.com/2014/03/13/programming-with-finite-fields/
LOG2 = lambda x:log(x, 2)

def is_power2(num):

	'states if a number is a power of two'

	return num != 0 and ((num & (num - 1)) == 0)

def is_cyclotomic(p):
    if p.coef[0] == 1 and p.coef[-1] == 1 and set(p.coef[1:-1]) == set([0]):
        return True
    else:
        return False

class Polynomial():
    #  Polynomial class supose this variable is a list where position
    # i stores the coefficient of x^i
    coef = None
    irr_poly = None
    mod = None
    aftermult = False

    def __init__(self,*arg,**kargs):
        if not kargs.has_key("coef"):
            self.coef = list(arg)
        else:
            self.coef = list(kargs["coef"])

        if kargs.has_key("irr_poly"):
            assert isinstance(kargs["irr_poly"],Polynomial) or kargs["irr_poly"] is None
            self.irr_poly = kargs["irr_poly"]

        if kargs.has_key("mod") and type(kargs["mod"]) in (int,long):
            # Coeff mod
            self.mod = kargs["mod"]
            self.coef = [c % self.mod for c in self.coef]

        # Adjust coefficients
        while len(self.coef) > 0 and (self.coef[-1] if self.mod is None else self.coef[-1] % self.mod) == 0:
            self.coef.pop()
        # if len(self.coef) == 0:
        #     self.coef.append(0)

    def encode(self,a):
        assert type(a) == int
        # encodes a integer a to this polynomial
        enc = []
        n = a
        while n != 0:
            enc.append(n % 2)
            n = n / 2
        self.coef = enc

    def decode(self):
        # decodes this polynomial to a integer

        n = 0
        for index,value in enumerate(self.coef):
            n += value*pow(2,index)
        return n


    def __xgcd(self,a, b):
        # gcd(a,b) = a*n + b*m
        # 1 = a*m % mod b, logo m eh a inversa
        #import pdb;pdb.set_trace()
        q = self.mod

        ri = max(a,b)
        ri1 = min(a,b)
        ri2 = None
        qi = None
        qi1 = None

        mi = 1
        ni = 0
        mi1 = 0
        ni1 = 1
        mi2 = None
        mi2 = None
        while ri2 != 0 and ri >= ri1:
            qi = qi1

            qi1,ri2 = ri // ri1
            assert ri1*qi1+ri2 == ri

            mi2 = mi - qi1*mi1
            ni2 = ni - qi1*ni1

            ri,ri1 = ri1,ri2
            mi,mi1 = mi1,mi2
            ni,ni1 = ni1,ni2

        #import pdb;pdb.set_trace()

        return ri,mi,ni # ri eh o gcd, ni eh a inversa

    def inverse(self):
        #import pdb;pdb.set_trace()
        assert self.irr_poly is not None
        assert self <= self.irr_poly

        gcd,m,n = self.__xgcd(self,self.irr_poly)

        if gcd.deg() == 0:
            # There is an inverse
            # if not isinstance(n,Polynomial):
            #     n = Polynomial(coef=[n],irr_poly=self.irr_poly,mod=self.mod).coeffmod()
            # if n*self % self.irr_poly == Polynomial(-1):
            #     return (-n) % self.irr_poly
            # else:
            #     return (n) % self.irr_poly
            return n / gcd
        else:
            raise ThereIsNoInverseException("There is no inverse for %s in mod %s" % (self,self.irr_poly))

    def near(self):
        near = lambda x:int(x)+1 if (x-int(x)) >= 0.5 else int(x)
        return Polynomial(coef=[near(x) for x in self.coef],irr_poly=self.irr_poly)

    def set_dth_cyclotomic(self,d):
       	assert type(d) in (int,long)
	assert is_power2(d)
	self.set_coeff(0, 1)
	self.set_coeff(d/2, 1)
	
        #phi = [Polynomial(coef=[0],irr_poly=self.irr_poly,mod=self.mod)]*(d+1)

        #for i in range(1,d+1):
        #    t = Polynomial(coef=[1],irr_poly=self.irr_poly,mod=self.mod)

        #    for j in range(1,i):
        #        if i % j == 0:
        #            t = t*phi[j]

        #    phi[i] = Polynomial(coef=[-1]+[0]*(i-1)+[1],irr_poly=self.irr_poly,mod=self.mod) / t
        #self.coef = phi[d].coef
        #self.coeffmod()
        return True

    def deg(self):
        return len(self)-1

    def map(self,q):
        # []_q: reduction modulo q into the interval (-q/2,q/2]
        #return Polynomial(coef=[long(x % q) for x in self.coef],irr_poly=self.irr_poly,mod=self.mod).coeffmod()
        return Polynomial(coef=[long(x % q)-q/2 for x in self.coef],irr_poly=self.irr_poly,mod=self.mod).coeffmod()
        #return self % Polynomial(q)

    def is_zero(self):
        return Polynomial(0) == self

    def to_JSON(self):
        return json.dumps(self, default=lambda o: o.__dict__,sort_keys=True)

    def from_JSON(self,j):
        self.coef = j["coef"]
        if j.has_key("irr_poly"):
            self.irr_poly = j["irr_poly"]
        return self

    def copy(self):
        # return a copy of the instance
        copy = Polynomial()
        if self.coef is not None:
            copy.coef=list(self.coef)
        if self.irr_poly is not None:
            copy.irr_poly = self.irr_poly.copy()
        if self.mod is not None:
            copy.mod = self.mod

        return copy

    def coeffmod(self):
        copy = self.copy()
        if self.mod is not None:
            copy.coef = [x % self.mod for x in copy]

        return copy

    def norm(self):
        # infinity norm ||self|| = max_i {|a_i|}
        return max([abs(x) for x in self.coef])

    def set_coeff(self,index,value):
        if index >= len(self.coef):
            self.coef.extend([0]*(index-len(self.coef)))
            self.coef.append(value)
        else:
            self.coef[index] = value

        assert(self.coef[index] == value)# This should be commented
        return

    def __mul__(self,other):
        # self*other
        if self.is_zero():
            return Polynomial(coef=[],irr_poly=self.irr_poly,mod=self.mod)

        # Receives another polynomial or number a and multiplies by this one
        if isinstance(other,Polynomial):
            # 	Modular multiplication of two polynomials
            if other.is_zero():
                return Polynomial(coef=[],irr_poly=self.irr_poly,mod=self.mod)

            s = self.coef
            o = other.coef
            c = [0 for _ in range(len(s) + len(o))]

            for i,a in enumerate(s):
                for j,b in enumerate(o):
                    c[i+j] += long(a*b)

            C = Polynomial(coef=c,irr_poly=self.irr_poly,mod=self.mod).coeffmod()

            # if self.irr_poly is not None and C.deg() >= self.irr_poly.deg():
            #     C = C % self.irr_poly
            return C

        elif type(other) in (int, float,long):
            # Do another thing
            return Polynomial(coef=[coef*other for coef in self.coef],irr_poly=self.irr_poly,mod=self.mod).coeffmod()
        else:
            raise Exception("Unknown operation - %s and %s",type(self),type(other))
    def __rmul__(self,other):
        # other*self
        return self*other


    def __divmod__(self,divisor):
        #return Polynomial(coef=polydiv(self.coef,b.coef)[0].tolist(),irr_poly=self.irr_poly)

       quotient, remainder = Polynomial(mod=self.mod), self.copy()
       divisorDeg = divisor.deg()
       divisorLC = divisor.coef[-1]

       while remainder.deg() >= divisorDeg:
          monomialExponent = remainder.deg() - divisorDeg
          monomialZeros = [0 for _ in range(monomialExponent)]
          monomialDivisor = Polynomial(coef=monomialZeros + [remainder.coef[-1] / divisorLC],irr_poly=self.irr_poly,mod=self.mod)

          quotient += monomialDivisor
          remainder -= monomialDivisor * divisor

       return quotient.coeffmod(), remainder.coeffmod()

    def __div__(self, b):
        # self / b
        if type(b) in (int,long,float):
            return Polynomial(coef = [x / b for x in self.coef])
        else:
            return (self // b)[0]

    def __rdiv__(self,b):
        # b / self
        if type(b) in (int,long,float):
            return (Polynomial(coef=[b],mod=self.mod) // self)[0]
        else:
            return (b // self)[0]

    def __mod__(self, b):
        if type(b) in (int,long,float):
            copy = self.copy()
            copy.coef = [x % b for x in copy]
            return copy

        if is_cyclotomic(b):
            N = b.deg()
            c_a = Polynomial(coef=self.coef[:N])
            c_b = Polynomial(coef=self.coef[N:])
            return c_a - c_b
        else:
            return (self // b)[1]

    def __floordiv__(self,other):
        dividend = self.copy()# Copy

        if type(other) in (int,long,float):
            divisor = self.copy()
            divisor.coef = [other]
        else:
            divisor = other.copy()# Copy
        div = long(divisor.coef[-1])

        if len(dividend) >= len(divisor):
            pass
        else:
            return Polynomial(coef=[0],irr_poly=self.irr_poly,mod=self.mod),Polynomial(coef=dividend,irr_poly=self.irr_poly,mod=self.mod)

        # Copy dividend's attributes and set coeffs to zero
        quot = dividend.copy()
        quot.coef = []

        mult = dividend.copy()
        mult.coef = []

        q = self.mod
        inverse = lambda y,x: (y*pow(x,q-2,q)) % q if q is not None else (float(y) / x)
        #import pdb;pdb.set_trace()

        diff_degree = dividend.deg() - divisor.deg()
        while diff_degree >= 0:
            #Get the next coefficient of the quotient.

            # Solves (-x+1) // 3
            #if mult % 1 != 0:
            #   break
            mult.coef = [0]*diff_degree + [inverse(dividend.coef[-1],div)]

            #quot = mult + quot
            quot += mult

            #Subtract mult * divisor from dividend, but don't bother if mult == 0
            #Note that when i==0, mult!=0; so quot is automatically normalized.
            if mult != 0:
                #d = [(mult * u) % q for u in divisor]
                #dividend = [(u - v) % q for u, v in zip(dividend, d)]
                d = mult*divisor
                dividend -= d
            diff_degree = dividend.deg() - divisor.deg()
            #dividend.pop()
            #divisor.pop(0)

        #return Polynomial(coef=quot,irr_poly=self.irr_poly,mod=self.mod).coeffmod(),Polynomial(coef=dividend,irr_poly=self.irr_poly,mod=self.mod).coeffmod() # Quotient,Remainder
        return quot, dividend

    def __add__(self,x):
        if isinstance(x,Polynomial):
            if x.is_zero():
                return self
            a = self.coef
            b = x.coef

            return Polynomial(coef=[sum(x) for x in itertools.izip_longest(a, b, fillvalue=0)],irr_poly=self.irr_poly,mod=self.mod).coeffmod()
        elif type(x) in (int,long,float):
            if x == 0:
                return self

            a = self.coef
            b = [x]

            return Polynomial(coef=[sum(x) for x in itertools.izip_longest(a, b, fillvalue=0)],irr_poly=self.irr_poly,mod=self.mod).coeffmod()
        else:
            raise Exception("Unknown operation - %s and %s",type(self),type(other))

    def __radd__(self,x):
        if isinstance(x,Polynomial):
            if self.is_zero():
                return x
            elif x.is_zero():
                return self
            else:
                b = x.coef
        elif type(x) in (int,long,float):
            if self.is_zero():
                return Polynomial(coef=[x],irr_poly=self.irr_poly,mod=self.mod).coeffmod()
            elif x == 0:
                return self
            else:
                b = [x]
        else:
            raise Exception("Ops! I don't know how to add this: %s",x)

        # Receives another polynomial  a and sum to this one
        a = self.coef
        return Polynomial(coef=[sum(x) for x in itertools.izip_longest(a, b, fillvalue=0)],irr_poly=self.irr_poly,mod=self.mod).coeffmod()

    def __sub__(self,x):
        # self - x
        if isinstance(x,Polynomial):
            if self.is_zero():
                return x
            elif x.is_zero():
                return self
            else:
                b = [-c for c in x.coef]
        elif type(x) in (int,long,float):
            if self.is_zero():
                return Polynomial(coef=[-x],irr_poly=self.irr_poly,mod=self.mod).coeffmod()
            elif x == 0:
                return self
            else:
                b = [-x]
        else:
            raise Exception("Invalid operation for %s and Polynomial" % type)

        a = self.coef
        return Polynomial(coef=[sum(x) for x in itertools.izip_longest(a, b, fillvalue=0)],irr_poly=self.irr_poly,mod=self.mod).coeffmod()

    def __rsub__(self,x):
        # x - self
        if isinstance(x,Polynomial):
            if self.is_zero():
                return x
            elif x.is_zero():
                return self
            else:
                b = x.coef
        elif type(x) in (int,long,float):
            if self.is_zero():
                return Polynomial(coef=[x],irr_poly=self.irr_poly,mod=self.mod).coeffmod()
            elif x == 0:
                return -self
            else:
                b = [x]
        else:
            raise Exception("Ops! I don't know how to add this: %s",b)

        a = [-c for c in self.coef]
        return Polynomial(coef=[sum(c) for c in itertools.izip_longest(a, b, fillvalue=0)],irr_poly=self.irr_poly,mod=self.mod).coeffmod()

    # Shift n bits to the least significant bit
    def __rshift__(self,n):
    	assert type(n) == int
        assert n >= 0
        v = list(self.coef)

        for count in range(n):
            v.pop(0) # Faz o shift no sentido do lsb

            v = v + [0]
    	return Polynomial(coef=v,irr_poly=self.irr_poly,mod=self.mod).coeffmod()

    # Shift n bits to the most significant bit
    def __lshift__(self,n):
    	assert type(n) == int
        assert n >= 0
    	v = list(self.coef)
        for count in range(n):
            #v.pop() # Faz o shift no sentido do msb

            v = [0] + v
        return Polynomial(coef=v,irr_poly=self.irr_poly,mod=self.mod).coeffmod()

    def __str__(self):
        return str(self.coef)

    def __repr__(self):
        output = ""
        for i in range(len(self.coef)):
            if self.coef[i] != 0:
                output +="%d*x^%d + " %(self.coef[i],i)
        return output[:-3]

    def __len__(self):
        if type(self.coef) in (list,tuple) and self != Polynomial(0):
            return len(self.coef)
        else:
            return 0

    def __abs__(self):
        return len(self)

    def __neg__(self):
        return Polynomial(coef=[-a for a in self],irr_poly=self.irr_poly,mod=self.mod).coeffmod()

    def __iter__(self):
        return iter(self.coef)

    def __eq__(self,other):
        if isinstance(other,Polynomial):
            b = other.copy()
        elif type(other) in (list,tuple):
            b = Polynomial(coef=other)
        elif type(other) in (int,long,float):
            b = Polynomial(coef=[other])
        else:
            return False

        #return repr(self) == repr(b)
        if self.mod is None:
            return self.coef == b.coef
        else:
            return [x % self.mod for x in self.coef] == [x % self.mod for x in b.coef]

    def __ne__(self,other):
        return not self.__eq__(other)

    # <
    def __lt__(self,b):
        return self.deg() < b.deg()
    # <=
    def __le__(self,b):
        return self.deg() <= b.deg()
    # >=
    def __ge__(self,b):
        return self.deg() >= b.deg()
    # >
    def __gt__(self,b):
        return self.deg() > b.deg()

    def __pow__(self,n):
        if n == 0:
            return Polynomial(1)
        if n < 0:
            return functools.reduce(operator.mul,[self]*n).inverse()
        return functools.reduce(operator.mul,[self]*n)

    def __iter__(self):
        return iter(self.coef)
