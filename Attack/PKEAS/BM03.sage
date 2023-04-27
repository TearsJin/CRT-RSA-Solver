'''

PKEAS-BM03

'''
class BM03Attack:
    def __init__(self,N=None,e=None,d0=None,dpbits=None,unknowbits=None) -> None:
        if None in (N,e,d0,dpbits,unknowbits):
            print("[+] Error, missing variable!")
            return
    
        self.N = N
        self.e = e
        self.d0 = d0
        self.alpha = float(log(e,N))
        self.dpbits = dpbits
        self.unknowbits = unknowbits

    def attack(self):
        N,e,d0,Alpha,dpbits,unknowbits = self.N,self.e,self.d0,self.alpha,self.dpbits,self.unknowbits
        kbits = dpbits - unknowbits
        M = 2 ^ kbits
        E = inverse_mod(e * M , N)
        PR.<x> = PolynomialRing(Zmod(N))

        for k in range(1,e):
            f = E * (e * d0 + k - 1) + x
            print(f.small_roots(X = 2 ^ unknowbits))
# N = 123718923712936178236172831212381902312371892837122111111111111111
# e = 65537
# d0 = 123213
# alpha = float(log(e,N))
# delta = 0.1
# beta = 0.5
# kbits = Integer(N.bits() * delta ) # known bits
# M = 2 ^ kbits # dp = d1 * M + d0

# # LSBs



# # MSBs
# PR.<x> = PolynomialRing(Zmod(N))
# f = e * d0 * M - Integer(N^(0.25)) + x
# f.small_roots(X = Integer(N ^ alpha))

# f = e * d0 * M + Integer(N^(0.25)) + x
# f.small_roots(X = Integer(N ^ alpha))




