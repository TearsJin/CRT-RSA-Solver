class TK15PKESAttack:
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

    def attack(self,m = None, eta = None, tau = None):
        N,e,d0,alpha,dpbits,unknowbits = self.N,self.e,self.d0,self.alpha,self.dpbits,self.unknowbits
        delta = unknowbits / N.nbits()
        beta = 1 / 2
        kbits = dpbits - unknowbits
        M = 2 ^ kbits
        if eta == None:
            eta = (1 - 2 * delta) / 2
            print("[+] Set eta =",eta)
        if tau == None:
            tau = (sqrt(1 - 4 * delta) - 2 * delta) / 2
            print("[+] Set tau =",tau)

        PR.<x,y,z,z_> = PolynomialRing(ZZ)
        c = 1 - e * d0 * M
        f = c + e * x + y * (z - 1)
        bounds = [int(N ^ delta),int(N ^ (alpha + beta)),int(N ^ 0.5),int(N ^ 0.5)]
        
        Monomials = f.monomials()
        Coefficients = f.coefficients()
        a = Integer(Coefficients[-1])

        for i in range(n):
            while gcd(bounds[i],a) != 1:
                bounds[i] -= 1

        # define W
        factors = [monomial(*bounds) for monomial in Monomials]
        W = max([i * j for i,j in zip(Coefficients,factors)])
        
        
        # extend z