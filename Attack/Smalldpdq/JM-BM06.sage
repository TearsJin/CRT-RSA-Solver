from time import time

class JM_BM06Attack():
    def __init__(self,N=None,e=None,delta=None) -> None:
        if None in (N,e,delta):
            print("[+] Error, missing variable!")
            return

        self.N = N
        self.e = e
        self.alpha = float(log(e,N))
        self.delta = delta

    def attack(self,m = 2,t = None, detail = False):
        N,e,Alpha,Delta = self.N,self.e,self.alpha,self.delta
        bounds = [int(N ^ Delta),int(N ^ Delta),int(N^(Alpha+Delta - 0.5)),int(N^(Alpha+Delta - 0.5))]
        if t == None:
            tau = 0.381788
            t = int(tau * m)
            print("[+] Set t =",t)
        n = len(bounds)

        # define polynomial
        PR.<x1,x2,x3,x4> = PolynomialRing(ZZ)
        f = e ^ 2 * x1 * x2 + e * x1 * x4 - e * x1 + e * x2 *x3 - e * x2 - (N - 1) * x3 * x4 - x3 - x4 + 1
        vars = [x1,x2,x3,x4]

        Monomials = f.monomials()
        Coefficients = f.coefficients()
        a = Integer(Coefficients[-1])

        for i in range(n):
            while gcd(bounds[i],a) != 1:
                bounds[i] -= 1

        # define W
        factors = [monomial(*bounds) for monomial in Monomials]
        W = max([i * j for i,j in zip(Coefficients,factors)])

        # JM extend, define S and M
        if Alpha > 0.5:
            # extend x1,x2
            S = []
            for j1 in range(t + 1):
                for j2 in range(t + 1):
                    for monomial in (f ^ (m - 1)).monomials():
                        new_monomial = monomial * x1 ^ j1 * x2 ^ j2
                        if new_monomial not in S:
                            S.append(new_monomial)
            # define lj
            L = [0 for _ in range(n)]
            for j in range(n):
                Max = 0
                for exp in sum(S).exponents():
                    if exp[j] > Max:
                        Max = exp[j]
                L[j] = Max
            # define R
            R =  W * prod([bounds[j] ^ L[j] for j in range(n)])

        # define M
        M = []
        for monomial in S:
            for monomial_ in (monomial * f).monomials():
                if monomial_ not in M:
                    M.append(monomial_) 

        g = lambda monomial,f:monomial*f*prod([(bounds[j]^(L[j]-monomial.exponents()[0][j])) for j in range(n)])
        g_ = lambda monomial:monomial * R

        ShiftPolys = Sequence([],f.parent())
        for monomial in S:
            ShiftPolys.append(g(monomial,f))
        for monomial in M:
            if monomial not in S:
                ShiftPolys.append(g_(monomial))


        B, monomials = ShiftPolys.coefficient_matrix()
        monomials = vector(monomials)
        nn = len(monomials)
        print("[+] Successfully defined polynomial.There are {} polynomials in total.".format(str(nn)))
        # 获得单项式的顺序
        monomialOrder = []
        for Poly in ShiftPolys:
            for x in Poly.monomials():
                if x not in monomialOrder:
                    monomialOrder.append(x)
        # 将上界乘进系数矩阵
        factors = [monomial(*bounds) for monomial in monomials]
        for i, factor in enumerate(factors):
            B.rescale_col(i, factor)
        # LLL
        print("[+] Start LLL...")
        LLLTime = 0
        S = time()
        B = B.dense_matrix().LLL()
        LLLTime = time()-S
        print("[+] LLL Done! The time taken is",LLLTime,"s")
        B = B.change_ring(QQ)
        for i, factor in enumerate(factors):
            B.rescale_col(i, 1/factor)
        print("[+] Start finding the root...")
        FindRootTime = 0
        S = time()
        # 求根
        H = Sequence([], f.parent().change_ring(QQ))
        for h in filter(None, B*monomials):
            H.append(h)
            I = H.ideal()
            if I.dimension() == -1:
                H.pop()
            elif I.dimension() == 0:
                roots = []
                for root in I.variety(ring=ZZ):
                    root = tuple(root[var] for var in [x1,x2,x3,x4])
                    roots.append(root)
                for root in roots:
                    if 0 not in root:
                        print("[+] Successfully found the roots:",root)
                        FindRootTime = time()-S
                        print("[+] The time taken is",FindRootTime,"s")