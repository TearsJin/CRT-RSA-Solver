class TLP19Attack():
    def __init__(self,N=None,e=None,delta=None) -> None:
        if None in (N,e,delta):
            print("[+] Error, missing variable!")
            return

        self.N = N
        self.e = e
        self.alpha = float(log(e,N))
        self.delta = delta

    def attack(self,m = 3, tau = None,detail = False):  
        N,e,Alpha,Delta = self.N,self.e,self.alpha,self.delta
        PR.<xp1,xq1,xp2,xq2,yp,yq> = PolynomialRing(ZZ)
        vars = (xp1,xq1,xp2,xq2,yp,yq)   
        bounds = [int(N ^ ((Alpha + Delta) - 1/2)),int(N ^ ((Alpha + Delta) - 1/2)),int(N ^ ((Alpha + Delta) - 1/2)),int(N ^ ((Alpha + Delta) - 1/2)),int(N ^ 0.5),int(N ^ 0.5)]

        if tau == None:
            tau = 1 - 2 * Delta
            print("[+] Set tau =",tau)

        # define function and polynomials
        f1 = N + xp1 * (N - yp)
        f2 = 1 + xp2 * (yq - 1)
        h = (N - 1) * xp1 * xp2 + xp1 + N * xp2
     
        g1 = lambda i1,i2,j1,j2,u:xp1 ^ j1 * xp2 ^ j2 * yq ^ floor((i1 + i2) / 2) * f1 ^ i1 * f2 ^ i2 * h ^ u * e ^ (m - i1 - i2 -u)
        gp = lambda i1,i2,j1:yq ^ (floor((i1 + i2) / 2) - j1) * f1 ^ i1 * f2 ^ i2  * e ^ (m - i1 - i2)
        gq = lambda i1,i2,j1:yq ^ (floor((i1 + i2) / 2) + j1) * f1 ^ i1 * f2 ^ i2  * e ^ (m - i1 - i2)


        # define Index 
        MHALF = floor(m / 2)
        Ix = [(i1,i2,0,0,u) for i1 in range(MHALF + 1) for i2 in range(MHALF + 1) for u in range(min(MHALF - i1 ,MHALF - i2) + 1)]
        Ix += [(i1,i2,1,0,u) for i1 in range(MHALF) for i2 in range(1,MHALF + 1) for u in range(min(MHALF - i1 - 1,MHALF - i2 ) + 1)]
        Ix += [(i1,0,j1,0,u) for i1 in range(MHALF + 1) for j1 in range(1,MHALF - i1 + 1) for u in range(MHALF - i1 - j1 + 1)]
        Ix += [(0,i2,0,j2,u) for i2 in range(MHALF + 1) for j2 in range(1,MHALF - i2 + 1) for u in range(MHALF - i2 - j2 + 1)]

        Iyp = [(i1,i2,j1) for i1 in range(MHALF + 1) for i2 in range(MHALF + 1) for j1 in range(1,floor(tau * (i1 + i2))-ceil((i1 + i2) / 2) + 1)]
        Iyq = [(i1,i2,j1) for i1 in range(MHALF + 1) for i2 in range(MHALF + 1) for j1 in range(1,floor(tau * (i1 + i2))-floor((i1 + i2) / 2) + 1)]


        ShiftPolys = Sequence([],h.parent())
        for _ in Ix:
            poly = g1(*_)
            ShiftPolys.append(self.transform(poly,vars))
        for _ in Iyp:
            poly = gp(*_)
            ShiftPolys.append(self.transform(poly,vars))
        for _ in Iyq:
            poly = gq(*_)
            ShiftPolys.append(self.transform(poly,vars))

        monomialOrder = []
        for Poly in ShiftPolys:
            for x in Poly.monomials():
                if x not in monomialOrder:
                    monomialOrder.append(x)
            
        for i in range(len(ShiftPolys)):
            Poly = ShiftPolys[i]
            # 除去e后如果非1 -1，则带有N
            while Poly % e == 0:
                Poly = Poly // e
            # 如果有N，则肯定不存在e ^ m，所以直接算就好了
            while abs(self.leadCoe(Poly,monomialOrder)) > 10000:
                Poly = (Poly * inverse_mod(N,e ^ m)) % e ^ m
                ShiftPolys[i] = ((ShiftPolys[i] * inverse_mod(N,e ^ m)) % e ^ m)
                


        B, monomials = ShiftPolys.coefficient_matrix()
        monomials = vector(monomials)
        nn = len(monomials)
        print("[+] Successfully defined polynomial.There are {} polynomials in total.".format(len(ShiftPolys)))

        print(len(ShiftPolys),nn)
        print(monomials)
        factors = [monomial(*bounds) for monomial in monomials]
        for i, factor in enumerate(factors):
            B.rescale_col(i, factor)
        # LLL
        print("[+] Start LLL...")
        LLLTime = 0
        S = time()
        B = B.dense_matrix().LLL()[-len(ShiftPolys):]
        print(len(B))
        LLLTime = time()-S
        print("[+] LLL Done! The time taken is",LLLTime,"s")
        B = B.change_ring(QQ)
        for i, factor in enumerate(factors):
            B.rescale_col(i, 1/factor)
        RES = []
        for i in filter(None, B * vector(monomials)):
            while 'yq' in str(i):
                i = self.transform(i * yp,vars)
            RES.append(i)
        # print(RES)
        # print("[+] Start finding the root...")
        # PR.<xp,yp> = PolynomialRing(QQ)
        # H = Sequence([], PR)
        # for h in RES:
        #     H.append(h)
        #     I = H.ideal()
        #     if I.dimension() == -1:
        #         H.pop()
        #     elif I.dimension() == 0:
        #         roots = []
        #         for root in I.variety(ring=ZZ):
        #             root = tuple(root[var] for var in [xp1,xp2,yp])
        #             roots.append(root)
        #         for root in roots:
        #             if root != (1,0) and root != (-1,0):
        #                 print("[+] Successfully found the roots:",root)
        #                 FindRootTime = time()-S
        #                 print("[+] The time taken is",FindRootTime,"s")
        #                 if detail:
        #                     print("[+] ------------------------- TLP19Attack ------------------------- ")
        #                     print("[+] Condition: Alpha = {}, Beta = {}, Delta = {}".format(str(Alpha),str(Beta),str(Delta)))
        #                     print("[+] Parameters: m = {}, LAMBDA = {},tau = {}".format(str(m),str(LAMBDA),str(tau)))
        #                     print("[+] The Lattice like this:")
        #                     self.display(ShiftPolys,monomialOrder)
        #                     print("[+] LLL time: {}s, Finding roots time: {}s, total time: {}s".format(str(LLLTime),str(FindRootTime),str(LLLTime + FindRootTime)))
        #                 return root
        #         return roots

    def transform(self,poly,vars):
        (xp1,xq1,xp2,xq2,yp,yq) = vars
        vars = [xp1,xq1,xp2,xq2,yp,yq]
        Exponents = poly.exponents()
        Coefficients = poly.coefficients()
        newPoly = 0
        for exp,coe in zip(Exponents,Coefficients):
            Newexp = [i for i in exp]
            POW = 0
            monomials = 0
            # replace yp * yq to N
            while Newexp[-1] * Newexp[-2] >= 1:
                Newexp[-1] -= 1
                Newexp[-2] -= 1
                POW += 1
            monomials += self.prod([i^j for i,j in zip(vars,Newexp)]) * (self.N ^ POW) * coe
            # replace xq1 if the powers of yp are nonagative, by xq1 = xp1 + 1
            # replace xp2 if the powers of yq are positive, by xp2 = xq2 - 1
            if Newexp[-2] >= 1:
                monomials = monomials(xq1 = xp1 + 1)
                monomials = monomials(xq2 = xp2 + 1)
            else:
                monomials = monomials(xp1 = xq1 - 1)
                monomials = monomials(xp2 = xq2 - 1)
            newPoly += monomials
        return newPoly
    
    def leadCoe(self,Poly,monomialOrder):
        monomials = Poly.monomials()
        coes = Poly.coefficients()
        List = sorted([(monomials.index(i),monomialOrder.index(i),i,j) for i,j in zip(monomials,coes)],key=lambda x:x[1],reverse=True)
        return List[0][3]

    def prod(self,LIST):
        PROD = 1
        for i in LIST:
            PROD *= i
        return PROD

    def display(self,ShiftPolys,monomialOrder):
        for Poly in ShiftPolys:
            ZERO = ['0'] * len(monomialOrder)
            monomials = Poly.monomials()
            for monomial in monomials:
                ZERO[monomialOrder.index(monomial)] = 'X'
            print(' '.join(ZERO))