class BM06Attack:
    def __init__(self,N=None,e=None,beta=None,delta=None) -> None:
        if None in (N,e,beta,delta):
            print("[+] Error, missing variable!")
            return
    
        self.N = N
        self.e = e
        self.alpha = float(log(e,N))
        self.beta = beta
        self.delta = delta

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

    def transform(self,poly,vars):
        (xp,yq,yp) = vars
        vars = [xp,yq,yp]
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
            newPoly += monomials
        return newPoly

    def attack(self, m = 3, s = None, t = None, detail = False):
        N,e,Alpha,Beta,Delta = self.N,self.e,self.alpha,self.beta,self.delta
        bounds = [int(N ^ ((Alpha + Beta + Delta) - 1)),int(N ^ (1 - Beta)),int(N ^ Beta)]
        if s == None:
            s = int((1 - Beta - Delta)/(2 * (1 - Beta)) * m)
            print("[+] Set s =",s)
        if t == None:
            t = int(((1 - Beta) ^ 2 - Delta) / (2 * Beta * (1 - Beta)) * m)
            print("[+] Set t =",t)


        PR.<xp,yq,yp> = PolynomialRing(ZZ)  
        fp = N + xp * (N - yp)
        vars = (xp,yq,yp)
        g = lambda i,j,s: xp ^ j * fp ^ i * e ^ (m - i) * yq ^ s
        h = lambda i,j,s: yp ^ j * fp ^ i * e ^ (m - i) * yq ^ s

        Ix = [(i,j,s) for i in range(0,m + 1) for j in range(0,m - i + 1)]
        Iy = [(i,j,s) for i in range(0,m + 1) for j in range(1,t + 1)] 
        ShiftPolys = Sequence([],fp.parent())
        for _ in Ix:
            poly = g(*_)
            ShiftPolys.append(self.transform(poly,vars))
        for _ in Iy:
            poly = h(*_)
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
            if abs(self.leadCoe(Poly,monomialOrder)) != 1:
                ShiftPolys[i] = ((ShiftPolys[i] * inverse_mod(N,e ^ m)) % e ^ m).change_ring(ZZ)
        B, monomials = ShiftPolys.coefficient_matrix()
        monomials = vector(monomials)
        nn = len(monomials)
        print("[+] Successfully defined polynomial.There are {} polynomials in total.".format(str(nn)))

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

        RES = []
        for i in filter(None, B * vector(monomials)):
            while 'yq' in str(i):
                i = self.transform(i * yp,vars)
            RES.append(i)
        PR.<xp,yp> = PolynomialRing(QQ)
        H = Sequence([], PR)
        for h in RES:
            H.append(h)
            I = H.ideal()
            if I.dimension() == -1:
                H.pop()
            elif I.dimension() == 0:
                roots = []
                for root in I.variety(ring=ZZ):
                    root = tuple(root[var] for var in [xp,yp])
                    roots.append(root)
                for root in roots:
                    if root != (1,0) and root != (-1,0):
                        print("[+] Successfully found the roots:",root)
                        FindRootTime = time()-S
                        print("[+] The time taken is",FindRootTime,"s")
                        if detail:
                            print("[+] ------------------------- BM06Attack ------------------------- ")
                            print("[+] Condition: Alpha = {}, Beta = {}, Delta = {}".format(str(Alpha),str(Beta),str(Delta)))
                            print("[+] Parameters: m = {}, t = {},s = {}".format(str(m),str(t),str(s)))
                            print("[+] The Lattice like this:")
                            self.display(ShiftPolys,monomialOrder)
                            print("[+] LLL time: {}s, Finding roots time: {}s, total time: {}s".format(str(LLLTime),str(FindRootTime),str(LLLTime + FindRootTime)))
                        return root
                return roots
    def display(self,ShiftPolys,monomialOrder):
        for Poly in ShiftPolys:
            ZERO = ['0'] * len(monomialOrder)
            monomials = Poly.monomials()
            for monomial in monomials:
                ZERO[monomialOrder.index(monomial)] = 'X'
            print(' '.join(ZERO))
