from time import time
class TK16PKESAttack:
    def __init__(self) -> None:
        # if None in (N,e,d0,dpbits,unknowbits):
        #     print("[+] Error, missing variable!")
        #     return
    
        # self.N = N
        # self.e = e
        # self.d0 = d0
        # self.alpha = float(log(e,N))
        # self.dpbits = dpbits
        # self.unknowbits = unknowbits
        pass

    def attack(self):
        # alpha = 3 / 8 beta = 1 / 2
        N = 113309068152935568893932620926523748561925916580328833761833709492836549360889407465934366065781314829245822754800772426000185646006813432371474526931829549338363892243205465109470135176086384087827228420777587664921704371825877672283627971299025154129511721048983330050794358456268791642813976689918889182193
        # alpha,e,d0 = self.alpha,self.e,self.d0
        delta = 0.1
        M = 2 ^ int(1024 * delta)
        e = 508
        d0 = 10
        alpha = 3 / 8   
        m = 7
        eta = 2 / 5
        s = floor(eta * m)
        print(eta,s)

        PR.<x,y,z> = PolynomialRing(ZZ)
        f = 1 - e * d0 + x * (y - 1)
        bounds = [int(alpha),int(N ^ 0.5),int(N ^ 0.5)]
        g = lambda i,j: (e * M) ^ (m - i) * x ^ j * z ^ s * f ^ i
        g_= lambda i,j: (e * M) ^ (m - i) * y ^ j * z ^ s * f ^ i

        Ix = [(i,j) for i in range(m + 1) for j in range(m - i + 1)]
        Ix_= [(i,j) for i in range(m + 1) for j in range(1, floor(s - 2 * delta * i) + 1)]


        ShiftPolys = Sequence([],f.parent())
        for _ in Ix:
            poly = g(*_)
            ShiftPolys.append(poly)
        for _ in Ix_:
            poly = g_(*_)
            ShiftPolys.append(poly)
        # 获得系数矩阵
        B, monomials = ShiftPolys.coefficient_matrix()
        monomials = vector(monomials)
        nn = len(monomials)
        print(nn,len(ShiftPolys))
        print("[+] Successfully defined polynomial.There are {} polynomials in total.".format(str(nn)))
        # # 获得单项式的顺序
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
        # B = B.change_ring(QQ)
        # for i, factor in enumerate(factors):
        #     B.rescale_col(i, 1/factor)
        # print("[+] Start finding the root...")