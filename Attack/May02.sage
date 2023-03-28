from Crypto.Util.number import *
from random import getrandbits,randrange
from time import time


class May02Attack:
    def __init__(self,N=None,e=None,beta=None,delta=None) -> None:
        if None in (N,e,beta,delta):
            print("[+] Error, missing variable!")
            return
    
        self.N = N
        self.e = e
        self.alpha = float(log(e,N))
        self.beta = beta
        self.delta = delta

    def tolatex(self,string = '', old = [], new = [],type = ''):
        string = "\\begin{" + type +"}\n" + string + "\n\\end{" + type + "}" 
        for x,y in zip(old,new):
            string = string.replace(x,y)
        while '&&' in string:
            string = string.replace("&&",'&')
        # 有无可以直接把e N这两个系数的位置直接调好的办法？
        return string.replace("[&"," ").replace("]","\\\\").replace("*"," ")

    def showMatrix(self,m = 3,t = 2,output = False ,filename = './Matrix.data'):
        PR.<x,y,N,e> = PolynomialRing(ZZ)
        f = x * (N - y) - N

        g = lambda i,j: x ^ j * f ^ i * e ^ (m - i)
        h = lambda i,j: y ^ j * f ^ i * e ^ (m - i)
        Ix = [(i,j) for i in range(0,m + 1) for j in range(0,m - i + 1)]
        Iy = [(i,j) for i in range(0,m + 1) for j in range(1,t+1)] 
        ShiftPolys = Sequence([],f.parent())
        for _ in Ix:
            poly = g(*_)
            ShiftPolys.append(poly)
        for _ in Iy:
            poly = h(*_)
            ShiftPolys.append(poly)

        monomialOrder = []
        for Poly in ShiftPolys:
            for x in Poly.monomials():
                if x(e=1,N=1) not in monomialOrder:
                    monomialOrder.append(x(e=1,N=1))
        nn = len(ShiftPolys)
        M = [[0 for _ in range(nn)] for _ in range(nn)]
        for i in range(nn):
            for monomial,b in zip(ShiftPolys[i].monomials(),ShiftPolys[i].coefficients()):
                M[i][monomialOrder.index(monomial(e=1,N=1))] = b * monomial
        print(str(Matrix(PR,M)))
        if output:
            with open(filename,'w') as f:
                f.write(str(Matrix(PR,M)) + '\n')
                f.write(self.tolatex(str(Matrix(PR,M)),[' ','x','y'],['&','X','Y'],"vmatrix"))

    def attack(self,m = 3,t = 2,detail = False):
        N,e,Alpha,Beta,Delta = self.N,self.e,self.alpha,self.beta,self.delta
        bounds = [Integer(int(N ^ ((Alpha + Beta + Delta) - 1))),Integer(int(N ^ Beta))]

        # 定义多项式
        PR.<y,z> = PolynomialRing(ZZ)
        f = y * (N - z) - N
        g = lambda i,j: y ^ j * f ^ i * e ^ (m - i)
        h = lambda i,j: z ^ j * f ^ i * e ^ (m - i)
        Ix = [(i,j) for i in range(0,m + 1) for j in range(0,m - i + 1)]
        Iy = [(i,j) for i in range(0,m + 1) for j in range(1,t+1)] 
        ShiftPolys = Sequence([],f.parent())
        for _ in Ix:
            poly = g(*_)
            ShiftPolys.append(poly)
        for _ in Iy:
            poly = h(*_)
            ShiftPolys.append(poly)
        # 获得系数矩阵
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
                    root = tuple(root[var] for var in [y,z])
                    roots.append(root)
                for root in roots:
                    if root != (1,0):
                        print("[+] Successfully found the roots:",root)
                        FindRootTime = time()-S
                        print("[+] The time taken is",FindRootTime,"s")
                        if detail:
                            print("[+] ------------------------- May03Attack ------------------------- ")
                            print("[+] Condition: Alpha = {}, Beta = {}, Delta = {}".format(str(Alpha),str(Beta),str(Delta)))
                            print("[+] Parameters: m = {}, t = {}".format(str(m),str(t)))
                            print("[+] The Lattice like this:")
                            self.display(ShiftPolys,monomialOrder)
                            print("[+] LLL time: {}s, Finding roots time: {}s, total time: {}s".format(str(LLLTime),str(FindRootTime),str(LLLTime + FindRootTime)))
                        return root
                print("[+] Can't not find the roots")
                return roots
    def display(self,ShiftPolys,monomialOrder):
        for Poly in ShiftPolys:
            ZERO = ['0'] * len(monomialOrder)
            monomials = Poly.monomials()
            for monomial in monomials:
                ZERO[monomialOrder.index(monomial)] = 'X'
            print(' '.join(ZERO))
