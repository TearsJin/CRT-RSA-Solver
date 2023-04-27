from math import prod
from time import time

N = 96418199744999808107808846459563245497752175518619867318081906772690766133904963002214789660537415496783249098583656299436686012061837292949321273850541814547052967090968681816906678221128983864333369509124282269679657826985653407375918917942608586289464846864623752502019942648264657027523012156001908149531

P = lambda f,i:(f^i).monomials() if i >= 0 else []
g = lambda monomial,l,k,f,m:monomial * f ^ k / l ^ k  * N ^ (m-k)
M = lambda f,l,m,k: list(filter(lambda x:x!=-1,[Xs if Xs/(l^k) in P(f,m-k) else -1 for Xs in P(f,m)]))


bounds = [1,1]

def issub(P1,P2):
    for Xs in P1:
        if Xs not in P2:
            return False
    return True


m = 2
# monic

PR.<x,y> = PolynomialRing(Zmod(N))
f = x ^ 2 * y + x * y ^ 2 + 1 
f = f * inverse_mod(Integer(f.coefficients()[0]),N)
f = f.change_ring(ZZ)
PR.<x,y> = PolynomialRing(ZZ)
l = f.monomials()[0]


for i in range(m+1):
    if not issub(P(f,i),P(f,m)):
        print("need change M")
        def M(f,l,m,k):
            res = []
            conditions1 = []
            for j in range(1,m+1):
                for Xs in P(f,j):
                    conditions1.append(Xs)
            for j in range(1,m+1-k):
                for Xs in P(f,j):
                    conditions2.append(Xs)
            for Xs in conditions1:
                if Xs / (l^k) in conditions2:
                    res.append(Xs)
            return res    

ShiftPolys = Sequence([],f.parent())

for k in range(m+1):
    Mk1 = M(f,l,m,k+1)
    Mk = M(f,l,m,k)
    for monomial in Mk:
        if monomial not in Mk1:
            ShiftPolys.append(g(monomial,l,k,f,m))

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