

# This file was *autogenerated* from the file ./Attack/JM.sage
from sage.all_cmdline import *   # import sage library

_sage_const_96418199744999808107808846459563245497752175518619867318081906772690766133904963002214789660537415496783249098583656299436686012061837292949321273850541814547052967090968681816906678221128983864333369509124282269679657826985653407375918917942608586289464846864623752502019942648264657027523012156001908149531 = Integer(96418199744999808107808846459563245497752175518619867318081906772690766133904963002214789660537415496783249098583656299436686012061837292949321273850541814547052967090968681816906678221128983864333369509124282269679657826985653407375918917942608586289464846864623752502019942648264657027523012156001908149531); _sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_2 = Integer(2)
from math import prod
from time import time

N = _sage_const_96418199744999808107808846459563245497752175518619867318081906772690766133904963002214789660537415496783249098583656299436686012061837292949321273850541814547052967090968681816906678221128983864333369509124282269679657826985653407375918917942608586289464846864623752502019942648264657027523012156001908149531 

P = lambda f,i:(f**i).monomials() if i >= _sage_const_0  else []
g = lambda monomial,l,k,f,m:monomial * f ** k / l ** k  * N ** (m-k)
M = lambda f,l,m,k: list(filter(lambda x:x!=-_sage_const_1 ,[Xs if Xs/(l**k) in P(f,m-k) else -_sage_const_1  for Xs in P(f,m)]))



bounds = [_sage_const_1 ,_sage_const_1 ]

def issub(P1,P2):
    for Xs in P1:
        if Xs not in P2:
            return False
    return True


m = _sage_const_2 
# monic

PR = PolynomialRing(Zmod(N), names=('x', 'y',)); (x, y,) = PR._first_ngens(2)
f = x ** _sage_const_2  * y + x * y ** _sage_const_2  + _sage_const_1  
f = f * inverse_mod(Integer(f.coefficients()[_sage_const_0 ]),N)
f = f.change_ring(ZZ)
PR = PolynomialRing(ZZ, names=('x', 'y',)); (x, y,) = PR._first_ngens(2)
l = f.monomials()[_sage_const_0 ]

ShiftPolys = Sequence([],f.parent())
for i in range(m+_sage_const_1 ):
    if not issub(P(f,i),P(f,m)):
        print("need change M")
        def M(f,l,m,k):
            res = []
            conditions1 = []
            for j in range(_sage_const_1 ,m+_sage_const_1 ):
                for Xs in P(f,j):
                    conditions1.append(Xs)
            for j in range(_sage_const_1 ,m+_sage_const_1 -k):
                for Xs in P(f,j):
                    conditions2.append(Xs)
            for Xs in conditions1:
                if Xs / (l**k) in conditions2:
                    res.append(Xs)
            return res    
            
for k in range(m+_sage_const_1 ):
    Mk1 = M(f,l,m,k+_sage_const_1 )
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
LLLTime = _sage_const_0 
S = time()
B = B.dense_matrix().LLL()
LLLTime = time()-S
print("[+] LLL Done! The time taken is",LLLTime,"s")
B = B.change_ring(QQ)
for i, factor in enumerate(factors):
    B.rescale_col(i, _sage_const_1 /factor)
print("[+] Start finding the root...")
FindRootTime = _sage_const_0 
S = time()
# 求根
H = Sequence([], f.parent().change_ring(QQ))
for h in filter(None, B*monomials):
    H.append(h)
    I = H.ideal()
    if I.dimension() == -_sage_const_1 :
        H.pop()
    elif I.dimension() == _sage_const_0 :
        roots = []
        for root in I.variety(ring=ZZ):
            root = tuple(root[var] for var in [y,z])
            roots.append(root)
        for root in roots:
            if root != (_sage_const_1 ,_sage_const_0 ):
                print("[+] Successfully found the roots:",root)
                FindRootTime = time()-S
                print("[+] The time taken is",FindRootTime,"s")
