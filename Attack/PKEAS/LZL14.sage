from time import time

'''
    LZL14 - LSBs
    dp < N ^ beta
    e < N ^ alpha
    d1 < N ^ delta

    dp = d1 * M + d0
'''

N = Integer(123718923712936178236172831212381902312371892837122111111111111111)
e = 65537
d0 = 123213
alpha = float(log(e,N))
delta = 0.1
beta = 0.5
kbits = int(N.nbits() * delta ) # known bits
M = 2 ^ kbits # dp = d1 * M + d0
bounds = [N^(beta + alpha - 0.5),N ^ 0.5,N ^ 0.5]
# tau = 0.5 - 2 * delta 
# sigma = 0.5 + delta 
tau = 2 / 3
sigma = 2 / 3
m = 7
print(tau,sigma)
t = int(tau * m)
s = int(sigma * m)

PR.<x,y,z> = PolynomialRing(ZZ)
f = x * (y - 1) - e * d0 + 1 
g = lambda i,j:(e * M) ^ (m - i) * x ^ j * z ^ s * f ^ i
h = lambda i,j:(e * M) ^ (m - i) * y ^ j * z ^ s * f ^ i
Ix = [(i,j) for i in range(m + 1) for j in range(m - i + 1)]
Iy = [(i,j) for i in range(m + 1) for j in range(1, t + 1)]

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
print("[+] Successfully defined polynomial.There are {} polynomials in total.".format(len(ShiftPolys)))
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
            root = tuple(root[var] for var in [x,y,z])
            roots.append(root)
        for root in roots:
            if root != (1,0):
                print("[+] Successfully found the roots:",root)
                FindRootTime = time()-S
                print("[+] The time taken is",FindRootTime,"s")


alpha = 0.2
print((5 - 2 * sqrt(1 + 6 * (alpha + 1/2))) / 6)