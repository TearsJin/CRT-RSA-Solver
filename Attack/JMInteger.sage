from time import time


bounds = [2,2]
m = 3
t = 0
n = len(bounds)

PR.<x,y> = PolynomialRing(ZZ)
f = x ^ 2 * y + x * y ^ 2 + 2
vars = [x,y]
Monomials = f.monomials() 
Coefficients = f.coefficients()
a = Integer(Coefficients[-1])

for i in range(n):
    while gcd(bounds[i],a) != 1:
        bounds[i] +=1 

# define dj
D = [0 for _ in range(len(vars))]
for j in range(len(D)):
    Max = 0
    for exp in f.exponents():
        if exp[j] > Max:
            Max = exp[j]
    D[j] = Max
factors = [monomial(*bounds) for monomial in Monomials]
W = max([i * j for i,j in zip(Coefficients,factors)])
R = W * prod([bounds[j] ^ (D[j]*(m-1)) for j in range(n)])
f = f * inverse_mod(a,R)
S = (f ^ (m - 1)).monomials()
M = [monomial * f for monomial in S]

L = [D[j] * (m - 1)for j in range(n)]
'''
extend 
S = 
M = 
L = [0 for _ in range(n)]
for j in range(n):
    Max = 0
    for exp in sum(S).exponents():
        if exp[j] > Max:
            Max = exp[j]
    L[j] = Max
R =  W * prod([bounds[j] ^ (L[j]) for j in range(n)])
'''
g = lambda monomial,f:monomial*f*prod([bounds[j]^(L[j]-monomial.exponents()[0][j]) for j in range(n)])
g_ = lambda monomial,f:monomial * R

ShiftPolys = Sequence([],f.parent())
for monomial in S:
    ShiftPolys.append(g(monomial,f))
for monomial in M:
    if monomial not in S:
        ShiftPolys.append(g_(monomial,f))

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