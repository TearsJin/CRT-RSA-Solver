from random import randint
from time import time
# generate 
n = 3
Bounds = [2^10 for _ in range(n)]
X = [randint(0,Bounds[i]) for i in range(n)]
c = X[0] * X[1] + X[0] * X[2] + X[1] * X[2]

# solve
PR.<x1,x2,x3> = PolynomialRing(ZZ)
f = x1 * x2 + x1 * x3 + x2 * x3 - c
m = 5
t = 2
Monomials = f.monomials()
Coefficients = f.coefficients()
a = Integer(Coefficients[-1])

for i in range(n):
    while gcd(Bounds[i],a) != 1:
        Bounds[i] -= 1

# define dj,lj
D = [0 for _ in range(n)]
for j in range(n):
    Max = 0
    for exp in f.exponents():
        if exp[j] > Max:
            Max = exp[j]
    D[j] = Max
L = [D[j] * (m - 1) for j in range(n)]
# define W,R
factors = [monomial(*Bounds) for monomial in Monomials]
W = max([i * j for i,j in zip(Coefficients,factors)])
R = W * prod([Bounds[j] ^ (D[j] * (m - 1)) for j in range(n)])
f = f * inverse_mod(a,R)

S = []
for monomial in (f ^ (m - 1)).monomials():
    for j in range(t):
        S.append(x1 ^ j * monomial)
M = []
for monomial in S:
    for monomial_ in (monomial * f).monomials():
        if monomial_ not in M:
            M.append(monomial_) 

g = lambda monomial,f:monomial * f * prod([Bounds[j]^(L[j]-monomial.exponents()[0][j]) for j in range(n)])
g_ = lambda monomial:monomial * R

ShiftPolys = Sequence([],f.parent())
for monomial in S:
    ShiftPolys.append(g(monomial,f))
for monomial in M:
    if monomial not in S:
        ShiftPolys.append(g_(monomial))

B, monomials = ShiftPolys.coefficient_matrix()
monomials = vector(monomials)
nn = len(ShiftPolys)
print(nn,len(monomials),len(S),len(M))
print("[+] Successfully defined polynomial.There are {} polynomials in total.".format(str(nn)))






# rescale the vectors
factors = [monomial(*Bounds) for monomial in monomials]
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
            root = tuple(root[var] for var in [x1,x2,x3])
            roots.append(root)
        for root in roots:
            if root != (1,0):
                print("[+] Successfully found the roots:",root)
                FindRootTime = time()-S
                print("[+] The time taken is",FindRootTime,"s")