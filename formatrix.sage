PR.<x,y,z,e,N> = PolynomialRing(ZZ)
m = 1
t = 2
s = 1
filename = 'BMMatrix'

f = x * (N - y) - N 

g = lambda i,j:e ^ (m - i) * x ^ j * z ^ s * f ^ i
h = lambda i,j:e ^ (m - i) * y ^ j * z ^ s * f ^ i

Ig = [(i,j) for i in range(m + 1) for j in range(m - i + 1)]
Ih = [(i,j) for i in range(m + 1) for j in range(1, t + 1)]

ShiftPolys = []

for I in Ig:
    ShiftPolys.append(g(*I))

for I in Ih:
    ShiftPolys.append(h(*I))


monomialOrder = []
for Poly in ShiftPolys:
    for x in Poly.monomials():
        if x(e=1,N=1) not in monomialOrder:
            monomialOrder.append(x(e=1,N=1))

with open(filename,'w') as f:
    nn = len(ShiftPolys)
    M = [[0 for _ in range(nn)] for _ in range(nn)]
    for i in range(nn):
        for monomial,b in zip(ShiftPolys[i].monomials(),ShiftPolys[i].coefficients()):
            M[i][monomialOrder.index(monomial(e=1,N=1))] = b * monomial

    f.write(str(Matrix(PR,M)) + '\n')

