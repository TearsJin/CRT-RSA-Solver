from time import time
m = 3

# Generate test data set
scale_factor = 1
p_size = int(512 * scale_factor)
l_size = int(170 * scale_factor)
d_size = int(220 * scale_factor)
M = 2^l_size

N, e = (113309068152935568893932620926523748561925916580328833761833709492836549360889407465934366065781314829245822754800772426000185646006813432371474526931829549338363892243205465109470135176086384087827228420777587664921704371825877672283627971299025154129511721048983330050794358456268791642813976689918889182193, 113309068152935568893932620926523748561925916580328833761833709492836549360889407465934366065781314829245822754800772426000185646006813432371474526931829549338363892243205465109470135176086384087827228420777587664921704371825877672283627971299025154129511721048983330050794358456268791642813976689918182742121)
dp_tilde, dq_tilde = (5013415024346389, 4333469053087705)

def nearest_below(x):
	c = floor(x)
	if not (c < x):
		c -= 1
	return c

############# PARAMETERS #################
delta = float((d_size - l_size) / int(N).bit_length())
beta = float((d_size) / int(N).bit_length())
alpha = float((int(e).bit_length()) / int(N).bit_length())
sigma = .859
tau = float(max(1/2, 1 - 2*beta))

X = int(e*2^(d_size - p_size))
Y = int(2^p_size)

print(f'beta: {float(beta)}')
print(f'(beta - delta)/beta: {float((beta - delta)/beta)}')
print(f'sigma: {sigma}')
print(f'delta: {delta}')
print(f'tau: {tau}')
P.<xp, xq, yp, yq, zp, zq> = PolynomialRing(ZZ, 6, order = 'lex')
index_map = {g:i for i,g in enumerate(P.gens())}

######################## HELPER DATA #########################
import itertools
M_sigma = [(a,b,c) for a, c, b in itertools.product(range(m + 1), range(m + 1), range(nearest_below(2*sigma*m) + 1))]

M_1 = [(a,b,c) for c in range(m + 1) for a in range(c + 1) for b in range(c - a + 1)]
M_2 = [(a,b,c) for c in range(m + 1) for a in range(c + 1, m + 1) for b in range(a - c)]
M_3 = [(a,b,c) for a in range(m + 1) for c in range(m + 1) for b in range(a + c + 1) if ((a,b,c) not in M_1 and (a,b,c) not in M_2 and (a + b + c) % 2 == 0)]
M_4 = [(a,b,c) for a in range(m + 1) for c in range(m + 1) for b in range(a + c + 1) if ((a,b,c) not in M_1 and (a,b,c) not in M_2 and (a,b,c) not in M_3)]
MM = [(a,b,c) for a in range(m + 1) for c in range(m + 1) for b in range(a + c + 1)]

M_tilde = [(a,b,c) for a, c, b in itertools.product(range(m + 1), range(m + 1), range(2*m + 1))]

# Sort M_sigma
M_tilde.sort(key = lambda p: xp^p[0]*yp^p[1]*zp^p[2])
M_sigma.sort(key = lambda p: xp^p[0]*yp^p[1]*zp^p[2])
MM.sort(key = lambda p: xp^p[0]*yp^p[1]*zp^p[2])

def E_f(a,b,c):
	if (a,b,c) in M_1:
		return 0
	elif (a,b,c) in M_2:
		return b
	elif (a,b,c) in M_3:
		return (a + b - c)//2
	elif (a,b,c) in M_4:
		return (a + b - c + 1)//2
	else:
		return a


def E_g(a, b, c):
	if (a,b,c) in M_1:
		return b
	elif (a,b,c) in M_2:
		return 0
	elif (a,b,c) in M_3:
		return (-a + b + c)//2
	elif (a,b,c) in M_4:
		return (-a + b + c - 1)//2
	else:
		return c

def E_h(a, b, c):
	if (a,b,c) in M_1:
		return a
	elif (a,b,c) in M_2:
		return c
	elif (a,b,c) in M_3:
		return (a - b + c)//2
	elif (a,b,c) in M_4:
		return (a - b + c - 1)//2
	else:
		return 0

def E_x(a, b, c):
	if (a,b,c) in M_1:
		return 0
	elif (a,b,c) in M_2:
		return a - b - c
	elif (a,b,c) in M_3:
		return 0
	elif (a,b,c) in M_4:
		return 0
	else:
		return 0

def E_z(a, b, c):
	if (a,b,c) in M_1:
		return -a - b + c
	elif (a,b,c) in M_2:
		return 0
	elif (a,b,c) in M_3:
		return 0
	elif (a,b,c) in M_4:
		return 1
	else:
		return 0

def trans(F):
	ypi = index_map[yp]
	yqi = index_map[yq]
	xpi = index_map[xp]
	zpi = index_map[zp]
	xqi = index_map[xq]
	zqi = index_map[zq]

	F = P(F)

	# Replace all instances of yp*yq by N
	new_dict = {}
	for t, v in F.dict().items():
		num = min(t[ypi], t[yqi])
		new_t = list(t)
		new_t[ypi] -= num
		new_t[yqi] -= num
		new_t = tuple(new_t)

		if new_t not in new_dict:
			new_dict[new_t] = 0

		new_dict[new_t] += int(N)^num * int(v)
	F = P(new_dict)

	# Step 2
	F_ = P(0)
	for t, v in F.dict().items():
		if t[ypi] != 0:
			F_ += P({t: v})
			continue

		new_t = list(t)
		xp_pow = t[xpi]
		zp_pow = t[zpi]
		new_t[xpi] = 0
		new_t[zpi] = 0
		new_t = tuple(new_t)

		mm = P({new_t: int(v)})
		F_ += mm*(xq + 1)^xp_pow*(zq - 1)^zp_pow
	F = F_

	# Step 3
	F_ = P(0)
	for t, v in F.dict().items():
		if t[ypi] == 0:
			F_ += P({t: int(v)})
			continue

		new_t = list(t)
		xq_pow = t[xqi]
		zq_pow = t[zqi]
		new_t[xqi] = 0
		new_t[zqi] = 0
		new_t = tuple(new_t)

		mm = P({new_t: v})
		F_ += mm*(xp - 1)^xq_pow*(zp + 1)^zq_pow
	F = F_

	return F

def lambda_abc(a, b, c):
	if b % 2 == 0:
		return xq^a*yq^(b//2)*zq^c
	else:
		return xp^a*yp^((b + 1)//2)*zp^c

# Highest monomial in the (zp, xp, yp) ordering
def highest_monomial(p):
	monomials = p.monomials()
	monomials.sort(key = lambda g: g(xp, xp, yp, yp, zp, zp))

	return monomials[-1]

# Rescale a poly to make the determinant of lattice smaller
def rescale_poly(poly):
	g = gcd(N - 1, e*M)

	monomials = poly.monomials()
	highest = highest_monomial(poly)

	d = int(poly[highest])
	t = next(iter(highest.dict().keys()))

	Xpow = t[index_map[xp]] + t[index_map[xq]] + t[index_map[zp]] + t[index_map[zq]]
	Ypow = t[index_map[yp]] + t[index_map[yq]]
	xy_pow = X^Xpow * Y^Ypow
	assert d % xy_pow == 0
	d //= xy_pow

	E4 = 0
	while d % N == 0:
		E4 += 1
		d //= N
	E5 = 0
	while d % (N - 1) == 0:
		E5 += 1
		d //= N - 1

	new_d = int(d) * int(g)^E5 * int(xy_pow)

	multiplier = pow(int(N), E4, (e*M)^(2*m))*pow(int((N - 1)//g), E5, (e*M)^(2*m))
	multiplier = pow(int(multiplier), -1, (e*M)^(2*m))

	p = P(0)
	for mm in monomials:
		if mm == highest:
			p += new_d * mm
		else:
			p += int(poly[mm]) * int(multiplier) * mm

	return p


################ PKE SHIFT POLYS #################
f_tilde = xp*yp - xq - e*dp_tilde
g_tilde = yp*zp - N*zq + e*dq_tilde*yp
h_tilde = N*xp*zq - xq*zp - e^2*dp_tilde*dq_tilde - e*dp_tilde*zp - e*dq_tilde*xq

def p_tilde(a, b, c):
	res = f_tilde^E_f(a, b, c)*g_tilde^E_g(a, b, c)*h_tilde^E_h(a, b, c)*xp^E_x(a, b, c)*zp^E_z(a, b, c)*(e*M)^(2*m - E_f(a, b, c) - E_g(a, b, c) - E_h(a, b, c))
	return P(res)

def pke_row(a, b, c):
	if (a,b,c) in MM:
		return trans(p_tilde(a, b, c)*yq^(b//2))(X*xp, X*xq, Y*yp, Y*yq, X*zp, X*zq)
	else:
		if b % 2 == 0:
			return trans(p_tilde(a, b, c)*yq^((a + c)//2)*yq^((b - a - c + 1)//2))(X*xp, X*xq, Y*yp, Y*yq, X*zp, X*zq)
		else:
			return trans(p_tilde(a, b, c)*yq^((a + c)//2)*yp^((b - a - c + 1)//2))(X*xp, X*xq, Y*yp, Y*yq, X*zp, X*zq)

################# TLP SHIFT POLYS #####################
f = M*(xp*yp - xq)
g = M*(yp*zp - N*zq)
h = M*(N*xp*zq - xq*zp)

def p(a, b, c):
	res = f^E_f(a, b, c)*g^E_g(a, b, c)*h^E_h(a, b, c)*xp^E_x(a, b, c)*zp^E_z(a, b, c)*(e*M)^(2*m - E_f(a, b, c) - E_g(a, b, c) - E_h(a, b, c))
	return P(res)

def p_ast(a, b, c, i, y):
	return trans(p(a, b, c)*yq^(b//2)*y^i)(X*xp, X*xq, Y*yp, Y*yq, X*zp, X*zq)

def tlp_row(a, b, c):
	return trans(p(a, b, c)*yq^(b//2))(X*xp, X*xq, Y*yp, Y*yq, X*zp, X*zq)

############ FETCH SHIFT POLYS ##################
# PKE polys
PKE_polys = []
for t in M_sigma:
	PKE_polys.append(pke_row(*t))

# TLP polys
TLP_polys = []
for t in MM:
	TLP_polys.append(tlp_row(*t))
	a, b, c = t
	if b == a + c:
		for i in range(1, floor(tau*b) - b//2 + 1):
			TLP_polys.append(p_ast(a, b, c, i, yq))
		for i in range(1, floor(tau*b) - (b + 1)//2 + 1):
			TLP_polys.append(p_ast(a, b, c, i, yp))

# Filter out un-needed tlp polys
for i in range(len(TLP_polys) - 1, -1, -1):
	p = TLP_polys[i]
	
	mm = highest_monomial(p)
	mm = next(iter(mm.dict().keys()))
	if mm[index_map[yq]] + mm[index_map[yp]] <= sigma*m:
		TLP_polys.pop(i)


############## COMBINE BOTH #################
shift_polys = PKE_polys + TLP_polys
# Rescale them for LLL
shift_polys = list(map(rescale_poly, shift_polys))

# Sort the monomials
monomials = set()
for p in shift_polys:
	monomials |= set(p.monomials())
monomials = list(monomials)
monomials.sort(key = lambda q: q(xp, xp, yp, yp, zp, zp))

# Build the coeff matrix
B = [[0 for _ in range(len(monomials))] for __ in range(len(shift_polys))]
for i, p in enumerate(shift_polys):
	for j, mm in enumerate(monomials):
		B[i][j] = p[mm]

B = matrix(ZZ, B)

############### APPLY COPPERSMITH TECHNIQUE TO FINISH ################
print(f'dim: {B.nrows()}')
S = time()
B = B.LLL()
print(f'Finished LLL {time() - S}')
B = B.change_ring(QQ)