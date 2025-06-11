#!/usr/bin/env sage
import sympy
x = polygen(QQ)

k.<a> = NumberField(x**10 + x**9 + x**8 + x**7 + x**6 + x**5 + x**4 + x**3 + x**2 + x +1, 'a')
print("Our cyclotomic field: ")
print(k)

U = k.S_unit_group(S=a)
print("Unit ring: ")
print(U)

print("Orders of units: ")
print([u.multiplicative_order() for u in U.gens()])

print("Generators: ")
print(U.gens_values())

B = []
for j in IntegerRange(0, 10):
	B += [[i**j % 11 for i in list(IntegerRange(1, 11))]]

var('b, c, a, d, e, f, g, h, i, j, k, l, w')
var('b2, c2, a2, d2, e2, f2, g2, h2, k2, l2')

B = Matrix(GF(11), B)

print("Coeff.'s of b, c, a, d, e, f, g, h, i, j, k, l(in a1, a2, ....): ")
print(B)
print("vector [a1, a2, a3, ..., a10]: ")

# Coefficients of 4 fundamental units
a = [-1, -1, -1, -1, -1, -1, -1, -1, -1, 0]
b = [-1, -1, -1, 0, -1, -1, -1, -1, -1, -1]
c = [0, 0, 0, 0, -1, 0, 0, -1, 0, 0]
d = [0, 0, 0, 0, -1, 0, -1, 0, 0, 0]

print("-------------------------------------------------")
print(B * vector(a)) # eta1
print(B * vector(b)) # eta2
print(B * vector(c)) # eta3
print(B * vector(d)) # eta4


print("-------------------------------------------------")
reset()
var('a, b, c, d, e, f, g, k, a0, b0, d0, e0, f0, g0, k0, eta1, eta2, eta3, eta4, m, n, o, p, alpha')
'''a = {eta1**m: m*2**(m - 1), eta2**n: n*2**(n+1), eta3**o: 5*o*3**(o - 1), eta4**p: p*2**(p+1), alpha: a0}
b = {eta1**m: 2**m, eta2**n: 2**n, eta3**o: 3**o, eta4**p: 2**p, alpha: b0}
d = {eta1**m: m * 2**(m - 1), eta2**n: n * 8 * 2**(n - 1), eta3**o: o * 9 * 3**(o - 1), eta4**p: p * 8 * 2**(p - 1), alpha: d0}
e = {eta1**m: 2**(m - 2)*(3*m**2 - m), eta2**n: 2**(n - 2) * (4*n**2 + 6*n), eta3**o: 3**(o - 2)*(9*o**2 - 2 * o), eta4**p: 2**(p - 2) * (4 * p**2 + 5*p), alpha: e0}
f = {eta1**m: 2**(m - 2) * (-m**2 + 3*m), eta2**n: 2**(n - 2) * (n**2 + 8*n), eta3**o: 3**(o - 2) * (1 - o) * o, eta4**p: 2**(p - 2) * (p - 1) * p, alpha: f0}
g = {eta1**m: 2**(m - 3) * (4*m**3 + 6*m**2 - 6*m), eta2**n: 2**(n - 3) * (3*n**3 + 10*n**2 + n), eta3**o: 3**(o - 3) * (5*o**3 + 7*o**2 + o), eta4**p: 2**p * (-p**3 + 2*p**2 + 2*p), alpha: g0}
k = {eta1**m: -5*2^(m - 4)*(m + 1)*(m - 1)*(m - 2)*m + 6*2^(m - 3)*(m - 1)*(m - 2)*m - 2*2^(m - 2)*(m - 1)*m + 2^(m - 1)*m, eta2**n: -1280*2^(n - 4)*(n + 1)*(n - 1)*(n - 2)*n + 1360*2^(n - 3)*(n - 1)*(n - 2)*n - 50*2^(n - 2)*(n - 1)*n + 4*2^(n - 1)*n, eta3**o: -3125*3^(o - 4)*(o + 1)*(o - 1)*(o - 2)*o + 2175*3^(o - 3)*(o - 1)*(o - 2)*o - 178*3^(o - 2)*(o - 1)*o + 5*3^(o - 1)*o, eta4**p: -1280*2^(p - 4)*(p + 1)*(p - 1)*(p - 2)*p + 1440*2^(p - 3)*(p - 1)*(p - 2)*p + 80*2^(p - 2)*(p - 1)*p + 4*2^(p - 1)*p, alpha: k0}
'''
# Fundamental units coefficients in basis: zeta, zeta^2, zeta^3, ...,  zeta^10


B = []
for j in IntegerRange(0, 10):
	B += [[i**j % 11 for i in list(IntegerRange(1, 11))]]

B = Matrix(GF(11), B)
print(B)
'''
e1 = [-1, -1, -1, -1, -1, -1, -1, -1, -1, 0]
e2 = [-1, -1, -1, 0, -1, -1, -1, -1, -1, -1]
e3 = [0, 0, 0, 0, -1, 0, 0, -1, 0, 0]
e4 = [0, 0, 0, 0, -1, 0, -1, 0, 0, 0]

e1 = [0, -1, -1, -1, -1, -1, -1, -1, -1, -1]
e2 = [-1, 0, -1, -1, -1, -1, -1, -1, -1, -1]
e3 = [0, 0, -1, -1, -1, -1, -1, -1, -1, -1]
e4 = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0]'''

# Coeff.'s are given in the basis \zeta, \zeta^2, ..., \zeta^10
e1 = [0, 0, 0, 0, 1, 1, 0, 0, 0, 0]
e2 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 1]
e3 = [0, -1, -1, -1, -1, -1, -1, -1, -1, 0]
e4 = [0, 0, 1, 0, 0, 0, 0, 1, 0, 0]

alpha_init = vector([-3, -4, 1, 0, 0, 0, 0, 0, 0, 0]) # b0 != 0, c0 == 0

print("Initial alpha for the computation of primary associate:")
alpha_init = B * alpha_init
print(alpha_init)

def fundam_powers_abc(B, e1, e2, e3, e4):
	var('a_ b_ e_ d_ g_ k_ f_ m n o p')

	a = b = d = e = f = g = k = {alpha: 0, eta1**m:0, eta2**n:0, eta3**o:0, eta4**p:0}
	a_power = n * a_ * b_ ** (n - 1)
	b_power = b_**n
	d_power = n * d_ * b_**(n - 1)
	e_power = n * e_ * b_**(n - 1) + 3 * n * (n - 1) * a_**2 * b_**(n - 2)
	f_power = n * f_ * b_** (n - 1) - n * (n - 1) * a_ * d_ * b_ ** (n - 2)
	g_power = n * g_ * b_**(n - 1) - n * (n - 1) * d_**2 * b_**(n - 2) + 4 * n * (n - 1) * a_ * e_ * b_**(n -2) + 4 * n * (n - 1) * (n - 2) * a_**3 * b_** (n - 3)
	k_power = 2 * n * (n - 1) * e_**2 * b_**(n - 2) + n* (n - 1) * d_ * f_ * b_**(n - 2) + 5 * n * (n - 1) * (n - 2) * a_ * d_**2 * b_**(n - 3) - 5 * n * (n - 1) * a_ * g_ * b_**(n - 2) + n * (n - 1) * (n - 2) * a_**2 * e_ * b_**(n - 3) - 5 * n * (n - 1) * (n -2) * (n - 3) * a_**4 * b_**(n - 4) + n * k_ * b_**(n - 1)

	# the evaluation of a, b, d, e, f, ... at the 4 fundamental units
	all_eval_e = vector(ZZ, B * vector(e1))
	all_eval_e2 = vector(ZZ, B * vector(e2)) # eta2
	all_eval_e3 = vector(ZZ, B * vector(e3)) # eta3
	all_eval_e4 = vector(ZZ, B * vector(e4)) # eta4
	'''print(all_eval_e)
	print(all_eval_e2) # eta2
	print(all_eval_e3) # eta3
	print(all_eval_e4) # eta4 '''

	# We will create the dictionaries for the evaluations at \alpha, eta1^m, eta2^n, eta3^o, eta4^p
	# for each a, b, d, ..., which will help us compute a(\alpha * eta1^m * eta2^n * eta3^o * eta4^p), ...

	# Here is the initialization
	lst = [('a', ({eta3^o: 0, eta4^p: 0, alpha: 0, eta1^m: 0, eta2^n: 0}, a_*b_^(n - 1)*n)), ('b', ({eta3^o: 0, eta4^p: 0, alpha: 0, eta1^m: 0, eta2^n: 0}, b_^n)), ('e', ({eta3^o: 0, eta4^p: 0, alpha: 0, eta1^m: 0, eta2^n: 0}, 3*a_^2*b_^(n - 2)*(n - 1)*n + b_^(n - 1)*e_*n)), ('d', ({eta3^o: 0, eta4^p: 0, alpha: 0, eta1^m: 0, eta2^n: 0}, b_^(n - 1)*d_*n)), ('g', ({eta3^o: 0, eta4^p: 0, alpha: 0, eta1^m: 0, eta2^n: 0}, 4*a_^3*b_^(n - 3)*(n + 1)*(n - 1)*n - b_^(n - 2)*d_^2*(n - 1)*n + 4*a_*b_^(n - 2)*e_*(n - 1)*n + b_^(n - 1)*g_*n)), ('f', ({eta3^o: 0, eta4^p: 0, alpha: 0, eta1^m: 0, eta2^n: 0}, -a_*b_^(n - 2)*d_*(n - 1)*n + b_^(n - 1)*f_*n)), ('k', ({eta3^o: 0, eta4^p: 0, alpha: 0, eta1^m: 0, eta2^n: 0}, 5*a_^4*b_^(n - 4)*(n + 1)*(n - 1)*(n - 2)*n + 5*a_*b_^(n - 3)*d_^2*(n - 1)*(n - 2)*n + a_^2*b_^(n - 3)*e_*(n - 1)*(n - 2)*n + 2*b_^(n - 2)*e_^2*(n - 1)*n + b_^(n - 2)*d_*f_*(n - 1)*n - 5*a_*b_^(n - 2)*g_*(n - 1)*n + b_^(n - 1)*k_*n))]

	# Looping through the list, at each item(corresponding to a dictionary),
	# we substitute the values computed in the vectors above into the formulas
	for i in range(len(lst)):
		# each item in the list is a 2-tuple of the form:
		# - the first element is the name of current function being considered
		# - the second element is a 2-tuple where first coord. is our desired dictionary,
		# the second coord. is our formula corresponding to the function;
		# e.g., the formula of a(y^n) = a_*b_^(n - 1)*n, where a_ = a(y), b(y)

		# The evaluation at \alpha is denoted _name_ + '0'
		lst[i][1][0][alpha] = lst[i][0]+"0"

		# The substitutions
		lst[i][1][0][eta1**m] = SR(lst[i][1][1].substitute(n == m, a_ == all_eval_e[2], b_ == all_eval_e[0], d_ == all_eval_e[3], e_ == all_eval_e[4], f_ == all_eval_e[5], g_ == all_eval_e[6], k_ == all_eval_e[8]))
		lst[i][1][0][eta2**n] = SR(lst[i][1][1].substitute(a_ == all_eval_e2[2], b_ == all_eval_e2[0], d_ == all_eval_e2[3], e_ == all_eval_e2[4], f_ == all_eval_e2[5], g_ == all_eval_e2[6], k_ == all_eval_e2[8]))
		lst[i][1][0][eta3**o] = SR(lst[i][1][1].substitute(n == o, a_ == all_eval_e3[2], b_ == all_eval_e3[0], d_ == all_eval_e3[3], e_ == all_eval_e3[4], f_ == all_eval_e3[5], g_ == all_eval_e3[6], k_ == all_eval_e3[8]))
		lst[i][1][0][eta4**p] = SR(lst[i][1][1].substitute(n == p, a_ == all_eval_e4[2], b_ == all_eval_e4[0], d_ == all_eval_e4[3], e_ == all_eval_e4[4], f_ == all_eval_e4[5], g_ == all_eval_e4[6], k_ == all_eval_e4[8]))
		#print(lst[i])
		#print(lst)

	return lst

def change_basis(old):
	# changing coordinates from basis \zeta, \zeta^2, ..., \zeta^10
	# to 1, \zeta, \zeta^2, ...\zeta^9
	new = [0] + [old[i] for i in range(0, 9)]
	for i in range(10):
		new[i] = (new[i] - old[9])
	return new

def change_basis_to_1_10(old):
	new = [old[1]- old[0], old[2]- old[0], old[3]- old[0], old[4]- old[0], old[5]- old[0], old[6]- old[0], old[7]- old[0], old[8]- old[0], old[9] - old[0], - old[0]]

	#new = matrix(QQ, new)
	new = vector(new)
	return new


# eps = +(-) \zeta^(alpha) * \eta1^(beta1) * \eta2^(beta2) * \eta3^(beta3) * \eta4^(beta4) * (1 - \zeta)^10
# find a unit epsilon by brute forcing beta, alpha and the sign
# epsilon satisfies that: 11 = \eps \times \omg^10
def find_ep_unit(e1, e2, e3, e4):
	var('zeta')
	k = CyclotomicField(11)
	unit_r = k.unit_group() # Z*[\zeta]
	ring = k.maximal_order() # Z[\zeta]
	gen_d = unit_r.gens_dict() # dictationary of generating values for unit_r,
							# in order: {-\zeta, [4 fundamental units]}

	eta1_p = unit_r([0] + e1)
	eta1_p = ring(eta1_p)
	print(eta1_p)
	eta2_p = unit_r([0] + e2)
	eta2_p = ring(eta2_p)
	print(eta2_p)
	eta3_p = unit_r([0] + e3)
	eta3_p = ring(eta3_p)
	print(eta3_p)
	eta4_p = unit_r([0] + e4)
	eta4_p = ring(eta4_p)
	print(eta4_p)
	dummy = unit_r([0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0])
	dummy = ring(dummy)
	print(dummy)
	omg = ring([1, -1])
	print(omg)
	# Check correctness: 11 * \omg^(-10) \in unit_r
	t = ring(11) * ring([1, -1])^(-10)
	assert t.is_unit()

	zet = ring([0, 1])


	alpha = beta1 = beta2 = beta3 = beta4 = dummyp = 1
	sign = [1, -1]

	for s in sign:
		for alpha in range(11):
			for beta1 in range(11):
				for beta2 in range(11):
					for beta3 in range(11):
						for beta4 in range(11):
							for dummyp in range(1, 11):
								eps = s * zet^alpha * eta1_p^(beta1) * eta2_p^(beta2) * eta3_p^(beta3) * eta4_p^(beta4) * dummy^(dummyp) * omg^10

								if ring(11) == ring(eps):
									print("hehe!!!!!")
									print(eps)
									print(s, alpha, beta1, beta2, beta3, beta4, dummyp)
									return (s, alpha, beta1, beta2, beta3, beta4, dummyp)
								else:
									print(eps)
									print("ooops!!!!!")

print("-------------------------------------")
print("-------------------------------------")
#find_ep_unit(e1, e2, e3, e4) --> (-1, 6, 4, 2, 4, 0, 2) for sign and the powers of zeta, \eta1, \eta2, \eta3, \eta4 and \dummy = \zeta^4 + \zeta^7
# Double-check
k_f = CyclotomicField(11)
unit_r = k_f.unit_group() # Z*[\zeta]
ring = k_f.maximal_order() # Z[\zeta]
eta1_p = unit_r([0] + e1)
eta1_p = ring(eta1_p)
print(eta1_p)
eta2_p = unit_r([0] + e2)
eta2_p = ring(eta2_p)
print(eta2_p)
eta3_p = unit_r([0] + e3)
eta3_p = ring(eta3_p)
print(eta3_p)
eta4_p = unit_r([0] + e4)
eta4_p = ring(eta4_p)
print(eta4_p)
dummy = unit_r([0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0])
dummy = ring(dummy)
print(dummy)
omg = ring([1, -1])
print(omg)
# Check correctness: 11 * \omg^(-10) \in unit_r
t = ring(11) * ring([1, -1])^(-10)
assert t.is_unit()

zet = ring([0, 1])

eps = (-1) * zet^6 * eta1_p^4 * eta2_p^2 * eta3_p^4 * eta4_p^0 * dummy^2 * omg^10
print("Checking... \\eps == 11")
assert eps == ring(11)

print("-------------------------------------")
print("-------------------------------------")

'''dict_a = {}
dict_b = {}
dict_d = {}
dict_e = {}
dict_f = {}
dict_g = {}
dict_k = {}'''


def look_up(dict_, fundam_power):
	if bool (fundam_power == eta1**m):
		res = dict_[eta1**m]
	elif bool(fundam_power == eta2**n):
		res = dict_[eta2**n]
	elif bool(fundam_power == eta3**o):
		res = dict_[eta3**o]
	elif bool(fundam_power == eta4**p):
		res = dict_[eta4**p]
	else:
		res = dict_[alpha]

	return res

def compute_a(l, alph, beta):
	if bool(len(alph) == 1):
		a_1 = look_up(l['dict_a'], alph[0])
	else:
		mid = len(alph) // 2
		a_1 = compute_a(l, alph[:mid], alph[mid:])

	if bool(len(beta) == 1):
		a_2 = look_up(l['dict_a'], beta[0])
	else:
		mid = len(beta) // 2
		a_2 = compute_a(l, beta[:mid], beta[mid:])

	if bool(len(alph) == 1):
		b_1 = look_up(l['dict_b'], alph[0])
	else:
		mid = len(alph) // 2
		b_1 = compute_b(l, alph[:mid], alph[mid:])

	if bool(len(beta) == 1):
		b_2 = look_up(l['dict_b'], beta[0])
	else:
		mid = len(beta) // 2
		b_2 = compute_b(l, beta[:mid], beta[mid:])

	# Since we're computing using symbolic expression
	# we need to define a 'symbolic' function for the returning result,
	# otherwise some term might be treated as integers and leading to error
	var('a1, a2, b1, b2')
	e = a1 * b2 + a2 * b1
	h(a1, a2, b1, b2) = e

	return h(a_1, a_2, b_1, b_2)

def compute_b(l, alph, beta):
	if bool(len(alph) == 1):
		b_1 = look_up(l['dict_b'], alph[0])
	else:
		mid = len(alph) // 2
		b_1 = compute_b(l, alph[:mid], alph[mid:])

	if bool(len(beta) == 1):
		b_2 = look_up(l['dict_b'], beta[0])
	else:
		mid = len(beta) // 2
		b_2 = compute_b(l, beta[:mid], beta[mid:])
	var('b1, b2')
	e = b1 * b2
	h(b1, b2) = e
	return h(b_1, b_2)

def compute_d(l, alph, beta):
	if bool(len(alph) == 1):
		d_1 = look_up(l['dict_d'], alph[0])
	else:
		mid = len(alph) // 2
		d_1 = compute_d(l, alph[:mid], alph[mid:])

	if bool(len(beta) == 1):
		d_2 = look_up(l['dict_d'], beta[0])
	else:
		mid = len(beta) // 2
		d_2 = compute_d(l, beta[:mid], beta[mid:])

	if bool(len(alph) == 1):
		b_1 = look_up(l['dict_b'], alph[0])
	else:
		mid = len(alph) // 2
		b_1 = compute_b(l, alph[:mid], alph[mid:])

	if bool(len(beta) == 1):
		b_2 = look_up(l['dict_b'], beta[0])
	else:
		mid = len(beta) // 2
		b_2 = compute_b(l, beta[:mid], beta[mid:])

	var('d1, d2, b1, b2')
	e = d1 * b2 + d2 * b1
	h(d1, d2, b1, b2) = e

	return h(d_1, d_2, b_1, b_2) # d_1 * b_2 + d_2 * b_1

def compute_e(l, alph, beta):
	if bool(len(alph) == 1):
		a_1 = look_up(l['dict_a'], alph[0])
		e_1 = look_up(l['dict_e'], alph[0])
		b_1 = look_up(l['dict_b'], alph[0])
	else:
		mid = len(alph) // 2
		a_1 = compute_a(l, alph[:mid], alph[mid:])
		e_1 = compute_e(l, alph[:mid], alph[mid:])
		b_1 = compute_b(l, alph[:mid], alph[mid:])

	if bool(len(beta) == 1):
		a_2 = look_up(l['dict_a'], beta[0])
		e_2 = look_up(l['dict_e'], beta[0])
		b_2 = look_up(l['dict_b'], beta[0])
	else:
		mid = len(beta) // 2
		a_2 = compute_a(l, beta[:mid], beta[mid:])
		e_2 = compute_e(l, beta[:mid], beta[mid:])
		b_2 = compute_b(l, beta[:mid], beta[mid:])

	var('a1, a2, b1, b2, e1, e2')
	e_ = 6 * a1 * a2 + b1 * e2 + b2 * e1

	h(a1, a2, b1, b2, e1, e2) = e_

	return h(a_1, a_2, b_1, b_2, e_1, e_2)

def compute_f(l, alph, beta):
	if bool(len(alph) == 1):
		a_1 = look_up(l['dict_a'], alph[0])
		f_1 = look_up(l['dict_f'], alph[0])
		b_1 = look_up(l['dict_b'], alph[0])
		d_1 = look_up(l['dict_d'], alph[0])
	else:
		mid = len(alph) // 2
		a_1 = compute_a(l, alph[:mid], alph[mid:])
		f_1 = compute_f(l, alph[:mid], alph[mid:])
		b_1 = compute_b(l, alph[:mid], alph[mid:])
		d_1 = compute_d(l, alph[:mid], alph[mid:])

	if bool(len(beta) == 1):
		a_2 = look_up(l['dict_a'], beta[0])
		f_2 = look_up(l['dict_f'], beta[0])
		b_2 = look_up(l['dict_b'], beta[0])
		d_2 = look_up(l['dict_d'], beta[0])
	else:
		mid = len(beta) // 2
		a_2 = compute_a(l, beta[:mid], beta[mid:])
		f_2 = compute_f(l, beta[:mid], beta[mid:])
		b_2 = compute_b(l, beta[:mid], beta[mid:])
		d_2 = compute_d(l, beta[:mid], beta[mid:])

	var('a1, a2, b1, b2, d1, d2, f1, f2')
	e = -a2 * d1 - a1 * d2 + b1 * f2 + b2 * f1
	h(a1, a2, b1, b2, d1, d2, f1, f2) = e

	return h(a_1, a_2, b_1, b_2, d_1, d_2, f_1, f_2)

def compute_g(l, alph, beta):
	if bool(len(alph) == 1):
		a_1 = look_up(l['dict_a'], alph[0])
		e_1 = look_up(l['dict_e'], alph[0])
		b_1 = look_up(l['dict_b'], alph[0])
		d_1 = look_up(l['dict_d'], alph[0])
		g_1 = look_up(l['dict_g'], alph[0])
	else:
		mid = len(alph) // 2
		a_1 = compute_a(l, alph[:mid], alph[mid:])
		e_1 = compute_e(l, alph[:mid], alph[mid:])
		b_1 = compute_b(l, alph[:mid], alph[mid:])
		d_1 = compute_d(l, alph[:mid], alph[mid:])
		g_1 = compute_g(l, alph[:mid], alph[mid:])

	if bool(len(beta) == 1):
		a_2 = look_up(l['dict_a'], beta[0])
		e_2 = look_up(l['dict_e'], beta[0])
		b_2 = look_up(l['dict_b'], beta[0])
		d_2 = look_up(l['dict_d'], beta[0])
		g_2 = look_up(l['dict_g'], beta[0])

	else:
		mid = len(beta) // 2
		a_2 = compute_a(l, beta[:mid], beta[mid:])
		e_2 = compute_e(l, beta[:mid], beta[mid:])
		b_2 = compute_b(l, beta[:mid], beta[mid:])
		d_2 = compute_d(l, beta[:mid], beta[mid:])
		g_2 = compute_g(l, beta[:mid], beta[mid:])

	var('a1, a2, b1, b2, d1, d2, e1, e2, g1, g2')
	e_ = 9*d1*d2 + 4*a1*e2 + 4*a2*e1 + b1*g2 + b2*g1
	h(a1, a2, b1, b2, d1, d2, e1, e2, g1, g2) = e_

	return h(a_1, a_2, b_1, b_2, d_1, d_2, e_1, e_2, g_1, g_2)
	#9*d_1*d_2 + 4*a_1*e_2 + 4*a_2*e_1 + b_1*g_2 + b_2*g_1

def compute_k(l, alph, beta):
	if bool(len(alph) == 1):
		a_1 = look_up(l['dict_a'], alph[0])
		e_1 = look_up(l['dict_e'], alph[0])
		b_1 = look_up(l['dict_b'], alph[0])
		d_1 = look_up(l['dict_d'], alph[0])
		g_1 = look_up(l['dict_g'], alph[0])
		f_1 = look_up(l['dict_f'], alph[0])
		k_1 = look_up(l['dict_k'], alph[0])
	else:
		mid = len(alph) // 2
		a_1 = compute_a(l, alph[:mid], alph[mid:])
		e_1 = compute_e(l, alph[:mid], alph[mid:])
		b_1 = compute_b(l, alph[:mid], alph[mid:])
		d_1 = compute_d(l, alph[:mid], alph[mid:])
		g_1 = compute_g(l, alph[:mid], alph[mid:])
		f_1 = compute_f(l, alph[:mid], alph[mid:])
		k_1 = compute_k(l, alph[:mid], alph[mid:])


	if bool(len(beta) == 1):
		a_2 = look_up(l['dict_a'], beta[0])
		e_2 = look_up(l['dict_e'], beta[0])
		b_2 = look_up(l['dict_b'], beta[0])
		d_2 = look_up(l['dict_d'], beta[0])
		g_2 = look_up(l['dict_g'], beta[0])
		f_2 = look_up(l['dict_f'], beta[0])
		k_2 = look_up(l['dict_k'], beta[0])

	else:
		mid = len(beta) // 2
		a_2 = compute_a(l, beta[:mid], beta[mid:])
		e_2 = compute_e(l, beta[:mid], beta[mid:])
		b_2 = compute_b(l, beta[:mid], beta[mid:])
		d_2 = compute_d(l, beta[:mid], beta[mid:])
		g_2 = compute_g(l, beta[:mid], beta[mid:])
		f_2 = compute_f(l, beta[:mid], beta[mid:])
		k_2 = compute_k(l, beta[:mid], beta[mid:])

	var('a1, a2, b1, b2, d1, d2, f1, f2, e1, e2, g1, g2, k1, k2')
	e_ = 4*e1*e2 + d2*f1 + d1*f2 - 5*a2*g1 - 5*a1*g2 + b1*k2 + b2*k1
	h(a1, a2, b1, b2, d1, d2, e1, e2, f1, f2, g1, g2, k1, k2) = e_

	return h(a_1, a_2, b_1, b_2, d_1, d_2, e_1, e_2, f_1, f_2, g_1, g_2, k_1, k_2)
	# 4*e_1*e_2 + d_2*f_1 + d_1*f_2 - 5*a_2*g_1 - 5*a_1*g_2 + b_1*k_2 + b_2*k_1

print("Compiled successfully!")

def check_prim(alpha):
	# input alpha is given as an element in the ring
	print("[+] Entering function check_prim() ...")
	alpha_conj = alpha.conjugate()

	omg = ring([1, -1]) # we use the same maximal_order in Q(\zeta11)

	v_alpha = vector(matrix(GF(11), ring.coordinates(alpha))) # this is coord.'s in basis 1, \zeta, ..., \zeta^9
	# convert to coordinates in \zeta, \zeta^2, ..., \zeta^10
	v_alpha = vector([v_alpha[1]- v_alpha[0], v_alpha[2]- v_alpha[0], v_alpha[3]- v_alpha[0], v_alpha[4]- v_alpha[0], v_alpha[5]- v_alpha[0], v_alpha[6]- v_alpha[0], v_alpha[7]- v_alpha[0], v_alpha[8]- v_alpha[0], v_alpha[9] - v_alpha[0], - v_alpha[0]])
	#we compute b, c, a, ..., l for alpha
	sum_pwrs_alpha = B * v_alpha
	print("[+] We compute b, c, a, ..., l of input: {}".format(sum_pwrs_alpha))

	v_prod_alpha_connj = vector(matrix(GF(11), ring.coordinates(alpha*alpha_conj)))
	# similar treatment
	v_prod_alpha_connj = vector([v_prod_alpha_connj[1]- v_prod_alpha_connj[0], v_prod_alpha_connj[2]- v_prod_alpha_connj[0], v_prod_alpha_connj[3]- v_prod_alpha_connj[0], v_prod_alpha_connj[4]- v_prod_alpha_connj[0], v_prod_alpha_connj[5]- v_prod_alpha_connj[0], v_prod_alpha_connj[6]- v_prod_alpha_connj[0], v_prod_alpha_connj[7]- v_prod_alpha_connj[0], v_prod_alpha_connj[8]- v_prod_alpha_connj[0], v_prod_alpha_connj[9] - v_prod_alpha_connj[0], - v_prod_alpha_connj[0]])
	#we compute b, c, a, ..., l for alpha*alpha_conj
	sum_pwrs_prod_conj_alpha = B * v_prod_alpha_connj
	print("[+] We compute b, c, a, ..., l for input*input_conj: {}".format(sum_pwrs_prod_conj_alpha))

	print("[+] Checking \\input != 0 (mod \\omg) ...")
	b1 = vector(sum_pwrs_alpha)[0] % 11 != 0 # \alpha != 0 (mod \omg) <=> b(\alpha) != 0 (mod 11)
	print("[+] Checking \\input == b(\\alpha) (mod \\omg^2) ...")
	b2 = vector(sum_pwrs_alpha)[1] % 11 == 0 # \alpha == b(\alpha) (mod \omg^2) <=> c(\alpha) == 0 (mod 11)
	print("[+] Checking \\input*\\input_conj == b(\\input)^2 (mod 11) ...")
	b3 = vector(sum_pwrs_prod_conj_alpha)[1] % 11 == 0 # aux. check to make sure that \alpha*\alpha_conj == b(\alpha*\alpha_conj) (mod \omg^2)
	b4 = vector(sum_pwrs_prod_conj_alpha)[0] % 11 == vector(sum_pwrs_alpha)[0]^2 % 11
	return (b1 and b2 and b3 and b4)


def compute_prim_assoc(alpha_input, e1, e2, e3, e4):
	# init the field, its maximal order and fundamental units
	k_f = CyclotomicField(11)
	unit_r = k_f.unit_group() # Z*[\zeta]
	ring = k_f.maximal_order() # Z[\zeta]
	zet = ring([0, 1])

	'''eta1_p = unit_r([0] + e1)
	eta1_p = ring(eta1_p)'''
	eta1_p = zet^5 + zet^6
	assert eta1_p.norm() == 1
	'''eta2_p = unit_r([0] + e2)
	eta2_p = ring(eta2_p)'''
	eta2_p = zet^10 + zet
	assert eta2_p.norm() == 1
	'''eta3_p = unit_r([0] + e3)
	eta3_p = ring(eta3_p)'''
	eta3_p = 1 + zet + zet^10
	assert eta3_p.norm() == 1
	'''eta4_p = unit_r([0] + e4)
	eta4_p = ring(eta4_p)'''
	eta4_p = zet^3 + zet^8
	assert eta4_p.norm() == 1

	omg = ring([1, -1])

	# Check correctness: 11 * \omg^(-10) \in unit_r
	t = ring(11) * ring([1, -1])^(-10)
	assert t.is_unit()

	# compute b, c, a, ... of \alpha
	# we need to treat the type of alpha, since we'll input the result from change_basis_to_1_10
	alpha_input = matrix(ZZ, alpha_input)
	alpha_input = vector(alpha_input)

	alpha_init = vector(alpha_input)
	alpha_init = B * alpha_init

	'''print "[+] alpha: {}".format(alpha_input)
	print "[+] Before nomalization: {}".format(alpha_init)
	print "[+] Checking b(\\alpha) != 0 mod 11 ..."'''
	assert alpha_init[0] % 11 != 0
	# compute the appropriate power for factor \zeta
	k_power = (-1 * 1/alpha_init[0] * alpha_init[1]) % 11
	#print "[+] k_power: {}".format(k_power)

	# do not change the semantic of 'alpha_init', it is used later in the computation
	alpha_in_ring_elem = ring(change_basis(alpha_input)) # b0 != 0
	'''if check_prim(alpha_in_ring_elem):
		return [0, 0, 0, 0, 0, alpha_in_ring_elem]'''
	alpha_normlzed = zet^k_power * alpha_in_ring_elem

	# get the coordinates of normalized alpha, then put them in basis \zeta, \zeta^2, .., \zeta^10
	v_alpha_normlzed = vector(matrix(GF(11), ring.coordinates(alpha_normlzed)))
	v_alpha_normlzed = vector([v_alpha_normlzed[1]- v_alpha_normlzed[0], v_alpha_normlzed[2]- v_alpha_normlzed[0], v_alpha_normlzed[3]- v_alpha_normlzed[0], v_alpha_normlzed[4]- v_alpha_normlzed[0], v_alpha_normlzed[5]- v_alpha_normlzed[0], v_alpha_normlzed[6]- v_alpha_normlzed[0], v_alpha_normlzed[7]- v_alpha_normlzed[0], v_alpha_normlzed[8]- v_alpha_normlzed[0], v_alpha_normlzed[9] - v_alpha_normlzed[0], - v_alpha_normlzed[0]])

	alpha_init = B * v_alpha_normlzed
	assert alpha_init[1] % 11 == 0


	# compute the evaluations of b, c, a, ... at product of powers of fundamental units
	dicts = fundam_powers_abc(B, e1, e2, e3, e4)
	dict_a = [item for item in dicts if item[0] == 'a'][0][1][0]
	dict_a[alpha] = SR(dict_a[alpha])

	dict_b = [item for item in dicts if item[0] == 'b'][0][1][0]
	dict_b[alpha] = SR(dict_b[alpha])

	dict_d = [item for item in dicts if item[0] == 'd'][0][1][0]
	dict_d[alpha] = SR(dict_d[alpha])

	dict_e = [item for item in dicts if item[0] == 'e'][0][1][0]
	dict_e[alpha] = SR(dict_e[alpha])


	dict_f = [item for item in dicts if item[0] == 'f'][0][1][0]
	dict_f[alpha] = SR(dict_f[alpha])

	dict_g = [item for item in dicts if item[0] == 'g'][0][1][0]
	dict_g[alpha] = SR(dict_g[alpha])

	dict_k = [item for item in dicts if item[0] == 'k'][0][1][0]
	dict_k[alpha] = SR(dict_k[alpha])

	# normalize \alpha to \alpha' s.t c(\alpha') == 0 mod 11, b(\alpha') = b(\alpha)
	alph = [alpha, eta1**m, eta2**n]
	beta = [eta3**o, eta4**p]

	l_dict = {'dict_a':dict_a, 'dict_b':dict_b, 'dict_d':dict_d, 'dict_e':dict_e, 'dict_f':dict_f, 'dict_g':dict_g, 'dict_k':dict_k}

	a_fundams = simplify(compute_a(l_dict, alph, beta).expand()).combine()
	b_fundams = simplify(compute_b(l_dict, alph, beta).expand()).combine()
	d_fundams = simplify(compute_d(l_dict, alph, beta).expand()).combine()
	e_fundams = simplify(compute_e(l_dict, alph, beta).expand()).combine()
	f_fundams = simplify(compute_f(l_dict, alph, beta).expand()).combine()
	g_fundams = simplify(compute_g(l_dict, alph, beta).expand()).combine()
	k_fundams = simplify(compute_k(l_dict, alph, beta).expand()).combine()

	# derive 4 equations for primary conditions
	eq1 = a_fundams
	eq2 = e_fundams
	eq3 = d_fundams**2 + b_fundams * g_fundams
	eq4 = d_fundams * f_fundams - b_fundams * k_fundams

	# b, c, a,... for 4 fundamental units
	vect_pwrs_e1 = B * vector(e1)
	vect_pwrs_e2 = B * vector(e2)
	vect_pwrs_e3 = B * vector(e3)
	vect_pwrs_e4 = B * vector(e4)

	# Tedious tranformation to vectors in QQ, for symbolic computation
	vect_pwrs_e1 = matrix(QQ, vect_pwrs_e1)
	vect_pwrs_e1 = vector(vect_pwrs_e1)
	vect_pwrs_e2 = matrix(QQ, vect_pwrs_e2)
	vect_pwrs_e2 = vector(vect_pwrs_e2)
	vect_pwrs_e3 = matrix(QQ, vect_pwrs_e3)
	vect_pwrs_e3 = vector(vect_pwrs_e3)
	vect_pwrs_e4 = matrix(QQ, vect_pwrs_e4)
	vect_pwrs_e4 = vector(vect_pwrs_e4)

	# note that only b() of each fundam. unit gives the basis for exponents
	# Remove the exponents(which must be > 0 in the field GF(11) since 11 is a prime)
	expo_factors = vect_pwrs_e1[0]^m * vect_pwrs_e2[0]^n * vect_pwrs_e3[0]^o * vect_pwrs_e4[0]^p
	eq1 = simplify((eq1 / (expo_factors)).expand())
	eq2 = simplify((eq2 / (expo_factors)).expand())
	eq3 = simplify((eq3 / (expo_factors)**2).expand()) # eq3 is of degree 2 in b, c, a ...
	eq4 = simplify((eq4 / (expo_factors)**2).expand()) # eq4 is of degree 2 in b, c, a ...

	m_s = solve(SR(eq1) == 0, [m])[0]
	m_s = m_s.right()

	eq2 = SR(eq2).substitute(m = m_s)
	n_s = solve(SR(eq2) == 0, [n])[0]
	n_s = n_s.right()

	'''print "[+] m from eq1: {}".format(m_s)
	print "[+] n from eq2 and m's subs.: {}".format(n_s)'''

	eq3 = (SR(eq3).substitute(m = m_s)).substitute(n = n_s)
	eq4 = (SR(eq4).substitute(m = m_s)).substitute(n = n_s)
	eq3 = simplify(eq3.expand())
	eq4 = simplify(eq4.expand())

	# mult to convert SR to poly(GF(11)), note that these are equations and b0 != 0
	e3_denom = eq3.denominator()
	e4_denom = eq4.denominator()
	eq3 = (eq3 * e3_denom).polynomial(GF(11))
	eq4 = (eq4 * e4_denom).polynomial(GF(11))

	'''print "[+] eq3 and eq4 only in variables o, p: "
	print "[-] eq3: {}".format(eq3)
	print "[-] eq4: {}".format(eq4)'''

	o_prim = 0
	p_prim = 0
	#print "[+] Normalized alpha : {}".format(alpha_init)
	for o_t in range(11):
		for p_t in range(11):
			e3_eval = mod(eq3(b0 = alpha_init[0], a0 = alpha_init[2], d0=alpha_init[3], e0=alpha_init[4], g0=alpha_init[6], o=o_t, p=p_t), 11)

			e4_eval = mod(eq4(b0 = alpha_init[0], a0 = alpha_init[2], d0=alpha_init[3], e0=alpha_init[4], f0=alpha_init[5], g0=alpha_init[6], k0=alpha_init[8], o=o_t, p=p_t), 11)
			if (e3_eval == mod(0, 11)) and (e4_eval == mod(0, 11)):
				o_prim = o_t
				p_prim = p_t
				#print "[+] Found o, p!!"
	alpha_init = matrix(QQ, alpha_init)
	alpha_init = vector(alpha_init) # convert back to vector over QQ, for the substitution

	n_prim = n_s.substitute(b0 = alpha_init[0], a0 = alpha_init[2], d0=alpha_init[3], e0=alpha_init[4], f0=alpha_init[5], g0=alpha_init[6], k0=alpha_init[8], o=o_prim, p=p_prim)

	m_prim = m_s.substitute(b0 = alpha_init[0], a0 = alpha_init[2], d0=alpha_init[3], e0=alpha_init[4], f0=alpha_init[5], g0=alpha_init[6], k0=alpha_init[8], o=o_prim, p=p_prim, n=n_prim)

	# after susbtitution, they're still rationals
	# need to convert back to modulo 11
	n_prim = Rational(n_prim) % 11
	m_prim = Rational(m_prim) % 11

	alpha_prim_assoc = eta1_p**m_prim * eta2_p**n_prim * eta3_p**o_prim * eta4_p**p_prim * alpha_normlzed

	print("[+] Checking correctness...")
	assert check_prim(alpha_prim_assoc)
	print("[+] Checking that \\alpha and \\alph_prim_assoc are associates ..")
	eps = alpha_in_ring_elem / alpha_prim_assoc
	assert eps.is_unit()

	print("[+] Done")
	return [k_power, m_prim, n_prim, o_prim, p_prim, alpha_prim_assoc]



# --------------------- Alternative implementation ---------------------
def compute_prim_assoc_alt(alpha_input, e1, e2, e3, e4):
	# init the field, its maximal order and fundamental units
	k_f = CyclotomicField(11)
	unit_r = k_f.unit_group() # Z*[\zeta]
	ring = k_f.maximal_order() # Z[\zeta]
	zet = ring([0, 1])

	'''eta1_p = unit_r([0] + e1)
	eta1_p = ring(eta1_p)'''
	eta1_p = zet^5 + zet^6
	assert eta1_p.norm() == 1
	'''eta2_p = unit_r([0] + e2)
	eta2_p = ring(eta2_p)'''
	eta2_p = zet^10 + zet
	assert eta2_p.norm() == 1
	'''eta3_p = unit_r([0] + e3)
	eta3_p = ring(eta3_p)'''
	eta3_p = 1 + zet + zet^10
	assert eta3_p.norm() == 1
	'''eta4_p = unit_r([0] + e4)
	eta4_p = ring(eta4_p)'''
	eta4_p = zet^3 + zet^8
	assert eta4_p.norm() == 1

	omg = ring([1, -1])

	# Check correctness: 11 * \omg^(-10) \in unit_r
	t = ring(11) * ring([1, -1])^(-10)
	assert t.is_unit()

	# compute b, c, a, ... of \alpha
	# we need to treat the type of alpha, since we'll input the result from change_basis_to_1_10
	alpha_input = matrix(ZZ, alpha_input)
	alpha_input = vector(alpha_input)

	alpha_init = vector(alpha_input)
	alpha_init = B * alpha_init

	assert alpha_init[0] % 11 != 0
	e0 = -alpha_init[1] / alpha_init[0] % 11

	# do not change the semantic of 'alpha_init', it is used later in the computation
	alpha_in_ring_elem = ring(change_basis(alpha_input)) # b0 != 0
	'''if check_prim(alpha_in_ring_elem):
		return [0, 0, 0, 0, 0, alpha_in_ring_elem]'''
	alpha_normlzed = zet^e0 * alpha_in_ring_elem

	# deterministic computation
	# get the coordinates of normalized alpha, then put them in basis \zeta, \zeta^2, .., \zeta^10
	v_alpha_normlzed = vector(matrix(GF(11), ring.coordinates(alpha_normlzed)))
	v_alpha_normlzed = vector([v_alpha_normlzed[1]- v_alpha_normlzed[0], v_alpha_normlzed[2]- v_alpha_normlzed[0], v_alpha_normlzed[3]- v_alpha_normlzed[0], v_alpha_normlzed[4]- v_alpha_normlzed[0], v_alpha_normlzed[5]- v_alpha_normlzed[0], v_alpha_normlzed[6]- v_alpha_normlzed[0], v_alpha_normlzed[7]- v_alpha_normlzed[0], v_alpha_normlzed[8]- v_alpha_normlzed[0], v_alpha_normlzed[9] - v_alpha_normlzed[0], - v_alpha_normlzed[0]])

	alpha_init = B * v_alpha_normlzed
	assert alpha_init[1] % 11 == 0
	w1 = -alpha_init[2] / alpha_init[0] % 11
	w2 = (-alpha_init[4] / alpha_init[0] + 3*(alpha_init[2] / alpha_init[0])^2) % 11
	w3 = (-alpha_init[6] / alpha_init[0] + 4*alpha_init[4]*alpha_init[2] / alpha_init[0]^2 - (alpha_init[3] / alpha_init[0])^2 + 3*(alpha_init[2]/alpha_init[0])^3) % 11
	w4 = (-alpha_init[8] / alpha_init[0] - 5*alpha_init[6]*alpha_init[2]/alpha_init[0]^2 + alpha_init[5]*alpha_init[3]/alpha_init[0]^2 + 2*(alpha_init[4]/alpha_init[0])^2 - 2*(alpha_init[4]/alpha_init[0])*(alpha_init[2]/alpha_init[0])^2 + (alpha_init[3]/alpha_init[0])^2*alpha_init[2]/alpha_init[0] + 3*(alpha_init[2]/alpha_init[0])^4) % 11

	e1 = (-2*w1 -4*w2+5*w3 -3*w4) % 11
	e2 = (4*w1 + 3*w2 - 5*w3 + 2*w4) % 11
	e3 = (w1 + 5*w2 - 3*w3 - 3*w4) % 11
	e4 = (-3*w1 - w2 + 4*w3 + w4) % 11

	alpha_prim_assoc = eta1_p**e1 * eta2_p**e2 * eta3_p**e3 * eta4_p**e4 * alpha_normlzed

	print("[+] Checking correctness...")
	assert check_prim(alpha_prim_assoc)
	print("[+] Checking that \\alpha and \\alph_prim_assoc are associates ..")
	eps = alpha_in_ring_elem / alpha_prim_assoc
	assert eps.is_unit()

	print("[+] Done")
	return [e0, e1, e2, e3, e4, alpha_prim_assoc]

# --------------------- Alternative implementation ---------------------

print("[+] Testing function ...")
# Coeff.'s are given in the basis \zeta, \zeta^2, ..., \zeta^10
e1 = [0, 0, 0, 0, 1, 1, 0, 0, 0, 0]
e2 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 1]
e3 = [0, -1, -1, -1, -1, -1, -1, -1, -1, 0]
e4 = [0, 0, 1, 0, 0, 0, 0, 1, 0, 0]

test = [ZZ.random_element(1, 2**10) for i in range(10)]
print(compute_prim_assoc_alt(test, e1, e2, e3, e4))
print("-------------------------------------")
'''
eq1: -5*b0*m - 3*b0*n - 5*b0*o + 4*b0*p + a0
eq2: -2*b0*m^2 + 2*b0*m*n + 5*b0*n^2 - 4*b0*m*o + 2*b0*n*o - 2*b0*o^2 + b0*m*p + 5*b0*n*p + b0*o*p + 4*b0*p^2 + 3*a0*m - 3*b0*m + 4*a0*n + 2*b0*n + 3*a0*o + 3*b0*o + 2*a0*p + 2*b0*p + e0
eq3: -5*b0^2*m^3 + 2*b0^2*m^2*n - b0^2*m*n^2 + 2*b0^2*n^3 - 4*b0^2*m^2*o + 4*b0^2*m*n*o - b0^2*n^2*o - 4*b0^2*m*o^2 + 2*b0^2*n*o^2 - 5*b0^2*o^3 + b0^2*m^2*p - b0^2*m*n*p + 3*b0^2*n^2*p + 2*b0^2*m*o*p - b0^2*n*o*p + b0^2*o^2*p - 3*b0^2*m*p^2 - 4*b0^2*n*p^2 - 3*b0^2*o*p^2 + 3*b0^2*p^3 + 3*a0*b0*m^2 + b0^2*m^2 - 3*a0*b0*m*n - 4*b0^2*m*n - 2*a0*b0*n^2 + 4*b0^2*n^2 - 5*a0*b0*m*o - 3*a0*b0*n*o + b0^2*n*o + 3*a0*b0*o^2 + 2*b0^2*o^2 + 4*a0*b0*m*p - 2*a0*b0*n*p - 3*b0^2*n*p + 4*a0*b0*o*p - 3*b0^2*o*p + 5*a0*b0*p^2 - 3*b0^2*p^2 - a0*b0*m + 2*b0^2*m + 2*b0*e0*m - 3*a0*b0*n - 3*b0^2*n - b0*e0*n + a0*b0*o - b0^2*o + 2*b0*e0*o - 3*a0*b0*p - 3*b0^2*p + 5*b0*e0*p + d0^2 + b0*g0
eq4: -b0^2*m^4 - 2*b0^2*m^3*n - 4*b0^2*m^2*n^2 + 5*b0^2*m*n^3 + 2*b0^2*n^4 + 4*b0^2*m^3*o + 5*b0^2*m^2*n*o + 3*b0^2*m*n^2*o + 5*b0^2*n^3*o - 5*b0^2*m^2*o^2 + 5*b0^2*m*n*o^2 - 4*b0^2*n^2*o^2 + 4*b0^2*m*o^3 - 2*b0^2*n*o^3 - b0^2*o^4 - b0^2*m^3*p - 4*b0^2*m^2*n*p + 2*b0^2*m*n^2*p - 4*b0^2*n^3*p - 3*b0^2*m^2*o*p + 3*b0^2*m*n*o*p + 2*b0^2*n^2*o*p - 3*b0^2*m*o^2*p - 4*b0^2*n*o^2*p - b0^2*o^3*p - b0^2*m^2*p^2 + b0^2*m*n*p^2 - 3*b0^2*n^2*p^2 - 2*b0^2*m*o*p^2 + b0^2*n*o*p^2 - b0^2*o^2*p^2 + 2*b0^2*m*p^3 - b0^2*n*p^3 + 2*b0^2*o*p^3 - 4*b0^2*p^4 - 3*a0*b0*m^3 - 3*b0^2*m^3 - a0*b0*m^2*n - 3*b0^2*m^2*n - 5*a0*b0*m*n^2 - 3*b0^2*m*n^2 - a0*b0*n^3 - 5*b0^2*n^3 + 2*a0*b0*m^2*o + 5*b0^2*m^2*o - 2*a0*b0*m*n*o + 2*b0^2*m*n*o - 5*a0*b0*n^2*o - 2*b0^2*n^2*o + 2*a0*b0*m*o^2 - 2*b0^2*m*o^2 - a0*b0*n*o^2 + 4*b0^2*n*o^2 - 3*a0*b0*o^3 + b0^2*o^3 + 5*a0*b0*m^2*p - 5*b0^2*m^2*p - 5*a0*b0*m*n*p - 2*b0^2*m*n*p + 4*a0*b0*n^2*p - 4*b0^2*n^2*p - a0*b0*m*o*p + 5*b0^2*m*o*p - 5*a0*b0*n*o*p - b0^2*n*o*p + 5*a0*b0*o^2*p + 4*b0^2*o^2*p - 4*a0*b0*m*p^2 + 3*b0^2*m*p^2 + 2*a0*b0*n*p^2 - 4*a0*b0*o*p^2 - 5*b0^2*o*p^2 + 4*a0*b0*p^3 + 4*b0^2*p^3 - 2*a0*b0*m^2 + 3*b0*d0*m^2 - 3*b0*e0*m^2 - 4*a0*b0*m*n + 3*b0^2*m*n - b0*d0*m*n + 3*b0*e0*m*n + 3*a0*b0*n^2 - 3*b0*d0*n^2 + 2*b0*e0*n^2 - 3*a0*b0*m*o - 5*b0*d0*m*o + 5*b0*e0*m*o - a0*b0*n*o - b0*d0*n*o + 3*b0*e0*n*o + 3*a0*b0*o^2 - 2*b0^2*o^2 + 3*b0*d0*o^2 - 3*b0*e0*o^2 - 4*a0*b0*m*p - 4*b0^2*m*p - 5*b0*d0*m*p - 4*b0*e0*m*p - a0*b0*n*p + 2*b0^2*n*p + 2*b0*d0*n*p + 2*b0*e0*n*p + 3*a0*b0*o*p - 2*b0^2*o*p - 5*b0*d0*o*p - 4*b0*e0*o*p - 3*a0*b0*p^2 - 5*b0^2*p^2 - b0*d0*p^2 - 5*b0*e0*p^2 - a0*b0*m + b0^2*m - 5*a0*d0*m + 5*d0^2*m + b0*e0*m - 3*b0*g0*m - 4*a0*b0*n - 2*b0^2*n + a0*d0*n + 3*d0^2*n + 3*b0*e0*n - 4*b0*g0*n - 5*a0*b0*o + 3*b0^2*o - 5*a0*d0*o + 5*d0^2*o - b0*e0*o - 3*b0*g0*o - 4*a0*b0*p + 4*b0^2*p - 3*a0*d0*p - 4*d0^2*p + 3*b0*e0*p - 2*b0*g0*p + d0*f0 - b0*k0
'''


import sympy

def ith_der (alpha, i):
	#print "The {}th deriv. of {}".format(i, alpha)
	var('v e x')

	variab = vector([1, sympy.E**v, sympy.E**(v*2), sympy.E**(v*3), sympy.E**(v*4), sympy.E**(v*5), sympy.E**(v*6), sympy.E**(v*7), sympy.E**(v*8), sympy.E**(v*9)])
	#alpha = vector([-alpha[9]] + [alpha[j-1]-alpha[9] for j in range(1, 10)])
	alpha = vector(change_basis(alpha))

	#print ("-------------------------------------")
	p_reduced = (derivative(log(variab.dot_product(alpha)), v, i))
	#print "p_red: {}".format(p_reduced)
	p_reduced = (derivative(log(variab.dot_product(alpha)), v, i)).subs(v == 0)
	#print "p_red(evaluated at v == 0): {}".format(p_reduced)

	a = int(p_reduced.numerator()) % 11
	b = int(p_reduced.denominator()) % 11
	inv_b = (GF(11)(b))^(-1)
	#print "= {}".format(a*inv_b)
	#print ("-------------------------------------")

	return a*inv_b

def pi_deriv(j): # compute the ith derivative at 1 of pi(t) = -a_{10} - a_{10}*\sum_{i=1}^{9}*t^i + \sum_{i=1}^{9}*a_i*t^i
	var('v t a1 a2 a3 a4 a5 a6 a7 a8 a9 a10')
	# bugs here !!! Not-usable
	t = exp(v)
	pi = -a10 - a10*(t + t^2 + t^3 + t^4 + t^5 +t^6 + t^7 + t^8 + t^9 + t^10) + (a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5 + a6*t^6 + a7*t^7 + a8*t^8 + a9*t^9)
	'''f_0 = (pi.derivative(v, 0)).substitute(v==0)
	f_1 = (pi.derivative(v, 1)).substitute(v==0)
	f_2 = (pi.derivative(v, 2)).substitute(v==0)
	f_3 = (pi.derivative(v, 3)).substitute(v==0)
	f_4 = (pi.derivative(v, 4)).substitute(v==0)
	f_5 = (pi.derivative(v, 5)).substitute(v==0)
	f_6 = (pi.derivative(v, 6)).substitute(v==0)
	f_7 = (pi.derivative(v, 7)).substitute(v==0)
	f_8 = (pi.derivative(v, 8)).substitute(v==0)
	f_9 = (pi.derivative(v, 9)).substitute(v==0)
	f_10 = (pi.derivative(v, 10)).substitute(v==0)
	f_11 = (pi.derivative(v, 11)).substitute(v==0)'''

	for i in range(12):
		print((pi.derivative(v, i)).substitute(v==0)).polynomial(GF(121))



	'''for i in range(j - 1):
		pi_reduced = simplify((pi_reduced.derivative(v)).expand())
	pi_reduced = simplify((pi_reduced.substitute(v==0)).expand())'''

var('v')
#pi_deriv(2) # try to find F(), F`(), F``(), ...
f = function('f')
g(v) = log(f(exp(v)))
#print (g.derivative(v, 11)).substitute(v==0)

'''
a = [0, -1, -1, -1, -1, -1, -1, -1, -1, -1]
b = [-1, 0, -1, -1, -1, -1, -1, -1, -1, -1]
c = [0, 0, -1, -1, -1, -1, -1, -1, -1, -1]
d = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0]

a = [-1, -1, -1, -1, -1, -1, -1, -1, -1, 0]
b = [-1, -1, -1, 0, -1, -1, -1, -1, -1, -1]
c = [0, 0, 0, 0, -1, 0, 0, -1, 0, 0]
d = [0, 0, 0, 0, -1, 0, -1, 0, 0, 0]'''

# All initialization for the computations below starts here ...
e1 = [0, 0, 0, 0, 1, 1, 0, 0, 0, 0]
e2 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 1]
e3 = [0, -1, -1, -1, -1, -1, -1, -1, -1, 0]
e4 = [0, 0, 1, 0, 0, 0, 0, 1, 0, 0]
e5 = [0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0] # aux. element, ~ 'dummy'
k_f = CyclotomicField(11)
unit_r = k_f.unit_group() # Z*[\zeta]
ring = k_f.maximal_order() # Z[\zeta]
spec_gamma = k_f(1)


def compute_D3579():
	var('a b c d e f g h k l v')
	# given \alpha, we have a, b, c, d, ..., l as above
	# moreover, since \alpha is primary, a == e == c == 0
	# Taylor expansion of ln(A(v)), deg == 9
	f_ = l*v^9/factorial(9)+k*v^8/factorial(8)+h*v^7/factorial(7)+g*v^6/factorial(6)+f*v^5/factorial(5)+e*v^4/24+d*v^3/6+a*v^2/2+c*v+b
	fun(v) = log(v)
	D3 = diff(fun(f_), v, 3).substitute(v ==0, a==0, e==0, c==0)
	D5 = diff(fun(f_), v, 5).substitute(v ==0, a==0, e==0, c==0)
	D7 = diff(fun(f_), v, 7).substitute(v ==0, a==0, e==0, c==0)
	D9 = diff(fun(f_), v, 9).substitute(v ==0, a==0, e==0, c==0)

	return D3, D5, D7, D9

def index_fundams(u):
	var('a b c d e f g h k l v')

	D3, D5, D7, D9 = compute_D3579() # c(fundams) == 0 mod 11


	# D1 = ith_der(a, 1) == 0
	D2 = ith_der(u, 2)
	D4 = ith_der(u, 4)
	D6 = ith_der(u, 6)
	D8 = ith_der(u, 8)

	e_ = int(D2) * D9 + int(D4) * D7 + int(D6) * D5 + int(D8) * D3
	h_(a, b, c, d, e, f, g, h, k, l) = e_

	return h_(a, b, c, d, e, f, g, h, k, l)


def norm(alph): # computing the norm of \alph in Z[\zeta]
	# in basis: \zeta, \zeta^2, ..., \zeta^10
	coord_alph = alph# since the norm formula is expressed in basis \zeta, \zeta^2,... ,\zeta^10

	zet = ring([0, 1])
	basis = [zet^i for i in range(1, 11)]
	res = ring([1])

	for j in range(1,11):
		res = res * (vector(coord_alph).dot_product(vector([z^j for z in basis])))

	return res

def trace(alph):
	# input as element in the cyclotomic fields
	zet = ring([0, 1])
	basis = [zet^i for i in range(1, 11)]
	res = ring([0])
	get_coords = zet.coordinates_in_terms_of_powers()
	coord_alph = change_basis_to_1_10(get_coords(alph))


	for j in range(1,11):
		res = res + (vector(coord_alph).dot_product(vector([z^j for z in basis])))

	return res

def trace_w_conj(x):
	# given x as an element in Q[\zeta]
	# compute T(x*x.conj)
	x_conj = x.conjugate()
	prod = x * x_conj
	return trace(prod)

def approx_by_cyclo_int(coord_x):
	# Given x \in Q[\zeta], expressed as a vector of coordinates
	# find y \in Z[zeta] s.t. N(x - y) < 1
	print("[+] Entering approximation ...")
	pwrs = [k_f.gen()**i for i in range(11)] # = [1, \zeta, \zeta^2, ..., \zeta^10]
	indices = [i for i in range(11)] # = [0, 1, 2, ..., 10]
	z11 = k_f.gen() # zeta11
	x = vector(coord_x).dot_product(vector(pwrs)) # now x is an element in Q[\zeta]
	print(x)
	get_coords = z11.coordinates_in_terms_of_powers()
	#print "[-] c_x = {}".format(coord_x)
	#print "[-] x = {}".format(x)

	# chck trivial case: if x \in Z[\zeta]
	if x.is_integral():
		y_ = x
		return y_
	else:
		coord_y = [floor(crd_x) for crd_x in coord_x]
		tmp = zip(coord_x, coord_y)
		coord_z = [t[0] - t[1] for t in tmp]

		y_ = vector(coord_y).dot_product(vector(pwrs))
		z_ = vector(coord_z).dot_product(vector(pwrs))
		# now y, z are two elements in Z[\zeta], Q[\zeta], resp.
		#print "[-] y = {}".format(y_)
		#print "[-] z = {}".format(z_)

		assert len(coord_z) == len(coord_y) == 11
		#we sort coord_z in non-descending order
		coord_z_with_ind = zip(coord_z, indices) # each element has form (crd, power of corr. zeta^i)

		sorted_coord_z_with_ind = sorted(coord_z_with_ind, key=lambda x:x[0], reverse=False)
		# sort the newly created list, using the first element, in non-descending order
		lower_bnd = (11**2 - 1) / 12 # from Lenstra's method
		print(k_f(lower_bnd))
		u = trace_w_conj(z_)
		print(u)

		while QQ(u) > QQ(lower_bnd):
			y_ = y_ - z11^sorted_coord_z_with_ind[0][1]
			z_ = z_ + z11^sorted_coord_z_with_ind[0][1]


			t = sorted_coord_z_with_ind[0]
			# we circulate z0 = z1, z1 = z2, ..., z(9) = z(10), z(10) = t + 1
			sorted_coord_z_with_ind = [(s_pair[0], s_pair[1]) for s_pair in sorted_coord_z_with_ind[1:]]
			assert len(sorted_coord_z_with_ind) == 10
			sorted_coord_z_with_ind += [(t[0] + 1, t[1])]
			#print "[-] circulated srted_crd_z: {}".format(sorted_coord_z_with_ind)
			u = trace_w_conj(z_)

	#print "[-] out while ..."
	#print "[-] y: {}".format(y_)
	#print "[-] x - y: {}".format((x - y_))
	norm_x_subs_y = QQ(norm(change_basis_to_1_10(get_coords(x - y_))))

	if norm_x_subs_y == 1:
		spec_gamma = y_
		raise Exception("Norm gamma == 1")
	assert norm_x_subs_y < QQ(1)

	print("[-] Done")
	return y_

def euclid_division(alph, beta):
	# alph, beta \in Z[\zeta]
	# expressed in basis \zeta, \zeta^2, .., \zeta^10
	print("[+] Entering euclid_division ...")

	pwrs = [k_f.gen()**i for i in range(1, 11)] # = [1, \zeta, \zeta^2, ..., \zeta^9]
	indices = [i for i in range(1, 11)] # = [0, 1, 2, ..., 9]
	z11 = k_f.gen() # zeta11
	get_coords = z11.coordinates_in_terms_of_powers()

	alp = vector(alph).dot_product(vector(pwrs))
	bet = vector(beta).dot_product(vector(pwrs))

	assert bet != ring([0])
	assert alp.is_integral() and bet.is_integral()
	x = alp / bet

	coord_x = list(get_coords(x)) + [0]# we append a coordinate for \zet^10, for the execution of
								# approximation function

	'''
	try:
		y = approx_by_cyclo_int(coord_x)
	except Exception as e:
		y = spec_gamma'''
	y = approx_by_cyclo_int(coord_x)
	rho = alp - y * bet

	norm_rho = QQ(norm(change_basis_to_1_10(get_coords(rho))))
	norm_bet = QQ(norm(beta))

	print("[-] Checking result has norm strictly smaller ...")
	assert norm_rho < norm_bet
	return rho

alph = [0, -1, -1, -1, -1, -1, -1, -1, -1, -1]
print("[+] Testing norm...")
print("[+] alph = {}".format(alph))
print("[+] norm(alph)= {}".format(norm(alph)))
assert norm(alph) == ring(change_basis(alph)).absolute_norm() # double check
print("-------------------------------------")

zet = ring([0, 1])
tst_pr = ring([1, 0, 0, 1, 0, 0, 0, 1, 0, 0])
assert ring(change_basis([-1, -1, 0, -1, -1, -1, 0, -1, -1, -1])) == tst_pr

print("[+] Initialize a primary prime ...")
alpha_init = vector([-1, -1, 0, -1, -1, -1, 0, -1, -1, -1])
alpha_init = B * alpha_init
print("[+] alpha: {}".format([-1, -1, 0, -1, -1, -1, 0, -1, -1, -1]))
print("[+] Checking b(\\alpha) != 0 mod 11 ...")
assert alpha_init[0] % 11 != 0
# compute the appropriate power for factor \zeta
print("[+] Normalizing ... ")
k_power = (-1 * 1/alpha_init[0] * alpha_init[1]) % 11
print("[+] k_power: {}".format(k_power))
pr = tst_pr * zet^k_power
coord_pr = ring.coordinates(pr) # get coord.'s in 1, \zeta, ..., \zeta^9
#change to basis \zeta, ..., \zeta^10
coord_pr = [coord_pr[1]- coord_pr[0], coord_pr[2]- coord_pr[0], coord_pr[3]- coord_pr[0], coord_pr[4]- coord_pr[0], coord_pr[5]- coord_pr[0], coord_pr[6]- coord_pr[0], coord_pr[7]- coord_pr[0], coord_pr[8]- coord_pr[0], coord_pr[9] - coord_pr[0], - coord_pr[0]]

print("[+] Check primary: {}".format(check_prim(pr))) # pr is now a primary prime
print("[+] pr = {}".format(pr))
print("[+] Done")
print("-------------------------------------")

x = [QQ.random_element(2*10, 2*10) for i in range(11)]
print("[+] Testing approximation ...")
print("[+] x: {}".format(x))
print("[+] Approx. x in Z[\\zeta]: {}".format(approx_by_cyclo_int(x)))

# in basis \zeta, \zeta^2, ..., \zeta^10
alp = [ZZ.random_element() for i in range(10)]
bet = [ZZ.random_element() for i in range(10)]
print("[+] Checking \\beta != 0 ...")
assert ring(change_basis(bet)) != ring([0])

print("[+] Testing Euclidean div. ...")
print("[+] alp: {}".format(alp))
print("[+] bet: {}".format(bet))
print("[+] \\rho s.t \\alp \\equiv \\rho (mod \\bet): {}".format(euclid_division(alp, bet)))
print("-------------------------------------")

# get a prime from k_f by:
# choosing a fractional ideal from k_f.primes_of_degree_one_list(3)
# let f_idl be that ideal, then call f_idl.gen(1)
omg = ring([1, -1])
zet = ring([0, 1])
get_coords = zet.coordinates_in_terms_of_powers()

B = []
for j in IntegerRange(0, 10):
	B += [[i**j % 11 for i in list(IntegerRange(1, 11))]]

B = Matrix(GF(11), B)
print(B * vector(e1))
print(~B)

e1 = [0, 0, 0, 0, 1, 1, 0, 0, 0, 0]
e2 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 1]
e3 = get_coords(zet^(-1) * (1 + zet + zet^2))
#e3 = change_basis_to_1_10([0, 0, -1, -1, -1, -1, -1, -1, -1, -1])
e3 = change_basis_to_1_10(e3)
e3 = matrix(ZZ, e3)
e3 = vector(e3)
e4 = [0, 0, 1, 0, 0, 0, 0, 1, 0, 0]

print("[+] All the indices(w.r.t pr) of: ")
print("[-] e1 = {}".format(e1))
print("[-] e2 = {}".format(e2))
print("[-] e3 = {}".format(e3))
print("[-] e4 = {}".format(e4))
#print "[-] zeta = {}".format(zet)
print("[-] \\omg = {}".format(omg))
'''print "e1 : ----"
ind_e1 = simplify(index_fundams(e1).expand())
print "e2 : ----"
ind_e2 = simplify(index_fundams(e2).expand())
print "e3 : ----"
ind_e3 = simplify(index_fundams(e3).expand())
print "e4 : ----"
ind_e4 = simplify(index_fundams(e4).expand())
print "ind_omg : ----"
ind_omg = simplify(index_omg(coord_pr, e1, e2, e3, e4, e5).expand())
print "[+] (In that order)\n {}\n {}\n {}\n {}\n {}\n".format(ind_e1, ind_e2, ind_e3, ind_e4, ind_omg.substitute(a==0, e==0, c==0))
'''
var('a b c d e f g h k l ind_zet')
# c.f. paper
ind_e1 = d/b + 3*f/b + 4*h/b + 4*g*d/b^2 + 3*l/b
ind_e2 = 3*d/b + 5*f/b + 9*h/b + 5*g*d/b^2 + l/b
ind_e3 = 2*d/b + 9*f/b + 3*h/b + 7*g*d/b^2 + 8*l/b
ind_e4 = 4*d/b + 4*f/b + 3*h/b + g*d/b^2 + 9*l/b


var('ind_zet')
#eps = (-1) * zet^6 * eta1_p^4 * eta2_p^2 * eta3_p^4 * eta4_p^0 * dummy^2 * omg^10
#dummy = \eta3^-1 * \eta4^-1
def index_omg(pi, e1, e2, e3, e4, e5):
	# computing ind_\pi(\omg)
	# pi is expressed in basis \zeta, ..., \zeta^10
	var('a b c c_ d e f g h k l v ind_11_')
	'''f_ = c*v^11/factorial(11) + l*v^9/factorial(9)+k*v^8/factorial(8)+h*v^7/factorial(7)+g*v^6/factorial(6)+f*v^5/factorial(5)+e*v^4/24+d*v^3/6+a*v^2/2+c*v+b

	fun(v) = log(v)
	ind_e1 = index_fundams(e1)
	ind_e2 = index_fundams(e2)
	ind_e3 = index_fundams(e3)
	ind_e4 = index_fundams(e4)
	ind_e5 = index_fundams(e5)'''

	coord_pi = pi


	ind_e1 = d/b + 3*f/b + 4*h/b + 4*g*d/b^2 + 3*l/b
	ind_e2 = 3*d/b + 5*f/b + 9*h/b + 5*g*d/b^2 + l/b
	ind_e3 = 2*d/b + 9*f/b + 3*h/b + 7*g*d/b^2 + 8*l/b
	ind_e4 = 4*d/b + 4*f/b + 3*h/b + g*d/b^2 + 9*l/b

	#f_11 = vector([i^11 for i in range(1, 11)]).dot_product(coord_pi) / (sum(coord_pi) - 11 * coord_pi[-1])
	f_1 = vector([i for i in range(1, 11)]).dot_product(coord_pi) / (vector([1 for i in range(1, 11)]).dot_product(coord_pi)) % 121

	#ind_11 = 2*f*g / b^2 + (f_11) / 11


	# multiply by b^3 then divide by b^3 to put the result in poly(GF(11))
	#ind_omg = (0 + 6*ind_zet + 4*ind_e1 + 2*ind_e2 + 2*ind_e3 + (-2)*ind_e4 - ind_11)
	ind_omg = 5*g*d/b^2 - 2*g*f/b^2 - f_1 / 11 -5 - 5*ind_zet
	# need to devide ind_11_ by 11, and convert ind_omg to a polynomial in GF(11)

	return ind_omg

def compute_resid_symbol(alpha_, beta_):
	print("[+] Checking gcd(\\alpha_ , \\beta) is unit ... ")
	assert (k_f.gcd(alpha_, beta_)).is_unit()
	assert (beta_.trace()) % 11 != 0

	e1 = [0, 0, 0, 0, 1, 1, 0, 0, 0, 0]
	e2 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 1]
	e3 = get_coords(zet^(-1) * (1 + zet + zet^2))
	#e3 = change_basis_to_1_10([0, 0, -1, -1, -1, -1, -1, -1, -1, -1])
	e3 = change_basis_to_1_10(e3)
	e3 = matrix(ZZ, e3)
	e3 = vector(e3)
	e4 = [0, 0, 1, 0, 0, 0, 0, 1, 0, 0]
	# compute s s.t (alpha/beta)_11 = \zeta^s
	pwrs = [k_f.gen()**i for i in range(10)]

	s = 0
	# these coordinates are in basis 1, \zeta, ..., \zeta^9
	coord_alp = get_coords(alpha_)
	coord_bet = get_coords(beta_)


	# change to basis \zeta, \zeta^2, ..., \zeta^10
	coord_alp = change_basis_to_1_10(coord_alp)
	norm_alp = QQ(norm(coord_alp))

	coord_bet = change_basis_to_1_10(coord_bet)

	# recompute ind_omg w.r.t \beta
	#ind_omg = simplify(index_omg(coord_bet, e1, e2, e3, e4, e5).expand())

	print("[+] Computing prim_assoc of beta ...")
	bet_prim = compute_prim_assoc_alt(coord_bet, e1, e2, e3, e4)[-1]

	# change to basis \zeta, \zeta^2, ..., \zeta^10
	coord_bet_prim = change_basis_to_1_10(get_coords(bet_prim))
	# coerce to a vetor over ZZ
	coord_bet_prim = matrix(ZZ, coord_bet_prim)
	coord_bet_prim = vector(coord_bet_prim)
	norm_bet_prim = ZZ(norm(coord_bet_prim))
	assert norm_bet_prim == norm(coord_bet)

	while norm_alp > QQ(1):
		old_norm = norm_alp

		sum_pwr_bet_prim = B * vector(coord_bet_prim) # compute a, b, c, ... of bet_prim
		c_subs = 0
		#c_subs = vector([i for i in range(1, 11)]).dot_product(coord_bet_prim) % 121


		print("[+] compute a, b, c, ... of bet_prim")#.format(sum_pwr_bet_prim)
		#find \gamma s.t. alp == gamma mod bet_prime
		print("[+] Finding gamma s.t. alph \\equiv gamma mod bet_prim ...")
		gamma = euclid_division(coord_alp, coord_bet_prim) #
		print("[+] gamma ...")#.format(gamma)
		coord_gamma = change_basis_to_1_10(get_coords(gamma))
		print(coord_gamma)
		coord_gamma = matrix(ZZ, coord_gamma)
		coord_gamma = vector(coord_gamma)

		sum_pwr_gamm = B * vector(coord_gamma) #

		i = 0
		t_gamm = ZZ(trace(ring(change_basis(coord_gamma))))
		while t_gamm % 11 == 0: # remove factor \\omg
			i += 1
			t_gamm = t_gamm / 11
			gamma = gamma / omg

		print("i = {}".format(i))

		#now b(gamma) != 0, we can find its prim. assoc
		coord_gamma_curr = change_basis_to_1_10(get_coords(gamma))
		coord_gamma_curr = matrix(ZZ, coord_gamma_curr)
		coord_gamma_curr = vector(coord_gamma_curr)

		print("[+] Computing prim_assoc of gamma ...")
		gamma_prim = compute_prim_assoc_alt(coord_gamma_curr, e1, e2, e3, e4)
		print("[+] gramma_prim = {}".format(gamma_prim))
		print("n_gmm_prim = {}".format(gamma_prim[-1].norm()))

		r = IntegerModRing(11)
		f_1 = vector([i for i in range(1, 11)]).dot_product(coord_bet_prim) / (vector([1 for i in range(1, 11)]).dot_product(coord_bet_prim)) % 121
		norm_bet_prim = ZZ(norm(coord_bet_prim))

		index_zet_subs = ZZ((norm_bet_prim - 1)/11)	# by multipliccativity of index
		ind_omg = 5*g*d/b^2 - 2*g*f/b^2 - (f_1 / 11) -5 - 5*ind_zet
		#ind_omg = simplify(index_omg(coord_bet_prim, e1, e2, e3, e4, e5).expand()) # recompute ind_omg


		# add multiple of omg's index
		num = (ind_omg.substitute(ind_zet=index_zet_subs, b=ZZ(sum_pwr_bet_prim[0]), c=c_subs, a=ZZ(sum_pwr_bet_prim[2]), d=ZZ(sum_pwr_bet_prim[3]), e=ZZ(sum_pwr_bet_prim[4]), f=ZZ(sum_pwr_bet_prim[5]), g=ZZ(sum_pwr_bet_prim[6]), h=ZZ(sum_pwr_bet_prim[7]), k=ZZ(sum_pwr_bet_prim[8]), l=ZZ(sum_pwr_bet_prim[9]))).numerator()
		denom = (ind_omg.substitute(ind_zet=index_zet_subs, b=ZZ(sum_pwr_bet_prim[0]), c=c_subs, a=ZZ(sum_pwr_bet_prim[2]), d=ZZ(sum_pwr_bet_prim[3]), e=ZZ(sum_pwr_bet_prim[4]), f=ZZ(sum_pwr_bet_prim[5]), g=ZZ(sum_pwr_bet_prim[6]), h=ZZ(sum_pwr_bet_prim[7]), k=ZZ(sum_pwr_bet_prim[8]), l=ZZ(sum_pwr_bet_prim[9]))).denominator()

		print("ind_omg: {}".format(ind_omg.substitute(ind_zet=index_zet_subs, b=ZZ(sum_pwr_bet_prim[0]), c=c_subs, a=ZZ(sum_pwr_bet_prim[2]), d=ZZ(sum_pwr_bet_prim[3]), e=ZZ(sum_pwr_bet_prim[4]), f=ZZ(sum_pwr_bet_prim[5]), g=ZZ(sum_pwr_bet_prim[6]), h=ZZ(sum_pwr_bet_prim[7]), k=ZZ(sum_pwr_bet_prim[8]), l=ZZ(sum_pwr_bet_prim[9]))))

		s = s + i * (r(num) / r(denom))
		# add multiple of zeta index
		s = s - r(gamma_prim[0]*index_zet_subs)#ind_zet)

		# add multiple of e1, e2, e3, e4 indexes
		tmp1 = (ind_e1.substitute(b=ZZ(sum_pwr_bet_prim[0]), c=c_subs, a=ZZ(sum_pwr_bet_prim[2]), d=ZZ(sum_pwr_bet_prim[3]), e=ZZ(sum_pwr_bet_prim[4]), f=ZZ(sum_pwr_bet_prim[5]), g=ZZ(sum_pwr_bet_prim[6]), h=ZZ(sum_pwr_bet_prim[7]), k=ZZ(sum_pwr_bet_prim[8]), l=ZZ(sum_pwr_bet_prim[9])))
		tmp1 = (r(tmp1.numerator()) * r(tmp1.denominator())^-1)
		s = s - r(gamma_prim[1]*tmp1)

		tmp2 = (ind_e2.substitute(b=ZZ(sum_pwr_bet_prim[0]), c=c_subs, a=ZZ(sum_pwr_bet_prim[2]), d=ZZ(sum_pwr_bet_prim[3]), e=ZZ(sum_pwr_bet_prim[4]), f=ZZ(sum_pwr_bet_prim[5]), g=ZZ(sum_pwr_bet_prim[6]), h=ZZ(sum_pwr_bet_prim[7]), k=ZZ(sum_pwr_bet_prim[8]), l=ZZ(sum_pwr_bet_prim[9])))
		tmp2 = (r(tmp2.numerator()) * r(tmp2.denominator())^-1)

		s = s - r(gamma_prim[2]*tmp2)

		tmp3 = (ind_e3.substitute(b=ZZ(sum_pwr_bet_prim[0]), c=c_subs, a=ZZ(sum_pwr_bet_prim[2]), d=ZZ(sum_pwr_bet_prim[3]), e=ZZ(sum_pwr_bet_prim[4]), f=ZZ(sum_pwr_bet_prim[5]), g=ZZ(sum_pwr_bet_prim[6]), h=ZZ(sum_pwr_bet_prim[7]), k=ZZ(sum_pwr_bet_prim[8]), l=ZZ(sum_pwr_bet_prim[9])))
		tmp3 = (r(tmp3.numerator()) * r(tmp3.denominator())^-1)

		s = s - r(gamma_prim[3]*tmp3)

		tmp4 = (ind_e4.substitute(b=ZZ(sum_pwr_bet_prim[0]), c=c_subs, a=ZZ(sum_pwr_bet_prim[2]), d=ZZ(sum_pwr_bet_prim[3]), e=ZZ(sum_pwr_bet_prim[4]), f=ZZ(sum_pwr_bet_prim[5]), g=ZZ(sum_pwr_bet_prim[6]), h=ZZ(sum_pwr_bet_prim[7]), k=ZZ(sum_pwr_bet_prim[8]), l=ZZ(sum_pwr_bet_prim[9])))
		tmp4 = (r(tmp4.numerator()) * r(tmp4.denominator())^-1)
		s = s - r(gamma_prim[4]*tmp4)

		print("[+] (So far) s = {}".format(s))

		if gamma_prim[-1].norm() == 1:
			break # norm(gamm_prim) == 1 --> (gamm_prim/bet) = 1 --> stop


		# these coordinates are in basis 1, \zeta, ..., \zeta^9
		coord_alp = coord_bet_prim #get_coords(bet_prim)
		coord_bet_prim = get_coords(gamma_prim[5]) # gamma-prim is currently a list returned from compute_prim_assoc

		# change to basis \zeta, \zeta^2, ..., \zeta^10
		coord_alp = change_basis_to_1_10(coord_alp)
		coord_alp = matrix(ZZ, coord_alp)
		coord_alp = vector(coord_alp)
		norm_alp = QQ(norm(coord_alp))

		coord_bet_prim = change_basis_to_1_10(coord_bet_prim)
		coord_bet_prim = matrix(ZZ, coord_bet_prim)
		coord_bet_prim = vector(coord_bet_prim)

		assert norm_alp < old_norm # otherwise, we need the norm decreasing strictly

	return s

'''alp = [ZZ.random_element() for i in range(10)]
bet = [ZZ.random_element() for i in range(10)]'''
alp = [-11 for i in range(10)]
bet1 = [18, 61, -43, -109, 24, 83, 178, 158, -23, -23]
bet2 = [-549, -3496, -5532, -4132, -1123, 213, -1919, -4870, -5216, -2659]
bet3 = [-482657, -162903, 585409, 887367, 389923, -325331, -422142, 212683, 836933, 720762]

s = compute_resid_symbol(ring(change_basis(alp)), ring(change_basis(bet3)))
print("[+] (alp/bet)_11 = \\zeta^{}".format(s))
