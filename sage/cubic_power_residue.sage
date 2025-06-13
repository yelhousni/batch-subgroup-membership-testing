def quo_rem_R(a, b, R):
    """
     Quotient and remainder in Eisenstein ring of integers R
    :param a: element in R
    :param b: element in R
    :return q, r: where a = q*b + r
    """
    num = a * b.conjugate()
    denum = b.norm()
    q1, r1 = Integer(num[0]).quo_rem(Integer(denum))
    q2, r2 = Integer(num[1]).quo_rem(Integer(denum))
    return R(q1+j*q2), R(r1+j*r2)/b.conjugate()

def cubic_symbol(alpha, beta, R):
    """
     cubic power residue symbol
    :param alpha: element in K
    :param beta: element in K
    :return {-1, 0, 1}
    """
    a = alpha[0]
    b = alpha[1]
    c = beta[0]
    d = beta[1]
    if alpha == R(1) or alpha == R(-1) or beta == R(1) or beta == R(-1):
        return 1
    q, r = quo_rem_R(alpha, beta, R)
    gamma = alpha - q*beta
    e, f = gamma[0], gamma[1]
    if gamma == R(0):
        return 0
    m = 0
    while ((2*e-f)%3 == 0) and ((e+f)%3 == 0):
        m+=1
        e_ = e
        e = (2*e-f)/3
        f = (e_+f)/3
    gamma = e+j*f
    if f%3 == 0:
        n = 0
    elif (-e)%3 == 0:
        n = 1
        gamma = (f-e) + j*(-e)
    elif (e-f)%3 == 0:
        n = 2
        gamma = (-f) + j*(e-f)
    return j^((-m * (c^2-1) + n * (c^2-c*d-1))%9 // 3) * cubic_symbol(beta, gamma, R)

# Test for BLS12-381
p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
Fp = GF(p)
L.<x,y> = PolynomialRing(ZZ, ['x', 'y'])
C=Conic(x^2 - x*y + y^2 - p, x, y)
x,y,_ = C.rational_point()
assert(x%3 == 0)
R.<j> = EisensteinIntegers()
beta = R(y+j*x)
assert(beta.norm() == p)
Fz.<z> = Fp[]

for i in range(20):
    a = Fp.random_element()
    s = cubic_symbol(R(a), beta, R)
    if len((z^3-a).roots()) > 0:
        assert(s[0] == 1)
    else:
        assert(s[0] == 0 or s[0] == -1)
