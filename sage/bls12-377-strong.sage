QQx.<x> = QQ[]
rx = x**4 - x**2 + 1
px = (x-1)**2*rx/3 + x
u=-9322841860747558913
p=Integer(px(u))
r=Integer(rx(u))
Fp=GF(p)
E1=EllipticCurve(Fp, [0, 1])
assert(E1.order()%r == 0)

h1=E1.order()//r
e1=(h1//3).sqrt()
e=3*e1
assert((e1//2).is_prime())
pi= 2**4 * 3**2 * 5 * 7 * 11 * 13
pi_=gcd(e,pi)
assert(pi_ == 2*3)

P3=E1(0,1)
assert(3*P3 == E1(0))
for i in range(10):
    Q = h1*E1.random_point()
    assert(r*Q == E1(0))
    assert(P3.tate_pairing(Q, 3, 1) == 1)
    # f_{3,P3}(Q) = y_Q-1
    assert((Q[1]-1)^((p-1)//3) == 1)

P2=E1(-1,0)
assert(2*P2 == E1(0))
for i in range(10):
    Q = h1*E1.random_point()
    assert(r*Q == E1(0))
    assert(P2.tate_pairing(Q, 2, 1) == 1)
    # f_{2,P2}(Q) = x_Q+1
    assert((Q[0]+1)^((p-1)//2) == 1)
