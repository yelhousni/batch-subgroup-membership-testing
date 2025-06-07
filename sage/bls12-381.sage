QQx.<x> = QQ[]
rx = x**4 - x**2 + 1
px = (x-1)**2*rx/3 + x
u=-15132376222941642752
p=Integer(px(u))
r=Integer(rx(u))
Fp=GF(p)
E1=EllipticCurve(Fp, [0, 4])
assert(E1.order()%r == 0)

h1=E1.order()//r
e1=(h1//3).sqrt()
e=3*e1
pi= 2**4 * 3**2 * 5 * 7 * 11 * 13
pi_=gcd(e,pi)
assert(pi_ == 3*11)

P3=E1(0,2)
assert(3*P3 == E1(0))
for i in range(10):
    Q = h1*E1.random_point()
    assert(r*Q == E1(0))
    assert(P3.tate_pairing(Q, 3, 1) == 1)
    # f_{3,P3}(Q) = y_Q-2
    assert((Q[1]-2)^((p-1)//3) == 1)

P=E1(4, Fp(4^3+4).sqrt())
P11=(P.order()//11)*P
assert(11*P11 == E1(0))
for i in range(10):
    Q = h1*E1.random_point()
    assert(r*Q == E1(0))
    assert(P11.tate_pairing(Q, 11, 1) == 1)
# f_{11,P} = (l_{P,P}^4 * (l_{4P,P} * l_{2P,2P})^2 * l_{5P,5P}) /
#   		 (v_{2P}^4 * (v_{5P} * v_{4P})^2)
Q2 = 2*P11
Q4 = 2*Q2
Q5 = Q4+P11
a_PP = 3*P11[0]^2/(2*P11[1])
Fxy.<x,y> = Fp[]
l_PP = (y-P11[1])-a_PP*(x-P11[0])
a_4PP = (Q4[1]-P11[1])/(Q4[0]-P11[0])
l_4PP = (y-Q4[1])-a_4PP*(x-Q4[0])
a_2P2P = 3*Q2[0]^2/(2*Q2[1])
l_2P2P = (y-Q2[1])-a_2P2P*(x-Q2[0])
a_5P5P = 3*Q5[0]^2/(2*Q5[1])
l_5P5P = (y-Q5[1])-a_5P5P*(x-Q5[0])
v_2P = x - Q2[0]
v_4P = x - Q4[0]
v_5P = x - Q5[0]

f1xy = (l_PP/v_2P)^4 * (l_4PP/v_5P * l_2P2P/v_4P)^2 * l_5P5P
f2xy = ((l_PP)^4 * (l_4PP * l_2P2P)^2 * l_5P5P) * (v_2P^4 * (v_4P * v_5P)^2)^10
for i in range(10):
    Q = h1*E1.random_point()
    assert(r*Q == E1(0))
    # f_{11,P11}(Q)
    assert(f1xy(Q[0],Q[1])^((p-1)//11) == 1)
    assert(f2xy(Q[0],Q[1])^((p-1)//11) == 1)
