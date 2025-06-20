from sage.all_cmdline import *   # import sage library
#
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from _utils import NAF, hamming_weight, NAF_hamming_weight
from sage.structure.proof.all import arithmetic

def search_bls12_curve():

    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    k = 12
    m = 24 
    n = 63 - m
    print("m = ", m)
    print("-------")

    j = 0
    if (j == 0):
        print("negative seeds: j = ", j)
    else:
        print("positive seeds: j = ", j)
    sign = (-1)**(1 - j)
    for z0 in range(598249620352, 2**(n+2)):
        z = ZZ(sign * z0 * 2**m - 1)
        e2 = 1 - z 
        if (e2 % 3) == 0: 
            r = z**4 - z**2 + 1
            re = z**4 - 3*z**2 + 3
            e1 = e2 // 3
            p = e1 * e2 * r + z 
            assert (e1 % 2) == 0
            e1Pr = e1 // 2
            if (p.nbits() < 384):
                if (p.is_pseudoprime() == True) & (r.is_pseudoprime() == True) & (abs(e1Pr).is_pseudoprime() == True) & (re.is_pseudoprime() == True):
                    h2 = (z**8 - 4*z**7 + 5*z**6 - 4*z**4 + 6*z**3 - 4*z**2 - 4*z + 13) // 9
                    hT = (z**20 - 8*z**19 + 25*z**18 - 32*z**17 - 8*z**16 + 76*z**15 - 93*z**14 + 36*z**13 + 51*z**12 - 112*z**11 
                        + 86*z**10 - 16*z**9 - 24*z**8 + 84*z**7 - 90*z**6 + 28*z**5 - 14*z**4 - 38*z**3 + 70*z**2 - 14*z + 73) // 81
                    if (h2.is_pseudoprime() == True):
                        vector_z = Integer(abs(z)).digits(2)
                        print("z = {:#x} {} bits".format(z, abs(z).nbits()))

                        assert( (r - 1) % 2**m == 0)
                        assert( (re - 1) % 2**m == 0)

                        if (hT.is_pseudoprime() == True):
                            vector_z = Integer(abs(z)).digits(2)
                            print("z = {:#x} {} bits".format(z, abs(z).nbits()))
                            print("r = {:#x} {} bits".format(r, r.nbits()))
                            print("re = {:#x} {} bits".format(re, re.nbits()))
                            print("p = {:#x} {} bits".format(p, p.nbits()))
                            print("h2 = {:#x} {} bits".format(h2, h2.nbits()))
                            print("hT = {:#x} {} bits".format(hT, hT.nbits()))

                            print("NAF(z) = ", NAF_hamming_weight(NAF(abs(z))))
                            print("HW(z) = ", hamming_weight(vector_z))

                            assert( (r - 1) % 2**m == 0)
                            assert( (re - 1) % 2**m == 0)

                            print("valuation(z + 1, 2) = ", (z + 1).valuation(2))
                            print("valuation(r - 1, 2) = ", (r - 1).valuation(2))
                            print("valuation(re - 1, 2) = ", (re - 1).valuation(2))

                            print("-----------------------")
            else: 
                break

if __name__ == "__main__":
    arithmetic(False)
    search_bls12_curve()