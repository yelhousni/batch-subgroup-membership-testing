from sage.all_cmdline import *   # import sage library
#
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.structure.proof.all import arithmetic

# functions for NAF and hamming weight
from _utils import NAF, hamming_weight, NAF_hamming_weight

def search_bls12_curve():
    v = 24 
    n = 63 - v
    print("v = ", v)
    print("-------")

    j = 0
    if (j == 0):
        print("negative seeds: j = ", j)
    else:
        print("positive seeds: j = ", j)
    sign = (-1)**(1 - j)
    for zPr in range(2**n, 2**(n+2)):
        # seeds of the form: z = 2^v * z' - 1 or z = -2^v * z' - 1
        # z + 1 has a large power of 2 as factor, at least 2^v
        # search range ensures that the size of the base field prime q is between 377- and 383-bits
        z = ZZ(sign * zPr * 2**v - 1)
        
        e = abs(z - 1)
        if (e % 3) == 0: 
            r = z**4 - z**2 + 1
            # r' = z^4 - 3z^2 + 3: prime order of embedded curve
            rPr = z**4 - 3*z**2 + 3
            # e1 satisfies h1 = 3*e1^2 where h1 is the cofactor of E/Fq
            e1 = e // 3
            assert (e1 % 2) == 0
            h1 = 3*e1**2
            # the base field prime: q = h1 * r + z
            q = h1 * r + z 
            # p = e1 / 2 must be prime
            p = e1 // 2
            if (q.nbits() < 384):
                if (q.is_pseudoprime() == True) & (r.is_pseudoprime() == True) & (abs(p).is_pseudoprime() == True) & (rPr.is_pseudoprime() == True):
                    # h2: cofactor of order of degree six twist, i.e. #E2(Fq^2) = h2 * r
                    h2 = (z**8 - 4*z**7 + 5*z**6 - 4*z**4 + 6*z**3 - 4*z**2 - 4*z + 13) // 9
                    # hT: cofactor of target group GT
                    hT = (z**20 - 8*z**19 + 25*z**18 - 32*z**17 - 8*z**16 + 76*z**15 - 93*z**14 + 36*z**13 + 51*z**12 - 112*z**11 
                       + 86*z**10 - 16*z**9 - 24*z**8 + 84*z**7 - 90*z**6 + 28*z**5 - 14*z**4 - 38*z**3 + 70*z**2 - 14*z + 73) // 81
                    if (h2.is_pseudoprime() == True):
                        # seeds which satisfy all conditions except primality of hT
                        vector_z = Integer(abs(z)).digits(2)
                        print("z = {:#x} {} bits".format(z, abs(z).nbits()))

                        # check that (r - 1) and (r' - 1) are both divisible by 2^v
                        assert( (r - 1) % 2**v == 0 )
                        assert( (rPr - 1) % 2**v == 0 )

                        if (hT.is_pseudoprime() == True):
                            # seeds which satisfy all conditions including primality of hT
                            vector_z = Integer(abs(z)).digits(2)
                            print("z = {:#x} {} bits".format(z, abs(z).nbits()))
                            print("r = {:#x} {} bits".format(r, r.nbits()))
                            print("rPr = {:#x} {} bits".format(rPr, rPr.nbits()))
                            print("q = {:#x} {} bits".format(q, q.nbits()))
                            print("h2 = {:#x} {} bits".format(h2, h2.nbits()))
                            print("hT = {:#x} {} bits".format(hT, hT.nbits()))

                            print("NAF(z) = ", NAF_hamming_weight(NAF(abs(z))))
                            print("HW(z) = ", hamming_weight(vector_z))

                            assert( (r - 1) % 2**v == 0)
                            assert( (rPr - 1) % 2**v == 0)

                            print("valuation(z + 1, 2) = ", (z + 1).valuation(2))
                            # valuation of r - 1 and r' - 1: v2*(r) and v2*(r')
                            print("valuation(r - 1, 2) = ", (r - 1).valuation(2))
                            print("valuation(rPr - 1, 2) = ", (rPr - 1).valuation(2))

                            print("-----------------------")
            else: 
                break

if __name__ == "__main__":
    arithmetic(False)
    search_bls12_curve()