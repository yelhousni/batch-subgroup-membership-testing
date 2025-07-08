import os
import time
from sage.all import Integer
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF


def NAF(x: int):
    """
    :param x: integer x
    :return: list corresponding to the NAF representation of x
    """
    t0 = time.time()

    naf_x = []
    xx = Integer(x)
    assert x >= 0
    while xx > 0:
        rr = xx % 4
        if rr == 3:
            rr = -1
        else:
            rr = rr % 2
        naf_x.append(rr)
        xx -= rr
        xx, rr = xx.quo_rem(2)
        assert rr == 0
    assert x == sum([r * 2 ** i for i, r in enumerate(naf_x)])

    return naf_x


def hamming_weight(bit_x: list) -> int:
    """
    :param bit_x: binary representartion of a positive integer x
    :return: Hamming weight of x
    """
    count = 0
    for i in range(len(bit_x)):
        if bit_x[i] == 1:
            count = count + 1

    return count


def NAF_hamming_weight(naf_x: list) -> int:
    """
    :param naf_x: NAF representation of a positive integer x
    :return: NAF Hamming weight of x
    """
    count = 0
    for i in range(0, len(naf_x)):
        if (naf_x[i] == 1) or (naf_x[i] == -1):
            count = count + 1

    return count
