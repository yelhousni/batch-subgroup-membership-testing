# Batch subgroup membership testing

Related code to the research paper "Batch subgroup membership testing on pairing-friendly curves".

*Authors: [Dimitri Koshelev](https://github.com/Dimitri-Koshelev), [Youssef El Housni](https://github.com/yelhousni) and [Georgios Fotiadis](https://github.com/gfotiadis).*

## Pre-requisites
Install Golang (https://go.dev/doc/install). This code was tested with the releases 1.23 and 1.22.

The code uses [gnark-crypto](https://github.com/ConsenSys/gnark-crypto) Golang library.

## Organization
- `bls12381/` contains the implementation of the new method for the BLS12-381 curve.
- `bls12377/` contains the implementation of the new method for the BLS12-377 curve.
- `bls12377-strong/` contains the implementation of a new BLS12-377 curve. It is G2-strong, GT-strong and optimal for batch SMT.
- `bls12376-strong/` contains the implementation of a new BLS12-376 curve. It is G2-strong and optimal for batch SMT.
- `sage/` contains SageMath scripts to validate formulas related to the elliptic curves and the Tate pairings. It also contains the scripts that were used to find the new BLS12 curves.
- `magma/` contains a MAGMA script to validate formulas and tables related to section 2.
