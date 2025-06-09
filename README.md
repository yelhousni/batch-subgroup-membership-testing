# Batch subgroup membership testing

Related code to the research paper "Batch subgroup membership testing on pairing-friendly curves".

*Authors: [Dimitri Koshelev](https://github.com/Dimitri-Koshelev) and [Youssef El Housni](https://github.com/yelhousni).*

## Pre-requisites
Install Golang (https://go.dev/doc/install). This code was tested with the last 2 major releases of Go (currently 1.23 and 1.22).

The code uses [gnark-crypto](https://github.com/ConsenSys/gnark-crypto) Golang library.

## Organization
- `bls12381/` contains the implementation of the new method for the BLS12-381 curve.
- `bls12377/` contains the implementation of the new method for the BLS12-377 curve.
- `bls12377-strong/` contains the implementation of a new BLS12-377 curve with optimal batch SMT.
- `sage/` contains SageMath scripts to validate formulas related to the Tate pairings.

## Benchmarks
On a Macbook Air M1

#### BLS12-381 (ℙ $\leq 2^{-64}$)
| Number of points | Naive method (ms) | This work (ms) | Speedup |
|------------------|-------------------|----------------|---------|
| $10^2$            | 5.055             | 3.955          | **1.277x**   |
| $10^3$            | 50.456            | 33.533         | **1.505x**   |
| $10^4$            | 508.432           | 325.157        | **1.564x**   |
| $10^5$            | 5976.587          | 3971.242       | **1.505x**   |
| $10^6$            | 50363.157         | 31377.842      | **1.605x**   |

#### BLS12-377 (ℙ $\leq 2^{-64}$)
| Number of points | Naive method (ms) | This work (ms) | Speedup |
|------------------|-------------------|----------------|---------|
| $10^2$            | 5.092             | 5.808          | **0.877x**   |
| $10^3$            | 50.657            | 28.114         | **1.802x**   |
| $10^4$            | 508.060           | 254.446        | **1.996x**   |
| $10^5$            | 5137.183          | 2452.670       | **2.094x**   |
| $10^6$            | 51130.979         | 24369.764      | **2.098x**   |

#### Strong BLS12-377 (ℙ $\leq 2^{-64}$)
| Number of points | Naive method (ms) | This work (ms) | Speedup   |
|------------------|-------------------|----------------|-----------|
| $10^2$            | 6.287           | 3.261            | **1.93x** |
| $10^3$            | 63.006          | 30.000           | **2.10x** |
| $10^4$            | 627.667         | 297.217          | **2.11x** |
| $10^5$            | 6339.666        | 2948.324         | **2.15x** |
| $10^6$            | 63122.227       | 29207.600        | **2.16x** |
