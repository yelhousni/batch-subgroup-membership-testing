package bls12377

import (
	"crypto/rand"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	curve "github.com/consensys/gnark-crypto/ecc/bls12-377"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fr"
)

// IsInSubGroupBatchNaive checks if a batch of points P_i are in G1.
// This is a naive method that checks each point individually using Scott test
// [Scott21].
//
// [Scott21]: https://eprint.iacr.org/2021/1130.pdf
func IsInSubGroupBatchNaive(points []curve.G1Affine) bool {
	for i := range points {
		if !points[i].IsInSubGroup() {
			return false
		}
	}
	return true
}

// IsInSubGroupBatch checks if a batch of points P_i are in G1.
// It generates random scalars s_i in the range [0, bound), performs
// n=rounds multi-scalar-multiplication Sj=∑[s_i]P_i of sizes N=len(points) and
// checks if Sj are on E[r] using Scott test [Scott21].
//
// [Scott21]: https://eprint.iacr.org/2021/1130.pdf
func IsInSubGroupBatch(points []curve.G1Affine, bound *big.Int, rounds int) bool {

	// Check Sj are on E[r]
	for i := 0; i < rounds; i++ {
		b, err := rand.Int(rand.Reader, bound)
		if err != nil {
			panic(err)
		}
		randoms := make([]fr.Element, len(points))
		for j := range randoms {
			randoms[j].SetBigInt(b)
		}
		var sum curve.G1Jac
		sum.MultiExp(points[:], randoms[:], ecc.MultiExpConfig{})
		if !sum.IsInSubGroup() {
			return false
		}
	}
	return true
}
