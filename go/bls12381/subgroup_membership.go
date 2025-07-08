package bls12381

import (
	"sync/atomic"

	curve "github.com/consensys/gnark-crypto/ecc/bls12-381"
	"github.com/yelhousni/batch-subgroup-membership/go/parallel"
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

func IsInSubGroupBatchNaiveParallel(points []curve.G1Affine) bool {
	var nbErrors int64
	parallel.Execute(len(points), func(start, end int) {
		for i := start; i < end; i++ {
			if !points[i].IsInSubGroup() {
				atomic.AddInt64(&nbErrors, 1)
				return
			}
		}
	})
	return nbErrors == 0
}

// IsInSubGroupBatch checks if a batch of points P_i are in G1.
// First, it checks that all points are on a larger torsion E[r*e'] using Tate
// pairings [Koshelev22].
// Second, it generates random scalars s_i in the range [0, bound), performs
// n=rounds multi-scalar-multiplication Sj=∑[s_i]P_i of sizes N=len(points) and
// checks if Sj are on E[r] using Scott test [Scott21].
//
// [Koshelev22]: https://eprint.iacr.org/2022/037.pdf
// [Scott21]: https://eprint.iacr.org/2021/1130.pdf
func IsInSubGroupBatch(points []curve.G1Affine, rounds int) bool {

	// 1. Check points are on E[r*e']
	for i := range points {
		// 1.1. Tate_{3,P3}(Q) = (y-2)^((p-1)/3) == 1, with P3 = (0,2).
		if !isFirstTateOne(points[i]) {
			return false
		}
		// 1.2. Tate_{11,P11}(Q) == 1
		if !isSecondTateOne(points[i]) {
			return false
		}
	}

	// 2. Check Sj are on E[r]
	const nbRounds = 5
	for i := 0; i < nbRounds; i++ {
		if !_msmCheck(points) {
			return false
		}
	}

	return true
}

func IsInSubGroupBatchParallel(points []curve.G1Affine, rounds int) bool {

	// 1. Check points are on E[r*e']
	var nbErrors int64
	parallel.Execute(len(points), func(start, end int) {
		for i := start; i < end; i++ {
			// 1.1. Tate_{3,P3}(Q) = (y-2)^((p-1)/3) == 1, with P3 = (0,2).
			if !isFirstTateOne(points[i]) {
				atomic.AddInt64(&nbErrors, 1)
				return
			}
			// 1.2. Tate_{11,P11}(Q) == 1
			if !isSecondTateOne(points[i]) {
				atomic.AddInt64(&nbErrors, 1)
				return
			}
		}
	})
	if nbErrors > 0 {
		return false
	}

	// 2. Check Sj are on E[r]
	const nbRounds = 5
	parallel.Execute(nbRounds, func(start, end int) {
		for i := start; i < end; i++ {
			if !_msmCheck(points) {
				atomic.AddInt64(&nbErrors, 1)
				return
			}
		}
	})

	return nbErrors == 0
}
