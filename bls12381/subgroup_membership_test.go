package bls12381

import (
	"fmt"
	"math/big"
	"testing"

	curve "github.com/consensys/gnark-crypto/ecc/bls12-381"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fp"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fr"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/hash_to_curve"

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
)

// Let h be the cofactor of (E/𝔽p) and let e=3√(h/3).
// bound=10177 is the smallest prime divisor of e'=e/gcd(π,e)
// where π= 2⁴·3²·5·7·11·13.
// For a failure probability of 2⁻ᵝ we need to set rounds=⌈β/log2(bound)⌉.
// For example β=64 gives rounds=5 and β=128 gives rounds=10.
var bound = big.NewInt(10177)
var rounds = 5

func TestIsInSubGroupBatch(t *testing.T) {
	t.Parallel()
	parameters := gopter.DefaultTestParameters()
	if testing.Short() {
		parameters.MinSuccessfulTests = 1
	} else {
		parameters.MinSuccessfulTests = 100
	}

	properties := gopter.NewProperties(parameters)

	// number of points to test
	const nbSamples = 100

	properties.Property("[BLS12-381] IsInSubGroupBatchNaive test should pass", prop.ForAll(
		func(mixer fr.Element) bool {
			// mixer ensures that all the words of a frElement are set
			var sampleScalars [nbSamples]fr.Element

			for i := 1; i <= nbSamples; i++ {
				sampleScalars[i-1].SetUint64(uint64(i)).
					Mul(&sampleScalars[i-1], &mixer)
			}

			_, _, g, _ := curve.Generators()
			result := curve.BatchScalarMultiplicationG1(&g, sampleScalars[:])

			return IsInSubGroupBatchNaive(result)
		},
		GenFr(),
	))

	properties.Property("[BLS12-381] IsInSubGroupBatchNaive test should not pass", prop.ForAll(
		func(mixer fr.Element, a fp.Element) bool {
			// mixer ensures that all the words of a frElement are set
			var sampleScalars [nbSamples]fr.Element

			for i := 1; i <= nbSamples; i++ {
				sampleScalars[i-1].SetUint64(uint64(i)).
					Mul(&sampleScalars[i-1], &mixer)
			}

			// random points in G1
			_, _, g, _ := curve.Generators()
			result := curve.BatchScalarMultiplicationG1(&g, sampleScalars[:])
			// random points in the h-torsion
			h := fuzzCofactorOfG1(a)
			result[0].FromJacobian(&h)
			h = fuzzCofactorOfG1(a)
			result[nbSamples-1].FromJacobian(&h)

			return !IsInSubGroupBatchNaive(result)
		},
		GenFr(),
		GenFp(),
	))

	properties.Property("[BLS12-381] IsInSubGroupBatch test should pass with probability 1-1/2^64", prop.ForAll(
		func(mixer fr.Element) bool {
			// mixer ensures that all the words of a frElement are set
			var sampleScalars [nbSamples]fr.Element

			for i := 1; i <= nbSamples; i++ {
				sampleScalars[i-1].SetUint64(uint64(i)).
					Mul(&sampleScalars[i-1], &mixer)
			}

			// random points in G1
			_, _, g, _ := curve.Generators()
			result := curve.BatchScalarMultiplicationG1(&g, sampleScalars[:])

			return IsInSubGroupBatch(result, rounds)
		},
		GenFr(),
	))

	properties.Property("[BLS12-381] IsInSubGroupBatch test should not pass", prop.ForAll(
		func(mixer fr.Element, a fp.Element) bool {
			// mixer ensures that all the words of a frElement are set
			var sampleScalars [nbSamples]fr.Element

			for i := 1; i <= nbSamples; i++ {
				sampleScalars[i-1].SetUint64(uint64(i)).
					Mul(&sampleScalars[i-1], &mixer)
			}

			// random points in G1
			_, _, g, _ := curve.Generators()
			result := curve.BatchScalarMultiplicationG1(&g, sampleScalars[:])

			// random points in the h-torsion
			h := fuzzCofactorOfG1(a)
			result[0].FromJacobian(&h)
			h = fuzzCofactorOfG1(a)
			result[nbSamples-1].FromJacobian(&h)

			return !IsInSubGroupBatch(result, rounds)
		},
		GenFr(),
		GenFp(),
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestTatePairings(t *testing.T) {
	t.Parallel()
	parameters := gopter.DefaultTestParameters()
	parameters.MinSuccessfulTests = 1

	properties := gopter.NewProperties(parameters)

	properties.Property("[BLS12-381] Tate(P3,Q) should be 1", prop.ForAll(
		func(a fr.Element) bool {
			var s big.Int
			a.BigInt(&s)
			_, _, g, _ := curve.Generators()
			g.ScalarMultiplication(&g, &s)
			return isFirstTateOne(g)
		},
		GenFr(),
	))

	properties.Property("[BLS12-381] Tate(P11,Q) should be 1", prop.ForAll(
		func(a fr.Element) bool {
			var s big.Int
			a.BigInt(&s)
			_, _, g, _ := curve.Generators()
			g.ScalarMultiplication(&g, &s)
			return isSecondTateOne(g)
		},
		GenFr(),
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

// benches
func BenchmarkIsInSubGroupBatchNaiveShort(b *testing.B) {
	const nbSamples = 100

	// mixer ensures that all the words of a frElement are set
	var mixer fr.Element
	mixer.SetRandom()
	var sampleScalars [nbSamples]fr.Element

	for i := 1; i <= nbSamples; i++ {
		sampleScalars[i-1].SetUint64(uint64(i)).
			Mul(&sampleScalars[i-1], &mixer)
	}

	_, _, g, _ := curve.Generators()
	result := curve.BatchScalarMultiplicationG1(&g, sampleScalars[:])
	b.ResetTimer()
	for j := 0; j < b.N; j++ {
		IsInSubGroupBatchNaive(result[:])
	}
}

func BenchmarkIsInSubGroupBatchShort(b *testing.B) {
	const nbSamples = 100

	// mixer ensures that all the words of a frElement are set
	var mixer fr.Element
	mixer.SetRandom()
	var sampleScalars [nbSamples]fr.Element

	for i := 1; i <= nbSamples; i++ {
		sampleScalars[i-1].SetUint64(uint64(i)).
			Mul(&sampleScalars[i-1], &mixer)
	}

	_, _, g, _ := curve.Generators()
	result := curve.BatchScalarMultiplicationG1(&g, sampleScalars[:])
	b.ResetTimer()
	for j := 0; j < b.N; j++ {
		IsInSubGroupBatch(result[:], rounds)
	}
}

func BenchmarkIsInSubGroupBatchNaive(b *testing.B) {
	const nbSamples = 1000000

	// mixer ensures that all the words of a frElement are set
	var mixer fr.Element
	mixer.SetRandom()
	var sampleScalars [nbSamples]fr.Element

	for i := 1; i <= nbSamples; i++ {
		sampleScalars[i-1].SetUint64(uint64(i)).
			Mul(&sampleScalars[i-1], &mixer)
	}

	_, _, g, _ := curve.Generators()
	result := curve.BatchScalarMultiplicationG1(&g, sampleScalars[:])

	for _, i := range []uint64{10, 100, 1000, 10000, 100000, 1000000} {
		b.Run(fmt.Sprintf("%d points", i), func(b *testing.B) {
			b.ResetTimer()
			for j := 0; j < b.N; j++ {
				IsInSubGroupBatchNaive(result[:i])
			}
		})
	}
}

func BenchmarkIsInSubGroupBatch(b *testing.B) {
	const nbSamples = 1000000

	// mixer ensures that all the words of a frElement are set
	var mixer fr.Element
	mixer.SetRandom()
	var sampleScalars [nbSamples]fr.Element

	for i := 1; i <= nbSamples; i++ {
		sampleScalars[i-1].SetUint64(uint64(i)).
			Mul(&sampleScalars[i-1], &mixer)
	}

	_, _, g, _ := curve.Generators()
	result := curve.BatchScalarMultiplicationG1(&g, sampleScalars[:])

	for _, i := range []uint64{10, 100, 1000, 10000, 100000, 1000000} {
		b.Run(fmt.Sprintf("%d points", i), func(b *testing.B) {
			b.ResetTimer()
			for j := 0; j < b.N; j++ {
				IsInSubGroupBatch(result[:i], rounds)
			}
		})
	}
}

// utils

// phi sets p to ϕ(a) where ϕ: (x,y) → (w x,y),
// where w is a third root of unity.
func phi(q *curve.G1Jac) *curve.G1Jac {
	var p curve.G1Jac
	p.Set(q)
	var w fp.Element
	w.SetString("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436")
	p.X.Mul(&p.X, &w)
	return &p
}

func fuzzCofactorOfG1(f fp.Element) curve.G1Jac {
	var res, jac curve.G1Jac
	aff := curve.MapToCurve1(&f)
	hash_to_curve.G1Isogeny(&aff.X, &aff.Y)
	jac.FromAffine(&aff)
	// p+x²ϕ(p) = [r]p
	res = *phi(&jac)
	res = *mulBySeed(&res)
	res = *mulBySeed(&res)
	res.AddAssign(&jac)
	return res
}

// mulBySeed multiplies the point q by the seed xGen in Jacobian coordinates
// using an optimized addition chain.
func mulBySeed(q *curve.G1Jac) *curve.G1Jac {
	// Generated by github.com/mmcloughlin/addchain v0.4.0.
	// Operations: 63 doublings 5 additions

	var res curve.G1Jac
	res.Double(q)
	res.AddAssign(q)
	for i := 0; i < 2; i++ {
		res.Double(&res)
	}
	res.AddAssign(q)
	for i := 0; i < 3; i++ {
		res.Double(&res)
	}
	res.AddAssign(q)
	for i := 0; i < 9; i++ {
		res.Double(&res)
	}
	res.AddAssign(q)
	for i := 0; i < 32; i++ {
		res.Double(&res)
	}
	res.AddAssign(q)
	for i := 0; i < 16; i++ {
		res.Double(&res)
	}
	return &res
}

// GenFr generates an Fr element
func GenFr() gopter.Gen {
	return func(genParams *gopter.GenParameters) *gopter.GenResult {
		var elmt fr.Element
		elmt.MustSetRandom()

		return gopter.NewGenResult(elmt, gopter.NoShrinker)
	}
}

// GenFp generates an Fp element
func GenFp() gopter.Gen {
	return func(genParams *gopter.GenParameters) *gopter.GenResult {
		var elmt fp.Element
		elmt.MustSetRandom()

		return gopter.NewGenResult(elmt, gopter.NoShrinker)
	}
}
