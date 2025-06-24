package bls12377

import (
	"fmt"
	"testing"

	curve "github.com/consensys/gnark-crypto/ecc/bls12-377"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fp"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fr"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/hash_to_curve"

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
)

// For highly 2-adic curves the bound is always 2.
// For a failure probability of 2⁻ᵝ we need to set rounds=β.
// For example β=64 gives rounds=64 and β=128 gives rounds=128.
var rounds = 64

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

	properties.Property("[BLS12-377] IsInSubGroupBatchNaive test should pass", prop.ForAll(
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

			return IsInSubGroupBatchNaive(result)
		},
		GenFr(),
	))

	properties.Property("[BLS12-377] IsInSubGroupBatchNaive test should not pass", prop.ForAll(
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

	properties.Property("[BLS12-377] IsInSubGroupBatch test should pass with probability 1-1/2^64", prop.ForAll(
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

	properties.Property("[BLS12-377] IsInSubGroupBatch test should not pass", prop.ForAll(
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
			// random point in the h-torsion
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

func BenchmarkComparison(b *testing.B) {
	const (
		pow       = 22
		nbSamples = 1 << pow
	)
	var sampleScalars [nbSamples]fr.Element
	fillBenchScalars(sampleScalars[:])

	_, _, g, _ := curve.Generators()
	result := curve.BatchScalarMultiplicationG1(&g, sampleScalars[:])

	for i := 5; i <= pow; i++ {
		using := 1 << i
		b.Run(fmt.Sprintf("%d points-naive", using), func(b *testing.B) {
			b.ResetTimer()
			for j := 0; j < b.N; j++ {
				IsInSubGroupBatchNaive(result[:using])
			}
		})
		b.Run(fmt.Sprintf("%d points", using), func(b *testing.B) {
			b.ResetTimer()
			for j := 0; j < b.N; j++ {
				IsInSubGroupBatch(result[:using], rounds)
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
	w.SetString("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945")
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
	// Operations: 61 doublings 7 additions

	// Allocate Temporaries.
	var res, t0, t1 curve.G1Jac
	res.Double(q)
	res.AddAssign(q)
	res.Double(&res)
	res.AddAssign(q)
	t0.Double(&res)
	for i := 1; i < 2; i++ {
		t0.Double(&t0)
	}
	res.AddAssign(&t0)
	t1.Double(&res)
	t1.AddAssign(&res)
	t0.AddAssign(&t1)
	for i := 0; i < 10; i++ {
		t0.Double(&t0)
	}
	res.AddAssign(&t0)
	for i := 0; i < 46; i++ {
		res.Double(&res)
	}
	res.AddAssign(q)
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

func fillBenchScalars(sampleScalars []fr.Element) {
	// ensure every words of the scalars are filled
	for i := 0; i < len(sampleScalars); i++ {
		sampleScalars[i].MustSetRandom()
	}
}
