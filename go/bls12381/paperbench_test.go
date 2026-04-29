package bls12381

import (
	"fmt"
	"testing"

	curve "github.com/consensys/gnark-crypto/ecc/bls12-381"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fr"
)

var paperBenchSizes = [...]int{32, 128, 512, 2048, 8192, 32768, 131072, 524288, 2097152}

func BenchmarkPaperComparison(b *testing.B) {
	const nbSamples = 1 << 21
	var sampleScalars [nbSamples]fr.Element
	fillBenchScalars(sampleScalars[:])

	_, _, g, _ := curve.Generators()
	result := curve.BatchScalarMultiplicationG1(&g, sampleScalars[:])

	for _, using := range paperBenchSizes {
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
		b.Run(fmt.Sprintf("%d points-step2", using), func(b *testing.B) {
			b.ResetTimer()
			for j := 0; j < b.N; j++ {
				_msmCheck(result[:using])
			}
		})
	}
}
