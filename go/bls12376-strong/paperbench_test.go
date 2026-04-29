package bls12376strong

import (
	"fmt"
	"testing"

	"github.com/yelhousni/batch-subgroup-membership/go/bls12376-strong/fr"
)

var paperBenchSizes = [...]int{32, 128, 512, 2048, 8192, 32768, 131072, 524288, 2097152}

func BenchmarkPaperComparison(b *testing.B) {
	const nbSamples = 1 << 21
	var sampleScalars [nbSamples]fr.Element
	fillBenchScalars(sampleScalars[:])

	result := BatchScalarMultiplicationG1(&g1GenAff, sampleScalars[:])

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
	}
}
