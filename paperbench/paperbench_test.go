package paperbench

import (
	"math/big"
	"testing"

	curve381 "github.com/consensys/gnark-crypto/ecc/bls12-381"
	fr381 "github.com/consensys/gnark-crypto/ecc/bls12-381/fr"
	curve376 "github.com/yelhousni/batch-subgroup-membership/go/bls12376-strong"
	fr376 "github.com/yelhousni/batch-subgroup-membership/go/bls12376-strong/fr"
	curve377 "github.com/yelhousni/batch-subgroup-membership/go/bls12377-strong"
	fr377 "github.com/yelhousni/batch-subgroup-membership/go/bls12377-strong/fr"
)

func fill381(dst []fr381.Element) {
	var mixer fr381.Element
	mixer.SetString("7716837800905789770901243404444209691916730933998574719964609384059111546487")
	for i := 1; i <= len(dst); i++ {
		dst[i-1].SetUint64(uint64(i)).Mul(&dst[i-1], &mixer)
	}
}

func fill377(dst []fr377.Element) {
	var mixer fr377.Element
	mixer.SetString("7716837800905789770901243404444209691916730933998574719964609384059111546487")
	for i := 1; i <= len(dst); i++ {
		dst[i-1].SetUint64(uint64(i)).Mul(&dst[i-1], &mixer)
	}
}

func fill376(dst []fr376.Element) {
	var mixer fr376.Element
	mixer.SetString("7716837800905789770901243404444209691916730933998574719964609384059111546487")
	for i := 1; i <= len(dst); i++ {
		dst[i-1].SetUint64(uint64(i)).Mul(&dst[i-1], &mixer)
	}
}

func bench381G1Scalar(b *testing.B) {
	var scalar big.Int
	_, _, g1, _ := curve381.Generators()
	scalar.SetString("5243587517512619047944770508185965837690552500527637822603658699938581184513", 10)
	scalar.Add(&scalar, fr381.Modulus())
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		var out curve381.G1Affine
		out.ScalarMultiplication(&g1, &scalar)
	}
}

func bench381G2Scalar(b *testing.B) {
	var scalar big.Int
	_, _, _, g2 := curve381.Generators()
	scalar.SetString("5243587517512619047944770508185965837690552500527637822603658699938581184513", 10)
	scalar.Add(&scalar, fr381.Modulus())
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		var out curve381.G2Affine
		out.ScalarMultiplication(&g2, &scalar)
	}
}

func bench381G1MSM(b *testing.B, n int) {
	_, _, g1, _ := curve381.Generators()
	scalars := make([]fr381.Element, n)
	fill381(scalars)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = curve381.BatchScalarMultiplicationG1(&g1, scalars)
	}
}

func bench381G2MSM(b *testing.B, n int) {
	_, _, _, g2 := curve381.Generators()
	scalars := make([]fr381.Element, n)
	fill381(scalars)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = curve381.BatchScalarMultiplicationG2(&g2, scalars)
	}
}

func bench381Pairing(b *testing.B) {
	_, _, g1, g2 := curve381.Generators()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = curve381.Pair([]curve381.G1Affine{g1}, []curve381.G2Affine{g2})
	}
}

func bench381PairingProduct(b *testing.B, n int) {
	_, _, g1, g2 := curve381.Generators()
	p := make([]curve381.G1Affine, n)
	q := make([]curve381.G2Affine, n)
	for i := 0; i < n; i++ {
		p[i] = g1
		q[i] = g2
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = curve381.Pair(p, q)
	}
}

func bench377G1Scalar(b *testing.B) {
	var scalar big.Int
	_, _, g1, _ := curve377.Generators()
	scalar.SetString("5243587517512619047944770508185965837690552500527637822603658699938581184513", 10)
	scalar.Add(&scalar, fr377.Modulus())
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		var out curve377.G1Affine
		out.ScalarMultiplication(&g1, &scalar)
	}
}

func bench377G2Scalar(b *testing.B) {
	var scalar big.Int
	_, _, _, g2 := curve377.Generators()
	scalar.SetString("5243587517512619047944770508185965837690552500527637822603658699938581184513", 10)
	scalar.Add(&scalar, fr377.Modulus())
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		var out curve377.G2Affine
		out.ScalarMultiplication(&g2, &scalar)
	}
}

func bench377G1MSM(b *testing.B, n int) {
	_, _, g1, _ := curve377.Generators()
	scalars := make([]fr377.Element, n)
	fill377(scalars)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = curve377.BatchScalarMultiplicationG1(&g1, scalars)
	}
}

func bench377G2MSM(b *testing.B, n int) {
	_, _, _, g2 := curve377.Generators()
	scalars := make([]fr377.Element, n)
	fill377(scalars)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = curve377.BatchScalarMultiplicationG2(&g2, scalars)
	}
}

func bench377Pairing(b *testing.B) {
	_, _, g1, g2 := curve377.Generators()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = curve377.Pair([]curve377.G1Affine{g1}, []curve377.G2Affine{g2})
	}
}

func bench377PairingProduct(b *testing.B, n int) {
	_, _, g1, g2 := curve377.Generators()
	p := make([]curve377.G1Affine, n)
	q := make([]curve377.G2Affine, n)
	for i := 0; i < n; i++ {
		p[i] = g1
		q[i] = g2
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = curve377.Pair(p, q)
	}
}

func bench376G1Scalar(b *testing.B) {
	var scalar big.Int
	_, _, g1, _ := curve376.Generators()
	scalar.SetString("5243587517512619047944770508185965837690552500527637822603658699938581184513", 10)
	scalar.Add(&scalar, fr376.Modulus())
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		var out curve376.G1Affine
		out.ScalarMultiplication(&g1, &scalar)
	}
}

func bench376G2Scalar(b *testing.B) {
	var scalar big.Int
	_, _, _, g2 := curve376.Generators()
	scalar.SetString("5243587517512619047944770508185965837690552500527637822603658699938581184513", 10)
	scalar.Add(&scalar, fr376.Modulus())
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		var out curve376.G2Affine
		out.ScalarMultiplication(&g2, &scalar)
	}
}

func bench376G1MSM(b *testing.B, n int) {
	_, _, g1, _ := curve376.Generators()
	scalars := make([]fr376.Element, n)
	fill376(scalars)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = curve376.BatchScalarMultiplicationG1(&g1, scalars)
	}
}

func bench376G2MSM(b *testing.B, n int) {
	_, _, _, g2 := curve376.Generators()
	scalars := make([]fr376.Element, n)
	fill376(scalars)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = curve376.BatchScalarMultiplicationG2(&g2, scalars)
	}
}

func bench376Pairing(b *testing.B) {
	_, _, g1, g2 := curve376.Generators()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = curve376.Pair([]curve376.G1Affine{g1}, []curve376.G2Affine{g2})
	}
}

func bench376PairingProduct(b *testing.B, n int) {
	_, _, g1, g2 := curve376.Generators()
	p := make([]curve376.G1Affine, n)
	q := make([]curve376.G2Affine, n)
	for i := 0; i < n; i++ {
		p[i] = g1
		q[i] = g2
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = curve376.Pair(p, q)
	}
}

func BenchmarkBLS12381G1Scalar(b *testing.B)        { bench381G1Scalar(b) }
func BenchmarkBLS12381G2Scalar(b *testing.B)        { bench381G2Scalar(b) }
func BenchmarkBLS12381G1MSM32(b *testing.B)         { bench381G1MSM(b, 1<<5) }
func BenchmarkBLS12381G1MSM2097152(b *testing.B)    { bench381G1MSM(b, 1<<21) }
func BenchmarkBLS12381G2MSM32(b *testing.B)         { bench381G2MSM(b, 1<<5) }
func BenchmarkBLS12381G2MSM2097152(b *testing.B)    { bench381G2MSM(b, 1<<21) }
func BenchmarkBLS12381Pairing(b *testing.B)         { bench381Pairing(b) }
func BenchmarkBLS12381PairingProduct2(b *testing.B) { bench381PairingProduct(b, 2) }
func BenchmarkBLS12381PairingProduct4(b *testing.B) { bench381PairingProduct(b, 4) }
func BenchmarkBLS12381PairingProduct10(b *testing.B) {
	bench381PairingProduct(b, 10)
}

func BenchmarkBLS12377StrongG1Scalar(b *testing.B)        { bench377G1Scalar(b) }
func BenchmarkBLS12377StrongG2Scalar(b *testing.B)        { bench377G2Scalar(b) }
func BenchmarkBLS12377StrongG1MSM32(b *testing.B)         { bench377G1MSM(b, 1<<5) }
func BenchmarkBLS12377StrongG1MSM2097152(b *testing.B)    { bench377G1MSM(b, 1<<21) }
func BenchmarkBLS12377StrongG2MSM32(b *testing.B)         { bench377G2MSM(b, 1<<5) }
func BenchmarkBLS12377StrongG2MSM2097152(b *testing.B)    { bench377G2MSM(b, 1<<21) }
func BenchmarkBLS12377StrongPairing(b *testing.B)         { bench377Pairing(b) }
func BenchmarkBLS12377StrongPairingProduct2(b *testing.B) { bench377PairingProduct(b, 2) }
func BenchmarkBLS12377StrongPairingProduct4(b *testing.B) { bench377PairingProduct(b, 4) }
func BenchmarkBLS12377StrongPairingProduct10(b *testing.B) {
	bench377PairingProduct(b, 10)
}

func BenchmarkBLS12376StrongG1Scalar(b *testing.B)        { bench376G1Scalar(b) }
func BenchmarkBLS12376StrongG2Scalar(b *testing.B)        { bench376G2Scalar(b) }
func BenchmarkBLS12376StrongG1MSM32(b *testing.B)         { bench376G1MSM(b, 1<<5) }
func BenchmarkBLS12376StrongG1MSM2097152(b *testing.B)    { bench376G1MSM(b, 1<<21) }
func BenchmarkBLS12376StrongG2MSM32(b *testing.B)         { bench376G2MSM(b, 1<<5) }
func BenchmarkBLS12376StrongG2MSM2097152(b *testing.B)    { bench376G2MSM(b, 1<<21) }
func BenchmarkBLS12376StrongPairing(b *testing.B)         { bench376Pairing(b) }
func BenchmarkBLS12376StrongPairingProduct2(b *testing.B) { bench376PairingProduct(b, 2) }
func BenchmarkBLS12376StrongPairingProduct4(b *testing.B) { bench376PairingProduct(b, 4) }
func BenchmarkBLS12376StrongPairingProduct10(b *testing.B) {
	bench376PairingProduct(b, 10)
}
