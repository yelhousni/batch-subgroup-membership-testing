package bls12376strong

import (
	"math/big"
	"math/bits"
	"sync"

	"github.com/yelhousni/batch-subgroup-membership/go/bls12376-strong/fp"
)

// 3⁻¹ mod 2⁶⁴, used for exact division by 3 via multiplication.
const inv3mod264 = 0xAAAAAAAAAAAAAAAB

// signed256 represents a signed 256-bit integer using 4 uint64 words.
type signed256 struct {
	w0, w1, w2, w3 uint64
	neg            bool
}

func (s signed256) isZero() bool {
	return s.w0 == 0 && s.w1 == 0 && s.w2 == 0 && s.w3 == 0
}

func neg256(s signed256) signed256 {
	if s.isZero() {
		return s
	}
	return signed256{s.w0, s.w1, s.w2, s.w3, !s.neg}
}

func cmpAbs256(a, b signed256) int {
	if a.w3 != b.w3 {
		if a.w3 > b.w3 {
			return 1
		}
		return -1
	}
	if a.w2 != b.w2 {
		if a.w2 > b.w2 {
			return 1
		}
		return -1
	}
	if a.w1 != b.w1 {
		if a.w1 > b.w1 {
			return 1
		}
		return -1
	}
	if a.w0 != b.w0 {
		if a.w0 > b.w0 {
			return 1
		}
		return -1
	}
	return 0
}

func add256(a, b signed256) signed256 {
	if a.neg == b.neg {
		w0, c := bits.Add64(a.w0, b.w0, 0)
		w1, c := bits.Add64(a.w1, b.w1, c)
		w2, c := bits.Add64(a.w2, b.w2, c)
		w3, _ := bits.Add64(a.w3, b.w3, c)
		return signed256{w0, w1, w2, w3, a.neg}
	}
	if cmpAbs256(a, b) >= 0 {
		w0, bw := bits.Sub64(a.w0, b.w0, 0)
		w1, bw := bits.Sub64(a.w1, b.w1, bw)
		w2, bw := bits.Sub64(a.w2, b.w2, bw)
		w3, _ := bits.Sub64(a.w3, b.w3, bw)
		return signed256{w0, w1, w2, w3, a.neg}
	}
	w0, bw := bits.Sub64(b.w0, a.w0, 0)
	w1, bw := bits.Sub64(b.w1, a.w1, bw)
	w2, bw := bits.Sub64(b.w2, a.w2, bw)
	w3, _ := bits.Sub64(b.w3, a.w3, bw)
	return signed256{w0, w1, w2, w3, b.neg}
}

func sub256(a, b signed256) signed256 { return add256(a, neg256(b)) }

func mulSmall256(s signed256, k int64) signed256 {
	neg := s.neg
	if k < 0 {
		neg = !neg
		k = -k
	}
	if k == 0 {
		return signed256{}
	}
	uk := uint64(k)
	h0, l0 := bits.Mul64(s.w0, uk)
	h1, l1 := bits.Mul64(s.w1, uk)
	h2, l2 := bits.Mul64(s.w2, uk)
	_, l3 := bits.Mul64(s.w3, uk)

	w0 := l0
	w1, c := bits.Add64(l1, h0, 0)
	w2, c := bits.Add64(l2, h1, c)
	w3, _ := bits.Add64(l3, h2, c)

	return signed256{w0, w1, w2, w3, neg}
}

func s256ToFloat(s signed256) float64 {
	f := float64(s.w3)*0x1p192 + float64(s.w2)*0x1p128 + float64(s.w1)*0x1p64 + float64(s.w0)
	if s.neg {
		f = -f
	}
	return f
}

// mod3_256 returns s mod 3 in [0,2]. Since 2^64 ≡ 1 mod 3, sum all words mod 3.
func mod3_256(s signed256) uint64 {
	v := (s.w0%3 + s.w1%3 + s.w2%3 + s.w3%3) % 3
	if s.neg && v != 0 {
		v = 3 - v
	}
	return v
}

// mod9_256 returns s mod 9 in [0,8].
// 2^64 ≡ 7 mod 9, 2^128 ≡ 4 mod 9, 2^192 ≡ 1 mod 9.
func mod9_256(s signed256) uint64 {
	v := (s.w0%9 + 7*(s.w1%9) + 4*(s.w2%9) + s.w3%9) % 9
	if s.neg && v != 0 {
		v = 9 - v
	}
	return v
}

// divExact3_256 divides s by 3 exactly using multiplicative inverse.
func divExact3_256(s signed256) signed256 {
	q0 := s.w0 * inv3mod264
	c0, _ := bits.Mul64(q0, 3)

	sub1, borrow1 := bits.Sub64(s.w1, c0, 0)
	q1 := sub1 * inv3mod264
	c1, _ := bits.Mul64(q1, 3)
	c1 += borrow1

	sub2, borrow2 := bits.Sub64(s.w2, c1, 0)
	q2 := sub2 * inv3mod264
	c2, _ := bits.Mul64(q2, 3)
	c2 += borrow2

	sub3, _ := bits.Sub64(s.w3, c2, 0)
	q3 := sub3 * inv3mod264

	return signed256{q0, q1, q2, q3, s.neg}
}

func bigToS256(s *signed256, x *big.Int) {
	s.neg = x.Sign() < 0
	w := x.Bits()
	s.w0, s.w1, s.w2, s.w3 = 0, 0, 0, 0
	if len(w) > 0 {
		s.w0 = uint64(w[0])
	}
	if len(w) > 1 {
		s.w1 = uint64(w[1])
	}
	if len(w) > 2 {
		s.w2 = uint64(w[2])
	}
	if len(w) > 3 {
		s.w3 = uint64(w[3])
	}
}

// --- signed128 type and arithmetic (for later GCD iterations) ---

type signed128 struct {
	lo, hi uint64
	neg    bool
}

func (s signed128) isZero() bool { return s.lo == 0 && s.hi == 0 }

func neg128(s signed128) signed128 {
	if s.isZero() {
		return s
	}
	return signed128{s.lo, s.hi, !s.neg}
}

func cmpAbs128(a, b signed128) int {
	if a.hi != b.hi {
		if a.hi > b.hi {
			return 1
		}
		return -1
	}
	if a.lo != b.lo {
		if a.lo > b.lo {
			return 1
		}
		return -1
	}
	return 0
}

func add128(a, b signed128) signed128 {
	if a.neg == b.neg {
		lo, c := bits.Add64(a.lo, b.lo, 0)
		hi, _ := bits.Add64(a.hi, b.hi, c)
		return signed128{lo, hi, a.neg}
	}
	if cmpAbs128(a, b) >= 0 {
		lo, bw := bits.Sub64(a.lo, b.lo, 0)
		hi, _ := bits.Sub64(a.hi, b.hi, bw)
		return signed128{lo, hi, a.neg}
	}
	lo, bw := bits.Sub64(b.lo, a.lo, 0)
	hi, _ := bits.Sub64(b.hi, a.hi, bw)
	return signed128{lo, hi, b.neg}
}

func sub128(a, b signed128) signed128 { return add128(a, neg128(b)) }

func mulSmall128(s signed128, k int64) signed128 {
	neg := s.neg
	if k < 0 {
		neg = !neg
		k = -k
	}
	if k == 0 {
		return signed128{}
	}
	uk := uint64(k)
	hi, lo := bits.Mul64(s.lo, uk)
	hi += s.hi * uk
	return signed128{lo, hi, neg}
}

func s128ToFloat(s signed128) float64 {
	f := float64(s.hi)*0x1p64 + float64(s.lo)
	if s.neg {
		f = -f
	}
	return f
}

// mod3_128 returns s mod 3 in [0,2]. Since 2^64 ≡ 1 mod 3, sum words mod 3.
func mod3_128(s signed128) uint64 {
	v := (s.lo%3 + s.hi%3) % 3
	if s.neg && v != 0 {
		v = 3 - v
	}
	return v
}

// mod9_128 returns s mod 9 in [0,8]. 2^64 ≡ 7 mod 9.
func mod9_128(s signed128) uint64 {
	v := (s.lo%9 + 7*(s.hi%9)) % 9
	if s.neg && v != 0 {
		v = 9 - v
	}
	return v
}

// divExact3_128 divides s by 3 exactly using multiplicative inverse.
func divExact3_128(s signed128) signed128 {
	q0 := s.lo * inv3mod264
	c0, _ := bits.Mul64(q0, 3)
	sub1, _ := bits.Sub64(s.hi, c0, 0)
	q1 := sub1 * inv3mod264
	return signed128{q0, q1, s.neg}
}

func isRealUnit128(re, im signed128) bool {
	return im.isZero() && re.hi == 0 && re.lo <= 1
}

// Eisenstein arithmetic with signed128

func eisQuotient128(a0, a1, b0, b1 signed128) (int64, int64) {
	af0 := s128ToFloat(a0)
	af1 := s128ToFloat(a1)
	bf0 := s128ToFloat(b0)
	bf1 := s128ToFloat(b1)

	nB := bf0*bf0 + bf1*bf1 - bf0*bf1
	if nB == 0 {
		return 0, 0
	}

	numRe := af0*bf0 - af0*bf1 + af1*bf1
	numIm := af1*bf0 - af0*bf1

	return cubicRoundFloat(numRe / nB), cubicRoundFloat(numIm / nB)
}

func eisRem128(a0, a1, b0, b1 signed128) (signed128, signed128) {
	qr, qi := eisQuotient128(a0, a1, b0, b1)
	qb0 := sub128(mulSmall128(b0, qr), mulSmall128(b1, qi))
	qb1 := sub128(add128(mulSmall128(b1, qr), mulSmall128(b0, qi)), mulSmall128(b1, qi))
	return sub128(a0, qb0), sub128(a1, qb1)
}

func divBy1MinusOmega128(e, f *signed128) {
	twoEMinusF := sub128(mulSmall128(*e, 2), *f)
	sum := add128(*e, *f)
	*e = divExact3_128(twoEMinusF)
	*f = divExact3_128(sum)
}

func makePrimaryEis128(e, f *signed128) int {
	r0 := mod3_128(*e)
	// mod3(e-f) = (mod3(e) + 3 - mod3(f)) % 3
	r1 := (r0 + 3 - mod3_128(*f)) % 3

	if r0 == 0 {
		newE := sub128(*f, *e)
		newF := neg128(*e)
		*e = newE
		*f = newF
		return 1
	}
	if r1 == 0 {
		newE := neg128(*f)
		newF := sub128(*e, *f)
		*e = newE
		*f = newF
		return 2
	}
	return 0
}

func cubicCorrection128(b0, b1 signed128, m uint64, n int) int {
	b0m9 := mod9_128(b0)
	b1m9 := mod9_128(b1)
	b0sqM9 := (b0m9 * b0m9) % 9
	b0b1M9 := (b0m9 * b1m9) % 9
	termM := (1 + 9 - b0sqM9) % 9
	termN := (b0sqM9 + 18 - b0b1M9 - 1) % 9
	q0m9 := (m%9*termM + uint64(n)*termN) % 9

	switch q0m9 {
	case 0:
		return 0
	case 3:
		return 1
	case 6:
		return 2
	default:
		return -1
	}
}

// fitsIn128 returns true if all four signed256 values safely fit in signed128.
// We require components < 2^96 (w1 < 2^32) to leave headroom for small
// multiplications (by quotient components ≤ 4) inside the signed128 GCD loop.
func fitsIn128(a0, a1, b0, b1 signed256) bool {
	const maxW1 = uint64(1) << 32
	return a0.w2 == 0 && a0.w3 == 0 && a0.w1 < maxW1 &&
		a1.w2 == 0 && a1.w3 == 0 && a1.w1 < maxW1 &&
		b0.w2 == 0 && b0.w3 == 0 && b0.w1 < maxW1 &&
		b1.w2 == 0 && b1.w3 == 0 && b1.w1 < maxW1
}

func s256to128(s signed256) signed128 {
	return signed128{s.w0, s.w1, s.neg}
}

// --- Eisenstein arithmetic with signed256 ---

func cubicRoundFloat(x float64) int64 {
	if x >= 0 {
		return int64(x + 0.5)
	}
	return -int64(-x + 0.5)
}

// eisQuotient256 computes the nearest Eisenstein quotient of a/b using float64.
// In ℤ[ω]: q = round(a·conj(b) / N(b)) where conj(b₀+b₁ω) = (b₀-b₁)+(-b₁)ω
// and N(b₀+b₁ω) = b₀²+b₁²-b₀b₁.
func eisQuotient256(a0, a1, b0, b1 signed256) (int64, int64) {
	af0 := s256ToFloat(a0)
	af1 := s256ToFloat(a1)
	bf0 := s256ToFloat(b0)
	bf1 := s256ToFloat(b1)

	nB := bf0*bf0 + bf1*bf1 - bf0*bf1
	if nB == 0 {
		return 0, 0
	}

	// a·conj(b):
	// Re = a₀b₀ - a₀b₁ + a₁b₁
	// Im = a₁b₀ - a₀b₁
	numRe := af0*bf0 - af0*bf1 + af1*bf1
	numIm := af1*bf0 - af0*bf1

	return cubicRoundFloat(numRe / nB), cubicRoundFloat(numIm / nB)
}

// eisRem256 computes the Eisenstein remainder a mod b in ℤ[ω].
func eisRem256(a0, a1, b0, b1 signed256) (signed256, signed256) {
	qr, qi := eisQuotient256(a0, a1, b0, b1)
	// q·b = (qr+qi·ω)(b₀+b₁·ω) = (qr·b₀-qi·b₁) + (qr·b₁+qi·b₀-qi·b₁)·ω
	qb0 := sub256(mulSmall256(b0, qr), mulSmall256(b1, qi))
	qb1 := sub256(add256(mulSmall256(b1, qr), mulSmall256(b0, qi)), mulSmall256(b1, qi))
	return sub256(a0, qb0), sub256(a1, qb1)
}

// divBy1MinusOmega256 divides (e + f·ω) by (1-ω).
// (e + f·ω)/(1-ω) = ((2e-f)/3) + ((e+f)/3)·ω
// Valid only when e+f ≡ 0 mod 3.
func divBy1MinusOmega256(e, f *signed256) {
	twoEMinusF := sub256(mulSmall256(*e, 2), *f)
	sum := add256(*e, *f)
	*e = divExact3_256(twoEMinusF)
	*f = divExact3_256(sum)
}

// makePrimaryEis256 finds n (0 ≤ n < 3) such that (e+f·ω)·ω^{-n} is primary.
// Modifies e, f in place to the primary associate. Returns n.
func makePrimaryEis256(e, f *signed256) int {
	r0 := mod3_256(*e)
	// mod3(e-f) = (mod3(e) + 3 - mod3(f)) % 3
	r1 := (r0 + 3 - mod3_256(*f)) % 3

	if r0 == 0 {
		// n=1: multiply by ω² → (f-e) + (-e)·ω
		newE := sub256(*f, *e)
		newF := neg256(*e)
		*e = newE
		*f = newF
		return 1
	}
	if r1 == 0 {
		// n=2: multiply by ω → (-f) + (e-f)·ω
		newE := neg256(*f)
		newF := sub256(*e, *f)
		*e = newE
		*f = newF
		return 2
	}
	return 0
}

func isRealUnit256(re, im signed256) bool {
	return im.isZero() && re.w1 == 0 && re.w2 == 0 && re.w3 == 0 && re.w0 <= 1
}

// --- Phase 1: pooled big.Int scratch ---

type cubicScratch struct {
	xBI, numRe, numIm, qRe, qIm, t1, t2, e, f big.Int
}

var cubicPool = sync.Pool{
	New: func() any {
		sc := new(cubicScratch)
		for _, p := range []*big.Int{&sc.xBI, &sc.numRe, &sc.numIm, &sc.qRe, &sc.qIm, &sc.t1, &sc.t2, &sc.e, &sc.f} {
			p.SetBits(make([]big.Word, 8))
			p.SetUint64(0)
		}
		return sc
	},
}

// Precomputed constants for the Eisenstein prime β = a + b·ω with norm p.
var (
	cubBetaA0BI   big.Int   // a (negative)
	cubBetaA1BI   big.Int   // b (positive)
	cubBetaConjBI big.Int   // a - b (conjugate's real part)
	cubNormBI     big.Int   // N(β) = a² + b² - ab = p
	cubBeta256A0  signed256 // a as signed256
	cubBeta256A1  signed256 // b as signed256
	cubBiOne      = big.NewInt(1)
)

func init() {
	cubBetaA0BI.SetString("433386200905713772878563252435522861392305022889786430806", 10)
	cubBetaA1BI.SetString("216693100452856886439281626217761430700483883102737443499", 10)
	cubBetaConjBI.Sub(&cubBetaA0BI, &cubBetaA1BI)

	var t1, t2, t3 big.Int
	t1.Mul(&cubBetaA0BI, &cubBetaA0BI)
	t2.Mul(&cubBetaA1BI, &cubBetaA1BI)
	t3.Mul(&cubBetaA0BI, &cubBetaA1BI)
	cubNormBI.Add(&t1, &t2)
	cubNormBI.Sub(&cubNormBI, &t3)

	bigToS256(&cubBeta256A0, &cubBetaA0BI)
	bigToS256(&cubBeta256A1, &cubBetaA1BI)
}

// cubicRoundDiv computes z = round(a/b) for b > 0, using pre-allocated temps.
func cubicRoundDiv(z, aa, bb *big.Int, q, r *big.Int) {
	q.QuoRem(aa, bb, r)
	r.Abs(r)
	r.Lsh(r, 1)
	if r.Cmp(bb) > 0 {
		if aa.Sign() >= 0 {
			z.Add(q, cubBiOne)
		} else {
			z.Sub(q, cubBiOne)
		}
	} else {
		z.Set(q)
	}
}

// cubicCorrection computes the power-of-ω correction for one GCD step.
// Given current denominator (b₀, b₁) and the (1-ω)-valuation m and primary-adjustment n,
// returns the exponent of ω to add to the result (0, 1, or 2), or -1 on error.
func cubicCorrection(b0, b1 signed256, m uint64, n int) int {
	b0m9 := mod9_256(b0)
	b1m9 := mod9_256(b1)
	b0sqM9 := (b0m9 * b0m9) % 9
	b0b1M9 := (b0m9 * b1m9) % 9
	// ((1-ω)/β)₃ = ω^{(1-b₀²)/3}: termM = (1-b₀²) mod 9
	termM := (1 + 9 - b0sqM9) % 9
	// (ω/β)₃ = ω^{(b₀²-b₀b₁-1)/3}: termN = (b₀²-b₀b₁-1) mod 9
	termN := (b0sqM9 + 18 - b0b1M9 - 1) % 9
	q0m9 := (m%9*termM + uint64(n)*termN) % 9

	switch q0m9 {
	case 0:
		return 0
	case 3:
		return 1
	case 6:
		return 2
	default:
		return -1
	}
}

// CubicSymbolFast computes the cubic residue symbol of x modulo the BLS12-376
// Eisenstein prime using a two-phase approach:
//   - Phase 1: one big.Int Euclidean step to reduce from ~376 bits to ~192 bits
//   - Phase 2: Eisenstein GCD loop using fixed-width signed256/signed128 arithmetic
//
// Returns 0 (symbol=1, cubic residue), 1 (symbol=ω), or 2 (symbol=ω²).
func CubicSymbolFast(x fp.Element) uint8 {
	if x.IsZero() {
		return 0
	}

	// Phase 1: compute first remainder = (x, 0) mod β using big.Int
	sc := cubicPool.Get().(*cubicScratch)
	x.BigInt(&sc.xBI)

	// q = round(x·conj(β) / N(β))
	// x·conj(β) where x=(xBI,0), conj(β)=(a-b, -b):
	//   Re = xBI·(a-b),  Im = -xBI·b
	sc.numRe.Mul(&sc.xBI, &cubBetaConjBI)
	sc.numIm.Mul(&sc.xBI, &cubBetaA1BI)
	sc.numIm.Neg(&sc.numIm)

	cubicRoundDiv(&sc.qRe, &sc.numRe, &cubNormBI, &sc.t1, &sc.t2)
	cubicRoundDiv(&sc.qIm, &sc.numIm, &cubNormBI, &sc.t1, &sc.t2)

	// remainder = (x, 0) - q·β
	// q·β = (qRe·a - qIm·b) + (qRe·b + qIm·(a-b))·ω
	sc.t1.Mul(&sc.qRe, &cubBetaA0BI)
	sc.t2.Mul(&sc.qIm, &cubBetaA1BI)
	sc.e.Sub(&sc.t1, &sc.t2)
	sc.e.Sub(&sc.xBI, &sc.e)

	sc.t1.Mul(&sc.qRe, &cubBetaA1BI)
	sc.t2.Mul(&sc.qIm, &cubBetaConjBI)
	sc.f.Add(&sc.t1, &sc.t2)
	sc.f.Neg(&sc.f)

	if sc.e.Sign() == 0 && sc.f.Sign() == 0 {
		cubicPool.Put(sc)
		return 0
	}

	var e0, e1 signed256
	bigToS256(&e0, &sc.e)
	bigToS256(&e1, &sc.f)
	cubicPool.Put(sc)

	// Process first remainder: remove (1-ω) factors and make primary
	result := uint64(0)

	m := uint64(0)
	for (mod3_256(e0)+mod3_256(e1))%3 == 0 {
		divBy1MinusOmega256(&e0, &e1)
		m++
	}
	n := makePrimaryEis256(&e0, &e1)

	corr := cubicCorrection(cubBeta256A0, cubBeta256A1, m, n)
	if corr < 0 {
		return cubicSymbolFallback(x)
	}
	result = uint64(corr)

	// Swap: a = β_orig, b = processed first remainder
	a0, a1 := cubBeta256A0, cubBeta256A1
	b0, b1 := e0, e1

	// Phase 2a: signed256 Eisenstein GCD loop until components fit in 128 bits
	for iter := 0; ; iter++ {
		if iter > 300 {
			return cubicSymbolFallback(x)
		}

		if isRealUnit256(a0, a1) || isRealUnit256(b0, b1) {
			return uint8(result)
		}

		// Check if we can switch to the faster signed128 loop
		if fitsIn128(a0, a1, b0, b1) {
			return cubicGCD128(s256to128(a0), s256to128(a1), s256to128(b0), s256to128(b1), result, x)
		}

		e0, e1 = eisRem256(a0, a1, b0, b1)

		if e0.isZero() && e1.isZero() {
			return 0
		}

		m = 0
		for (mod3_256(e0)+mod3_256(e1))%3 == 0 {
			divBy1MinusOmega256(&e0, &e1)
			m++
		}

		n = makePrimaryEis256(&e0, &e1)

		corr = cubicCorrection(b0, b1, m, n)
		if corr < 0 {
			return cubicSymbolFallback(x)
		}
		result = (result + uint64(corr)) % 3

		a0, a1 = b0, b1
		b0, b1 = e0, e1
	}
}

// cubicGCD128 continues the Eisenstein GCD using signed128 arithmetic.
func cubicGCD128(a0, a1, b0, b1 signed128, result uint64, x fp.Element) uint8 {
	for iter := 0; ; iter++ {
		if iter > 200 {
			return cubicSymbolFallback(x)
		}

		if isRealUnit128(a0, a1) || isRealUnit128(b0, b1) {
			break
		}

		e0, e1 := eisRem128(a0, a1, b0, b1)

		if e0.isZero() && e1.isZero() {
			return 0
		}

		m := uint64(0)
		for (mod3_128(e0)+mod3_128(e1))%3 == 0 {
			divBy1MinusOmega128(&e0, &e1)
			m++
		}

		n := makePrimaryEis128(&e0, &e1)

		corr := cubicCorrection128(b0, b1, m, n)
		if corr < 0 {
			return cubicSymbolFallback(x)
		}
		result = (result + uint64(corr)) % 3

		a0, a1 = b0, b1
		b0, b1 = e0, e1
	}

	return uint8(result)
}

// cubicSymbolFallback uses the exponentiation-based cubic character.
func cubicSymbolFallback(x fp.Element) uint8 {
	sym := expByp3(&x)
	if sym.IsOne() {
		return 0
	}
	if sym.Equal(&thirdRootOneG1) {
		return 1
	}
	return 2
}

// IsCubicResidueFast checks whether x is a cubic residue mod p using the fast
// Eisenstein GCD algorithm.
func IsCubicResidueFast(x *fp.Element) bool {
	return CubicSymbolFast(*x) == 0
}
