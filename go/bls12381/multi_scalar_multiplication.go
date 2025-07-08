package bls12381

import (
	"crypto/rand"
	curve "github.com/consensys/gnark-crypto/ecc/bls12-381"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fp"
	"unsafe"
)

// g1JacExtended is a point in extended Jacobian coordinates (x=X/ZZ, y=Y/ZZZ, ZZ³=ZZZ²)
type g1JacExtended struct {
	X, Y, ZZ, ZZZ fp.Element
}

// Set sets p to a in extended Jacobian coordinates.
func (p *g1JacExtended) Set(q *g1JacExtended) *g1JacExtended {
	p.X, p.Y, p.ZZ, p.ZZZ = q.X, q.Y, q.ZZ, q.ZZZ
	return p
}

// SetInfinity sets p to the infinity point (1,1,0,0).
func (p *g1JacExtended) SetInfinity() *g1JacExtended {
	p.X.SetOne()
	p.Y.SetOne()
	p.ZZ = fp.Element{}
	p.ZZZ = fp.Element{}
	return p
}

// IsInfinity checks if the p is infinity, i.e. p.ZZ=0.
func (p *g1JacExtended) IsInfinity() bool {
	return p.ZZ.IsZero()
}

// unsafeFromJacExtended converts an extended Jacobian point, distinct from Infinity, to a Jacobian point.
func unsafeFromJacExtended(q *g1JacExtended) *curve.G1Jac {
	var p curve.G1Jac
	p.X.Square(&q.ZZ).Mul(&p.X, &q.X)
	p.Y.Square(&q.ZZZ).Mul(&p.Y, &q.Y)
	p.Z = q.ZZZ
	return &p
}

// add sets p to p+q in extended Jacobian coordinates.
//
// https://www.hyperelliptic.org/EFD/g1p/auto-shortw-xyzz.html#addition-add-2008-s
func (p *g1JacExtended) add(q *g1JacExtended) *g1JacExtended {
	//if q is infinity return p
	if q.ZZ.IsZero() {
		return p
	}
	// p is infinity, return q
	if p.ZZ.IsZero() {
		p.Set(q)
		return p
	}

	var A, B, U1, U2, S1, S2 fp.Element

	// p2: q, p1: p
	U2.Mul(&q.X, &p.ZZ)
	U1.Mul(&p.X, &q.ZZ)
	A.Sub(&U2, &U1)
	S2.Mul(&q.Y, &p.ZZZ)
	S1.Mul(&p.Y, &q.ZZZ)
	B.Sub(&S2, &S1)

	if A.IsZero() {
		if B.IsZero() {
			return p.double(q)

		}
		p.ZZ = fp.Element{}
		p.ZZZ = fp.Element{}
		return p
	}

	var P, R, PP, PPP, Q, V fp.Element
	P.Sub(&U2, &U1)
	R.Sub(&S2, &S1)
	PP.Square(&P)
	PPP.Mul(&P, &PP)
	Q.Mul(&U1, &PP)
	V.Mul(&S1, &PPP)

	p.X.Square(&R).
		Sub(&p.X, &PPP).
		Sub(&p.X, &Q).
		Sub(&p.X, &Q)
	p.Y.Sub(&Q, &p.X).
		Mul(&p.Y, &R).
		Sub(&p.Y, &V)
	p.ZZ.Mul(&p.ZZ, &q.ZZ).
		Mul(&p.ZZ, &PP)
	p.ZZZ.Mul(&p.ZZZ, &q.ZZZ).
		Mul(&p.ZZZ, &PPP)

	return p
}

// double sets p to [2]q in Jacobian extended coordinates.
//
// http://www.hyperelliptic.org/EFD/g1p/auto-shortw-xyzz.html#doubling-dbl-2008-s-1
// N.B.: since we consider any point on Z=0 as the point at infinity
// this doubling formula works for infinity points as well.
func (p *g1JacExtended) double(q *g1JacExtended) *g1JacExtended {
	var U, V, W, S, XX, M fp.Element

	U.Double(&q.Y)
	V.Square(&U)
	W.Mul(&U, &V)
	S.Mul(&q.X, &V)
	XX.Square(&q.X)
	M.Double(&XX).
		Add(&M, &XX) // -> + A, but A=0 here
	U.Mul(&W, &q.Y)

	p.X.Square(&M).
		Sub(&p.X, &S).
		Sub(&p.X, &S)
	p.Y.Sub(&S, &p.X).
		Mul(&p.Y, &M).
		Sub(&p.Y, &U)
	p.ZZ.Mul(&V, &q.ZZ)
	p.ZZZ.Mul(&W, &q.ZZZ)

	return p
}

// addMixed sets p to p+q in extended Jacobian coordinates, where a.ZZ=1.
//
// http://www.hyperelliptic.org/EFD/g1p/auto-shortw-xyzz.html#addition-madd-2008-s
func (p *g1JacExtended) addMixed(a *curve.G1Affine) *g1JacExtended {

	//if a is infinity return p
	if a.IsInfinity() {
		return p
	}
	// p is infinity, return a
	if p.ZZ.IsZero() {
		p.X = a.X
		p.Y = a.Y
		p.ZZ.SetOne()
		p.ZZZ.SetOne()
		return p
	}

	var P, R fp.Element

	// p2: a, p1: p
	P.Mul(&a.X, &p.ZZ)
	P.Sub(&P, &p.X)

	R.Mul(&a.Y, &p.ZZZ)
	R.Sub(&R, &p.Y)

	if P.IsZero() {
		if R.IsZero() {
			return p.doubleMixed(a)

		}
		p.ZZ = fp.Element{}
		p.ZZZ = fp.Element{}
		return p
	}

	var PP, PPP, Q, Q2, RR, X3, Y3 fp.Element

	PP.Square(&P)
	PPP.Mul(&P, &PP)
	Q.Mul(&p.X, &PP)
	RR.Square(&R)
	X3.Sub(&RR, &PPP)
	Q2.Double(&Q)
	p.X.Sub(&X3, &Q2)
	Y3.Sub(&Q, &p.X).Mul(&Y3, &R)
	R.Mul(&p.Y, &PPP)
	p.Y.Sub(&Y3, &R)
	p.ZZ.Mul(&p.ZZ, &PP)
	p.ZZZ.Mul(&p.ZZZ, &PPP)

	return p

}

// doubleMixed sets p to [2]a in Jacobian extended coordinates, where a.ZZ=1.
//
// http://www.hyperelliptic.org/EFD/g1p/auto-shortw-xyzz.html#doubling-dbl-2008-s-1
func (p *g1JacExtended) doubleMixed(a *curve.G1Affine) *g1JacExtended {

	var U, V, W, S, XX, M, S2, L fp.Element

	U.Double(&a.Y)
	V.Square(&U)
	W.Mul(&U, &V)
	S.Mul(&a.X, &V)
	XX.Square(&a.X)
	M.Double(&XX).
		Add(&M, &XX) // -> + A, but A=0 here
	S2.Double(&S)
	L.Mul(&W, &a.Y)

	p.X.Square(&M).
		Sub(&p.X, &S2)
	p.Y.Sub(&S, &p.X).
		Mul(&p.Y, &M).
		Sub(&p.Y, &L)
	p.ZZ.Set(&V)
	p.ZZZ.Set(&W)

	return p
}

// --- MSM ---
type bucketg1JacExtendedC6 [32]g1JacExtended

func _msmCheck(points []curve.G1Affine) bool {
	// const nbBitsBounds = 13
	const c = 6
	const nbChunks = 3 //(nbBitsBounds + c - 1) / c

	// for each chunk, spawn one go routine that'll loop through all the scalars in the
	// corresponding bit-window
	// note that buckets is an array allocated on the stack and this is critical for performance

	// each go routine sends its result in chChunks[i] channel
	chChunks := make([]chan g1JacExtended, nbChunks)
	for i := 0; i < len(chChunks); i++ {
		chChunks[i] = make(chan g1JacExtended, 1)
	}

	for j := int(nbChunks - 1); j >= 0; j-- {
		go processChunkG1Simplified[bucketg1JacExtendedC6](uint64(j), chChunks[j], c, points)
	}

	p := msmReduceChunkG1Affine(int(c), chChunks[:])

	return p.IsInSubGroup()
}

func processChunkG1Simplified[B bucketg1JacExtendedC6](chunk uint64,
	chRes chan<- g1JacExtended,
	c uint64,
	points []curve.G1Affine) {

	const windowSize = 1024
	var br [windowSize * 2]byte

	// interpret br as an array of uint16 of size windowSize/2
	randomScalars := (*[windowSize]uint16)(unsafe.Pointer(&br[0]))

	// we need a mask to get only the (c-1) lowest bits of each scalar
	mask := uint16((1 << (c - 1)) - 1)

	var buckets B
	for i := 0; i < len(buckets); i++ {
		buckets[i].SetInfinity()
	}

	// for each scalars, get the digit corresponding to the chunk we're processing.
	for i := range points {
		if i%windowSize == 0 {
			// fill the lowest c bits of each scalar with random bytes
			rand.Read(br[:]) // does not return an error, always fills br
		}
		digit := randomScalars[i%windowSize] & mask
		if digit == 0 {
			continue
		}
		buckets[digit-1].addMixed(&points[i])
	}

	// reduce buckets into total
	// total =  bucket[0] + 2*bucket[1] + 3*bucket[2] ... + n*bucket[n-1]

	var runningSum, total g1JacExtended
	runningSum.SetInfinity()
	total.SetInfinity()
	for k := len(buckets) - 1; k >= 0; k-- {
		if !buckets[k].IsInfinity() {
			runningSum.add(&buckets[k])
		}
		total.add(&runningSum)
	}

	chRes <- total
}

// msmReduceChunkG1Affine reduces the weighted sum of the buckets into the result of the multiExp
func msmReduceChunkG1Affine(c int, chChunks []chan g1JacExtended) *curve.G1Jac {
	var _p g1JacExtended
	totalj := <-chChunks[len(chChunks)-1]
	_p.Set(&totalj)
	for j := len(chChunks) - 2; j >= 0; j-- {
		for l := 0; l < c; l++ {
			_p.double(&_p)
		}
		totalj := <-chChunks[j]
		_p.add(&totalj)
	}

	return unsafeFromJacExtended(&_p)
}
