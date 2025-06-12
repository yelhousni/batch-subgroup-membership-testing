package bls12377

import (
	curve "github.com/consensys/gnark-crypto/ecc/bls12-377"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fp"
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

// fromJacExtended converts an extended Jacobian point to a Jacobian point.
func fromJacExtended(q *g1JacExtended) *curve.G1Jac {
	var p curve.G1Jac
	if q.ZZ.IsZero() {
		p.X.SetOne()
		p.Y.SetOne()
		return &p
	}
	p.X.Mul(&q.ZZ, &q.X).Mul(&p.X, &q.ZZ)
	p.Y.Mul(&q.ZZZ, &q.Y).Mul(&p.Y, &q.ZZZ)
	p.Z.Set(&q.ZZZ)
	return &p
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
