// Copyright 2020-2025 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

package fptower

import "github.com/yelhousni/batch-subgroup-membership/bls12376-strong/fp"

// Frobenius set z to Frobenius(x), return z
func (z *E12) Frobenius(x *E12) *E12 {
	/*
		// Algorithm 28 from https://eprint.iacr.org/2010/354.pdf (beware typos!)
		var t [6]E2

		// Frobenius acts on fp2 by conjugation
		t[0].Conjugate(&x.C0.B0)
		t[1].Conjugate(&x.C0.B1)
		t[2].Conjugate(&x.C0.B2)
		t[3].Conjugate(&x.C1.B0)
		t[4].Conjugate(&x.C1.B1)
		t[5].Conjugate(&x.C1.B2)

		t[1].MulByNonResidue1Power2(&t[1])
		t[2].MulByNonResidue1Power4(&t[2])
		t[3].MulByNonResidue1Power1(&t[3])
		t[4].MulByNonResidue1Power3(&t[4])
		t[5].MulByNonResidue1Power5(&t[5])

		z.C0.B0 = t[0]
		z.C0.B1 = t[1]
		z.C0.B2 = t[2]
		z.C1.B0 = t[3]
		z.C1.B1 = t[4]
		z.C1.B2 = t[5]
	*/

	return z.Exp(*x, fp.Modulus())
}

// FrobeniusSquare set z to Frobenius^2(x), and return z
func (z *E12) FrobeniusSquare(x *E12) *E12 {
	// Algorithm 29 from https://eprint.iacr.org/2010/354.pdf (beware typos!)
	var t [6]E2

	t[1].MulByNonResidue2Power2(&x.C0.B1)
	t[2].MulByNonResidue2Power4(&x.C0.B2)
	t[3].MulByNonResidue2Power1(&x.C1.B0)
	t[4].MulByNonResidue2Power3(&x.C1.B1)
	t[5].MulByNonResidue2Power5(&x.C1.B2)

	z.C0.B0 = x.C0.B0
	z.C0.B1 = t[1]
	z.C0.B2 = t[2]
	z.C1.B0 = t[3]
	z.C1.B1 = t[4]
	z.C1.B2 = t[5]

	return z
}

// MulByNonResidue1Power1 set z=x*(2,1)^(1*(p^1-1)/6) and return z
func (z *E2) MulByNonResidue1Power1(x *E2) *E2 {
	var b E2
	b.A0.SetString("29177675865650992594787457015448945452360868311852759961698790731039461202404454988570800873188746199730051751713")
	b.A1.SetString("116765314714690436723542560131291747950248794923099216398770044227380785464118002303740288921054076719620637736859")
	z.Mul(x, &b)
	return z
}

// MulByNonResidue1Power2 set z=x*(2,1)^(2*(p^1-1)/6) and return z
func (z *E2) MulByNonResidue1Power2(x *E2) *E2 {
	var b E2
	b.A0.SetString("140561165204752693968442986751009159918684358697519646318104245152543753703565053581517093024779089196165329841598")
	b.A1.SetString("12501110777296501362649498893840435564216723874639740444505037171775817995576077929995824260237754684131039407530")
	z.Mul(x, &b)
	return z
}

// MulByNonResidue1Power3 set z=x*(2,1)^(3*(p^1-1)/6) and return z
func (z *E2) MulByNonResidue1Power3(x *E2) *E2 {
	var b E2
	b.A0.SetString("75201584747538240500516978080273132192742755666327224708764357155094419395322684271120834934885605769625959889371")
	b.A1.SetString("9535470143460704912540493225097022982725398892628277965570340605126442295165261819974918101063502935169374636299")
	z.Mul(x, &b)
	return z
}

// MulByNonResidue1Power4 set z=x*(2,1)^(4*(p^1-1)/6) and return z
func (z *E2) MulByNonResidue1Power4(x *E2) *E2 {
	var b E2
	b.A0.SetString("65552494440426522975681974708270100852610477749656678440380944937857180175567487366783165104019873695771422617175")
	b.A1.SetString("121089828545086329871480407182982929530836308917828143949408395120098741430317462730021421762991875455523736207342")
	z.Mul(x, &b)
	return z
}

// MulByNonResidue1Power5 set z=x*(2,1)^(5*(p^1-1)/6) and return z
func (z *E2) MulByNonResidue1Power5(x *E2) *E2 {
	var b E2
	b.A0.SetString("97644110826275522386934201807896205628646150259488954007340862889685155355370443857087601114384841457871232815005")
	b.A1.SetString("74329664808463695693511827733761867325216837452802813672391973783213790041437903564173726002770420825766866016855")
	z.Mul(x, &b)
	return z
}

// MulByNonResidue2Power1 set z=x*(2,1)^(1*(p^2-1)/6) and return z
func (z *E2) MulByNonResidue2Power1(x *E2) *E2 {
	var b fp.Element
	b.SetString("48783980161284704423109040876777787467604488351383894715195417022243404367584773787115656839173")
	z.A0.Mul(&x.A0, &b)
	z.A1.Mul(&x.A1, &b)
	return z
}

// MulByNonResidue2Power2 set z=x*(2,1)^(2*(p^2-1)/6) and return z
func (z *E2) MulByNonResidue2Power2(x *E2) *E2 {
	var b fp.Element
	b.SetString("48783980161284704423109040876777787467604488351383894715195417022243404367584773787115656839172")
	z.A0.Mul(&x.A0, &b)
	z.A1.Mul(&x.A1, &b)
	return z
}

// MulByNonResidue2Power3 set z=x*(2,1)^(3*(p^2-1)/6) and return z
func (z *E2) MulByNonResidue2Power3(x *E2) *E2 {
	var b fp.Element
	b.SetString("140867699351615776088493462935449241402760112440026171451958373705062396495480106722266751768707708604082545142442")
	z.A0.Mul(&x.A0, &b)
	z.A1.Mul(&x.A1, &b)
	return z
}

// MulByNonResidue2Power4 set z=x*(2,1)^(4*(p^2-1)/6) and return z
func (z *E2) MulByNonResidue2Power4(x *E2) *E2 {
	var b fp.Element
	b.SetString("140867699351615776039709482774164536979651071563248383984353885353678501780284689700023347401122934816966888303270")
	z.A0.Mul(&x.A0, &b)
	z.A1.Mul(&x.A1, &b)
	return z
}

// MulByNonResidue2Power5 set z=x*(2,1)^(5*(p^2-1)/6) and return z
func (z *E2) MulByNonResidue2Power5(x *E2) *E2 {
	var b fp.Element
	b.SetString("140867699351615776039709482774164536979651071563248383984353885353678501780284689700023347401122934816966888303271")
	z.A0.Mul(&x.A0, &b)
	z.A1.Mul(&x.A1, &b)
	return z
}
