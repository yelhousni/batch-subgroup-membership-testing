// Package bls12376strong efficient elliptic curve, pairing and hash to curve implementation for bls12376-strong.
//
// bls12376-strong: A batch-SMT-friendly G2-strong and Gt-strong Barreto--Lynn--Scott curve with
//
//	embedding degree k=12
//	seed x₀=-0x78383f2600000001
//	𝔽r: r=0xc734c9428876f0384334fc1cb852a71933b8f2c9102a034f0707e4c00000001 (x₀⁴-x₀²+1)
//	𝔽p: p=0xea4ce6919c7e9da513d75a5627c38344bcd510d6c9940f5ed48752a49ef938f617216eb1f402881353681caaaaaaab ((x₀-1)² ⋅ r(x₀)/3+x₀)
//	(E/𝔽p): Y²=X³+1
//	(Et/𝔽p²): Y² = X³+(u+2) (D-type twist)
//	r ∣ #E(Fp) and r ∣ #Et(𝔽p²)
//
// Extension fields tower:
//
//	𝔽p²[u] = 𝔽p/u²+1
//	𝔽p⁶[v] = 𝔽p²/v³-u-2
//	𝔽p¹²[w] = 𝔽p⁶/w²-v
//
// optimal Ate loop size:
//
//	x₀
//
// Security: estimated 126-bit level following [https://eprint.iacr.org/2019/885.pdf]
// (r is 253 bits and p¹² is 4517 bits)
//
// # Warning
//
// This code has not been audited and is provided as-is. In particular, there is no security guarantees such as constant time implementation or side-channel attack resistance.
package bls12376strong

import (
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	"github.com/yelhousni/batch-subgroup-membership/bls12376-strong/fp"
	"github.com/yelhousni/batch-subgroup-membership/bls12376-strong/fr"
	"github.com/yelhousni/batch-subgroup-membership/bls12376-strong/internal/fptower"
)

// aCurveCoeff is the a coefficients of the curve Y²=X³+ax+b
var aCurveCoeff fp.Element
var bCurveCoeff fp.Element

// bTwistCurveCoeff b coeff of the twist (defined over 𝔽p²) curve
var bTwistCurveCoeff fptower.E2

// generators of the r-torsion group, resp. in ker(pi-id), ker(Tr)
var g1Gen G1Jac
var g2Gen G2Jac

var g1GenAff G1Affine
var g2GenAff G2Affine

// point at infinity
var g1Infinity G1Jac
var g2Infinity G2Jac

// optimal Ate loop counter
var LoopCounter [64]int8

// Parameters useful for the GLV scalar multiplication. The third roots define the
// endomorphisms ϕ₁ and ϕ₂ for <G1Affine> and <G2Affine>. lambda is such that <r, ϕ-λ> lies above
// <r> in the ring Z[ϕ]. More concretely it's the associated eigenvalue
// of ϕ₁ (resp ϕ₂) restricted to <G1Affine> (resp <G2Affine>)
// see https://www.cosic.esat.kuleuven.be/nessie/reports/phase2/GLV.pdf
var thirdRootOneG1 fp.Element
var thirdRootOneG2 fp.Element
var lambdaGLV big.Int

// glvBasis stores R-linearly independent vectors (a,b), (c,d)
// in ker((u,v) → u+vλ[r]), and their determinant
var glvBasis ecc.Lattice

// ψ o π o ψ⁻¹, where ψ:E → E' is the degree 6 iso defined over 𝔽p¹²
var endo struct {
	u fptower.E2
	v fptower.E2
}

// seed x₀ of the curve
var xGen big.Int

// expose the tower -- github.com/consensys/gnark uses it in a gnark circuit

// 𝔽p²
type E2 = fptower.E2

// 𝔽p⁶
type E6 = fptower.E6

// 𝔽p¹²
type E12 = fptower.E12

func init() {
	aCurveCoeff.SetUint64(0)
	bCurveCoeff.SetUint64(1)

	// M-twist
	bTwistCurveCoeff.A0.SetUint64(2)
	bTwistCurveCoeff.A1.SetUint64(1)

	g1Gen.X.SetString("61330030814471061323712620583080524012204268040824655382537602348181775880210169900398027460586532789718009095159")
	g1Gen.Y.SetString("52868851560239652006481474438593256759490795966216789134321594112084450134421472552309055168836178468410056480251")
	g1Gen.Z.SetOne()

	g2Gen.X.SetString("28988179907377812980831649841643169916204148000111384526394004226223791135945078543147889732228863982955130291721",
		"70871326616439600719313602985971642611927053780911795168783665427924556382691605055542283811201110134450518849666")
	g2Gen.Y.SetString("9553056756834044529237883456099116967066707965755647302212023049612458685085883427795013851521809030648063150703",
		"37395328089683907784108448508154777420462345444396929856831817173416965886766627854198375253760548226373570882881")
	g2Gen.Z.SetString("1",
		"0")

	g1GenAff.FromJacobian(&g1Gen)
	g2GenAff.FromJacobian(&g2Gen)

	// (X,Y,Z) = (1,1,0)
	g1Infinity.X.SetOne()
	g1Infinity.Y.SetOne()
	g2Infinity.X.SetOne()
	g2Infinity.Y.SetOne()

	thirdRootOneG1.SetString("140867699351615776039709482774164536979651071563248383984353885353678501780284689700023347401122934816966888303270") // x₀^5-3x₀^4+3x₀^3-x₀+1
	thirdRootOneG2.Square(&thirdRootOneG1)
	lambdaGLV.SetString("75043121753505027792636171519680053248", 10) //(x₀²-1)
	_r := fr.Modulus()
	ecc.PrecomputeLattice(_r, &lambdaGLV, &glvBasis)

	endo.u.A0.SetString("111500873050125781194417014601114806570209374443562897647051899312710190755662650302772863774467674177942918775363")
	endo.u.A1.SetString("68748812938746705868036438902226141839486122574411760541166443540320030400452424842532684433527940788298361761168")
	endo.v.A0.SetString("65666114604077535587976484855176109210017356773698946743194016549967977100157422451145916833822102834456585253072")
	endo.v.A1.SetString("9535470143460704912540493225097022982725398892628277965570340605126442295165261819974918101063502935169374636299")

	// -x₀
	xGen.SetString("8662743315688456193", 10)

	// 2-NAF decomposition of -x₀ little endian
	ecc.NafDecomposition(&xGen, LoopCounter[:])
}

// Generators return the generators of the r-torsion group, resp. in ker(pi-id), ker(Tr)
func Generators() (g1Jac G1Jac, g2Jac G2Jac, g1Aff G1Affine, g2Aff G2Affine) {
	g1Aff = g1GenAff
	g2Aff = g2GenAff
	g1Jac = g1Gen
	g2Jac = g2Gen
	return
}

// CurveCoefficients returns the a, b coefficients of the curve equation.
func CurveCoefficients() (a, b fp.Element) {
	return aCurveCoeff, bCurveCoeff
}
