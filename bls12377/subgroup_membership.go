package bls12377

import (
	"crypto/rand"
	"math/big"

	"github.com/consensys/gnark-crypto/ecc"
	curve "github.com/consensys/gnark-crypto/ecc/bls12-377"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fp"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fr"
)

// IsInSubGroupBatchNaive checks if a batch of points P_i are in G1.
// This is a naive method that checks each point individually using Scott test
// [Scott21].
//
// [Scott21]: https://eprint.iacr.org/2021/1130.pdf
func IsInSubGroupBatchNaive(points []curve.G1Affine) bool {
	for i := range points {
		if !points[i].IsInSubGroup() {
			return false
		}
	}
	return true
}

// IsInSubGroupBatch checks if a batch of points P_i are in G1.
// First, it checks that all points are on a larger torsion E[r*e'] using Tate
// pairings [Koshelev22].
// Second, it generates random scalars s_i in the range [0, bound), performs
// n=rounds multi-scalar-multiplication Sj=∑[s_i]P_i of sizes N=len(points) and
// checks if Sj are on E[r] using Scott test [Scott21].
//
// [Koshelev22]: https://eprint.iacr.org/2022/037.pdf
// [Scott21]: https://eprint.iacr.org/2021/1130.pdf
func IsInSubGroupBatch(points []curve.G1Affine, bound *big.Int, rounds int) bool {

	// 1. Check points are on E[r*e']
	for i := range points {
		// 1.1. Tate_{2^4,P16}(Q) == 1
		if !isFirstTateOne(points[i]) {
			return false
		}
		// 1.2. Tate_{3,P3}(Q) == 1
		if !isSecondTateOne(points[i]) {
			return false
		}
		// 1.3. Tate_{7,P7}(Q) == 1
		if !isThirdTateOne(points[i]) {
			return false
		}
		// 1.4. Tate_{13,P13}(Q) == 1
		if !isFourthTateOne(points[i]) {
			return false
		}
	}

	// 2. Check Sj are on E[r]
	for i := 0; i < rounds; i++ {
		b, err := rand.Int(rand.Reader, bound)
		if err != nil {
			panic(err)
		}
		randoms := make([]fr.Element, len(points))
		for j := range randoms {
			randoms[j].SetBigInt(b)
		}
		var sum curve.G1Jac
		sum.MultiExp(points[:], randoms[:], ecc.MultiExpConfig{})
		if !sum.IsInSubGroup() {
			return false
		}
	}
	return true
}

// isFirstTateOne checks that Tate_{16,P16}(Q) == 1
// where P16=(x,y) a point of order 16 on the curve.
// x = 0x62af02c766f6bdd9a611da0db5cbd52102751f8476d97903a05a3274fb0b498d99628c38327d58081e1e8107930a64
// y = 0x17ffe960fa1502725f36844f63a83638e2a7e321948984a95ba562cf719f6785139614a1c1de0d2da5211384f762bd4
func isFirstTateOne(point curve.G1Affine) bool {
	// f_{16,P} = (l_{P,P}^8 * (l_{2P,2P}^4 * l_{4P,4P})^2) /
	// 			  (v_{2P}^8 * v_{4P}^4 * v_{8P})

	var num, denom, tate, f1, f2, f3 fp.Element

	// l_{P,P}^8
	f1.Mul(&point.X, &lines16[0].a).Add(&f1, &point.Y).Add(&f1, &lines16[0].b)
	num.Square(&f1).Square(&num).Square(&num)
	// l_{2P,2P}^4
	f1.Mul(&point.X, &lines16[1].a).Add(&f1, &point.Y).Add(&f1, &lines16[1].b)
	f1.Square(&f1).Square(&f1)
	num.Mul(&num, &f1)
	// l_{4P,4P}^2
	f1.Mul(&point.X, &lines16[2].a).Add(&f1, &point.Y).Add(&f1, &lines16[2].b)
	f1.Square(&f1)
	num.Mul(&num, &f1)

	// v_{2P}^4
	f1.Add(&point.X, &lines16[3].b)
	denom.Square(&f1).Square(&denom).Square(&denom)
	// v_{4P}^4
	f1.Add(&point.X, &lines16[4].b)
	f1.Square(&f1).Square(&f1)
	denom.Mul(&denom, &f1)
	// v_{8P}
	f1.Add(&point.X, &lines16[5].b)
	denom.Mul(&denom, &f1)

	// denom^{-1} = denom^{15} inside the 16-th power residue symbol
	f1.Square(&denom)
	f2.Square(&f1)
	f3.Square(&f2)
	f3.Mul(&f3, &f2).Mul(&f3, &f1)
	denom.Mul(&denom, &f3)

	// tate = num * denom^{-1}
	tate.Mul(&num, &denom)

	// tate^((p-1)/16)
	tate = *expByp16(&tate)

	return tate.IsOne()
}

// isSecondTateOne checks that Tate_{3,P3}(Q) = (y-1)^((p-1)/3) == 1
// where P3 = (0,1) a point of order 3 on the curve.
func isSecondTateOne(point curve.G1Affine) bool {
	var tate, one fp.Element
	one.SetOne()
	tate.Sub(&point.Y, &one)
	tate = *expByp3(&tate)
	return tate.IsOne()
}

// isThirdTateOne checks that Tate_{7,P7}(Q) == 1
// where P7 = (x, y) a point of order 7 on the curve.
// x = 0x7881bc38aea0c08601e9b87c6e6dffca55f96ce64faa8a54fcc76146dee6a5468ef503cfeac9f4c169e308abc0746
// y = 0xcf189de61c5daaddcb28032fd26caae5b8b0e27bc3ba0d4b4b357ac3343fadc164b2581f32600d6be5df353018380b
func isThirdTateOne(point curve.G1Affine) bool {
	// f_{7,P} = ((l_{P,P} * l_{2P,P})^2 * l_{3P,3P}) /
	// 			  (v_{2P} * v_{3P})^2

	var num, denom, tate, f1, f2 fp.Element

	// l_{P,P}
	f1.Mul(&point.X, &lines7[0].a).Add(&f1, &point.Y).Add(&f1, &lines7[0].b)
	// l_{2P,P}
	f2.Mul(&point.X, &lines7[1].a).Add(&f2, &point.Y).Add(&f2, &lines7[1].b)
	// (l_{P,P} * l_{2P,P})^2
	num.Mul(&f1, &f2).Square(&num)
	// l_{3P,3P}
	f1.Mul(&point.X, &lines7[2].a).Add(&f1, &point.Y).Add(&f1, &lines7[2].b)
	num.Mul(&num, &f1)

	// v_{2P}
	denom.Add(&point.X, &lines7[3].b)
	// v_{3P}
	f1.Add(&point.X, &lines7[4].b)
	denom.Mul(&denom, &f1).Square(&denom)

	// denom^{-1} = denom^{6} inside the 7-th power residue symbol
	f1.Square(&denom)
	f2.Square(&f1)
	denom.Mul(&f2, &f1)

	// tate = num * denom^{-1}
	tate.Mul(&num, &denom)

	// tate^((p-1)/7)
	tate = *expByp7(&tate)

	return tate.IsOne()
}

// isFourthTateOne checks that Tate_{13,P13}(Q) == 1
// where P13 = (x, y) a point of order 13 on the curve.
// x = 0x0x93a9ac62f3fa678eb04330e6973412e364875495c12a889ff93652e223784c93b2e63c6440cf9ab49ef3bd91c6f0b3
// y = 0x0x197214129826b8308cc6c9eb538506f5e5b9955731f6b654958e1d4134fbb6792508aa63c152574937d81ab9574c4da
func isFourthTateOne(point curve.G1Affine) bool {
	// f_{13,P} = ((l_{P,P} * l_{2P,P})^4 * l_{3P,3P}^2 * l_{6P,6P}) /
	// 			  ((v_{2P} * v_{3P})^4 * v_{6P}^2)

	var num, denom, tate, f1, f2 fp.Element

	// l_{P,P}
	f1.Mul(&point.X, &lines13[0].a).Add(&f1, &point.Y).Add(&f1, &lines13[0].b)
	// l_{2P,P}
	f2.Mul(&point.X, &lines13[1].a).Add(&f2, &point.Y).Add(&f2, &lines13[1].b)
	// (l_{P,P} * l_{2P,P})^4
	num.Mul(&f1, &f2).Square(&num).Square(&num)
	// l_{3P,3P}^2
	f1.Mul(&point.X, &lines13[2].a).Add(&f1, &point.Y).Add(&f1, &lines13[2].b)
	f1.Square(&f1)
	num.Mul(&num, &f1)
	// l_{6P,6P}
	f1.Mul(&point.X, &lines13[3].a).Add(&f1, &point.Y).Add(&f1, &lines13[3].b)
	num.Mul(&num, &f1)

	// v_{2P}
	denom.Add(&point.X, &lines13[4].b)
	// v_{3P}
	f1.Add(&point.X, &lines13[5].b)
	denom.Mul(&denom, &f1).Square(&denom).Square(&denom)
	// v_{6P}
	f1.Add(&point.X, &lines13[6].b)
	f1.Square(&f1)
	denom.Mul(&denom, &f1)

	// denom^{-1} = denom^{12} inside the 13-th power residue symbol
	f1.Square(&denom).Square(&f1)
	f2.Square(&f1)
	denom.Mul(&f2, &f1)

	// tate = num * denom^{-1}
	tate.Mul(&num, &denom)

	// tate^((p-1)/13)
	tate = *expByp13(&tate)

	return tate.IsOne()
}

// --------------------

// line represents a line in the form y + ax + b = 0.
//
// A a line through P=(x1,y1) and Q=(x2,y2) has:
//
//	a = (y1 - y2) / (x2 - x1)
//	b = -y1 - a*x1
//
// A tangent line to P=(x1,y1) has:
//
//	a = -3*x1^2 / (2*y1)
//	b = -y1 - a*x1
//
// A vertical line through P=(x1,y1) has::
//
// a = 0
// b = -x1
type line struct {
	a, b fp.Element
}

var lines16 [6]line
var lines7 [5]line
var lines13 [7]line

func init() {
	// P = (
	// 	   0x62af02c766f6bdd9a611da0db5cbd52102751f8476d97903a05a3274fb0b498d99628c38327d58081e1e8107930a64,
	// 	   0x17ffe960fa1502725f36844f63a83638e2a7e321948984a95ba562cf719f6785139614a1c1de0d2da5211384f762bd4,
	// ) of order 16
	//
	// l_{P,P}
	lines16[0].a.SetString("186753577267350681523223042672752493139503988550545545075290854342660310908569018683427166604708127195493070987924")
	lines16[0].b.SetString("107330485521934419149968690552530385427291718866313009555984104959818466429125580230768172226069235556223352865086")
	// l_{2P,2P}
	lines16[1].a.SetString("175935236032103269343375946024283166710560620391198844633277396555567531758612875037785990437334766668107363736967")
	lines16[1].b.SetString("236848165819904464687783593859859323610430407448167350950413296931613834965794012549522186594903892882952053511479")
	// l_{4P,4P}
	lines16[2].a.SetString("112478608797383214241485104901425712751775004500176940507327430758489657972141195128249787755176905310338369153303")
	lines16[2].b.SetString("112478608797383214241485104901425712751775004500176940507327430758489657972141195128249787755176905310338369153303")
	// v_{2P}
	lines16[3].b.SetString("48230533005596826793855842116095834521794558834178910932061719547839888944383298695107214925898849215589302992700")
	// v_{4P}
	lines16[4].b.SetString("228097355113300204138531148905234651262148041026195375645000724271212049151994375092458297304264351187709081232385")
	// v_{8P}
	lines16[5].b.SetOne()

	// P = (
	// 	   0x7881bc38aea0c08601e9b87c6e6dffca55f96ce64faa8a54fcc76146dee6a5468ef503cfeac9f4c169e308abc0746,
	// 	   0xcf189de61c5daaddcb28032fd26caae5b8b0e27bc3ba0d4b4b357ac3343fadc164b2581f32600d6be5df353018380b,
	// ) of order 7
	//
	// l_{P,P}
	lines7[0].a.SetString("102548817759151958053854624226127297976516351221993975628189733611446056212576494970386236978637761344561627738139")
	lines7[0].b.SetString("17030867045621927057629273133228761122091739345682975116557206480246500156223480777129261588141939420596970482342")
	// l_{2P,P}
	lines7[1].a.SetZero()
	lines7[1].b.SetString("134152643881269286775097829598623639630244859176246873534346527217779962437931261386223067655864012941095227344886")
	// l_{3P,3P}
	lines7[2].a.SetString("41739981286577398547281522846729756657097943860142380000234932931732311129520705533178816526268710200851448556869")
	lines7[2].b.SetString("241633558967347166953023460561664772414301773409231685423327056186473968192117341997839626551431420703843350975835")
	// v_{2P}
	lines7[3].b.SetString("204744603551330385445560097166282514611054159489724542883675765492852590078480732567937474615634803114031397601760")
	// v_{3P}
	lines7[4].b.SetString("58448067266526559821346123515157347081128640428262535756251141576595431412326053326376932871054281490576697066855")

	// P = (
	// 	   0x0x93a9ac62f3fa678eb04330e6973412e364875495c12a889ff93652e223784c93b2e63c6440cf9ab49ef3bd91c6f0b3,
	// 	   0x0x197214129826b8308cc6c9eb538506f5e5b9955731f6b654958e1d4134fbb6792508aa63c152574937d81ab9574c4da,
	// ) of order 13
	//
	// l_{P,P}
	lines13[0].a.SetString("183579943793544804897476052857594258755927986483725121520421281946549956407388028820239156530887269727288254885281")
	lines13[0].b.SetString("97063117930000100298119690927243291363224138368031034609964133138462297929511222884901733844308712431010578723049")
	// l_{2P,P}
	lines13[1].a.SetString("180035543548967714696397569102579437592770312927595736230857439639951087408103141540164300303358323539791676553053")
	lines13[1].b.SetString("210126715578130270696858689901060712003078807911225921583698281160874216530440284717694059485909536402363871525298")
	// l_{3P,3P}
	lines13[2].a.SetString("4282054308984930878992226539979437885658694062438259447982317858786230548397814401985288490821580648342517728752")
	lines13[2].b.SetString("127318834423694657140856755906867202309831593130954153295642640165587379046580483351497220667165568133918610893638")
	// l_{6P,6P}
	lines13[3].a.SetString("22202681511834823051272872251736326029858955988451643352298670507253383954515418615323155224049137987765327070690")
	lines13[3].b.SetString("214072754304739335135165757692888637954951222476496042054929387894092676271856840132869780631338519011959933234102")
	// v_{2P}
	lines13[4].b.SetString("244726416345920407590548041236872672250773857846891814638438483340822716811062521184620839895813441502845271856410")
	// v_{3P}
	lines13[5].b.SetString("54566584751129659042206658331333639595618281273806595914801276316885696696917410056322751824942598554030653131398")
	// v_{6P}
	lines13[6].b.SetString("178759099517859104288140376745009522978072207792085573130942352457759231448263986053703503899198102036110549151122")

}

// expByp3 uses a short addition chain to compute x^p3 where p3=(p-1)/3 .
func expByp3(x *fp.Element) *fp.Element {
	// Operations: 370 squares 62 multiplies
	//
	// Generated by github.com/mmcloughlin/addchain v0.4.0.

	// Allocate Temporaries.
	var z = new(fp.Element)
	var (
		t0  = new(fp.Element)
		t1  = new(fp.Element)
		t2  = new(fp.Element)
		t3  = new(fp.Element)
		t4  = new(fp.Element)
		t5  = new(fp.Element)
		t6  = new(fp.Element)
		t7  = new(fp.Element)
		t8  = new(fp.Element)
		t9  = new(fp.Element)
		t10 = new(fp.Element)
		t11 = new(fp.Element)
		t12 = new(fp.Element)
		t13 = new(fp.Element)
		t14 = new(fp.Element)
		t15 = new(fp.Element)
		t16 = new(fp.Element)
		t17 = new(fp.Element)
	)

	// Step 1: t15 = x^0x2
	t15.Square(x)

	// Step 2: t11 = x^0x3
	t11.Mul(x, t15)

	// Step 3: t2 = x^0x4
	t2.Mul(x, t11)

	// Step 4: t8 = x^0x5
	t8.Mul(x, t2)

	// Step 5: t3 = x^0x6
	t3.Mul(x, t8)

	// Step 6: t12 = x^0x9
	t12.Mul(t11, t3)

	// Step 7: z = x^0xb
	z.Mul(t15, t12)

	// Step 8: t0 = x^0xc
	t0.Mul(x, z)

	// Step 9: t1 = x^0xd
	t1.Mul(x, t0)

	// Step 10: t5 = x^0xf
	t5.Mul(t15, t1)

	// Step 11: t0 = x^0x1b
	t0.Mul(t0, t5)

	// Step 12: t13 = x^0x1d
	t13.Mul(t15, t0)

	// Step 13: t16 = x^0x23
	t16.Mul(t3, t13)

	// Step 14: t7 = x^0x25
	t7.Mul(t15, t16)

	// Step 15: t6 = x^0x29
	t6.Mul(t2, t7)

	// Step 16: t14 = x^0x2f
	t14.Mul(t3, t6)

	// Step 17: t4 = x^0x31
	t4.Mul(t15, t14)

	// Step 18: t9 = x^0x35
	t9.Mul(t2, t4)

	// Step 19: t2 = x^0x39
	t2.Mul(t2, t9)

	// Step 20: t10 = x^0x3b
	t10.Mul(t15, t2)

	// Step 21: t3 = x^0x3d
	t3.Mul(t15, t10)

	// Step 22: t17 = x^0x46
	t17.Mul(t12, t3)

	// Step 26: t17 = x^0x460
	for s := 0; s < 4; s++ {
		t17.Square(t17)
	}

	// Step 27: t17 = x^0x47b
	t17.Mul(t0, t17)

	// Step 34: t17 = x^0x23d80
	for s := 0; s < 7; s++ {
		t17.Square(t17)
	}

	// Step 35: t16 = x^0x23da3
	t16.Mul(t16, t17)

	// Step 40: t16 = x^0x47b460
	for s := 0; s < 5; s++ {
		t16.Square(t16)
	}

	// Step 41: t16 = x^0x47b461
	t16.Mul(x, t16)

	// Step 53: t16 = x^0x47b461000
	for s := 0; s < 12; s++ {
		t16.Square(t16)
	}

	// Step 54: t15 = x^0x47b461002
	t15.Mul(t15, t16)

	// Step 55: t15 = x^0x47b46103f
	t15.Mul(t3, t15)

	// Step 62: t15 = x^0x23da3081f80
	for s := 0; s < 7; s++ {
		t15.Square(t15)
	}

	// Step 63: t15 = x^0x23da3081fb1
	t15.Mul(t4, t15)

	// Step 70: t15 = x^0x11ed1840fd880
	for s := 0; s < 7; s++ {
		t15.Square(t15)
	}

	// Step 71: t15 = x^0x11ed1840fd8b5
	t15.Mul(t9, t15)

	// Step 75: t15 = x^0x11ed1840fd8b50
	for s := 0; s < 4; s++ {
		t15.Square(t15)
	}

	// Step 76: t15 = x^0x11ed1840fd8b5f
	t15.Mul(t5, t15)

	// Step 85: t15 = x^0x23da3081fb16be00
	for s := 0; s < 9; s++ {
		t15.Square(t15)
	}

	// Step 86: t15 = x^0x23da3081fb16be3b
	t15.Mul(t10, t15)

	// Step 94: t15 = x^0x23da3081fb16be3b00
	for s := 0; s < 8; s++ {
		t15.Square(t15)
	}

	// Step 95: t14 = x^0x23da3081fb16be3b2f
	t14.Mul(t14, t15)

	// Step 101: t14 = x^0x8f68c207ec5af8ecbc0
	for s := 0; s < 6; s++ {
		t14.Square(t14)
	}

	// Step 102: t14 = x^0x8f68c207ec5af8ecbe5
	t14.Mul(t7, t14)

	// Step 108: t14 = x^0x23da3081fb16be3b2f940
	for s := 0; s < 6; s++ {
		t14.Square(t14)
	}

	// Step 109: t13 = x^0x23da3081fb16be3b2f95d
	t13.Mul(t13, t14)

	// Step 121: t13 = x^0x23da3081fb16be3b2f95d000
	for s := 0; s < 12; s++ {
		t13.Square(t13)
	}

	// Step 122: t12 = x^0x23da3081fb16be3b2f95d009
	t12.Mul(t12, t13)

	// Step 132: t12 = x^0x8f68c207ec5af8ecbe57402400
	for s := 0; s < 10; s++ {
		t12.Square(t12)
	}

	// Step 133: t12 = x^0x8f68c207ec5af8ecbe57402435
	t12.Mul(t9, t12)

	// Step 135: t12 = x^0x23da3081fb16be3b2f95d0090d4
	for s := 0; s < 2; s++ {
		t12.Square(t12)
	}

	// Step 136: t11 = x^0x23da3081fb16be3b2f95d0090d7
	t11.Mul(t11, t12)

	// Step 146: t11 = x^0x8f68c207ec5af8ecbe57402435c00
	for s := 0; s < 10; s++ {
		t11.Square(t11)
	}

	// Step 147: t11 = x^0x8f68c207ec5af8ecbe57402435c31
	t11.Mul(t4, t11)

	// Step 155: t11 = x^0x8f68c207ec5af8ecbe57402435c3100
	for s := 0; s < 8; s++ {
		t11.Square(t11)
	}

	// Step 156: t10 = x^0x8f68c207ec5af8ecbe57402435c313b
	t10.Mul(t10, t11)

	// Step 163: t10 = x^0x47b46103f62d7c765f2ba0121ae189d80
	for s := 0; s < 7; s++ {
		t10.Square(t10)
	}

	// Step 164: t10 = x^0x47b46103f62d7c765f2ba0121ae189d9b
	t10.Mul(t0, t10)

	// Step 173: t10 = x^0x8f68c207ec5af8ecbe57402435c313b3600
	for s := 0; s < 9; s++ {
		t10.Square(t10)
	}

	// Step 174: t10 = x^0x8f68c207ec5af8ecbe57402435c313b360f
	t10.Mul(t5, t10)

	// Step 182: t10 = x^0x8f68c207ec5af8ecbe57402435c313b360f00
	for s := 0; s < 8; s++ {
		t10.Square(t10)
	}

	// Step 183: t9 = x^0x8f68c207ec5af8ecbe57402435c313b360f35
	t9.Mul(t9, t10)

	// Step 187: t9 = x^0x8f68c207ec5af8ecbe57402435c313b360f350
	for s := 0; s < 4; s++ {
		t9.Square(t9)
	}

	// Step 188: t9 = x^0x8f68c207ec5af8ecbe57402435c313b360f351
	t9.Mul(x, t9)

	// Step 200: t9 = x^0x8f68c207ec5af8ecbe57402435c313b360f351000
	for s := 0; s < 12; s++ {
		t9.Square(t9)
	}

	// Step 201: t8 = x^0x8f68c207ec5af8ecbe57402435c313b360f351005
	t8.Mul(t8, t9)

	// Step 209: t8 = x^0x8f68c207ec5af8ecbe57402435c313b360f35100500
	for s := 0; s < 8; s++ {
		t8.Square(t8)
	}

	// Step 210: t8 = x^0x8f68c207ec5af8ecbe57402435c313b360f3510051b
	t8.Mul(t0, t8)

	// Step 219: t8 = x^0x11ed1840fd8b5f1d97cae80486b862766c1e6a200a3600
	for s := 0; s < 9; s++ {
		t8.Square(t8)
	}

	// Step 220: t7 = x^0x11ed1840fd8b5f1d97cae80486b862766c1e6a200a3625
	t7.Mul(t7, t8)

	// Step 226: t7 = x^0x47b46103f62d7c765f2ba0121ae189d9b079a88028d8940
	for s := 0; s < 6; s++ {
		t7.Square(t7)
	}

	// Step 227: t7 = x^0x47b46103f62d7c765f2ba0121ae189d9b079a88028d897d
	t7.Mul(t3, t7)

	// Step 233: t7 = x^0x11ed1840fd8b5f1d97cae80486b862766c1e6a200a3625f40
	for s := 0; s < 6; s++ {
		t7.Square(t7)
	}

	// Step 234: t7 = x^0x11ed1840fd8b5f1d97cae80486b862766c1e6a200a3625f69
	t7.Mul(t6, t7)

	// Step 240: t7 = x^0x47b46103f62d7c765f2ba0121ae189d9b079a88028d897da40
	for s := 0; s < 6; s++ {
		t7.Square(t7)
	}

	// Step 241: t7 = x^0x47b46103f62d7c765f2ba0121ae189d9b079a88028d897da7d
	t7.Mul(t3, t7)

	// Step 247: t7 = x^0x11ed1840fd8b5f1d97cae80486b862766c1e6a200a3625f69f40
	for s := 0; s < 6; s++ {
		t7.Square(t7)
	}

	// Step 248: t7 = x^0x11ed1840fd8b5f1d97cae80486b862766c1e6a200a3625f69f79
	t7.Mul(t2, t7)

	// Step 254: t7 = x^0x47b46103f62d7c765f2ba0121ae189d9b079a88028d897da7de40
	for s := 0; s < 6; s++ {
		t7.Square(t7)
	}

	// Step 255: t7 = x^0x47b46103f62d7c765f2ba0121ae189d9b079a88028d897da7de5b
	t7.Mul(t0, t7)

	// Step 263: t7 = x^0x47b46103f62d7c765f2ba0121ae189d9b079a88028d897da7de5b00
	for s := 0; s < 8; s++ {
		t7.Square(t7)
	}

	// Step 264: t6 = x^0x47b46103f62d7c765f2ba0121ae189d9b079a88028d897da7de5b29
	t6.Mul(t6, t7)

	// Step 268: t6 = x^0x47b46103f62d7c765f2ba0121ae189d9b079a88028d897da7de5b290
	for s := 0; s < 4; s++ {
		t6.Square(t6)
	}

	// Step 269: t5 = x^0x47b46103f62d7c765f2ba0121ae189d9b079a88028d897da7de5b29f
	t5.Mul(t5, t6)

	// Step 282: t5 = x^0x8f68c207ec5af8ecbe57402435c313b360f3510051b12fb4fbcb653e000
	for s := 0; s < 13; s++ {
		t5.Square(t5)
	}

	// Step 283: t4 = x^0x8f68c207ec5af8ecbe57402435c313b360f3510051b12fb4fbcb653e031
	t4.Mul(t4, t5)

	// Step 284: t4 = x^0x11ed1840fd8b5f1d97cae80486b862766c1e6a200a3625f69f796ca7c062
	t4.Square(t4)

	// Step 285: t4 = x^0x11ed1840fd8b5f1d97cae80486b862766c1e6a200a3625f69f796ca7c063
	t4.Mul(x, t4)

	// Step 307: t4 = x^0x47b46103f62d7c765f2ba0121ae189d9b079a88028d897da7de5b29f018c00000
	for s := 0; s < 22; s++ {
		t4.Square(t4)
	}

	// Step 308: t3 = x^0x47b46103f62d7c765f2ba0121ae189d9b079a88028d897da7de5b29f018c0003d
	t3.Mul(t3, t4)

	// Step 315: t3 = x^0x23da3081fb16be3b2f95d0090d70c4ecd83cd440146c4bed3ef2d94f80c60001e80
	for s := 0; s < 7; s++ {
		t3.Square(t3)
	}

	// Step 316: t2 = x^0x23da3081fb16be3b2f95d0090d70c4ecd83cd440146c4bed3ef2d94f80c60001eb9
	t2.Mul(t2, t3)

	// Step 320: t2 = x^0x23da3081fb16be3b2f95d0090d70c4ecd83cd440146c4bed3ef2d94f80c60001eb90
	for s := 0; s < 4; s++ {
		t2.Square(t2)
	}

	// Step 321: t1 = x^0x23da3081fb16be3b2f95d0090d70c4ecd83cd440146c4bed3ef2d94f80c60001eb9d
	t1.Mul(t1, t2)

	// Step 329: t1 = x^0x23da3081fb16be3b2f95d0090d70c4ecd83cd440146c4bed3ef2d94f80c60001eb9d00
	for s := 0; s < 8; s++ {
		t1.Square(t1)
	}

	// Step 330: t0 = x^0x23da3081fb16be3b2f95d0090d70c4ecd83cd440146c4bed3ef2d94f80c60001eb9d1b
	t0.Mul(t0, t1)

	// Step 336: t0 = x^0x8f68c207ec5af8ecbe57402435c313b360f3510051b12fb4fbcb653e03180007ae746c0
	for s := 0; s < 6; s++ {
		t0.Square(t0)
	}

	// Step 337: t0 = x^0x8f68c207ec5af8ecbe57402435c313b360f3510051b12fb4fbcb653e03180007ae746c1
	t0.Mul(x, t0)

	// Step 371: t0 = x^0x23da3081fb16be3b2f95d0090d70c4ecd83cd440146c4bed3ef2d94f80c60001eb9d1b0400000000
	for s := 0; s < 34; s++ {
		t0.Square(t0)
	}

	// Step 372: t0 = x^0x23da3081fb16be3b2f95d0090d70c4ecd83cd440146c4bed3ef2d94f80c60001eb9d1b040000000b
	t0.Mul(z, t0)

	// Step 379: t0 = x^0x11ed1840fd8b5f1d97cae80486b862766c1e6a200a3625f69f796ca7c0630000f5ce8d820000000580
	for s := 0; s < 7; s++ {
		t0.Square(t0)
	}

	// Step 380: z = x^0x11ed1840fd8b5f1d97cae80486b862766c1e6a200a3625f69f796ca7c0630000f5ce8d82000000058b
	z.Mul(z, t0)

	// Step 385: z = x^0x23da3081fb16be3b2f95d0090d70c4ecd83cd440146c4bed3ef2d94f80c60001eb9d1b040000000b160
	for s := 0; s < 5; s++ {
		z.Square(z)
	}

	// Step 386: z = x^0x23da3081fb16be3b2f95d0090d70c4ecd83cd440146c4bed3ef2d94f80c60001eb9d1b040000000b161
	z.Mul(x, z)

	// Step 432: z = x^0x8f68c207ec5af8ecbe57402435c313b360f3510051b12fb4fbcb653e03180007ae746c100000002c58400000000000
	for s := 0; s < 46; s++ {
		z.Square(z)
	}

	return z
}

// expByp16 uses a short addition chain to compute x^p16 where p16=(p-1)/16 .
func expByp16(x *fp.Element) *fp.Element {
	// Operations: 368 squares 62 multiplies
	//
	// Generated by github.com/mmcloughlin/addchain v0.4.0.

	// Allocate Temporaries.
	var z = new(fp.Element)
	var (
		t0  = new(fp.Element)
		t1  = new(fp.Element)
		t2  = new(fp.Element)
		t3  = new(fp.Element)
		t4  = new(fp.Element)
		t5  = new(fp.Element)
		t6  = new(fp.Element)
		t7  = new(fp.Element)
		t8  = new(fp.Element)
		t9  = new(fp.Element)
		t10 = new(fp.Element)
		t11 = new(fp.Element)
	)

	// Step 1: t6 = x^0x2
	t6.Square(x)

	// Step 2: t1 = x^0x3
	t1.Mul(x, t6)

	// Step 3: t5 = x^0x4
	t5.Mul(x, t1)

	// Step 4: t0 = x^0x5
	t0.Mul(x, t5)

	// Step 5: t9 = x^0x7
	t9.Mul(t6, t0)

	// Step 6: t4 = x^0x9
	t4.Mul(t6, t9)

	// Step 7: t3 = x^0xb
	t3.Mul(t6, t4)

	// Step 8: t8 = x^0xf
	t8.Mul(t5, t3)

	// Step 9: z = x^0x11
	z.Mul(t6, t8)

	// Step 10: t10 = x^0x13
	t10.Mul(t6, z)

	// Step 11: t2 = x^0x17
	t2.Mul(t5, t10)

	// Step 12: t7 = x^0x1b
	t7.Mul(t5, t2)

	// Step 13: t5 = x^0x1d
	t5.Mul(t6, t7)

	// Step 14: t6 = x^0x1f
	t6.Mul(t6, t5)

	// Step 15: t11 = x^0x34
	t11.Mul(t2, t5)

	// Step 17: t11 = x^0xd0
	for s := 0; s < 2; s++ {
		t11.Square(t11)
	}

	// Step 18: t11 = x^0xd7
	t11.Mul(t9, t11)

	// Step 26: t11 = x^0xd700
	for s := 0; s < 8; s++ {
		t11.Square(t11)
	}

	// Step 27: t11 = x^0xd71d
	t11.Mul(t5, t11)

	// Step 34: t11 = x^0x6b8e80
	for s := 0; s < 7; s++ {
		t11.Square(t11)
	}

	// Step 35: t11 = x^0x6b8e91
	t11.Mul(z, t11)

	// Step 36: t11 = x^0xd71d22
	t11.Square(t11)

	// Step 37: t11 = x^0xd71d23
	t11.Mul(x, t11)

	// Step 46: t11 = x^0x1ae3a4600
	for s := 0; s < 9; s++ {
		t11.Square(t11)
	}

	// Step 47: t11 = x^0x1ae3a4617
	t11.Mul(t2, t11)

	// Step 49: t11 = x^0x6b8e9185c
	for s := 0; s < 2; s++ {
		t11.Square(t11)
	}

	// Step 50: t11 = x^0x6b8e9185f
	t11.Mul(t1, t11)

	// Step 56: t11 = x^0x1ae3a4617c0
	for s := 0; s < 6; s++ {
		t11.Square(t11)
	}

	// Step 57: t11 = x^0x1ae3a4617c5
	t11.Mul(t0, t11)

	// Step 61: t11 = x^0x1ae3a4617c50
	for s := 0; s < 4; s++ {
		t11.Square(t11)
	}

	// Step 62: t11 = x^0x1ae3a4617c51
	t11.Mul(x, t11)

	// Step 71: t11 = x^0x35c748c2f8a200
	for s := 0; s < 9; s++ {
		t11.Square(t11)
	}

	// Step 72: t11 = x^0x35c748c2f8a21d
	t11.Mul(t5, t11)

	// Step 77: t11 = x^0x6b8e9185f1443a0
	for s := 0; s < 5; s++ {
		t11.Square(t11)
	}

	// Step 78: t11 = x^0x6b8e9185f1443ab
	t11.Mul(t3, t11)

	// Step 83: t11 = x^0xd71d230be2887560
	for s := 0; s < 5; s++ {
		t11.Square(t11)
	}

	// Step 84: t11 = x^0xd71d230be2887563
	t11.Mul(t1, t11)

	// Step 92: t11 = x^0xd71d230be288756300
	for s := 0; s < 8; s++ {
		t11.Square(t11)
	}

	// Step 93: t11 = x^0xd71d230be28875631d
	t11.Mul(t5, t11)

	// Step 94: t11 = x^0x1ae3a4617c510eac63a
	t11.Square(t11)

	// Step 95: t11 = x^0x1ae3a4617c510eac63b
	t11.Mul(x, t11)

	// Step 105: t11 = x^0x6b8e9185f1443ab18ec00
	for s := 0; s < 10; s++ {
		t11.Square(t11)
	}

	// Step 106: t11 = x^0x6b8e9185f1443ab18ec17
	t11.Mul(t2, t11)

	// Step 118: t11 = x^0x6b8e9185f1443ab18ec17000
	for s := 0; s < 12; s++ {
		t11.Square(t11)
	}

	// Step 119: t11 = x^0x6b8e9185f1443ab18ec1701b
	t11.Mul(t7, t11)

	// Step 124: t11 = x^0xd71d230be28875631d82e0360
	for s := 0; s < 5; s++ {
		t11.Square(t11)
	}

	// Step 125: t11 = x^0xd71d230be28875631d82e0365
	t11.Mul(t0, t11)

	// Step 132: t11 = x^0x6b8e9185f1443ab18ec1701b280
	for s := 0; s < 7; s++ {
		t11.Square(t11)
	}

	// Step 133: t11 = x^0x6b8e9185f1443ab18ec1701b285
	t11.Mul(t0, t11)

	// Step 139: t11 = x^0x1ae3a4617c510eac63b05c06ca140
	for s := 0; s < 6; s++ {
		t11.Square(t11)
	}

	// Step 140: t11 = x^0x1ae3a4617c510eac63b05c06ca149
	t11.Mul(t4, t11)

	// Step 147: t11 = x^0xd71d230be28875631d82e03650a480
	for s := 0; s < 7; s++ {
		t11.Square(t11)
	}

	// Step 148: t11 = x^0xd71d230be28875631d82e03650a49d
	t11.Mul(t5, t11)

	// Step 153: t11 = x^0x1ae3a4617c510eac63b05c06ca1493a0
	for s := 0; s < 5; s++ {
		t11.Square(t11)
	}

	// Step 154: t11 = x^0x1ae3a4617c510eac63b05c06ca1493b1
	t11.Mul(z, t11)

	// Step 157: t11 = x^0xd71d230be28875631d82e03650a49d88
	for s := 0; s < 3; s++ {
		t11.Square(t11)
	}

	// Step 158: t11 = x^0xd71d230be28875631d82e03650a49d8d
	t11.Mul(t0, t11)

	// Step 166: t11 = x^0xd71d230be28875631d82e03650a49d8d00
	for s := 0; s < 8; s++ {
		t11.Square(t11)
	}

	// Step 167: t11 = x^0xd71d230be28875631d82e03650a49d8d11
	t11.Mul(z, t11)

	// Step 173: t11 = x^0x35c748c2f8a21d58c760b80d942927634440
	for s := 0; s < 6; s++ {
		t11.Square(t11)
	}

	// Step 174: t11 = x^0x35c748c2f8a21d58c760b80d94292763445b
	t11.Mul(t7, t11)

	// Step 181: t11 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d80
	for s := 0; s < 7; s++ {
		t11.Square(t11)
	}

	// Step 182: t11 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f
	t11.Mul(t6, t11)

	// Step 186: t11 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f0
	for s := 0; s < 4; s++ {
		t11.Square(t11)
	}

	// Step 187: t11 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f3
	t11.Mul(t1, t11)

	// Step 199: t11 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f3000
	for s := 0; s < 12; s++ {
		t11.Square(t11)
	}

	// Step 200: t11 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f
	t11.Mul(t8, t11)

	// Step 204: t11 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f0
	for s := 0; s < 4; s++ {
		t11.Square(t11)
	}

	// Step 205: t11 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5
	t11.Mul(t0, t11)

	// Step 213: t11 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f500
	for s := 0; s < 8; s++ {
		t11.Square(t11)
	}

	// Step 214: t10 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f513
	t10.Mul(t10, t11)

	// Step 219: t10 = x^0x35c748c2f8a21d58c760b80d94292763445b3e601ea260
	for s := 0; s < 5; s++ {
		t10.Square(t10)
	}

	// Step 220: t10 = x^0x35c748c2f8a21d58c760b80d94292763445b3e601ea271
	t10.Mul(z, t10)

	// Step 223: t10 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f51388
	for s := 0; s < 3; s++ {
		t10.Square(t10)
	}

	// Step 224: t9 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f
	t9.Mul(t9, t10)

	// Step 231: t9 = x^0xd71d230be28875631d82e03650a49d8d116cf9807a89c780
	for s := 0; s < 7; s++ {
		t9.Square(t9)
	}

	// Step 232: t9 = x^0xd71d230be28875631d82e03650a49d8d116cf9807a89c78f
	t9.Mul(t8, t9)

	// Step 237: t9 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1e0
	for s := 0; s < 5; s++ {
		t9.Square(t9)
	}

	// Step 238: t8 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef
	t8.Mul(t8, t9)

	// Step 245: t8 = x^0xd71d230be28875631d82e03650a49d8d116cf9807a89c78f780
	for s := 0; s < 7; s++ {
		t8.Square(t8)
	}

	// Step 246: t7 = x^0xd71d230be28875631d82e03650a49d8d116cf9807a89c78f79b
	t7.Mul(t7, t8)

	// Step 254: t7 = x^0xd71d230be28875631d82e03650a49d8d116cf9807a89c78f79b00
	for s := 0; s < 8; s++ {
		t7.Square(t7)
	}

	// Step 255: t7 = x^0xd71d230be28875631d82e03650a49d8d116cf9807a89c78f79b11
	t7.Mul(z, t7)

	// Step 261: t7 = x^0x35c748c2f8a21d58c760b80d94292763445b3e601ea271e3de6c440
	for s := 0; s < 6; s++ {
		t7.Square(t7)
	}

	// Step 262: t6 = x^0x35c748c2f8a21d58c760b80d94292763445b3e601ea271e3de6c45f
	t6.Mul(t6, t7)

	// Step 268: t6 = x^0xd71d230be28875631d82e03650a49d8d116cf9807a89c78f79b117c0
	for s := 0; s < 6; s++ {
		t6.Square(t6)
	}

	// Step 269: t5 = x^0xd71d230be28875631d82e03650a49d8d116cf9807a89c78f79b117dd
	t5.Mul(t5, t6)

	// Step 278: t5 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba00
	for s := 0; s < 9; s++ {
		t5.Square(t5)
	}

	// Step 279: t5 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba09
	t5.Mul(t4, t5)

	// Step 284: t5 = x^0x35c748c2f8a21d58c760b80d94292763445b3e601ea271e3de6c45f74120
	for s := 0; s < 5; s++ {
		t5.Square(t5)
	}

	// Step 285: t4 = x^0x35c748c2f8a21d58c760b80d94292763445b3e601ea271e3de6c45f74129
	t4.Mul(t4, t5)

	// Step 304: t4 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba09480000
	for s := 0; s < 19; s++ {
		t4.Square(t4)
	}

	// Step 305: t4 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba09480017
	t4.Mul(t2, t4)

	// Step 313: t4 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba0948001700
	for s := 0; s < 8; s++ {
		t4.Square(t4)
	}

	// Step 314: t3 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba094800170b
	t3.Mul(t3, t4)

	// Step 320: t3 = x^0x6b8e9185f1443ab18ec1701b28524ec688b67cc03d44e3c7bcd88bee82520005c2c0
	for s := 0; s < 6; s++ {
		t3.Square(t3)
	}

	// Step 321: t2 = x^0x6b8e9185f1443ab18ec1701b28524ec688b67cc03d44e3c7bcd88bee82520005c2d7
	t2.Mul(t2, t3)

	// Step 325: t2 = x^0x6b8e9185f1443ab18ec1701b28524ec688b67cc03d44e3c7bcd88bee82520005c2d70
	for s := 0; s < 4; s++ {
		t2.Square(t2)
	}

	// Step 326: t2 = x^0x6b8e9185f1443ab18ec1701b28524ec688b67cc03d44e3c7bcd88bee82520005c2d75
	t2.Mul(t0, t2)

	// Step 330: t2 = x^0x6b8e9185f1443ab18ec1701b28524ec688b67cc03d44e3c7bcd88bee82520005c2d750
	for s := 0; s < 4; s++ {
		t2.Square(t2)
	}

	// Step 331: t2 = x^0x6b8e9185f1443ab18ec1701b28524ec688b67cc03d44e3c7bcd88bee82520005c2d751
	t2.Mul(x, t2)

	// Step 337: t2 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba094800170b5d440
	for s := 0; s < 6; s++ {
		t2.Square(t2)
	}

	// Step 338: t1 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba094800170b5d443
	t1.Mul(t1, t2)

	// Step 367: t1 = x^0x35c748c2f8a21d58c760b80d94292763445b3e601ea271e3de6c45f741290002e16ba8860000000
	for s := 0; s < 29; s++ {
		t1.Square(t1)
	}

	// Step 368: t1 = x^0x35c748c2f8a21d58c760b80d94292763445b3e601ea271e3de6c45f741290002e16ba8860000001
	t1.Mul(x, t1)

	// Step 375: t1 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba094800170b5d443000000080
	for s := 0; s < 7; s++ {
		t1.Square(t1)
	}

	// Step 376: t0 = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba094800170b5d443000000085
	t0.Mul(t0, t1)

	// Step 385: t0 = x^0x35c748c2f8a21d58c760b80d94292763445b3e601ea271e3de6c45f741290002e16ba88600000010a00
	for s := 0; s < 9; s++ {
		t0.Square(t0)
	}

	// Step 386: z = x^0x35c748c2f8a21d58c760b80d94292763445b3e601ea271e3de6c45f741290002e16ba88600000010a11
	z.Mul(z, t0)

	// Step 387: z = x^0x6b8e9185f1443ab18ec1701b28524ec688b67cc03d44e3c7bcd88bee82520005c2d7510c00000021422
	z.Square(z)

	// Step 388: z = x^0x6b8e9185f1443ab18ec1701b28524ec688b67cc03d44e3c7bcd88bee82520005c2d7510c00000021423
	z.Mul(x, z)

	// Step 430: z = x^0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba094800170b5d44300000008508c0000000000
	for s := 0; s < 42; s++ {
		z.Square(z)
	}

	return z
}

// expByp7 uses a short addition chain to compute x^p7 where p7=(p-1)/7 .
func expByp7(x *fp.Element) *fp.Element {
	// Operations: 370 squares 58 multiplies
	//
	// Generated by github.com/mmcloughlin/addchain v0.4.0.

	// Allocate Temporaries.
	var z = new(fp.Element)
	var (
		t0  = new(fp.Element)
		t1  = new(fp.Element)
		t2  = new(fp.Element)
		t3  = new(fp.Element)
		t4  = new(fp.Element)
		t5  = new(fp.Element)
		t6  = new(fp.Element)
		t7  = new(fp.Element)
		t8  = new(fp.Element)
		t9  = new(fp.Element)
		t10 = new(fp.Element)
		t11 = new(fp.Element)
		t12 = new(fp.Element)
		t13 = new(fp.Element)
	)

	// Step 1: t11 = x^0x2
	t11.Square(x)

	// Step 2: t5 = x^0x3
	t5.Mul(x, t11)

	// Step 3: t0 = x^0x4
	t0.Mul(x, t5)

	// Step 4: z = x^0x5
	z.Mul(x, t0)

	// Step 5: t12 = x^0x7
	t12.Mul(t11, z)

	// Step 6: t6 = x^0x9
	t6.Mul(t11, t12)

	// Step 7: t4 = x^0xd
	t4.Mul(t0, t6)

	// Step 8: t8 = x^0xf
	t8.Mul(t11, t4)

	// Step 9: t2 = x^0x11
	t2.Mul(t11, t8)

	// Step 10: t0 = x^0x13
	t0.Mul(t11, t2)

	// Step 11: t3 = x^0x15
	t3.Mul(t11, t0)

	// Step 12: t1 = x^0x17
	t1.Mul(t11, t3)

	// Step 13: t9 = x^0x19
	t9.Mul(t11, t1)

	// Step 14: t10 = x^0x1b
	t10.Mul(t11, t9)

	// Step 15: t7 = x^0x1d
	t7.Mul(t11, t10)

	// Step 16: t13 = x^0x1e
	t13.Mul(x, t7)

	// Step 21: t13 = x^0x3c0
	for s := 0; s < 5; s++ {
		t13.Square(t13)
	}

	// Step 22: t13 = x^0x3d7
	t13.Mul(t1, t13)

	// Step 25: t13 = x^0x1eb8
	for s := 0; s < 3; s++ {
		t13.Square(t13)
	}

	// Step 26: t13 = x^0x1ebb
	t13.Mul(t5, t13)

	// Step 34: t13 = x^0x1ebb00
	for s := 0; s < 8; s++ {
		t13.Square(t13)
	}

	// Step 35: t13 = x^0x1ebb05
	t13.Mul(z, t13)

	// Step 47: t13 = x^0x1ebb05000
	for s := 0; s < 12; s++ {
		t13.Square(t13)
	}

	// Step 48: t13 = x^0x1ebb0501b
	t13.Mul(t10, t13)

	// Step 55: t13 = x^0xf5d8280d80
	for s := 0; s < 7; s++ {
		t13.Square(t13)
	}

	// Step 56: t13 = x^0xf5d8280d95
	t13.Mul(t3, t13)

	// Step 63: t13 = x^0x7aec1406ca80
	for s := 0; s < 7; s++ {
		t13.Square(t13)
	}

	// Step 64: t13 = x^0x7aec1406ca97
	t13.Mul(t1, t13)

	// Step 71: t13 = x^0x3d760a03654b80
	for s := 0; s < 7; s++ {
		t13.Square(t13)
	}

	// Step 72: t13 = x^0x3d760a03654b8f
	t13.Mul(t8, t13)

	// Step 74: t13 = x^0xf5d8280d952e3c
	for s := 0; s < 2; s++ {
		t13.Square(t13)
	}

	// Step 75: t13 = x^0xf5d8280d952e3d
	t13.Mul(x, t13)

	// Step 84: t13 = x^0x1ebb0501b2a5c7a00
	for s := 0; s < 9; s++ {
		t13.Square(t13)
	}

	// Step 85: t13 = x^0x1ebb0501b2a5c7a07
	t13.Mul(t12, t13)

	// Step 92: t13 = x^0xf5d8280d952e3d0380
	for s := 0; s < 7; s++ {
		t13.Square(t13)
	}

	// Step 93: t13 = x^0xf5d8280d952e3d038f
	t13.Mul(t8, t13)

	// Step 97: t13 = x^0xf5d8280d952e3d038f0
	for s := 0; s < 4; s++ {
		t13.Square(t13)
	}

	// Step 98: t12 = x^0xf5d8280d952e3d038f7
	t12.Mul(t12, t13)

	// Step 102: t12 = x^0xf5d8280d952e3d038f70
	for s := 0; s < 4; s++ {
		t12.Square(t12)
	}

	// Step 103: t12 = x^0xf5d8280d952e3d038f71
	t12.Mul(x, t12)

	// Step 118: t12 = x^0x7aec1406ca971e81c7b88000
	for s := 0; s < 15; s++ {
		t12.Square(t12)
	}

	// Step 119: t12 = x^0x7aec1406ca971e81c7b8801d
	t12.Mul(t7, t12)

	// Step 120: t11 = x^0x7aec1406ca971e81c7b8801f
	t11.Mul(t11, t12)

	// Step 129: t11 = x^0xf5d8280d952e3d038f71003e00
	for s := 0; s < 9; s++ {
		t11.Square(t11)
	}

	// Step 130: t11 = x^0xf5d8280d952e3d038f71003e13
	t11.Mul(t0, t11)

	// Step 140: t11 = x^0x3d760a03654b8f40e3dc400f84c00
	for s := 0; s < 10; s++ {
		t11.Square(t11)
	}

	// Step 141: t11 = x^0x3d760a03654b8f40e3dc400f84c15
	t11.Mul(t3, t11)

	// Step 149: t11 = x^0x3d760a03654b8f40e3dc400f84c1500
	for s := 0; s < 8; s++ {
		t11.Square(t11)
	}

	// Step 150: t11 = x^0x3d760a03654b8f40e3dc400f84c1519
	t11.Mul(t9, t11)

	// Step 153: t11 = x^0x1ebb0501b2a5c7a071ee2007c260a8c8
	for s := 0; s < 3; s++ {
		t11.Square(t11)
	}

	// Step 154: t11 = x^0x1ebb0501b2a5c7a071ee2007c260a8cb
	t11.Mul(t5, t11)

	// Step 165: t11 = x^0xf5d8280d952e3d038f71003e1305465800
	for s := 0; s < 11; s++ {
		t11.Square(t11)
	}

	// Step 166: t11 = x^0xf5d8280d952e3d038f71003e1305465813
	t11.Mul(t0, t11)

	// Step 171: t11 = x^0x1ebb0501b2a5c7a071ee2007c260a8cb0260
	for s := 0; s < 5; s++ {
		t11.Square(t11)
	}

	// Step 172: t11 = x^0x1ebb0501b2a5c7a071ee2007c260a8cb027d
	t11.Mul(t7, t11)

	// Step 177: t11 = x^0x3d760a03654b8f40e3dc400f84c1519604fa0
	for s := 0; s < 5; s++ {
		t11.Square(t11)
	}

	// Step 178: t11 = x^0x3d760a03654b8f40e3dc400f84c1519604fa9
	t11.Mul(t6, t11)

	// Step 188: t11 = x^0xf5d8280d952e3d038f71003e1305465813ea400
	for s := 0; s < 10; s++ {
		t11.Square(t11)
	}

	// Step 189: t10 = x^0xf5d8280d952e3d038f71003e1305465813ea41b
	t10.Mul(t10, t11)

	// Step 195: t10 = x^0x3d760a03654b8f40e3dc400f84c1519604fa906c0
	for s := 0; s < 6; s++ {
		t10.Square(t10)
	}

	// Step 196: t10 = x^0x3d760a03654b8f40e3dc400f84c1519604fa906dd
	t10.Mul(t7, t10)

	// Step 201: t10 = x^0x7aec1406ca971e81c7b8801f0982a32c09f520dba0
	for s := 0; s < 5; s++ {
		t10.Square(t10)
	}

	// Step 202: t10 = x^0x7aec1406ca971e81c7b8801f0982a32c09f520dbb3
	t10.Mul(t0, t10)

	// Step 207: t10 = x^0xf5d8280d952e3d038f71003e1305465813ea41b7660
	for s := 0; s < 5; s++ {
		t10.Square(t10)
	}

	// Step 208: t10 = x^0xf5d8280d952e3d038f71003e1305465813ea41b7677
	t10.Mul(t1, t10)

	// Step 213: t10 = x^0x1ebb0501b2a5c7a071ee2007c260a8cb027d4836ecee0
	for s := 0; s < 5; s++ {
		t10.Square(t10)
	}

	// Step 214: t10 = x^0x1ebb0501b2a5c7a071ee2007c260a8cb027d4836ecef1
	t10.Mul(t2, t10)

	// Step 219: t10 = x^0x3d760a03654b8f40e3dc400f84c1519604fa906dd9de20
	for s := 0; s < 5; s++ {
		t10.Square(t10)
	}

	// Step 220: t9 = x^0x3d760a03654b8f40e3dc400f84c1519604fa906dd9de39
	t9.Mul(t9, t10)

	// Step 230: t9 = x^0xf5d8280d952e3d038f71003e1305465813ea41b76778e400
	for s := 0; s < 10; s++ {
		t9.Square(t9)
	}

	// Step 231: t9 = x^0xf5d8280d952e3d038f71003e1305465813ea41b76778e411
	t9.Mul(t2, t9)

	// Step 236: t9 = x^0x1ebb0501b2a5c7a071ee2007c260a8cb027d4836ecef1c8220
	for s := 0; s < 5; s++ {
		t9.Square(t9)
	}

	// Step 237: t9 = x^0x1ebb0501b2a5c7a071ee2007c260a8cb027d4836ecef1c8235
	t9.Mul(t3, t9)

	// Step 241: t9 = x^0x1ebb0501b2a5c7a071ee2007c260a8cb027d4836ecef1c82350
	for s := 0; s < 4; s++ {
		t9.Square(t9)
	}

	// Step 242: t8 = x^0x1ebb0501b2a5c7a071ee2007c260a8cb027d4836ecef1c8235f
	t8.Mul(t8, t9)

	// Step 247: t8 = x^0x3d760a03654b8f40e3dc400f84c1519604fa906dd9de39046be0
	for s := 0; s < 5; s++ {
		t8.Square(t8)
	}

	// Step 248: t8 = x^0x3d760a03654b8f40e3dc400f84c1519604fa906dd9de39046be9
	t8.Mul(t6, t8)

	// Step 254: t8 = x^0xf5d8280d952e3d038f71003e1305465813ea41b76778e411afa40
	for s := 0; s < 6; s++ {
		t8.Square(t8)
	}

	// Step 255: t7 = x^0xf5d8280d952e3d038f71003e1305465813ea41b76778e411afa5d
	t7.Mul(t7, t8)

	// Step 262: t7 = x^0x7aec1406ca971e81c7b8801f0982a32c09f520dbb3bc7208d7d2e80
	for s := 0; s < 7; s++ {
		t7.Square(t7)
	}

	// Step 263: t7 = x^0x7aec1406ca971e81c7b8801f0982a32c09f520dbb3bc7208d7d2e91
	t7.Mul(t2, t7)

	// Step 271: t7 = x^0x7aec1406ca971e81c7b8801f0982a32c09f520dbb3bc7208d7d2e9100
	for s := 0; s < 8; s++ {
		t7.Square(t7)
	}

	// Step 272: t6 = x^0x7aec1406ca971e81c7b8801f0982a32c09f520dbb3bc7208d7d2e9109
	t6.Mul(t6, t7)

	// Step 278: t6 = x^0x1ebb0501b2a5c7a071ee2007c260a8cb027d4836ecef1c8235f4ba44240
	for s := 0; s < 6; s++ {
		t6.Square(t6)
	}

	// Step 279: t6 = x^0x1ebb0501b2a5c7a071ee2007c260a8cb027d4836ecef1c8235f4ba44253
	t6.Mul(t0, t6)

	// Step 281: t6 = x^0x7aec1406ca971e81c7b8801f0982a32c09f520dbb3bc7208d7d2e91094c
	for s := 0; s < 2; s++ {
		t6.Square(t6)
	}

	// Step 282: t5 = x^0x7aec1406ca971e81c7b8801f0982a32c09f520dbb3bc7208d7d2e91094f
	t5.Mul(t5, t6)

	// Step 303: t5 = x^0xf5d8280d952e3d038f71003e1305465813ea41b76778e411afa5d22129e00000
	for s := 0; s < 21; s++ {
		t5.Square(t5)
	}

	// Step 304: t4 = x^0xf5d8280d952e3d038f71003e1305465813ea41b76778e411afa5d22129e0000d
	t4.Mul(t4, t5)

	// Step 311: t4 = x^0x7aec1406ca971e81c7b8801f0982a32c09f520dbb3bc7208d7d2e91094f0000680
	for s := 0; s < 7; s++ {
		t4.Square(t4)
	}

	// Step 312: t3 = x^0x7aec1406ca971e81c7b8801f0982a32c09f520dbb3bc7208d7d2e91094f0000695
	t3.Mul(t3, t4)

	// Step 317: t3 = x^0xf5d8280d952e3d038f71003e1305465813ea41b76778e411afa5d22129e0000d2a0
	for s := 0; s < 5; s++ {
		t3.Square(t3)
	}

	// Step 318: t2 = x^0xf5d8280d952e3d038f71003e1305465813ea41b76778e411afa5d22129e0000d2b1
	t2.Mul(t2, t3)

	// Step 327: t2 = x^0x1ebb0501b2a5c7a071ee2007c260a8cb027d4836ecef1c8235f4ba44253c0001a56200
	for s := 0; s < 9; s++ {
		t2.Square(t2)
	}

	// Step 328: t1 = x^0x1ebb0501b2a5c7a071ee2007c260a8cb027d4836ecef1c8235f4ba44253c0001a56217
	t1.Mul(t1, t2)

	// Step 333: t1 = x^0x3d760a03654b8f40e3dc400f84c1519604fa906dd9de39046be974884a7800034ac42e0
	for s := 0; s < 5; s++ {
		t1.Square(t1)
	}

	// Step 334: t1 = x^0x3d760a03654b8f40e3dc400f84c1519604fa906dd9de39046be974884a7800034ac42e5
	t1.Mul(z, t1)

	// Step 370: t1 = x^0x3d760a03654b8f40e3dc400f84c1519604fa906dd9de39046be974884a7800034ac42e5000000000
	for s := 0; s < 36; s++ {
		t1.Square(t1)
	}

	// Step 371: t0 = x^0x3d760a03654b8f40e3dc400f84c1519604fa906dd9de39046be974884a7800034ac42e5000000013
	t0.Mul(t0, t1)

	// Step 381: t0 = x^0xf5d8280d952e3d038f71003e1305465813ea41b76778e411afa5d22129e0000d2b10b9400000004c00
	for s := 0; s < 10; s++ {
		t0.Square(t0)
	}

	// Step 382: z = x^0xf5d8280d952e3d038f71003e1305465813ea41b76778e411afa5d22129e0000d2b10b9400000004c05
	z.Mul(z, t0)

	// Step 428: z = x^0x3d760a03654b8f40e3dc400f84c1519604fa906dd9de39046be974884a7800034ac42e500000001301400000000000
	for s := 0; s < 46; s++ {
		z.Square(z)
	}

	return z
}

// expByp13 uses a short addition chain to compute x^p13 where p13=(p-1)/13 .
func expByp13(x *fp.Element) *fp.Element {
	// Operations: 369 squares 62 multiplies
	//
	// Generated by github.com/mmcloughlin/addchain v0.4.0.

	// Allocate Temporaries.
	var z = new(fp.Element)
	var (
		t0  = new(fp.Element)
		t1  = new(fp.Element)
		t2  = new(fp.Element)
		t3  = new(fp.Element)
		t4  = new(fp.Element)
		t5  = new(fp.Element)
		t6  = new(fp.Element)
		t7  = new(fp.Element)
		t8  = new(fp.Element)
		t9  = new(fp.Element)
		t10 = new(fp.Element)
		t11 = new(fp.Element)
		t12 = new(fp.Element)
		t13 = new(fp.Element)
	)

	// Step 1: t0 = x^0x2
	t0.Square(x)

	// Step 2: t3 = x^0x3
	t3.Mul(x, t0)

	// Step 3: t1 = x^0x5
	t1.Mul(t0, t3)

	// Step 4: z = x^0x7
	z.Mul(t0, t1)

	// Step 5: t7 = x^0x9
	t7.Mul(t0, z)

	// Step 6: t12 = x^0xb
	t12.Mul(t0, t7)

	// Step 7: t5 = x^0xd
	t5.Mul(t0, t12)

	// Step 8: t2 = x^0xf
	t2.Mul(t0, t5)

	// Step 9: t6 = x^0x11
	t6.Mul(t0, t2)

	// Step 10: t10 = x^0x13
	t10.Mul(t0, t6)

	// Step 11: t9 = x^0x15
	t9.Mul(t0, t10)

	// Step 12: t4 = x^0x17
	t4.Mul(t0, t9)

	// Step 13: t8 = x^0x19
	t8.Mul(t0, t4)

	// Step 14: t11 = x^0x1b
	t11.Mul(t0, t8)

	// Step 15: t0 = x^0x1d
	t0.Mul(t0, t11)

	// Step 16: t13 = x^0x20
	t13.Mul(t3, t0)

	// Step 20: t13 = x^0x200
	for s := 0; s < 4; s++ {
		t13.Square(t13)
	}

	// Step 21: t13 = x^0x211
	t13.Mul(t6, t13)

	// Step 22: t13 = x^0x422
	t13.Square(t13)

	// Step 23: t13 = x^0x423
	t13.Mul(x, t13)

	// Step 32: t13 = x^0x84600
	for s := 0; s < 9; s++ {
		t13.Square(t13)
	}

	// Step 33: t12 = x^0x8460b
	t12.Mul(t12, t13)

	// Step 37: t12 = x^0x8460b0
	for s := 0; s < 4; s++ {
		t12.Square(t12)
	}

	// Step 38: t12 = x^0x8460b3
	t12.Mul(t3, t12)

	// Step 46: t12 = x^0x8460b300
	for s := 0; s < 8; s++ {
		t12.Square(t12)
	}

	// Step 47: t12 = x^0x8460b31b
	t12.Mul(t11, t12)

	// Step 56: t12 = x^0x108c1663600
	for s := 0; s < 9; s++ {
		t12.Square(t12)
	}

	// Step 57: t12 = x^0x108c1663603
	t12.Mul(t3, t12)

	// Step 64: t12 = x^0x8460b31b0180
	for s := 0; s < 7; s++ {
		t12.Square(t12)
	}

	// Step 65: t12 = x^0x8460b31b018f
	t12.Mul(t2, t12)

	// Step 73: t12 = x^0x8460b31b018f00
	for s := 0; s < 8; s++ {
		t12.Square(t12)
	}

	// Step 74: t12 = x^0x8460b31b018f0d
	t12.Mul(t5, t12)

	// Step 79: t12 = x^0x108c16636031e1a0
	for s := 0; s < 5; s++ {
		t12.Square(t12)
	}

	// Step 80: t12 = x^0x108c16636031e1a5
	t12.Mul(t1, t12)

	// Step 85: t12 = x^0x21182cc6c063c34a0
	for s := 0; s < 5; s++ {
		t12.Square(t12)
	}

	// Step 86: t12 = x^0x21182cc6c063c34a5
	t12.Mul(t1, t12)

	// Step 92: t12 = x^0x8460b31b018f0d2940
	for s := 0; s < 6; s++ {
		t12.Square(t12)
	}

	// Step 93: t12 = x^0x8460b31b018f0d294d
	t12.Mul(t5, t12)

	// Step 99: t12 = x^0x21182cc6c063c34a5340
	for s := 0; s < 6; s++ {
		t12.Square(t12)
	}

	// Step 100: t12 = x^0x21182cc6c063c34a534f
	t12.Mul(t2, t12)

	// Step 107: t12 = x^0x108c16636031e1a529a780
	for s := 0; s < 7; s++ {
		t12.Square(t12)
	}

	// Step 108: t12 = x^0x108c16636031e1a529a79b
	t12.Mul(t11, t12)

	// Step 116: t12 = x^0x108c16636031e1a529a79b00
	for s := 0; s < 8; s++ {
		t12.Square(t12)
	}

	// Step 117: t12 = x^0x108c16636031e1a529a79b17
	t12.Mul(t4, t12)

	// Step 122: t12 = x^0x21182cc6c063c34a534f362e0
	for s := 0; s < 5; s++ {
		t12.Square(t12)
	}

	// Step 123: t12 = x^0x21182cc6c063c34a534f362fb
	t12.Mul(t11, t12)

	// Step 128: t12 = x^0x4230598d80c78694a69e6c5f60
	for s := 0; s < 5; s++ {
		t12.Square(t12)
	}

	// Step 129: t12 = x^0x4230598d80c78694a69e6c5f7b
	t12.Mul(t11, t12)

	// Step 135: t12 = x^0x108c16636031e1a529a79b17dec0
	for s := 0; s < 6; s++ {
		t12.Square(t12)
	}

	// Step 136: t12 = x^0x108c16636031e1a529a79b17ded1
	t12.Mul(t6, t12)

	// Step 140: t12 = x^0x108c16636031e1a529a79b17ded10
	for s := 0; s < 4; s++ {
		t12.Square(t12)
	}

	// Step 141: t12 = x^0x108c16636031e1a529a79b17ded19
	t12.Mul(t7, t12)

	// Step 147: t12 = x^0x4230598d80c78694a69e6c5f7b4640
	for s := 0; s < 6; s++ {
		t12.Square(t12)
	}

	// Step 148: t12 = x^0x4230598d80c78694a69e6c5f7b4657
	t12.Mul(t4, t12)

	// Step 153: t12 = x^0x8460b31b018f0d294d3cd8bef68cae0
	for s := 0; s < 5; s++ {
		t12.Square(t12)
	}

	// Step 154: t11 = x^0x8460b31b018f0d294d3cd8bef68cafb
	t11.Mul(t11, t12)

	// Step 158: t11 = x^0x8460b31b018f0d294d3cd8bef68cafb0
	for s := 0; s < 4; s++ {
		t11.Square(t11)
	}

	// Step 159: t11 = x^0x8460b31b018f0d294d3cd8bef68cafb9
	t11.Mul(t7, t11)

	// Step 165: t11 = x^0x21182cc6c063c34a534f362fbda32bee40
	for s := 0; s < 6; s++ {
		t11.Square(t11)
	}

	// Step 166: t11 = x^0x21182cc6c063c34a534f362fbda32bee51
	t11.Mul(t6, t11)

	// Step 170: t11 = x^0x21182cc6c063c34a534f362fbda32bee510
	for s := 0; s < 4; s++ {
		t11.Square(t11)
	}

	// Step 171: t11 = x^0x21182cc6c063c34a534f362fbda32bee517
	t11.Mul(z, t11)

	// Step 178: t11 = x^0x108c16636031e1a529a79b17ded195f728b80
	for s := 0; s < 7; s++ {
		t11.Square(t11)
	}

	// Step 179: t11 = x^0x108c16636031e1a529a79b17ded195f728b99
	t11.Mul(t8, t11)

	// Step 183: t11 = x^0x108c16636031e1a529a79b17ded195f728b990
	for s := 0; s < 4; s++ {
		t11.Square(t11)
	}

	// Step 184: t11 = x^0x108c16636031e1a529a79b17ded195f728b99d
	t11.Mul(t5, t11)

	// Step 193: t11 = x^0x21182cc6c063c34a534f362fbda32bee51733a00
	for s := 0; s < 9; s++ {
		t11.Square(t11)
	}

	// Step 194: t10 = x^0x21182cc6c063c34a534f362fbda32bee51733a13
	t10.Mul(t10, t11)

	// Step 196: t10 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84c
	for s := 0; s < 2; s++ {
		t10.Square(t10)
	}

	// Step 197: t10 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f
	t10.Mul(t3, t10)

	// Step 201: t10 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f0
	for s := 0; s < 4; s++ {
		t10.Square(t10)
	}

	// Step 202: t10 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1
	t10.Mul(x, t10)

	// Step 212: t10 = x^0x21182cc6c063c34a534f362fbda32bee51733a13c400
	for s := 0; s < 10; s++ {
		t10.Square(t10)
	}

	// Step 213: t9 = x^0x21182cc6c063c34a534f362fbda32bee51733a13c415
	t9.Mul(t9, t10)

	// Step 220: t9 = x^0x108c16636031e1a529a79b17ded195f728b99d09e20a80
	for s := 0; s < 7; s++ {
		t9.Square(t9)
	}

	// Step 221: t9 = x^0x108c16636031e1a529a79b17ded195f728b99d09e20a99
	t9.Mul(t8, t9)

	// Step 228: t9 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1054c80
	for s := 0; s < 7; s++ {
		t9.Square(t9)
	}

	// Step 229: t8 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1054c99
	t8.Mul(t8, t9)

	// Step 235: t8 = x^0x21182cc6c063c34a534f362fbda32bee51733a13c41532640
	for s := 0; s < 6; s++ {
		t8.Square(t8)
	}

	// Step 236: t8 = x^0x21182cc6c063c34a534f362fbda32bee51733a13c4153264d
	t8.Mul(t5, t8)

	// Step 241: t8 = x^0x4230598d80c78694a69e6c5f7b4657dca2e67427882a64c9a0
	for s := 0; s < 5; s++ {
		t8.Square(t8)
	}

	// Step 242: t8 = x^0x4230598d80c78694a69e6c5f7b4657dca2e67427882a64c9af
	t8.Mul(t2, t8)

	// Step 247: t8 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1054c9935e0
	for s := 0; s < 5; s++ {
		t8.Square(t8)
	}

	// Step 248: t7 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1054c9935e9
	t7.Mul(t7, t8)

	// Step 254: t7 = x^0x21182cc6c063c34a534f362fbda32bee51733a13c4153264d7a40
	for s := 0; s < 6; s++ {
		t7.Square(t7)
	}

	// Step 255: t6 = x^0x21182cc6c063c34a534f362fbda32bee51733a13c4153264d7a51
	t6.Mul(t6, t7)

	// Step 261: t6 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1054c9935e9440
	for s := 0; s < 6; s++ {
		t6.Square(t6)
	}

	// Step 262: t6 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1054c9935e945d
	t6.Mul(t0, t6)

	// Step 268: t6 = x^0x21182cc6c063c34a534f362fbda32bee51733a13c4153264d7a51740
	for s := 0; s < 6; s++ {
		t6.Square(t6)
	}

	// Step 269: t6 = x^0x21182cc6c063c34a534f362fbda32bee51733a13c4153264d7a5175d
	t6.Mul(t0, t6)

	// Step 275: t6 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1054c9935e945d740
	for s := 0; s < 6; s++ {
		t6.Square(t6)
	}

	// Step 276: t6 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1054c9935e945d745
	t6.Mul(t1, t6)

	// Step 283: t6 = x^0x4230598d80c78694a69e6c5f7b4657dca2e67427882a64c9af4a2eba280
	for s := 0; s < 7; s++ {
		t6.Square(t6)
	}

	// Step 284: t5 = x^0x4230598d80c78694a69e6c5f7b4657dca2e67427882a64c9af4a2eba28d
	t5.Mul(t5, t6)

	// Step 305: t5 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1054c9935e945d7451a00000
	for s := 0; s < 21; s++ {
		t5.Square(t5)
	}

	// Step 306: t5 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1054c9935e945d7451a00007
	t5.Mul(z, t5)

	// Step 314: t5 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1054c9935e945d7451a0000700
	for s := 0; s < 8; s++ {
		t5.Square(t5)
	}

	// Step 315: t4 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1054c9935e945d7451a0000717
	t4.Mul(t4, t5)

	// Step 319: t4 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1054c9935e945d7451a00007170
	for s := 0; s < 4; s++ {
		t4.Square(t4)
	}

	// Step 320: t4 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1054c9935e945d7451a00007173
	t4.Mul(t3, t4)

	// Step 327: t4 = x^0x4230598d80c78694a69e6c5f7b4657dca2e67427882a64c9af4a2eba28d000038b980
	for s := 0; s < 7; s++ {
		t4.Square(t4)
	}

	// Step 328: t3 = x^0x4230598d80c78694a69e6c5f7b4657dca2e67427882a64c9af4a2eba28d000038b983
	t3.Mul(t3, t4)

	// Step 335: t3 = x^0x21182cc6c063c34a534f362fbda32bee51733a13c4153264d7a5175d14680001c5cc180
	for s := 0; s < 7; s++ {
		t3.Square(t3)
	}

	// Step 336: t2 = x^0x21182cc6c063c34a534f362fbda32bee51733a13c4153264d7a5175d14680001c5cc18f
	t2.Mul(t2, t3)

	// Step 371: t2 = x^0x108c16636031e1a529a79b17ded195f728b99d09e20a99326bd28bae8a340000e2e60c7800000000
	for s := 0; s < 35; s++ {
		t2.Square(t2)
	}

	// Step 372: t1 = x^0x108c16636031e1a529a79b17ded195f728b99d09e20a99326bd28bae8a340000e2e60c7800000005
	t1.Mul(t1, t2)

	// Step 380: t1 = x^0x108c16636031e1a529a79b17ded195f728b99d09e20a99326bd28bae8a340000e2e60c780000000500
	for s := 0; s < 8; s++ {
		t1.Square(t1)
	}

	// Step 381: t0 = x^0x108c16636031e1a529a79b17ded195f728b99d09e20a99326bd28bae8a340000e2e60c78000000051d
	t0.Mul(t0, t1)

	// Step 384: t0 = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1054c9935e945d7451a00007173063c000000028e8
	for s := 0; s < 3; s++ {
		t0.Square(t0)
	}

	// Step 385: z = x^0x8460b31b018f0d294d3cd8bef68cafb945cce84f1054c9935e945d7451a00007173063c000000028ef
	z.Mul(z, t0)

	// Step 431: z = x^0x21182cc6c063c34a534f362fbda32bee51733a13c4153264d7a5175d14680001c5cc18f00000000a3bc00000000000
	for s := 0; s < 46; s++ {
		z.Square(z)
	}

	return z
}
