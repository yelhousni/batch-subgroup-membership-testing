package bls12381

import (
	"math"
	"math/big"
	"unsafe"

	curve "github.com/consensys/gnark-crypto/ecc/bls12-381"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fp"
	"github.com/yelhousni/batch-subgroup-membership/go/eisenstein"
)

// isFirstTateOne checks that Tate_{3,P3}(Q) = (y-2)^((p-1)/3) == 1
// where P3 = (0,2) a point of order 3 on the
func isFirstTateOne(point curve.G1Affine) bool {
	var tate fp.Element
	tate.Sub(&point.Y, &two_p)
	tate = expByp3(tate)
	return tate.IsOne()
}

// isSecondTateOne checks that Tate_{11,P11}(Q) == Tate_{11,P'11}(Q) == 1
// where P11 = (x,y) and P'11 = (x',y') are points of order 11 on the curve
//
//	x  = 0x1147c19050b3c4b663a4ca29c4859eeb1ac05a91659009602e7443347ad659e9f838f4ed07337c4c6d3a48d612b4bb92
//	y  = 0x8d7c25237c7dcea6ea0c6c37053882c59cc0ee424b3545bb25116d53e383574063149edb438b959dd169d0e01b2d3bc
//	x' = 0xb9529a7b23788075a6c33c7b77b3dcf4da4f58af5310f32e739a6c653a5a8f7cf7f19a297bd6a8f3f19ea82cf9419
//	y' = 0x2ecc645926cbd45f215b3fa17df0d7a50e5814f9631c502f2b2c2457926089a452bd11bf89ee72baa1981f99f88acb2
func isSecondTateOne(point curve.G1Affine) bool {
	return tateP11(point, lines1).IsOne() && tateP11(point, lines2).IsOne()
}

func tateP11(point curve.G1Affine, lines [7]line) *fp.Element {

	// f_{11,P} = (l_{P,P}^4 * (l_{4P,P} * l_{2P,2P})^2 * l_{5P,5P}) /
	// 			  (v_{2P}^4 * (v_{5P} * v_{4P})^2)

	var num, denom, tate, f1, f2 fp.Element

	// l_{P,P}^4
	f1.Mul(&point.X, &lines[0].a).Add(&f1, &point.Y).Add(&f1, &lines[0].b)
	num.Square(&f1).Square(&num)
	// (l_{4P,P} * l_{2P,2P})^2
	f1.Mul(&point.X, &lines[1].a).Add(&f1, &point.Y).Add(&f1, &lines[1].b)
	f2.Mul(&point.X, &lines[2].a).Add(&f2, &point.Y).Add(&f2, &lines[2].b)
	f1.Mul(&f1, &f2).Square(&f1)
	num.Mul(&num, &f1)
	// l_{5P,5P}
	f1.Mul(&point.X, &lines[3].a).Add(&f1, &point.Y).Add(&f1, &lines[3].b)
	num.Mul(&num, &f1)

	// v_{2P}^4
	f1.Add(&point.X, &lines[4].b)
	denom.Square(&f1).Square(&denom)
	// (v_{5P} * v_{4P})^2
	f1.Add(&point.X, &lines[5].b)
	f2.Add(&point.X, &lines[6].b)
	f1.Mul(&f1, &f2).Square(&f1)
	denom.Mul(&denom, &f1)

	// denom^{-1} = denom^{10} inside the 11-th power residue symbol
	f1.Square(&denom)
	f2.Square(&f1).Square(&f2)
	denom.Mul(&f1, &f2)

	// tate = num * denom^{-1}
	tate.Mul(&num, &denom)

	// tate^((p-1)/11)
	tate = expByp11(tate)

	return &tate
}

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

var lines1, lines2 [7]line

var a, b big.Int
var beta eisenstein.ComplexNumber
var two_p fp.Element
var one, mone, three, mInt big.Int

func init() {
	// P = (
	// 	   0x1147c19050b3c4b663a4ca29c4859eeb1ac05a91659009602e7443347ad659e9f838f4ed07337c4c6d3a48d612b4bb92,
	// 	   0x8d7c25237c7dcea6ea0c6c37053882c59cc0ee424b3545bb25116d53e383574063149edb438b959dd169d0e01b2d3bc,
	// )
	// l_{P,P}
	lines1[0].a.SetString("3299267897751996220385221613595438297307729278788392950891735532675429610222287484309000239394909811727236076945185")
	lines1[0].b.SetString("18310294753454426021959467664309789431911031424237667094992897857206723725083858249069477932398292510179589984700")
	// l_{4P,P}
	lines1[1].a.SetString("1326060517996439830854349154180133688409720663931962928965182360010832188530135117396607306069977833911060775112646")
	lines1[1].b.SetString("1860711533062130732961381660170466745245304904508895668176066939612956152476857706132410954150204925316691548907270")
	// l_{2P,2P}
	lines1[2].a.SetString("2853179109878601234046782899570881274317753498202772373425894224814228501182619893673996877217795247222500060654119")
	lines1[2].b.SetString("3919148652435634127739626306550673271134953250084791759799938026680308413200580025545156518087656006387815930965981")
	// l_{5P,5P}
	lines1[3].a.SetString("76869151700841882652162232262323355256984607267248845269243875865905828905681810750474587342425543874725249474214")
	lines1[3].b.SetString("2712605722009629562960998908199761937467351166907817800849564303111095879681627219760162201521231845288626389023754")
	// v_{2P}
	lines1[4].b.SetString("149300628287719084363805442027015732343796298341825369398280665579062340373588762829663721363015489369012576683403")
	// v_{5P}
	lines1[5].b.SetString("3190598233284899258971997366336263489839614042767959121007866745814914033177431592603096075902487433824628934272481")
	// v_{4P}
	lines1[6].b.SetString("801510375849910381258432150154296010042040970731859091646614021822277212541925850988502284218759419518696898884221")

	// P = (
	//		0xb9529a7b23788075a6c33c7b77b3dcf4da4f58af5310f32e739a6c653a5a8f7cf7f19a297bd6a8f3f19ea82cf9419
	//		0x2ecc645926cbd45f215b3fa17df0d7a50e5814f9631c502f2b2c2457926089a452bd11bf89ee72baa1981f99f88acb2
	// )
	// l_{P,P}
	lines2[0].a.SetString("2235951733532706502589776802947624104221891830898748507534133878243037267233337743484969160310120281034890952723855")
	lines2[0].b.SetString("3984099260468212967395830358071594367124971788514770218237065238266824926765754006193618151196617371527714682575087")
	// l_{4P,P}
	lines2[1].a.SetString("3113526804415798250523539787369118304805432039590334316115598711790390730680464357060352593763557367287683707731426")
	lines2[1].b.SetString("2141698022159536660456408165565437411311577915430112217155991196511075498013980158310276674978810738721202723652517")
	// l_{2P,2P}
	lines2[2].a.SetString("2087932601249983316438079912808932212890646191989315900972772301900616622702132982227825023298901701929602218097069")
	lines2[2].b.SetString("83260902786033265678163519185230885421929569854216125532120109443723237290257838897531111041359657650078341593806")
	// l_{5P,5P}
	lines2[3].a.SetString("2231745861652326832984051517563825681412940931618559967380562501075717893398129749522302298367351947562599726222509")
	lines2[3].b.SetString("1289803833212037830456790917536142219089531653031190084482493833012935770809210644682525427607783818749267883536033")
	// v_{2P}
	lines2[4].b.SetString("3925304434066507134345633153253941624275432978844429686751215338284346393393611223652014992179037913992576863075846")
	// v_{5P}
	lines2[5].b.SetString("2168153076519464482853727572352404894820967148884886371318019694225522318905442037294956076553151994964407937624199")
	// v_{4P}
	lines2[6].b.SetString("2936782925110917657004772320034083346502209258972949212348266083811354308950004356358196405296552596411593186618131")

	// beta = a+ω*b with b primary and norm(beta)=p
	a.SetString("-1155048275357884106335086113613464118783412807316232579754", 10)
	b.SetString("1155048275357884106335086113613464118768280431093290937003", 10)
	beta.A0 = a
	beta.A1 = b

	// some constatns
	two_p.SetUint64(2)
	one.SetUint64(1)
	mone.Neg(&one)
	three.SetUint64(3)

}

// expByp11 uses a short addition chain to compute x^p11 where p11=(p-1)/11 .
func expByp11(x fp.Element) fp.Element {
	// Operations: 372 squares 77 multiplies
	//
	// Generated by github.com/mmcloughlin/addchain v0.4.0.

	// Allocate Temporaries.
	var z fp.Element
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
		t18 = new(fp.Element)
		t19 = new(fp.Element)
		t20 = new(fp.Element)
		t21 = new(fp.Element)
		t22 = new(fp.Element)
		t23 = new(fp.Element)
		t24 = new(fp.Element)
		t25 = new(fp.Element)
		t26 = new(fp.Element)
		t27 = new(fp.Element)
		t28 = new(fp.Element)
		t29 = new(fp.Element)
		t30 = new(fp.Element)
	)

	// Step 1: t0 = x^0x2
	t0.Square(&x)

	// Step 2: t9 = x^0x4
	t9.Square(t0)

	// Step 3: t1 = x^0x8
	t1.Square(t9)

	// Step 4: t16 = x^0xc
	t16.Mul(t9, t1)

	// Step 5: t15 = x^0xe
	t15.Mul(t0, t16)

	// Step 6: t5 = x^0x12
	t5.Mul(t9, t15)

	// Step 7: t29 = x^0x13
	t29.Mul(&x, t5)

	// Step 8: t2 = x^0x17
	t2.Mul(t9, t29)

	// Step 9: t11 = x^0x1a
	t11.Mul(t1, t5)

	// Step 10: t8 = x^0x1b
	t8.Mul(&x, t11)

	// Step 11: t6 = x^0x1d
	t6.Mul(t0, t8)

	// Step 12: t18 = x^0x1e
	t18.Mul(&x, t6)

	// Step 13: z = x^0x1f
	z.Mul(&x, t18)

	// Step 14: t24 = x^0x27
	t24.Mul(t1, &z)

	// Step 15: t12 = x^0x29
	t12.Mul(t0, t24)

	// Step 16: t28 = x^0x43
	t28.Mul(t11, t12)

	// Step 17: t21 = x^0x4f
	t21.Mul(t16, t28)

	// Step 18: t14 = x^0x51
	t14.Mul(t0, t21)

	// Step 19: t4 = x^0x55
	t4.Mul(t9, t14)

	// Step 20: t7 = x^0x59
	t7.Mul(t9, t4)

	// Step 21: t27 = x^0x5b
	t27.Mul(t0, t7)

	// Step 22: t3 = x^0x5d
	t3.Mul(t0, t27)

	// Step 23: t10 = x^0x65
	t10.Mul(t1, t3)

	// Step 24: t23 = x^0x69
	t23.Mul(t9, t10)

	// Step 25: t17 = x^0x6b
	t17.Mul(t0, t23)

	// Step 26: t5 = x^0x7d
	t5.Mul(t5, t17)

	// Step 27: t26 = x^0x81
	t26.Mul(t9, t5)

	// Step 28: t25 = x^0x83
	t25.Mul(t0, t26)

	// Step 29: t9 = x^0x89
	t9.Mul(t1, t26)

	// Step 30: t19 = x^0x91
	t19.Mul(t1, t9)

	// Step 31: t1 = x^0x97
	t1.Mul(t15, t9)

	// Step 32: t20 = x^0xa1
	t20.Mul(t18, t25)

	// Step 33: t22 = x^0xa7
	t22.Mul(t18, t9)

	// Step 34: t13 = x^0xb1
	t13.Mul(t11, t1)

	// Step 35: t11 = x^0xc5
	t11.Mul(t18, t22)

	// Step 36: t15 = x^0xd3
	t15.Mul(t15, t11)

	// Step 37: t18 = x^0xf1
	t18.Mul(t18, t15)

	// Step 38: t16 = x^0xfd
	t16.Mul(t16, t18)

	// Step 39: t0 = x^0xff
	t0.Mul(t0, t16)

	// Step 40: t30 = x^0x12e
	t30.Mul(t27, t15)

	// Step 45: t30 = x^0x25c0
	for s := 0; s < 5; s++ {
		t30.Square(t30)
	}

	// Step 46: t29 = x^0x25d3
	t29.Mul(t29, t30)

	// Step 59: t29 = x^0x4ba6000
	for s := 0; s < 13; s++ {
		t29.Square(t29)
	}

	// Step 60: t29 = x^0x4ba6059
	t29.Mul(t7, t29)

	// Step 69: t29 = x^0x974c0b200
	for s := 0; s < 9; s++ {
		t29.Square(t29)
	}

	// Step 70: t28 = x^0x974c0b243
	t28.Mul(t28, t29)

	// Step 76: t28 = x^0x25d302c90c0
	for s := 0; s < 6; s++ {
		t28.Square(t28)
	}

	// Step 77: t28 = x^0x25d302c90dd
	t28.Mul(t6, t28)

	// Step 88: t28 = x^0x12e9816486e800
	for s := 0; s < 11; s++ {
		t28.Square(t28)
	}

	// Step 89: t28 = x^0x12e9816486e8a7
	t28.Mul(t22, t28)

	// Step 96: t28 = x^0x974c0b243745380
	for s := 0; s < 7; s++ {
		t28.Square(t28)
	}

	// Step 97: t27 = x^0x974c0b2437453db
	t27.Mul(t27, t28)

	// Step 110: t27 = x^0x12e9816486e8a7b6000
	for s := 0; s < 13; s++ {
		t27.Square(t27)
	}

	// Step 111: t26 = x^0x12e9816486e8a7b6081
	t26.Mul(t26, t27)

	// Step 120: t26 = x^0x25d302c90dd14f6c10200
	for s := 0; s < 9; s++ {
		t26.Square(t26)
	}

	// Step 121: t25 = x^0x25d302c90dd14f6c10283
	t25.Mul(t25, t26)

	// Step 127: t25 = x^0x974c0b2437453db040a0c0
	for s := 0; s < 6; s++ {
		t25.Square(t25)
	}

	// Step 128: t24 = x^0x974c0b2437453db040a0e7
	t24.Mul(t24, t25)

	// Step 139: t24 = x^0x4ba605921ba29ed8205073800
	for s := 0; s < 11; s++ {
		t24.Square(t24)
	}

	// Step 140: t23 = x^0x4ba605921ba29ed8205073869
	t23.Mul(t23, t24)

	// Step 149: t23 = x^0x974c0b2437453db040a0e70d200
	for s := 0; s < 9; s++ {
		t23.Square(t23)
	}

	// Step 150: t22 = x^0x974c0b2437453db040a0e70d2a7
	t22.Mul(t22, t23)

	// Step 159: t22 = x^0x12e9816486e8a7b608141ce1a54e00
	for s := 0; s < 9; s++ {
		t22.Square(t22)
	}

	// Step 160: t21 = x^0x12e9816486e8a7b608141ce1a54e4f
	t21.Mul(t21, t22)

	// Step 170: t21 = x^0x4ba605921ba29ed82050738695393c00
	for s := 0; s < 10; s++ {
		t21.Square(t21)
	}

	// Step 171: t20 = x^0x4ba605921ba29ed82050738695393ca1
	t20.Mul(t20, t21)

	// Step 181: t20 = x^0x12e9816486e8a7b608141ce1a54e4f28400
	for s := 0; s < 10; s++ {
		t20.Square(t20)
	}

	// Step 182: t19 = x^0x12e9816486e8a7b608141ce1a54e4f28491
	t19.Mul(t19, t20)

	// Step 194: t19 = x^0x12e9816486e8a7b608141ce1a54e4f28491000
	for s := 0; s < 12; s++ {
		t19.Square(t19)
	}

	// Step 195: t18 = x^0x12e9816486e8a7b608141ce1a54e4f284910f1
	t18.Mul(t18, t19)

	// Step 205: t18 = x^0x4ba605921ba29ed82050738695393ca12443c400
	for s := 0; s < 10; s++ {
		t18.Square(t18)
	}

	// Step 206: t17 = x^0x4ba605921ba29ed82050738695393ca12443c46b
	t17.Mul(t17, t18)

	// Step 215: t17 = x^0x974c0b2437453db040a0e70d2a727942488788d600
	for s := 0; s < 9; s++ {
		t17.Square(t17)
	}

	// Step 216: t16 = x^0x974c0b2437453db040a0e70d2a727942488788d6fd
	t16.Mul(t16, t17)

	// Step 226: t16 = x^0x25d302c90dd14f6c102839c34a9c9e509221e235bf400
	for s := 0; s < 10; s++ {
		t16.Square(t16)
	}

	// Step 227: t15 = x^0x25d302c90dd14f6c102839c34a9c9e509221e235bf4d3
	t15.Mul(t15, t16)

	// Step 236: t15 = x^0x4ba605921ba29ed82050738695393ca12443c46b7e9a600
	for s := 0; s < 9; s++ {
		t15.Square(t15)
	}

	// Step 237: t14 = x^0x4ba605921ba29ed82050738695393ca12443c46b7e9a651
	t14.Mul(t14, t15)

	// Step 246: t14 = x^0x974c0b2437453db040a0e70d2a727942488788d6fd34ca200
	for s := 0; s < 9; s++ {
		t14.Square(t14)
	}

	// Step 247: t13 = x^0x974c0b2437453db040a0e70d2a727942488788d6fd34ca2b1
	t13.Mul(t13, t14)

	// Step 255: t13 = x^0x974c0b2437453db040a0e70d2a727942488788d6fd34ca2b100
	for s := 0; s < 8; s++ {
		t13.Square(t13)
	}

	// Step 256: t12 = x^0x974c0b2437453db040a0e70d2a727942488788d6fd34ca2b129
	t12.Mul(t12, t13)

	// Step 266: t12 = x^0x25d302c90dd14f6c102839c34a9c9e509221e235bf4d328ac4a400
	for s := 0; s < 10; s++ {
		t12.Square(t12)
	}

	// Step 267: t12 = x^0x25d302c90dd14f6c102839c34a9c9e509221e235bf4d328ac4a41b
	t12.Mul(t8, t12)

	// Step 278: t12 = x^0x12e9816486e8a7b608141ce1a54e4f284910f11adfa6994562520d800
	for s := 0; s < 11; s++ {
		t12.Square(t12)
	}

	// Step 279: t11 = x^0x12e9816486e8a7b608141ce1a54e4f284910f11adfa6994562520d8c5
	t11.Mul(t11, t12)

	// Step 287: t11 = x^0x12e9816486e8a7b608141ce1a54e4f284910f11adfa6994562520d8c500
	for s := 0; s < 8; s++ {
		t11.Square(t11)
	}

	// Step 288: t10 = x^0x12e9816486e8a7b608141ce1a54e4f284910f11adfa6994562520d8c565
	t10.Mul(t10, t11)

	// Step 298: t10 = x^0x4ba605921ba29ed82050738695393ca12443c46b7e9a65158948363159400
	for s := 0; s < 10; s++ {
		t10.Square(t10)
	}

	// Step 299: t9 = x^0x4ba605921ba29ed82050738695393ca12443c46b7e9a65158948363159489
	t9.Mul(t9, t10)

	// Step 304: t9 = x^0x974c0b2437453db040a0e70d2a727942488788d6fd34ca2b12906c62b29120
	for s := 0; s < 5; s++ {
		t9.Square(t9)
	}

	// Step 305: t8 = x^0x974c0b2437453db040a0e70d2a727942488788d6fd34ca2b12906c62b2913b
	t8.Mul(t8, t9)

	// Step 320: t8 = x^0x4ba605921ba29ed82050738695393ca12443c46b7e9a65158948363159489d8000
	for s := 0; s < 15; s++ {
		t8.Square(t8)
	}

	// Step 321: t7 = x^0x4ba605921ba29ed82050738695393ca12443c46b7e9a65158948363159489d8059
	t7.Mul(t7, t8)

	// Step 328: t7 = x^0x25d302c90dd14f6c102839c34a9c9e509221e235bf4d328ac4a41b18aca44ec02c80
	for s := 0; s < 7; s++ {
		t7.Square(t7)
	}

	// Step 329: t6 = x^0x25d302c90dd14f6c102839c34a9c9e509221e235bf4d328ac4a41b18aca44ec02c9d
	t6.Mul(t6, t7)

	// Step 339: t6 = x^0x974c0b2437453db040a0e70d2a727942488788d6fd34ca2b12906c62b2913b00b27400
	for s := 0; s < 10; s++ {
		t6.Square(t6)
	}

	// Step 340: t6 = x^0x974c0b2437453db040a0e70d2a727942488788d6fd34ca2b12906c62b2913b00b2745d
	t6.Mul(t3, t6)

	// Step 351: t6 = x^0x4ba605921ba29ed82050738695393ca12443c46b7e9a65158948363159489d80593a2e800
	for s := 0; s < 11; s++ {
		t6.Square(t6)
	}

	// Step 352: t5 = x^0x4ba605921ba29ed82050738695393ca12443c46b7e9a65158948363159489d80593a2e87d
	t5.Mul(t5, t6)

	// Step 360: t5 = x^0x4ba605921ba29ed82050738695393ca12443c46b7e9a65158948363159489d80593a2e87d00
	for s := 0; s < 8; s++ {
		t5.Square(t5)
	}

	// Step 361: t4 = x^0x4ba605921ba29ed82050738695393ca12443c46b7e9a65158948363159489d80593a2e87d55
	t4.Mul(t4, t5)

	// Step 371: t4 = x^0x12e9816486e8a7b608141ce1a54e4f284910f11adfa6994562520d8c56522760164e8ba1f55400
	for s := 0; s < 10; s++ {
		t4.Square(t4)
	}

	// Step 372: t3 = x^0x12e9816486e8a7b608141ce1a54e4f284910f11adfa6994562520d8c56522760164e8ba1f5545d
	t3.Mul(t3, t4)

	// Step 380: t3 = x^0x12e9816486e8a7b608141ce1a54e4f284910f11adfa6994562520d8c56522760164e8ba1f5545d00
	for s := 0; s < 8; s++ {
		t3.Square(t3)
	}

	// Step 381: t2 = x^0x12e9816486e8a7b608141ce1a54e4f284910f11adfa6994562520d8c56522760164e8ba1f5545d17
	t2.Mul(t2, t3)

	// Step 392: t2 = x^0x974c0b2437453db040a0e70d2a727942488788d6fd34ca2b12906c62b2913b00b2745d0faaa2e8b800
	for s := 0; s < 11; s++ {
		t2.Square(t2)
	}

	// Step 393: t1 = x^0x974c0b2437453db040a0e70d2a727942488788d6fd34ca2b12906c62b2913b00b2745d0faaa2e8b897
	t1.Mul(t1, t2)

	// Step 403: t1 = x^0x25d302c90dd14f6c102839c34a9c9e509221e235bf4d328ac4a41b18aca44ec02c9d1743eaa8ba2e25c00
	for s := 0; s < 10; s++ {
		t1.Square(t1)
	}

	// Step 404: t1 = x^0x25d302c90dd14f6c102839c34a9c9e509221e235bf4d328ac4a41b18aca44ec02c9d1743eaa8ba2e25cff
	t1.Mul(t0, t1)

	// Step 412: t1 = x^0x25d302c90dd14f6c102839c34a9c9e509221e235bf4d328ac4a41b18aca44ec02c9d1743eaa8ba2e25cff00
	for s := 0; s < 8; s++ {
		t1.Square(t1)
	}

	// Step 413: t1 = x^0x25d302c90dd14f6c102839c34a9c9e509221e235bf4d328ac4a41b18aca44ec02c9d1743eaa8ba2e25cffff
	t1.Mul(t0, t1)

	// Step 421: t1 = x^0x25d302c90dd14f6c102839c34a9c9e509221e235bf4d328ac4a41b18aca44ec02c9d1743eaa8ba2e25cffff00
	for s := 0; s < 8; s++ {
		t1.Square(t1)
	}

	// Step 422: t1 = x^0x25d302c90dd14f6c102839c34a9c9e509221e235bf4d328ac4a41b18aca44ec02c9d1743eaa8ba2e25cffffff
	t1.Mul(t0, t1)

	// Step 430: t1 = x^0x25d302c90dd14f6c102839c34a9c9e509221e235bf4d328ac4a41b18aca44ec02c9d1743eaa8ba2e25cffffff00
	for s := 0; s < 8; s++ {
		t1.Square(t1)
	}

	// Step 431: t0 = x^0x25d302c90dd14f6c102839c34a9c9e509221e235bf4d328ac4a41b18aca44ec02c9d1743eaa8ba2e25cffffffff
	t0.Mul(t0, t1)

	// Step 436: t0 = x^0x4ba605921ba29ed82050738695393ca12443c46b7e9a65158948363159489d80593a2e87d551745c4b9fffffffe0
	for s := 0; s < 5; s++ {
		t0.Square(t0)
	}

	// Step 437: t0 = x^0x4ba605921ba29ed82050738695393ca12443c46b7e9a65158948363159489d80593a2e87d551745c4b9fffffffff
	t0.Mul(&z, t0)

	// Step 447: t0 = x^0x12e9816486e8a7b608141ce1a54e4f284910f11adfa6994562520d8c56522760164e8ba1f5545d1712e7ffffffffc00
	for s := 0; s < 10; s++ {
		t0.Square(t0)
	}

	// Step 448: z = x^0x12e9816486e8a7b608141ce1a54e4f284910f11adfa6994562520d8c56522760164e8ba1f5545d1712e7ffffffffc1f
	z.Mul(&z, t0)

	// Step 449: z = x^0x25d302c90dd14f6c102839c34a9c9e509221e235bf4d328ac4a41b18aca44ec02c9d1743eaa8ba2e25cfffffffff83e
	z.Square(&z)

	return z
}

func IsCubicResidue(x *fp.Element) bool {
	return CubicSymbol(*x).A0.Cmp(&one) == 0
}

func CubicSymbol(x fp.Element) *eisenstein.ComplexNumber {
	// α = x + ω * 0
	var alpha eisenstein.ComplexNumber
	x.BigInt(&alpha.A0)
	alpha.A1.SetUint64(0)
	// norm(β) = p, with β[1] = 0 mod 3 and β[0] ≠ 0
	var _beta eisenstein.ComplexNumber
	_beta.Set(&beta)
	return cubicSymbol(&alpha, &_beta)
}

func cubicSymbol(alpha, beta *eisenstein.ComplexNumber) *eisenstein.ComplexNumber {
	var result, gamma eisenstein.ComplexNumber
	var q0 big.Int
	result.SetOne()

	for {
		// Base cases: if α = ±1 or β = ±1 return result
		if (alpha.A1.Sign() == 0 && (alpha.A0.Cmp(&one) == 0 || alpha.A0.Cmp(&mone) == 0)) ||
			(beta.A1.Sign() == 0 && (beta.A0.Cmp(&one) == 0 || beta.A0.Cmp(&mone) == 0)) {
			return &result
		}

		// q = ⌊α/β⌉
		// γ = α - q * β
		gamma.Quo(alpha, beta)
		gamma.Mul(&gamma, beta)
		gamma.Sub(alpha, &gamma)

		// If γ = 0 return 0
		if gamma.A0.Sign() == 0 && gamma.A1.Sign() == 0 {
			return &eisenstein.ComplexNumber{}
		}

		// Remove ramified factors:
		// Compute m≥0 s.t. γ' = γ / (1-ω)^m is not divisible by 1-ω
		// i.e. γ'[0]+γ'[1] ≠ 0 mod 3.
		//
		// We use:
		// 		(γ / (1-ω))[0] = (2γ[0] - γ[1]) / 3
		// 		(γ / (1-ω))[1] = (γ[0] + γ[1]) / 3
		m := uint64(0)
		alpha.A1.Add(&gamma.A0, &gamma.A1)
		r1 := mod3(&alpha.A1)

		// the result is integral iff γ[0] + γ[1] = 0 mod 3
		for r1 == 0 {
			alpha.A0.Add(&gamma.A0, &gamma.A0).
				Sub(&alpha.A0, &gamma.A1)
			gamma.A0.Quo(&alpha.A0, &three)
			gamma.A1.Quo(&alpha.A1, &three)

			m++

			alpha.A1.Add(&gamma.A0, &gamma.A1)
			r1 = mod3(&alpha.A1)
		}

		// Make primary:
		// Find n with 0 ≤ n < 3 so that γ'/ω^n is primary.
		//
		// Division by ω is a multiplication by ω^2 so the result is one of:
		// 		γ[0] + ω * γ[1]
		// 		γ[1] - γ[0] - ω * γ[0]
		// 		-γ[1] - ω * (γ[0] - γ[1])
		n := 0
		alpha.A1.Sub(&gamma.A0, &gamma.A1)
		r0 := mod3(&gamma.A0)
		r1 = mod3(&alpha.A1)
		if r0 == 0 {
			n = 1
			gamma.A1.Neg(&gamma.A0)
			gamma.A0.Neg(&alpha.A1)
		} else if r1 == 0 {
			n = 2
			gamma.A0.Neg(&gamma.A1)
			gamma.A1.Set(&alpha.A1)
		}

		// Compute ω^exp, where
		// 		exp = ( n * (β[0]^2 − β[0]*β[1] − 1) + m * (1 - β[0]^2) ) / 3

		// m can be arbitrary but we specialize cases for m = 0, 1, 2.
		alpha.A1.Mul(&beta.A0, &beta.A0)
		switch m {
		case 0:
			q0.SetUint64(0)
		case 1:
			q0.Sub(&one, &alpha.A1)
		case 2:
			q0.Sub(&one, &alpha.A1).Add(&q0, &q0)
		default:
			mInt.SetUint64(m)
			q0.Sub(&one, &alpha.A1).Mul(&q0, &mInt)
		}

		// n can only be {0,1,2}
		switch n {
		case 0:
			// nada
		case 1:
			alpha.A0.Mul(&beta.A0, &beta.A1)
			increment(&alpha.A0)
			alpha.A0.Sub(&alpha.A1, &alpha.A0)
			q0.Add(&alpha.A0, &q0)
		case 2:
			alpha.A0.Mul(&beta.A0, &beta.A1)
			increment(&alpha.A0)
			alpha.A0.Sub(&alpha.A1, &alpha.A0).Add(&alpha.A0, &alpha.A0)
			q0.Add(&alpha.A0, &q0)
		default:
			panic("unexpected value for n in cubic symbol computation")
		}

		// Multiply the result by one of 1, ω, ω^2.
		//
		// We avoid dividing exp by 3 and switch on the value of exp mod 9.
		r0 = mod9(&q0)
		switch r0 {
		case 0:
			// muliply the result by 1
		case 3:
			// muliply the result by ω
			alpha.A0.Neg(&result.A1)
			result.A1.Add(&result.A0, &alpha.A0)
			result.A0.Set(&alpha.A0)
		case 6:
			// muliply the result by ω^2
			alpha.A1.Neg(&result.A0)
			result.A0.Add(&result.A1, &alpha.A1)
			result.A1.Set(&alpha.A1)
		default:
			panic("unexpected value in cubic symbol computation")
		}

		// Swap for the next iteration
		alpha.Set(beta)
		beta.Set(&gamma)
	}
}

func mod3(z *big.Int) uint64 {
	limbs := z.Bits()
	var sumOfLimbsMod3 uint64
	for _, limb := range limbs {
		sumOfLimbsMod3 = (sumOfLimbsMod3 + (uint64(limb) % 3)) % 3
	}
	return sumOfLimbsMod3
}

func mod9(z *big.Int) uint64 {
	limbs := z.Bits()
	var sumOfLimbsMod9 uint64

	// Determine the effective power of 2 for each limb based on Word size
	// big.Word is an unsigned integer type for a machine word.
	// We need 2^WordSize % 9
	var wordPow2Mod9 uint64

	switch unsafe.Sizeof(big.Word(0)) {
	case 4: // uint32
		// 2^32 % 9
		// 32 = 6*5 + 2
		// 2^32 = (2^6)^5 * 2^2
		// (2^6 % 9)^5 * 2^2 % 9
		// (1)^5 * 4 % 9 = 4 % 9 = 4
		wordPow2Mod9 = 4
	case 8: // uint64
		// 2^64 % 9
		// 64 = 6*10 + 4
		// 2^64 = (2^6)^10 * 2^4
		// (2^6 % 9)^10 * 2^4 % 9
		// (1)^10 * 16 % 9 = 16 % 9 = 7
		wordPow2Mod9 = 7
	default:
		// Fallback for unexpected big.Word size, though highly unlikely
		// This path would indicate a very unusual or non-standard Go compilation.
		// In a real-world scenario, you might want to panic or use z.Mod directly.
		return z.Mod(z, big.NewInt(9)).Uint64()
	}

	currentMultiplier := uint64(1) // Represents (2^W)^i % 9 for the current limb

	for _, limb := range limbs {
		// Each limb represents limb_val * (2^W)^i
		// The overall sum is sum (limb_val_i * (2^W)^i)
		// We need to sum (limb_val_i % 9 * ((2^W)^i % 9)) % 9

		limbValMod9 := uint64(limb) % 9
		term := (limbValMod9 * currentMultiplier) % 9
		sumOfLimbsMod9 = (sumOfLimbsMod9 + term) % 9

		// Update multiplier for the next limb: next_multiplier = current_multiplier * (2^W % 9)
		currentMultiplier = (currentMultiplier * wordPow2Mod9) % 9
	}

	// Handle negative numbers: math/big's Rem returns a result with the same sign as z.
	// For modulo, we typically want a non-negative result.
	if z.Sign() < 0 && sumOfLimbsMod9 != 0 {
		return 9 - sumOfLimbsMod9
	}
	return sumOfLimbsMod9
}

func increment(z *big.Int) {
	if z.Sign() > 0 {
		zBits := z.Bits()
		if zBits[0] < math.MaxUint64 {
			zBits[0] = big.Word(uint64(zBits[0]) + 1)
			return
		}
	}
	z.Add(z, &one)
}

// expByp3 uses a short addition chain to compute x^p3 where p3=(p-1)/3 .
func expByp3(x fp.Element) fp.Element {
	// Operations: 375 squares 79 multiplies
	//
	// Generated by github.com/mmcloughlin/addchain v0.4.0.

	// Allocate Temporaries.
	var z fp.Element
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
	)

	// Step 1: t0 = x^0x2
	t0.Square(&x)

	// Step 2: t6 = x^0x3
	t6.Mul(&x, t0)

	// Step 3: t3 = x^0x5
	t3.Mul(t0, t6)

	// Step 4: z = x^0x7
	z.Mul(t0, t3)

	// Step 5: t13 = x^0x9
	t13.Mul(t0, &z)

	// Step 6: t10 = x^0xb
	t10.Mul(t0, t13)

	// Step 7: t4 = x^0xd
	t4.Mul(t0, t10)

	// Step 8: t7 = x^0xf
	t7.Mul(t0, t4)

	// Step 9: t8 = x^0x11
	t8.Mul(t0, t7)

	// Step 10: t2 = x^0x13
	t2.Mul(t0, t8)

	// Step 11: t1 = x^0x15
	t1.Mul(t0, t2)

	// Step 12: t5 = x^0x17
	t5.Mul(t0, t1)

	// Step 13: t14 = x^0x19
	t14.Mul(t0, t5)

	// Step 14: t15 = x^0x1b
	t15.Mul(t0, t14)

	// Step 15: t11 = x^0x1d
	t11.Mul(t0, t15)

	// Step 16: t9 = x^0x1f
	t9.Mul(t0, t11)

	// Step 17: t0 = x^0x3e
	t0.Square(t9)

	// Step 18: t12 = x^0x3f
	t12.Mul(&x, t0)

	// Step 20: t0 = x^0xfc
	t0.Square(t12)
	for s := 1; s < 2; s++ {
		t0.Square(t0)
	}

	// Step 21: t0 = x^0xff
	t0.Mul(t6, t0)

	// Step 22: t16 = x^0x110
	t16.Mul(t8, t0)

	// Step 24: t16 = x^0x440
	for s := 0; s < 2; s++ {
		t16.Square(t16)
	}

	// Step 25: t16 = x^0x455
	t16.Mul(t1, t16)

	// Step 26: t16 = x^0x8aa
	t16.Square(t16)

	// Step 27: t16 = x^0x8ab
	t16.Mul(&x, t16)

	// Step 33: t16 = x^0x22ac0
	for s := 0; s < 6; s++ {
		t16.Square(t16)
	}

	// Step 34: t16 = x^0x22ac1
	t16.Mul(&x, t16)

	// Step 41: t16 = x^0x1156080
	for s := 0; s < 7; s++ {
		t16.Square(t16)
	}

	// Step 42: t16 = x^0x11560bf
	t16.Mul(t12, t16)

	// Step 50: t16 = x^0x11560bf00
	for s := 0; s < 8; s++ {
		t16.Square(t16)
	}

	// Step 51: t16 = x^0x11560bf17
	t16.Mul(t5, t16)

	// Step 56: t16 = x^0x22ac17e2e0
	for s := 0; s < 5; s++ {
		t16.Square(t16)
	}

	// Step 57: t16 = x^0x22ac17e2f7
	t16.Mul(t5, t16)

	// Step 63: t16 = x^0x8ab05f8bdc0
	for s := 0; s < 6; s++ {
		t16.Square(t16)
	}

	// Step 64: t16 = x^0x8ab05f8bdd5
	t16.Mul(t1, t16)

	// Step 70: t16 = x^0x22ac17e2f7540
	for s := 0; s < 6; s++ {
		t16.Square(t16)
	}

	// Step 71: t16 = x^0x22ac17e2f7553
	t16.Mul(t2, t16)

	// Step 78: t16 = x^0x11560bf17baa980
	for s := 0; s < 7; s++ {
		t16.Square(t16)
	}

	// Step 79: t15 = x^0x11560bf17baa99b
	t15.Mul(t15, t16)

	// Step 81: t15 = x^0x45582fc5eeaa66c
	for s := 0; s < 2; s++ {
		t15.Square(t15)
	}

	// Step 82: t15 = x^0x45582fc5eeaa66f
	t15.Mul(t6, t15)

	// Step 91: t15 = x^0x8ab05f8bdd54cde00
	for s := 0; s < 9; s++ {
		t15.Square(t15)
	}

	// Step 92: t14 = x^0x8ab05f8bdd54cde19
	t14.Mul(t14, t15)

	// Step 100: t14 = x^0x8ab05f8bdd54cde1900
	for s := 0; s < 8; s++ {
		t14.Square(t14)
	}

	// Step 101: t13 = x^0x8ab05f8bdd54cde1909
	t13.Mul(t13, t14)

	// Step 105: t13 = x^0x8ab05f8bdd54cde19090
	for s := 0; s < 4; s++ {
		t13.Square(t13)
	}

	// Step 106: t13 = x^0x8ab05f8bdd54cde19093
	t13.Mul(t6, t13)

	// Step 113: t13 = x^0x45582fc5eeaa66f0c84980
	for s := 0; s < 7; s++ {
		t13.Square(t13)
	}

	// Step 114: t12 = x^0x45582fc5eeaa66f0c849bf
	t12.Mul(t12, t13)

	// Step 121: t12 = x^0x22ac17e2f75533786424df80
	for s := 0; s < 7; s++ {
		t12.Square(t12)
	}

	// Step 122: t11 = x^0x22ac17e2f75533786424df9d
	t11.Mul(t11, t12)

	// Step 127: t11 = x^0x45582fc5eeaa66f0c849bf3a0
	for s := 0; s < 5; s++ {
		t11.Square(t11)
	}

	// Step 128: t11 = x^0x45582fc5eeaa66f0c849bf3b5
	t11.Mul(t1, t11)

	// Step 131: t11 = x^0x22ac17e2f75533786424df9da8
	for s := 0; s < 3; s++ {
		t11.Square(t11)
	}

	// Step 132: t11 = x^0x22ac17e2f75533786424df9daf
	t11.Mul(&z, t11)

	// Step 141: t11 = x^0x45582fc5eeaa66f0c849bf3b5e00
	for s := 0; s < 9; s++ {
		t11.Square(t11)
	}

	// Step 142: t11 = x^0x45582fc5eeaa66f0c849bf3b5e1f
	t11.Mul(t9, t11)

	// Step 149: t11 = x^0x22ac17e2f75533786424df9daf0f80
	for s := 0; s < 7; s++ {
		t11.Square(t11)
	}

	// Step 150: t11 = x^0x22ac17e2f75533786424df9daf0f91
	t11.Mul(t8, t11)

	// Step 158: t11 = x^0x22ac17e2f75533786424df9daf0f9100
	for s := 0; s < 8; s++ {
		t11.Square(t11)
	}

	// Step 159: t11 = x^0x22ac17e2f75533786424df9daf0f911f
	t11.Mul(t9, t11)

	// Step 163: t11 = x^0x22ac17e2f75533786424df9daf0f911f0
	for s := 0; s < 4; s++ {
		t11.Square(t11)
	}

	// Step 164: t11 = x^0x22ac17e2f75533786424df9daf0f911f3
	t11.Mul(t6, t11)

	// Step 173: t11 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e600
	for s := 0; s < 9; s++ {
		t11.Square(t11)
	}

	// Step 174: t11 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613
	t11.Mul(t2, t11)

	// Step 177: t11 = x^0x22ac17e2f75533786424df9daf0f911f3098
	for s := 0; s < 3; s++ {
		t11.Square(t11)
	}

	// Step 178: t11 = x^0x22ac17e2f75533786424df9daf0f911f309f
	t11.Mul(&z, t11)

	// Step 186: t11 = x^0x22ac17e2f75533786424df9daf0f911f309f00
	for s := 0; s < 8; s++ {
		t11.Square(t11)
	}

	// Step 187: t11 = x^0x22ac17e2f75533786424df9daf0f911f309f0f
	t11.Mul(t7, t11)

	// Step 192: t11 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1e0
	for s := 0; s < 5; s++ {
		t11.Square(t11)
	}

	// Step 193: t10 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb
	t10.Mul(t10, t11)

	// Step 199: t10 = x^0x11560bf17baa99bc32126fced787c88f984f87ac0
	for s := 0; s < 6; s++ {
		t10.Square(t10)
	}

	// Step 200: t9 = x^0x11560bf17baa99bc32126fced787c88f984f87adf
	t9.Mul(t9, t10)

	// Step 205: t9 = x^0x22ac17e2f75533786424df9daf0f911f309f0f5be0
	for s := 0; s < 5; s++ {
		t9.Square(t9)
	}

	// Step 206: t9 = x^0x22ac17e2f75533786424df9daf0f911f309f0f5bef
	t9.Mul(t7, t9)

	// Step 212: t9 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbc0
	for s := 0; s < 6; s++ {
		t9.Square(t9)
	}

	// Step 213: t9 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7
	t9.Mul(t5, t9)

	// Step 220: t9 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb80
	for s := 0; s < 7; s++ {
		t9.Square(t9)
	}

	// Step 221: t9 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb83
	t9.Mul(t6, t9)

	// Step 232: t9 = x^0x22ac17e2f75533786424df9daf0f911f309f0f5bef5c1800
	for s := 0; s < 11; s++ {
		t9.Square(t9)
	}

	// Step 233: t9 = x^0x22ac17e2f75533786424df9daf0f911f309f0f5bef5c18ff
	t9.Mul(t0, t9)

	// Step 239: t9 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fc0
	for s := 0; s < 6; s++ {
		t9.Square(t9)
	}

	// Step 240: t9 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd
	t9.Mul(t4, t9)

	// Step 244: t9 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd0
	for s := 0; s < 4; s++ {
		t9.Square(t9)
	}

	// Step 245: t9 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd1
	t9.Mul(&x, t9)

	// Step 255: t9 = x^0x22ac17e2f75533786424df9daf0f911f309f0f5bef5c18ff34400
	for s := 0; s < 10; s++ {
		t9.Square(t9)
	}

	// Step 256: t9 = x^0x22ac17e2f75533786424df9daf0f911f309f0f5bef5c18ff34411
	t9.Mul(t8, t9)

	// Step 261: t9 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb831fe688220
	for s := 0; s < 5; s++ {
		t9.Square(t9)
	}

	// Step 262: t8 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb831fe688231
	t8.Mul(t8, t9)

	// Step 267: t8 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104620
	for s := 0; s < 5; s++ {
		t8.Square(t8)
	}

	// Step 268: t8 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635
	t8.Mul(t1, t8)

	// Step 271: t8 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb831fe688231a8
	for s := 0; s < 3; s++ {
		t8.Square(t8)
	}

	// Step 272: t8 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb831fe688231ad
	t8.Mul(t3, t8)

	// Step 278: t8 = x^0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b40
	for s := 0; s < 6; s++ {
		t8.Square(t8)
	}

	// Step 279: t7 = x^0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f
	t7.Mul(t7, t8)

	// Step 282: t7 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a78
	for s := 0; s < 3; s++ {
		t7.Square(t7)
	}

	// Step 283: t7 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a79
	t7.Mul(&x, t7)

	// Step 291: t7 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a7900
	for s := 0; s < 8; s++ {
		t7.Square(t7)
	}

	// Step 292: t7 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a7905
	t7.Mul(t3, t7)

	// Step 295: t7 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb831fe688231ad3c828
	for s := 0; s < 3; s++ {
		t7.Square(t7)
	}

	// Step 296: t7 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb831fe688231ad3c829
	t7.Mul(&x, t7)

	// Step 303: t7 = x^0x22ac17e2f75533786424df9daf0f911f309f0f5bef5c18ff344118d69e41480
	for s := 0; s < 7; s++ {
		t7.Square(t7)
	}

	// Step 304: t6 = x^0x22ac17e2f75533786424df9daf0f911f309f0f5bef5c18ff344118d69e41483
	t6.Mul(t6, t7)

	// Step 313: t6 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb831fe688231ad3c8290600
	for s := 0; s < 9; s++ {
		t6.Square(t6)
	}

	// Step 314: t6 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb831fe688231ad3c8290605
	t6.Mul(t3, t6)

	// Step 320: t6 = x^0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a418140
	for s := 0; s < 6; s++ {
		t6.Square(t6)
	}

	// Step 321: t6 = x^0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a418147
	t6.Mul(&z, t6)

	// Step 328: t6 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a380
	for s := 0; s < 7; s++ {
		t6.Square(t6)
	}

	// Step 329: t6 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395
	t6.Mul(t1, t6)

	// Step 335: t6 = x^0x22ac17e2f75533786424df9daf0f911f309f0f5bef5c18ff344118d69e41483028e540
	for s := 0; s < 6; s++ {
		t6.Square(t6)
	}

	// Step 336: t6 = x^0x22ac17e2f75533786424df9daf0f911f309f0f5bef5c18ff344118d69e41483028e555
	t6.Mul(t1, t6)

	// Step 340: t6 = x^0x22ac17e2f75533786424df9daf0f911f309f0f5bef5c18ff344118d69e41483028e5550
	for s := 0; s < 4; s++ {
		t6.Square(t6)
	}

	// Step 341: t6 = x^0x22ac17e2f75533786424df9daf0f911f309f0f5bef5c18ff344118d69e41483028e5555
	t6.Mul(t3, t6)

	// Step 346: t6 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb831fe688231ad3c82906051caaaa0
	for s := 0; s < 5; s++ {
		t6.Square(t6)
	}

	// Step 347: t6 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb831fe688231ad3c82906051caaaa7
	t6.Mul(&z, t6)

	// Step 354: t6 = x^0x22ac17e2f75533786424df9daf0f911f309f0f5bef5c18ff344118d69e41483028e5555380
	for s := 0; s < 7; s++ {
		t6.Square(t6)
	}

	// Step 355: t5 = x^0x22ac17e2f75533786424df9daf0f911f309f0f5bef5c18ff344118d69e41483028e5555397
	t5.Mul(t5, t6)

	// Step 362: t5 = x^0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb80
	for s := 0; s < 7; s++ {
		t5.Square(t5)
	}

	// Step 363: t4 = x^0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d
	t4.Mul(t4, t5)

	// Step 369: t4 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb831fe688231ad3c82906051caaaa72e340
	for s := 0; s < 6; s++ {
		t4.Square(t4)
	}

	// Step 370: t4 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb831fe688231ad3c82906051caaaa72e355
	t4.Mul(t1, t4)

	// Step 376: t4 = x^0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d540
	for s := 0; s < 6; s++ {
		t4.Square(t4)
	}

	// Step 377: t4 = x^0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555
	t4.Mul(t1, t4)

	// Step 381: t4 = x^0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d5550
	for s := 0; s < 4; s++ {
		t4.Square(t4)
	}

	// Step 382: t3 = x^0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d5555
	t3.Mul(t3, t4)

	// Step 389: t3 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa80
	for s := 0; s < 7; s++ {
		t3.Square(t3)
	}

	// Step 390: t2 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa93
	t2.Mul(t2, t3)

	// Step 396: t2 = x^0x22ac17e2f75533786424df9daf0f911f309f0f5bef5c18ff344118d69e41483028e55553971aaaaa4c0
	for s := 0; s < 6; s++ {
		t2.Square(t2)
	}

	// Step 397: t1 = x^0x22ac17e2f75533786424df9daf0f911f309f0f5bef5c18ff344118d69e41483028e55553971aaaaa4d5
	t1.Mul(t1, t2)

	// Step 407: t1 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa935400
	for s := 0; s < 10; s++ {
		t1.Square(t1)
	}

	// Step 408: t1 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa9354ff
	t1.Mul(t0, t1)

	// Step 416: t1 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa9354ff00
	for s := 0; s < 8; s++ {
		t1.Square(t1)
	}

	// Step 417: t1 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa9354ffff
	t1.Mul(t0, t1)

	// Step 425: t1 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa9354ffff00
	for s := 0; s < 8; s++ {
		t1.Square(t1)
	}

	// Step 426: t1 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa9354ffffff
	t1.Mul(t0, t1)

	// Step 434: t1 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa9354ffffff00
	for s := 0; s < 8; s++ {
		t1.Square(t1)
	}

	// Step 435: t0 = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa9354ffffffff
	t0.Mul(t0, t1)

	// Step 438: t0 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb831fe688231ad3c82906051caaaa72e3555549aa7fffffff8
	for s := 0; s < 3; s++ {
		t0.Square(t0)
	}

	// Step 439: t0 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb831fe688231ad3c82906051caaaa72e3555549aa7ffffffff
	t0.Mul(&z, t0)

	// Step 445: t0 = x^0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc0
	for s := 0; s < 6; s++ {
		t0.Square(t0)
	}

	// Step 446: t0 = x^0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc7
	t0.Mul(&z, t0)

	// Step 452: t0 = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb831fe688231ad3c82906051caaaa72e3555549aa7ffffffff1c0
	for s := 0; s < 6; s++ {
		t0.Square(t0)
	}

	// Step 453: z = x^0x45582fc5eeaa66f0c849bf3b5e1f223e613e1eb7deb831fe688231ad3c82906051caaaa72e3555549aa7ffffffff1c7
	z.Mul(&z, t0)

	// Step 454: z = x^0x8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa9354ffffffffe38e
	z.Square(&z)

	return z
}
