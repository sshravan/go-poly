// +build !bignum_pure,!bignum_hol256,!bignum_kilic,!bignum_hbls

package ff

import (
	"unsafe"

	gmcl "github.com/alinush/go-mcl"
)

func init() {
	gmcl.InitFromString("bls12-381")
	initGlobals()
	ClearG1(&ZERO_G1)
	initG1G2()
}

type Fr gmcl.Fr

func SetFr(dst *Fr, v string) {
	if err := (*gmcl.Fr)(dst).SetString(v, 10); err != nil {
		panic(err)
	}
}

// FrFrom32 mutates the fr num. The value v is little-endian 32-bytes.
func FrFrom32(dst *Fr, v [32]byte) {
	(*gmcl.Fr)(dst).SetLittleEndian(v[:])
}

// FrTo32 serializes a fr number to 32 bytes. Encoded little-endian.
func FrTo32(src *Fr) (v [32]byte) {
	b := (*gmcl.Fr)(src).Serialize()
	last := len(b) - 1
	// reverse endianness, Herumi outputs big-endian bytes
	for i := 0; i < 16; i++ {
		b[i], b[last-i] = b[last-i], b[i]
	}
	copy(v[:], b)
	return
}

func CopyFr(dst *Fr, v *Fr) {
	*dst = *v
}

func AsFr(dst *Fr, i uint64) {
	(*gmcl.Fr)(dst).SetInt64(int64(i))
}

func FrStr(b *Fr) string {
	if b == nil {
		return "<nil>"
	}
	return (*gmcl.Fr)(b).GetString(10)
}

func EqualOne(v *Fr) bool {
	return (*gmcl.Fr)(v).IsOne()
}

func EqualZero(v *Fr) bool {
	return (*gmcl.Fr)(v).IsZero()
}

func EqualFr(a *Fr, b *Fr) bool {
	return (*gmcl.Fr)(a).IsEqual((*gmcl.Fr)(b))
}

func RandomFr() *Fr {
	var out gmcl.Fr
	out.Random()
	return (*Fr)(&out)
}

func SubModFr(dst *Fr, a, b *Fr) {
	gmcl.FrSub((*gmcl.Fr)(dst), (*gmcl.Fr)(a), (*gmcl.Fr)(b))
}

func AddModFr(dst *Fr, a, b *Fr) {
	gmcl.FrAdd((*gmcl.Fr)(dst), (*gmcl.Fr)(a), (*gmcl.Fr)(b))
}

func DivModFr(dst *Fr, a, b *Fr) {
	gmcl.FrDiv((*gmcl.Fr)(dst), (*gmcl.Fr)(a), (*gmcl.Fr)(b))
}

func MulModFr(dst *Fr, a, b *Fr) {
	gmcl.FrMul((*gmcl.Fr)(dst), (*gmcl.Fr)(a), (*gmcl.Fr)(b))
}

func InvModFr(dst *Fr, v *Fr) {
	gmcl.FrInv((*gmcl.Fr)(dst), (*gmcl.Fr)(v))
}

func NegModFr(dst *Fr, v *Fr) {
	gmcl.FrNeg((*gmcl.Fr)(dst), (*gmcl.Fr)(v))
}

//func SqrModFr(dst *Fr, v *Fr) {
//	gmcl.FrSqr((*gmcl.Fr)(dst), (*gmcl.Fr)(v))
//}

func EvalPolyAt(dst *Fr, p []Fr, x *Fr) {
	if err := gmcl.FrEvaluatePolynomial(
		(*gmcl.Fr)(dst),
		*(*[]gmcl.Fr)(unsafe.Pointer(&p)),
		(*gmcl.Fr)(x),
	); err != nil {
		panic(err) // TODO: why does the herumi API return an error? When coefficients are empty?
	}
}

func IntAsFr(dst *Fr, i int64) {
	(*gmcl.Fr)(dst).SetInt64(i)
}

func FromInt64Vec(in []int64) []Fr {
	n := len(in)
	dst := make([]Fr, n, n)
	for i := 0; i < n; i++ {
		(*gmcl.Fr)(&dst[i]).SetInt64(in[i])
	}
	return dst
}

func MulVecFr(a, b []Fr) []Fr {

	n := len(a)
	if n == len(b) && n > 0 {
		result := make([]Fr, n, n)
		for i := 0; i < n; i++ {
			MulModFr(&result[i], &a[i], &b[i])
		}
		return result
	}
	result := make([]Fr, 0)
	return result
}
