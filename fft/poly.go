package fft

import (
	"fmt"
	"math/bits"

	"github.com/sshravan/go-poly/ff"
)

// Returns true if polynomial A is a zero polynomial.
func IsPolyZero(a []ff.Fr) bool {

	n := len(a)
	if n == 0 {
		panic("IsPolyZero Error")
	}
	var flag bool
	flag = true
	for i := 0; i < n && flag == true; i++ {
		flag = flag && ff.EqualZero(&a[i])
	}
	return flag
}

//  Removes extraneous zero entries from in vector representation of polynomial.
//  Example - Degree-4 Polynomial: [0, 1, 2, 3, 4, 0, 0, 0, 0] -> [0, 1, 2, 3, 4]
//  Note: Simplest condensed form is a zero polynomial of vector form: [0]
func PolyCondense(a []ff.Fr) []ff.Fr {
	n := len(a)
	if n == 0 {
		panic("PolyCondense Error")
	}

	i := n
	for i > 1 {
		if ff.EqualZero(&a[i-1]) != true {
			break
		}
		i--
	}
	return a[:i]
}

// Computes the standard polynomial addition, polynomial A + polynomial B, and stores result in polynomial C.
func PolyAdd(a []ff.Fr, b []ff.Fr) []ff.Fr {

	if IsPolyZero(a) {
		return PolyCondense(b)
	}

	if IsPolyZero(b) {
		return PolyCondense(a)
	}

	aLen := len(a)
	bLen := len(b)
	n := ff.Max(aLen, bLen)
	c := make([]ff.Fr, n, n)

	for i := 0; i < n; i++ {
		if i < aLen {
			ff.AddModFr(&c[i], &c[i], &a[i])
		}
		if i < bLen {
			ff.AddModFr(&c[i], &c[i], &b[i])
		}
	}
	c = PolyCondense(c)
	return c
}

// Computes the standard polynomial subtraction, polynomial A - polynomial B, and stores result in polynomial C.
func PolySub(a []ff.Fr, b []ff.Fr) []ff.Fr {

	if IsPolyZero(b) {
		return a
	}

	aLen := len(a)
	bLen := len(b)
	n := ff.Max(aLen, bLen)
	c := make([]ff.Fr, n, n)

	for i := 0; i < n; i++ {
		if i < aLen {
			ff.AddModFr(&c[i], &c[i], &a[i])
		}
		if i < bLen {
			ff.SubModFr(&c[i], &c[i], &b[i])
		}
	}
	c = PolyCondense(c)
	return c
}

// Compute a(x) * b(x)
func PolyMul(a []ff.Fr, b []ff.Fr) []ff.Fr {
	if IsPolyZero(a) || IsPolyZero(b) {
		return []ff.Fr{ff.ZERO}
	}

	aLen := len(a)
	bLen := len(b)
	if aLen == bLen && aLen == 1 {
		c := make([]ff.Fr, 1, 1)
		ff.MulModFr(&c[0], &a[0], &b[0])
		return c
	}
	n := uint64(2 * ff.Max(aLen, bLen))
	n = nextPowOf2(n)

	var padding []ff.Fr

	padding = make([]ff.Fr, n-uint64(aLen), n-uint64(aLen))
	a = append(a, padding...)

	padding = make([]ff.Fr, n-uint64(bLen), n-uint64(bLen))
	b = append(b, padding...)

	l := uint8(bits.Len64(n)) - 1 // n = 8 => 3 or 4?
	fs := NewFFTSettings(l)

	evalsA, _ := fs.FFT(a, false)
	evalsB, _ := fs.FFT(b, false)

	res, _ := fs.FFT(ff.MulVecFr(evalsA, evalsB), true)
	res = res[:(aLen + bLen - 1)]
	res = PolyCondense(res)
	return res
}

// Builds the polynomial from its roots
// (x - a_1)(x - a_2)(x - a_3)(x - a_4)(x - a_5)(1)(1)(1)
// Need not be a power of two
func PolyTree(a []ff.Fr) []ff.Fr {

	n := uint64(len(a))
	aLen := n
	n = nextPowOf2(n)

	var padding []ff.Fr

	padding = make([]ff.Fr, n-aLen, n-aLen)
	a = append(a, padding...)

	l := uint8(bits.Len64(n)) - 1

	var M [][]ff.Fr

	M = make([][]ff.Fr, n, n)
	for j := uint64(0); j < n; j++ {
		if j < aLen {
			M[j] = make([]ff.Fr, 2)
			ff.NegModFr(&M[j][0], &a[j])
			ff.IntAsFr(&M[j][1], 1)
		} else {
			M[j] = []ff.Fr{ff.ONE}
		}
	}

	var x []ff.Fr
	var y []ff.Fr
	var index int64
	index = 0
	for i := uint8(1); i <= l; i++ {
		L := uint64(1) << (l - i)
		m := make([][]ff.Fr, L, L)
		for j := uint64(0); j < L; j++ {
			x = M[index]
			index++
			y = M[index]
			index++
			m[j] = PolyMul(x, y)
		}
		index = 0
		M = m
	}
	return PolyCondense(M[0])
}

// Invert the divisor, then multiply
func polyFactorDiv(dst *ff.Fr, a *ff.Fr, b *ff.Fr) {
	// TODO: use divmod instead.
	var tmp ff.Fr
	ff.InvModFr(&tmp, b)
	ff.MulModFr(dst, &tmp, a)
}

// Long polynomial division for two polynomials in coefficient form
// Need to check divide by zero
func PolyLongDiv(A []ff.Fr, B []ff.Fr) []ff.Fr {
	if IsPolyZero(B) == true {
		panic("PolyLongDiv: Cannot divide by zero polynomial.")
	}
	a := make([]ff.Fr, len(A), len(A))
	for i := 0; i < len(a); i++ {
		ff.CopyFr(&a[i], &A[i])
	}
	aPos := len(a) - 1
	bPos := len(B) - 1
	diff := aPos - bPos
	out := make([]ff.Fr, diff+1, diff+1)
	for diff >= 0 {
		quot := &out[diff]
		polyFactorDiv(quot, &a[aPos], &B[bPos])
		var tmp, tmp2 ff.Fr
		for i := bPos; i >= 0; i-- {
			// In steps: a[diff + i] -= b[i] * quot
			// tmp =  b[i] * quot
			ff.MulModFr(&tmp, quot, &B[i])
			// tmp2 = a[diff + i] - tmp
			ff.SubModFr(&tmp2, &a[diff+i], &tmp)
			// a[diff + i] = tmp2
			ff.CopyFr(&a[diff+i], &tmp2)
		}
		aPos -= 1
		diff -= 1
	}
	return out
}

// Computes q(x) and r(x) s.t. a(x) = q(x) * b(x) + r(x)
func PolyDiv(A []ff.Fr, B []ff.Fr) ([]ff.Fr, []ff.Fr) {
	if IsPolyZero(B) == true {
		panic("PolyDiv: Cannot divide by zero polynomial.")
	}

	if len(B) > len(A) {
		panic("PolyDiv: Deg(B) should be <= Ded(A)")
	}

	a := make([]ff.Fr, len(A), len(A))
	for i := 0; i < len(a); i++ {
		ff.CopyFr(&a[i], &A[i])
	}
	aPos := len(a) - 1
	bPos := len(B) - 1
	diff := aPos - bPos
	out := make([]ff.Fr, diff+1, diff+1)
	for diff >= 0 {
		quot := &out[diff]
		polyFactorDiv(quot, &a[aPos], &B[bPos])
		var tmp, tmp2 ff.Fr
		for i := bPos; i >= 0; i-- {
			// In steps: a[diff + i] -= b[i] * quot
			// tmp =  b[i] * quot
			ff.MulModFr(&tmp, quot, &B[i])
			// tmp2 = a[diff + i] - tmp
			ff.SubModFr(&tmp2, &a[diff+i], &tmp)
			// a[diff + i] = tmp2
			ff.CopyFr(&a[diff+i], &tmp2)
		}
		aPos -= 1
		diff -= 1
	}
	a = PolyCondense(a)
	return out, a
}

// Given M(x) computes  M'(x)
// For first derivative and multi-point evaluation:
// This [paper](https://arxiv.org/pdf/2002.10304.pdf) claims there is a faster way to do this
// Page 10 under interpolation
// Not sure how it is different from doing subproduct tree on M'(x)
func PolyDifferentiate(a []ff.Fr) []ff.Fr {
	n := int64(len(a))
	if n == 0 {
		panic("PolyDifferentiate: Input is empty")
	}
	if n == 1 {
		return make([]ff.Fr, 1)
	}
	c := make([]ff.Fr, n, n)
	var temp ff.Fr
	for i := int64(1); i < n; i++ {
		ff.IntAsFr(&temp, i)
		ff.MulModFr(&c[i], &a[i], &temp)
	}
	return c[1:]
}

// Extended GCG: Computes u(x) and v(x) s.t. u(x) * a(x) + v(x) * b(x) = g(x)
func xGCD(a []ff.Fr, b []ff.Fr) (g []ff.Fr, u []ff.Fr, v []ff.Fr) {
	return xGCD2(a, b)
}

// Computes Extended GCD using pseudocode **#1** here:
// https://en.wikipedia.org/w/index.php?title=Extended_Euclidean_algorithm&oldid=1003613686
// a * u + b * v = g
func xGCD1(a []ff.Fr, b []ff.Fr) (g []ff.Fr, u []ff.Fr, v []ff.Fr) {

	if len(b) > len(a) {
		g, v, u := xGCD1(b, a)
		return g, u, v
	}

	old_r, r := a, b
	old_s, s := []ff.Fr{ff.ONE}, []ff.Fr{ff.ZERO}
	old_t, t := []ff.Fr{ff.ZERO}, []ff.Fr{ff.ONE}

	for IsPolyZero(r) == false {
		quotient, remainder := PolyDiv(old_r, r)
		old_r, r = r, remainder
		old_s, s = s, PolySub(old_s, PolyMul(quotient, s))
		old_t, t = t, PolySub(old_t, PolyMul(quotient, t))
	}

	old_r = PolyCondense(old_r)
	old_s = PolyCondense(old_s)
	old_t = PolyCondense(old_t)
	return old_r, old_s, old_t
}

// Computes Extended GCD using pseudocode **#2** here:
// https://en.wikipedia.org/w/index.php?title=Extended_Euclidean_algorithm&oldid=1003613686
// a * u + b * v = g
func xGCD2(a []ff.Fr, b []ff.Fr) (g []ff.Fr, u []ff.Fr, v []ff.Fr) {

	if len(b) > len(a) {
		g, v, u := xGCD2(b, a)
		return g, u, v
	} else {
		s := []ff.Fr{ff.ZERO}
		old_s := []ff.Fr{ff.ONE}

		r := b
		old_r := a

		for IsPolyZero(r) == false {
			quotient, remainder := PolyDiv(old_r, r)
			old_r, r = r, remainder
			old_s, s = s, PolySub(old_s, PolyMul(quotient, s))
		}

		var bezout_t []ff.Fr
		if IsPolyZero(b) == false {
			bezout_t, _ = PolyDiv(PolySub(old_r, PolyMul(old_s, a)), b)
		} else {
			bezout_t = []ff.Fr{ff.ZERO}
		}
		return old_r, old_s, bezout_t
	}
}

// SubProdTree
// Needs to be power of two
// (x - a_1)(x - a_2)(x - a_3)(x - a_4)(x - a_5)(x - a_6)(x - a_7)(x - a_8)
func SubProductTree(a []ff.Fr) [][][]ff.Fr {

	n := uint64(len(a))

	if ff.IsPowerOfTwo(n) == false {
		panic("SubProductTree inputs needs to be power of two")
	}

	l := uint8(bits.Len64(n)) - 1
	// fmt.Println(l)
	var M [][][]ff.Fr
	M = make([][][]ff.Fr, l+1, l+1)
	M[0] = make([][]ff.Fr, n, n)
	for j := uint64(0); j < n; j++ {
		M[0][j] = make([]ff.Fr, 2)
		ff.NegModFr(&M[0][j][0], &a[j])
		ff.IntAsFr(&M[0][j][1], 1)
	}

	var x []ff.Fr
	var y []ff.Fr
	var index int64
	index = 0
	for i := uint8(1); i <= l; i++ {
		L := uint64(1) << (l - i)
		M[i] = make([][]ff.Fr, L, L)
		for j := uint64(0); j < L; j++ {
			x = M[i-1][index]
			index++
			y = M[i-1][index]
			index++
			M[i][j] = PolyMul(x, y)
		}
		index = 0
	}

	return M
}

// Fast multi-point evaluation using subproduct tree
// For n^ directly use EvalPolyAt $n$ times
func PolyMultiEvaluate(f []ff.Fr, M [][][]ff.Fr) []ff.Fr {
	n := int64(len(f))
	if n == 0 {
		panic("PolyDifferentiate: Input is empty")
	}
	if n == 1 {
		return f
	}

	k := len(M) - 1
	if !(1<<k >= n) {
		msg := fmt.Sprintf("Subproduct tree and the polynomial size did not match\n\t len(f): %d len(M): %d", n, 1<<k)
		panic(msg)
	}
	_, aL := PolyDiv(f, M[k-1][0])
	_, aR := PolyDiv(f, M[k-1][1])

	mL, mR := splitSubProdTree(M)

	l := PolyMultiEvaluate(aL, mL)
	r := PolyMultiEvaluate(aR, mR)

	return append(l, r...)
}

// Grossly assumes that it is a proper tree
// Given a SubProduct tree, divides into 2 sub-product trees
func splitSubProdTree(M [][][]ff.Fr) ([][][]ff.Fr, [][][]ff.Fr) {

	k := len(M) - 1
	L := make([][][]ff.Fr, k)
	R := make([][][]ff.Fr, k)

	for i := 0; i < k; i++ {
		n := len(M[i])
		L[i] = append(L[i], M[i][:(n/2)]...)
		R[i] = append(R[i], M[i][(n/2):]...)
	}

	return L, R

}
