// +build !bignum_pure,!bignum_hol256

package fft

import (
	"fmt"

	"github.com/sshravan/go-poly/ff"
)

func (fs *FFTSettings) simpleFTG1(vals []ff.G1Point, valsOffset uint64, valsStride uint64, rootsOfUnity []ff.Fr, rootsOfUnityStride uint64, out []ff.G1Point) {
	l := uint64(len(out))
	var v ff.G1Point
	var tmp ff.G1Point
	var last ff.G1Point
	for i := uint64(0); i < l; i++ {
		jv := &vals[valsOffset]
		r := &rootsOfUnity[0]
		ff.MulG1(&v, jv, r)
		ff.CopyG1(&last, &v)

		for j := uint64(1); j < l; j++ {
			jv := &vals[valsOffset+j*valsStride]
			r := &rootsOfUnity[((i*j)%l)*rootsOfUnityStride]
			ff.MulG1(&v, jv, r)
			ff.CopyG1(&tmp, &last)
			ff.AddG1(&last, &tmp, &v)
		}
		ff.CopyG1(&out[i], &last)
	}
}

func (fs *FFTSettings) _fftG1(vals []ff.G1Point, valsOffset uint64, valsStride uint64, rootsOfUnity []ff.Fr, rootsOfUnityStride uint64, out []ff.G1Point) {
	if len(out) <= 4 { // if the value count is small, run the unoptimized version instead. // TODO tune threshold. (can be different for G1)
		fs.simpleFTG1(vals, valsOffset, valsStride, rootsOfUnity, rootsOfUnityStride, out)
		return
	}

	half := uint64(len(out)) >> 1
	// L will be the left half of out
	fs._fftG1(vals, valsOffset, valsStride<<1, rootsOfUnity, rootsOfUnityStride<<1, out[:half])
	// R will be the right half of out
	fs._fftG1(vals, valsOffset+valsStride, valsStride<<1, rootsOfUnity, rootsOfUnityStride<<1, out[half:]) // just take even again

	var yTimesRoot ff.G1Point
	var x, y ff.G1Point
	for i := uint64(0); i < half; i++ {
		// temporary copies, so that writing to output doesn't conflict with input
		ff.CopyG1(&x, &out[i])
		ff.CopyG1(&y, &out[i+half])
		root := &rootsOfUnity[i*rootsOfUnityStride]
		ff.MulG1(&yTimesRoot, &y, root)
		ff.AddG1(&out[i], &x, &yTimesRoot)
		ff.SubG1(&out[i+half], &x, &yTimesRoot)
	}
}

func (fs *FFTSettings) FFTG1(vals []ff.G1Point, inv bool) ([]ff.G1Point, error) {
	n := uint64(len(vals))
	if n > fs.MaxWidth {
		return nil, fmt.Errorf("got %d values but only have %d roots of unity", n, fs.MaxWidth)
	}
	if !ff.IsPowerOfTwo(n) {
		return nil, fmt.Errorf("got %d values but not a power of two", n)
	}
	// We make a copy so we can mutate it during the work.
	valsCopy := make([]ff.G1Point, n, n)
	for i := 0; i < len(vals); i++ { // TODO: maybe optimize this away, and write back to original input array?
		ff.CopyG1(&valsCopy[i], &vals[i])
	}
	if inv {
		var invLen ff.Fr
		ff.AsFr(&invLen, n)
		ff.InvModFr(&invLen, &invLen)
		rootz := fs.ReverseRootsOfUnity[:fs.MaxWidth]
		stride := fs.MaxWidth / n

		out := make([]ff.G1Point, n, n)
		fs._fftG1(valsCopy, 0, 1, rootz, stride, out)
		var tmp ff.G1Point
		for i := 0; i < len(out); i++ {
			ff.MulG1(&tmp, &out[i], &invLen)
			ff.CopyG1(&out[i], &tmp)
		}
		return out, nil
	} else {
		out := make([]ff.G1Point, n, n)
		rootz := fs.ExpandedRootsOfUnity[:fs.MaxWidth]
		stride := fs.MaxWidth / n
		// Regular FFT
		fs._fftG1(valsCopy, 0, 1, rootz, stride, out)
		return out, nil
	}
}

// rearrange G1 elements in reverse bit order. Supports 2**31 max element count.
func ReverseBitOrderG1(values []ff.G1Point) {
	if len(values) > (1 << 31) {
		panic("list too large")
	}
	var tmp ff.G1Point
	reverseBitOrder(uint32(len(values)), func(i, j uint32) {
		ff.CopyG1(&tmp, &values[i])
		ff.CopyG1(&values[i], &values[j])
		ff.CopyG1(&values[j], &tmp)
	})
}
