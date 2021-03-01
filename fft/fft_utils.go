package fft

import (
	"fmt"

	"github.com/sshravan/go-poly/ff"
)

func (fs *FFTSettings) mulPolysWithFFT(a []ff.Fr, b []ff.Fr, rootsOfUnityStride uint64) []ff.Fr {
	size := fs.MaxWidth / rootsOfUnityStride
	aVals := make([]ff.Fr, size, size)
	bVals := make([]ff.Fr, size, size)
	for i := 0; i < len(a); i++ {
		aVals[i] = a[i]
	}
	for i := len(a); i < len(aVals); i++ {
		aVals[i] = ff.ZERO
	}
	for i := 0; i < len(b); i++ {
		bVals[i] = b[i]
	}
	for i := len(b); i < len(bVals); i++ {
		bVals[i] = ff.ZERO
	}
	rootz := fs.ExpandedRootsOfUnity[:fs.MaxWidth]
	// Get FFT of a and b
	x1 := make([]ff.Fr, len(aVals), len(aVals))
	fs._fft(aVals, 0, 1, rootz, rootsOfUnityStride, x1)

	x2 := make([]ff.Fr, len(bVals), len(bVals))
	fs._fft(bVals, 0, 1, rootz, rootsOfUnityStride, x2)

	// multiply the two. Hack: store results in x1
	var tmp ff.Fr
	for i := 0; i < len(x1); i++ {
		ff.CopyFr(&tmp, &x1[i])
		ff.MulModFr(&x1[i], &tmp, &x2[i])
	}
	revRootz := fs.ReverseRootsOfUnity[:fs.MaxWidth]

	out := make([]ff.Fr, len(x1), len(x1))
	// compute the FFT of the multiplied values.
	fs._fft(x1, 0, 1, revRootz, rootsOfUnityStride, out)
	return out
}

// unshift poly, in-place. Multiplies each coeff with 1/shift_factor**i
func (fs *FFTSettings) ShiftPoly(poly []ff.Fr) {
	var shiftFactor ff.Fr
	ff.AsFr(&shiftFactor, 5) // primitive root of unity
	var factorPower ff.Fr
	ff.CopyFr(&factorPower, &ff.ONE)
	var invFactor ff.Fr
	ff.InvModFr(&invFactor, &shiftFactor)
	var tmp ff.Fr
	for i := 0; i < len(poly); i++ {
		ff.CopyFr(&tmp, &poly[i])
		ff.MulModFr(&poly[i], &tmp, &factorPower)
		// TODO: pre-compute all these shift scalars
		ff.CopyFr(&tmp, &factorPower)
		ff.MulModFr(&factorPower, &tmp, &invFactor)
	}
}

// unshift poly, in-place. Multiplies each coeff with shift_factor**i
func (fs *FFTSettings) UnshiftPoly(poly []ff.Fr) {
	var shiftFactor ff.Fr
	ff.AsFr(&shiftFactor, 5) // primitive root of unity
	var factorPower ff.Fr
	ff.CopyFr(&factorPower, &ff.ONE)
	var tmp ff.Fr
	for i := 0; i < len(poly); i++ {
		ff.CopyFr(&tmp, &poly[i])
		ff.MulModFr(&poly[i], &tmp, &factorPower)
		// TODO: pre-compute all these shift scalars
		ff.CopyFr(&tmp, &factorPower)
		ff.MulModFr(&factorPower, &tmp, &shiftFactor)
	}
}

func (fs *FFTSettings) makeZeroPolyMulLeaf(dst []ff.Fr, indices []uint64, domainStride uint64) {
	if len(dst) < len(indices)+1 {
		panic(fmt.Sprintf("expected bigger destination length: %d, got: %d", len(indices)+1, len(dst)))
	}
	// zero out the unused slots
	for i := len(indices) + 1; i < len(dst); i++ {
		ff.CopyFr(&dst[i], &ff.ZERO)
	}
	ff.CopyFr(&dst[len(indices)], &ff.ONE)
	var negDi ff.Fr
	for i, v := range indices {
		ff.SubModFr(&negDi, &ff.ZERO, &fs.ExpandedRootsOfUnity[v*domainStride])
		ff.CopyFr(&dst[i], &negDi)
		if i > 0 {
			ff.AddModFr(&dst[i], &dst[i], &dst[i-1])
			for j := i - 1; j > 0; j-- {
				ff.MulModFr(&dst[j], &dst[j], &negDi)
				ff.AddModFr(&dst[j], &dst[j], &dst[j-1])
			}
			ff.MulModFr(&dst[0], &dst[0], &negDi)
		}
	}
}

func (fs *FFTSettings) reduceLeaves(scratch []ff.Fr, dst []ff.Fr, ps [][]ff.Fr) {
	n := uint64(len(dst))
	if !ff.IsPowerOfTwo(n) {
		panic("destination must be a power of two")
	}
	if len(ps) == 0 {
		panic("empty leaves")
	}
	if min := uint64(len(ps[0]) * len(ps)); min > n {
		panic(fmt.Sprintf("expected larger destination length: %d, got: %d", min, n))
	}
	if uint64(len(scratch)) < 3*n {
		panic("not enough scratch space")
	}
	// TODO: good to optimize, there's lots of padding
	pPadded := scratch[:n]
	prep := func(pi uint64) {
		p := ps[pi]
		for i := 0; i < len(p); i++ {
			ff.CopyFr(&pPadded[i], &p[i])
		}
		for i := uint64(len(p)); i < n; i++ {
			ff.CopyFr(&pPadded[i], &ff.ZERO)
		}
	}
	mulEvalPs := scratch[n : 2*n]
	pEval := scratch[2*n : 3*n]
	prep(0)
	if err := fs.InplaceFFT(pPadded, mulEvalPs, false); err != nil {
		panic(err)
	}
	for i := uint64(1); i < uint64(len(ps)); i++ {
		prep(i)
		if err := fs.InplaceFFT(pPadded, pEval, false); err != nil {
			panic(err)
		}
		for j := uint64(0); j < n; j++ {
			ff.MulModFr(&mulEvalPs[j], &mulEvalPs[j], &pEval[j])
		}
	}
	if err := fs.InplaceFFT(mulEvalPs, dst, true); err != nil {
		panic(err)
	}
	return
}
