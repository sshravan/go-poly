// +build !bignum_pure,!bignum_hol256

package fft

import (
	"fmt"
	"testing"

	"github.com/sshravan/go-poly/ff"
)

func benchFFTG1(scale uint8, b *testing.B) {
	fs := NewFFTSettings(scale)
	data := make([]ff.G1Point, fs.MaxWidth, fs.MaxWidth)
	for i := uint64(0); i < fs.MaxWidth; i++ {
		var tmpG1 ff.G1Point
		ff.CopyG1(&tmpG1, &ff.GenG1)
		ff.MulG1(&data[i], &tmpG1, ff.RandomFr())
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		out, err := fs.FFTG1(data, false)
		if err != nil {
			b.Fatal(err)
		}
		if len(out) != len(data) {
			panic("output len doesn't match input")
		}
	}
}

func BenchmarkFFTSettings_FFTG1(b *testing.B) {
	for scale := uint8(4); scale < 16; scale++ {
		b.Run(fmt.Sprintf("scale_%d", scale), func(b *testing.B) {
			benchFFTG1(scale, b)
		})
	}
}
