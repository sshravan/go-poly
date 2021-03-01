package fft

import (
	"fmt"
	"testing"

	"github.com/sshravan/go-poly/ff"
)

func BenchmarkPolyMul(b *testing.B) {

	for scale := uint8(10); scale < 20; scale++ {
		n := uint64(1) << scale
		A := make([]ff.Fr, n, n)
		B := make([]ff.Fr, n, n)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
			B[i] = *(ff.RandomFr())
		}
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			b.ResetTimer()
			_ = PolyMul(A, B)
		})
	}
}

func BenchmarkPolyTree(b *testing.B) {

	for scale := uint8(10); scale < 19; scale++ {
		n := uint64(1) << scale
		A := make([]ff.Fr, n, n)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
		}
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			b.ResetTimer()
			_ = PolyTree(A)
		})
	}
}

func BenchmarkPolyDiv(b *testing.B) {

	for scale := uint8(10); scale < 20; scale++ {
		n := uint64(1) << scale
		A := make([]ff.Fr, n, n)
		B := make([]ff.Fr, 2, 2)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
		}
		B[0] = *(ff.RandomFr())
		B[1] = *(ff.RandomFr())
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			b.ResetTimer()
			_, _ = PolyDiv(A, B)
		})
	}
}

func BenchmarkPolySubProductTree(b *testing.B) {

	for scale := uint8(10); scale < 15; scale++ {
		n := uint64(1) << scale
		A := make([]ff.Fr, n, n)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
		}
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			b.ResetTimer()
			_ = SubProductTree(A)
		})
	}
}

func BenchmarkPolyMultiEvaluate(b *testing.B) {

	for scale := uint8(10); scale < 15; scale++ {
		n := uint64(1) << scale
		A := make([]ff.Fr, n, n)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
		}
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			M := SubProductTree(A)
			b.ResetTimer()
			_ = PolyMultiEvaluate(A, M)
		})
	}
}

func BenchmarkPolyDifferentiate(b *testing.B) {

	for scale := uint8(10); scale < 20; scale++ {
		n := uint64(1) << scale
		data := make([]ff.Fr, n, n)
		for i := uint64(0); i < n; i++ {
			data[i] = *(ff.RandomFr())
		}
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			b.ResetTimer()
			_ = PolyDifferentiate(data)
		})
	}
}

func BenchmarkPolyXGCD1Balanced(b *testing.B) {

	for scale := uint8(10); scale < 12; scale++ {
		n := uint64(1) << scale
		A := make([]ff.Fr, n, n)
		B := make([]ff.Fr, n, n)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
			B[i] = *(ff.RandomFr())
		}
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			b.ResetTimer()
			_, _, _ = xGCD1(A, B)
		})
	}
}

func BenchmarkPolyXGCD1Lopsided(b *testing.B) {

	for scale := uint8(10); scale < 12; scale++ {
		n := uint64(1) << scale
		A := make([]ff.Fr, n, n)
		B := make([]ff.Fr, 2, 2)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
		}
		B[0] = *(ff.RandomFr())
		B[1] = *(ff.RandomFr())
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			b.ResetTimer()
			_, _, _ = xGCD1(A, B)
		})
	}
}

func BenchmarkPolyXGCD2Balanced(b *testing.B) {

	for scale := uint8(10); scale < 12; scale++ {
		n := uint64(1) << scale
		A := make([]ff.Fr, n, n)
		B := make([]ff.Fr, n, n)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
			B[i] = *(ff.RandomFr())
		}
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			b.ResetTimer()
			_, _, _ = xGCD2(A, B)
		})
	}
}

func BenchmarkPolyXGCD2Lopsided(b *testing.B) {

	for scale := uint8(10); scale < 12; scale++ {
		n := uint64(1) << scale
		A := make([]ff.Fr, n, n)
		B := make([]ff.Fr, 2, 2)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
		}
		B[0] = *(ff.RandomFr())
		B[1] = *(ff.RandomFr())
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			b.ResetTimer()
			_, _, _ = xGCD2(A, B)
		})
	}
}
