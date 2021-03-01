package fft

import (
	"fmt"
	"math/rand"
	"testing"

	"github.com/sshravan/go-poly/ff"
)

func TestFFTSettings_reduceLeaves(t *testing.T) {
	fs := NewFFTSettings(4)

	var fromTreeReduction []ff.Fr

	{
		// prepare some leaves
		leaves := [][]ff.Fr{make([]ff.Fr, 3), make([]ff.Fr, 3), make([]ff.Fr, 3), make([]ff.Fr, 3)}
		leafIndices := [][]uint64{{1, 3}, {7, 8}, {9, 10}, {12, 13}}
		for i := 0; i < 4; i++ {
			fs.makeZeroPolyMulLeaf(leaves[i], leafIndices[i], 1)
		}

		dst := make([]ff.Fr, 16, 16)
		scratch := make([]ff.Fr, 16*3, 16*3)
		fs.reduceLeaves(scratch, dst, leaves)
		fromTreeReduction = dst[:2*4+1]
	}

	var fromDirect []ff.Fr
	{
		dst := make([]ff.Fr, 9, 9)
		indices := []uint64{1, 3, 7, 8, 9, 10, 12, 13}
		fs.makeZeroPolyMulLeaf(dst, indices, 1)
		fromDirect = dst
	}

	if len(fromDirect) != len(fromTreeReduction) {
		t.Fatal("length mismatch")
	}
	for i := 0; i < len(fromDirect); i++ {
		a, b := &fromDirect[i], &fromTreeReduction[i]
		if !ff.EqualFr(a, b) {
			t.Errorf("zero poly coeff %d is different. direct: %s, tree: %s", i, ff.FrStr(a), ff.FrStr(b))
		}
	}
	//debug.DebugFrs("zero poly (tree reduction)", fromTreeReduction)
	//debug.DebugFrs("zero poly (direct slow)", fromDirect)
}

func TestFFTSettings_reduceLeaves_parametrized(t *testing.T) {
	ratios := []float64{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7}
	for scale := uint8(5); scale < 13; scale++ {
		t.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.T) {
			for i, ratio := range ratios {
				t.Run(fmt.Sprintf("ratio_%.3f", ratio), func(t *testing.T) {
					seed := int64(1000*int(scale) + i)
					testReduceLeaves(scale, ratio, seed, t)
				})
			}
		})
	}
}

func testReduceLeaves(scale uint8, missingRatio float64, seed int64, t *testing.T) {
	fs := NewFFTSettings(scale)
	rng := rand.New(rand.NewSource(seed))
	pointCount := uint64(1) << scale
	missingCount := uint64(int(float64(pointCount) * missingRatio))

	// select the missing points
	missing := make([]uint64, pointCount, pointCount)
	for i := uint64(0); i < pointCount; i++ {
		missing[i] = i
	}
	rng.Shuffle(int(pointCount), func(i, j int) {
		missing[i], missing[j] = missing[j], missing[i]
	})
	missing = missing[:missingCount]

	// build the leaves
	pointsPerLeaf := uint64(63)
	leafCount := (missingCount + pointsPerLeaf - 1) / pointsPerLeaf
	leaves := make([][]ff.Fr, leafCount, leafCount)
	for i := uint64(0); i < leafCount; i++ {
		start := i * pointsPerLeaf
		end := start + pointsPerLeaf
		if end > missingCount {
			end = missingCount
		}
		leafSize := end - start
		leaf := make([]ff.Fr, leafSize+1, leafSize+1)
		indices := make([]uint64, leafSize, leafSize)
		for j := uint64(0); j < leafSize; j++ {
			indices[j] = missing[i*pointsPerLeaf+j]
		}
		fs.makeZeroPolyMulLeaf(leaf, indices, 1)
		leaves[i] = leaf
	}

	var fromTreeReduction []ff.Fr

	{
		dst := make([]ff.Fr, pointCount, pointCount)
		scratch := make([]ff.Fr, pointCount*3, pointCount*3)
		fs.reduceLeaves(scratch, dst, leaves)
		fromTreeReduction = dst[:missingCount+1]
	}

	var fromDirect []ff.Fr
	{
		dst := make([]ff.Fr, missingCount+1, missingCount+1)
		fs.makeZeroPolyMulLeaf(dst, missing, fs.MaxWidth/pointCount)
		fromDirect = dst
	}

	if len(fromDirect) != len(fromTreeReduction) {
		t.Fatal("length mismatch")
	}
	for i := 0; i < len(fromDirect); i++ {
		a, b := &fromDirect[i], &fromTreeReduction[i]
		if !ff.EqualFr(a, b) {
			t.Errorf("zero poly coeff %d is different. direct: %s, tree: %s", i, ff.FrStr(a), ff.FrStr(b))
		}
	}
	//debug.DebugFrs("zero poly (tree reduction)", fromTreeReduction)
	//debug.DebugFrs("zero poly (direct slow)", fromDirect)
}
