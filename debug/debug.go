package debug

import (
	"fmt"
	"strings"

	"github.com/sshravan/go-poly/ff"
)

func DebugFrPtrs(msg string, values []*ff.Fr) {
	var out strings.Builder
	out.WriteString("---")
	out.WriteString(msg)
	out.WriteString("---\n")
	for i := range values {
		out.WriteString(fmt.Sprintf("#%4d: %s\n", i, ff.FrStr(values[i])))
	}
	fmt.Println(out.String())
}

func DebugFrs(msg string, values []ff.Fr) {
	fmt.Println("---------------------------")
	var out strings.Builder
	for i := range values {
		out.WriteString(fmt.Sprintf("%s %d: %s\n", msg, i, ff.FrStr(&values[i])))
	}
	fmt.Print(out.String())
}
