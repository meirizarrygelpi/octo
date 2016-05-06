// Copyright (c) 2016 Melvin Eloy Irizarry-Gelpí
// Licenced under the MIT License.

package octo

import (
	"fmt"
	"math"
	"strings"

	"github.com/meirizarrygelpi/quat"
)

var symbKlein = [8]string{"", "i", "j", "k", "s", "t", "u", "v"}

// A Klein represents a Klein octonion (also known as a split-octonion) as an
// ordered array of two pointers to quat.Hamilton values.
type Klein [2]*quat.Hamilton

// String.
func (z *Klein) String() string {
	v := make([]float64, 8)
	v[0], v[1] = real((z[0])[0]), imag((z[0])[0])
	v[2], v[3] = real((z[0])[1]), imag((z[0])[1])
	v[4], v[5] = real((z[1])[0]), imag((z[1])[0])
	v[6], v[7] = real((z[1])[1]), imag((z[1])[1])
	a := make([]string, 17)
	a[0] = "("
	a[1] = fmt.Sprintf("%g", v[0])
	i := 1
	for j := 2; j < 16; j = j + 2 {
		switch {
		case math.Signbit(v[i]):
			a[j] = fmt.Sprintf("%g", v[i])
		case math.IsInf(v[i], +1):
			a[j] = "+Inf"
		default:
			a[j] = fmt.Sprintf("+%g", v[i])
		}
		a[j+1] = symbKlein[i]
		i++
	}
	a[16] = ")"
	return strings.Join(a, "")
}

// Equals returns true if y and z are equal.
func (z *Klein) Equals(y *Klein) bool {
	if !z[0].Equals(y[0]) || !z[1].Equals(y[1]) {
		return false
	}
	return true
}

// Copy copies y onto z, and returns z.
func (z *Klein) Copy(y *Klein) *Klein {
	z[0] = new(quat.Hamilton).Copy(y[0])
	z[1] = new(quat.Hamilton).Copy(y[1])
	return z
}

// NewKlein returns a pointer to a Klein value made from eight given float64
// values.
func NewKlein(a, b, c, d, e, f, g, h float64) *Klein {
	z := new(Klein)
	z[0] = quat.NewHamilton(a, b, c, d)
	z[1] = quat.NewHamilton(e, f, g, h)
	return z
}

// IsInf returns true if any of the components of z are infinite.
func (z *Klein) IsInf() bool {
	if z[0].IsInf() || z[1].IsInf() {
		return true
	}
	return false
}

// KleinInf returns a pointer to a Klein octonionic infinity value.
func KleinInf(a, b, c, d, e, f, g, h int) *Klein {
	z := new(Klein)
	z[0] = quat.HamiltonInf(a, b, c, d)
	z[1] = quat.HamiltonInf(e, f, g, h)
	return z
}

// IsNaN returns true if any component of z is NaN and neither is an
// infinity.
func (z *Klein) IsNaN() bool {
	if z[0].IsInf() || z[1].IsInf() {
		return false
	}
	if z[0].IsNaN() || z[1].IsNaN() {
		return true
	}
	return false
}

// KleinNaN returns a pointer to a Klein octonionic NaN value.
func KleinNaN() *Klein {
	z := new(Klein)
	z[0] = quat.HamiltonNaN()
	z[1] = quat.HamiltonNaN()
	return z
}

// ScalR sets z equal to y scaled by a on the right, and returns z.
//
// This is a special case of Mul:
// 		ScalR(y, a) = Mul(y, Hamilton{a, 0})
func (z *Klein) ScalR(y *Klein, a *quat.Hamilton) *Klein {
	z[0] = new(quat.Hamilton).Mul(y[0], a)
	z[1] = new(quat.Hamilton).Mul(y[1], a)
	return z
}

// ScalL sets z equal to y scaled by a on the left, and returns z.
//
// This is a special case of Mul:
// 		ScalL(y, a) = Mul(Hamilton{a, 0}, y)
func (z *Klein) ScalL(a *quat.Hamilton, y *Klein) *Klein {
	z[0] = new(quat.Hamilton).Mul(a, y[0])
	z[1] = new(quat.Hamilton).Mul(a, y[1])
	return z
}

// Dil sets z equal to the dilation of y by a, and returns z.
//
// This is a special case of Mul:
// 		Dil(y, a) = Mul(y, Hamilton{quat.Hamilton{a, 0, 0, 0}, 0})
func (z *Klein) Dil(y *Klein, a float64) *Klein {
	z[0] = new(quat.Hamilton).Dil(y[0], a)
	z[1] = new(quat.Hamilton).Dil(y[1], a)
	return z
}

// Neg sets z equal to the negative of y, and returns z.
func (z *Klein) Neg(y *Klein) *Klein {
	return z.Dil(y, -1)
}

// Conj sets z equal to the conjugate of y, and returns z.
func (z *Klein) Conj(y *Klein) *Klein {
	z[0] = new(quat.Hamilton).Conj(y[0])
	z[1] = new(quat.Hamilton).Neg(y[1])
	return z
}

// Add sets z equal to the sum of x and y, and returns z.
func (z *Klein) Add(x, y *Klein) *Klein {
	z[0] = new(quat.Hamilton).Add(x[0], y[0])
	z[1] = new(quat.Hamilton).Add(x[1], y[1])
	return z
}

// Sub sets z equal to the difference of x and y, and returns z.
func (z *Klein) Sub(x, y *Klein) *Klein {
	z[0] = new(quat.Hamilton).Sub(x[0], y[0])
	z[1] = new(quat.Hamilton).Sub(x[1], y[1])
	return z
}

// Mul sets z equal to the noncommutative, nonassociative product of x and y,
// and returns z.
func (z *Klein) Mul(x, y *Klein) *Klein {
	p := new(Klein).Copy(x)
	q := new(Klein).Copy(y)
	z[0] = new(quat.Hamilton).Add(
		new(quat.Hamilton).Mul(p[0], q[0]),
		new(quat.Hamilton).Mul(new(quat.Hamilton).Conj(q[1]), p[1]),
	)
	z[1] = new(quat.Hamilton).Add(
		new(quat.Hamilton).Mul(q[1], p[0]),
		new(quat.Hamilton).Mul(p[1], q[0].Conj(q[0])),
	)
	return z
}

// Commutator sets z equal to the commutator of x and y, and returns z.
func (z *Klein) Commutator(x, y *Klein) *Klein {
	return z.Sub(new(Klein).Mul(x, y), new(Klein).Mul(y, x))
}

// Associator sets z equal to the associator of w, x, and y, and returns z.
func (z *Klein) Associator(w, x, y *Klein) *Klein {
	return z.Sub(
		new(Klein).Mul(new(Klein).Mul(w, x), y),
		new(Klein).Mul(w, new(Klein).Mul(x, y)),
	)
}

// Quad returns the non-negative quadrance of z.
func (z *Klein) Quad() float64 {
	a, b := z[0].Quad(), z[1].Quad()
	return a - b
}
