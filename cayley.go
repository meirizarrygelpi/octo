package octo

import (
	"fmt"
	"math"
	"strings"

	"github.com/meirizarrygelpi/quat"
)

var symbCayley = [8]string{"", "i", "j", "k", "m", "n", "p", "q"}

// A Cayley represents a Cayley octonion (i.e. a traditional octonion) as an
// ordered array of two pointers to quat.Hamilton values.
type Cayley [2]*quat.Hamilton

// String returns the string version of a Cayley value. If z corresponds to the
// Cayley octonion a + bi + cj + dk + em + fn + gp + hq, then the string is
// "(a+bi+cj+dk+em+fn+gp+hq)", similar to complex128 values.
func (z *Cayley) String() string {
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
		a[j+1] = symbCayley[i]
		i++
	}
	a[16] = ")"
	return strings.Join(a, "")
}

// Equals returns true if y and z are equal.
func (z *Cayley) Equals(y *Cayley) bool {
	if !z[0].Equals(y[0]) || !z[1].Equals(y[1]) {
		return false
	}
	return true
}

// Copy copies y onto z, and returns z.
func (z *Cayley) Copy(y *Cayley) *Cayley {
	z[0] = new(quat.Hamilton).Copy(y[0])
	z[1] = new(quat.Hamilton).Copy(y[1])
	return z
}

// NewCayley returns a pointer to a Cayley value made from eight given float64
// values.
func NewCayley(a, b, c, d, e, f, g, h float64) *Cayley {
	z := new(Cayley)
	z[0] = quat.NewHamilton(a, b, c, d)
	z[1] = quat.NewHamilton(e, f, g, h)
	return z
}

// IsInf returns true if any of the components of z are infinite.
func (z *Cayley) IsInf() bool {
	if z[0].IsInf() || z[1].IsInf() {
		return true
	}
	return false
}

// CayleyInf returns a pointer to a Cayley octonionic infinity value.
func CayleyInf(a, b, c, d, e, f, g, h int) *Cayley {
	z := new(Cayley)
	z[0] = quat.HamiltonInf(a, b, c, d)
	z[1] = quat.HamiltonInf(e, f, g, h)
	return z
}

// IsNaN returns true if any component of z is NaN and neither is an
// infinity.
func (z *Cayley) IsNaN() bool {
	if z[0].IsInf() || z[1].IsInf() {
		return false
	}
	if z[0].IsNaN() || z[1].IsNaN() {
		return true
	}
	return false
}

// CayleyNaN returns a pointer to a Cayley octonionic NaN value.
func CayleyNaN() *Cayley {
	z := new(Cayley)
	z[0] = quat.HamiltonNaN()
	z[1] = quat.HamiltonNaN()
	return z
}

// ScalR sets z equal to y scaled by a on the right, and returns z.
//
// This is a special case of Mul:
// 		ScalR(y, a) = Mul(y, Hamilton{a, 0})
func (z *Cayley) ScalR(y *Cayley, a *quat.Hamilton) *Cayley {
	z[0] = new(quat.Hamilton).Mul(y[0], a)
	z[1] = new(quat.Hamilton).Mul(y[1], a)
	return z
}

// ScalL sets z equal to y scaled by a on the left, and returns z.
//
// This is a special case of Mul:
// 		ScalL(y, a) = Mul(Hamilton{a, 0}, y)
func (z *Cayley) ScalL(a *quat.Hamilton, y *Cayley) *Cayley {
	z[0] = new(quat.Hamilton).Mul(a, y[0])
	z[1] = new(quat.Hamilton).Mul(a, y[1])
	return z
}

// Dil sets z equal to the dilation of y by a, and returns z.
//
// This is a special case of Mul:
// 		Dil(y, a) = Mul(y, Hamilton{quat.Hamilton{a, 0, 0, 0}, 0})
func (z *Cayley) Dil(y *Cayley, a float64) *Cayley {
	z[0] = new(quat.Hamilton).Dil(y[0], a)
	z[1] = new(quat.Hamilton).Dil(y[1], a)
	return z
}

// Neg sets z equal to the negative of y, and returns z.
func (z *Cayley) Neg(y *Cayley) *Cayley {
	return z.Dil(y, -1)
}

// Conj sets z equal to the conjugate of y, and returns z.
func (z *Cayley) Conj(y *Cayley) *Cayley {
	z[0] = new(quat.Hamilton).Conj(y[0])
	z[1] = new(quat.Hamilton).Neg(y[1])
	return z
}

// Add sets z equal to the sum of x and y, and returns z.
func (z *Cayley) Add(x, y *Cayley) *Cayley {
	z[0] = new(quat.Hamilton).Add(x[0], y[0])
	z[1] = new(quat.Hamilton).Add(x[1], y[1])
	return z
}

// Sub sets z equal to the difference of x and y, and returns z.
func (z *Cayley) Sub(x, y *Cayley) *Cayley {
	z[0] = new(quat.Hamilton).Sub(x[0], y[0])
	z[1] = new(quat.Hamilton).Sub(x[1], y[1])
	return z
}

// Mul sets z equal to the noncommutative, nonassociative product of x and y,
// and returns z.
func (z *Cayley) Mul(x, y *Cayley) *Cayley {
	p := new(Cayley).Copy(x)
	q := new(Cayley).Copy(y)
	z[0] = new(quat.Hamilton).Sub(
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
func (z *Cayley) Commutator(x, y *Cayley) *Cayley {
	return z.Sub(new(Cayley).Mul(x, y), new(Cayley).Mul(y, x))
}

// Associator sets z equal to the associator of w, x, and y, and returns z.
func (z *Cayley) Associator(w, x, y *Cayley) *Cayley {
	return z.Sub(
		new(Cayley).Mul(new(Cayley).Mul(w, x), y),
		new(Cayley).Mul(w, new(Cayley).Mul(x, y)),
	)
}

// Quad returns the non-negative quadrance of z.
func (z *Cayley) Quad() float64 {
	a, b := z[0].Quad(), z[1].Quad()
	return a + b
}

// Inv sets z equal to the inverse of y, and returns z. If y is zero, then Inv
// panics.
func (z *Cayley) Inv(y *Cayley) *Cayley {
	if y.Equals(&Cayley{&quat.Hamilton{0, 0}, &quat.Hamilton{0, 0}}) {
		panic("inverse of zero")
	}
	return z.Dil(new(Cayley).Conj(y), 1/y.Quad())
}

// Quo sets z equal to the quotient of x and y, and returns z. If y is zero,
// then Quo panics.
func (z *Cayley) Quo(x, y *Cayley) *Cayley {
	if y.Equals(&Cayley{&quat.Hamilton{0, 0}, &quat.Hamilton{0, 0}}) {
		panic("denominator is zero")
	}
	return z.Dil(new(Cayley).Mul(x, new(Cayley).Conj(y)), 1/y.Quad())
}

// RectCayley returns a Cayley value made from given curvilinear
// coordinates.
func RectCayley(r, θ1, θ2, θ3, θ4, θ5, θ6, θ7 float64) *Cayley {
	z := new(Cayley)
	// z[0] = r * math.Cos(θ1)
	// z[1] = r * math.Sin(θ1) * math.Cos(θ2)
	// z[2] = r * math.Sin(θ1) * math.Sin(θ2) * math.Cos(θ3)
	// z[3] = r * math.Sin(θ1) * math.Sin(θ2) * math.Sin(θ3) * math.Cos(θ4)
	// z[4] = r * math.Sin(θ1) * math.Sin(θ2) * math.Sin(θ3) * math.Sin(θ4) *
	// 	math.Cos(θ5)
	// z[5] = r * math.Sin(θ1) * math.Sin(θ2) * math.Sin(θ3) * math.Sin(θ4) *
	// 	math.Sin(θ5) * math.Cos(θ6)
	// z[6] = r * math.Sin(θ1) * math.Sin(θ2) * math.Sin(θ3) * math.Sin(θ4) *
	// 	math.Sin(θ5) * math.Sin(θ6) * math.Cos(θ7)
	// z[7] = r * math.Sin(θ1) * math.Sin(θ2) * math.Sin(θ3) * math.Sin(θ4) *
	// 	math.Sin(θ5) * math.Sin(θ6) * math.Sin(θ7)
	return z
}

// Curv returns the curvilinear coordinates of a Cayley value.
func (z *Cayley) Curv() (r, θ1, θ2, θ3, θ4, θ5, θ6, θ7 float64) {
	// h76 := math.Hypot(z[7], z[6])
	// h75 := math.Hypot(h76, z[5])
	// h74 := math.Hypot(h75, z[4])
	// h73 := math.Hypot(h74, z[3])
	// h72 := math.Hypot(h73, z[2])
	// h71 := math.Hypot(h72, z[1])
	// r = math.Sqrt(z.Quad())
	// θ1 = math.Atan(h71 / z[0])
	// θ2 = math.Atan(h72 / z[1])
	// θ3 = math.Atan(h73 / z[2])
	// θ4 = math.Atan(h74 / z[3])
	// θ5 = math.Atan(h75 / z[4])
	// θ6 = math.Atan(h76 / z[5])
	// if z[7] < 0 {
	// 	θ7 = math.Atan(z[7] / z[6])
	// 	return
	// }
	// θ7 = math.Pi + math.Atan(z[7]/z[6])
	return
}
