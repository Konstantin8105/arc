package arc_test

import (
	"bytes"
	"fmt"
	"math"
	"os"
)

// # The arc length function solves the 2nd order equation w.r.t. ddl
// # Returns two values as ddl1 and ddl2
//
//	𝜓         is psi
//	𝐪         is iq
//	Δu        is da, dao
//	δū        is dab
//	δut       is dat
//	Δλ        is dl
//	δu        is dat, dda
//	δu1       is dda1
//	δu2       is dda2
//	δλ        is ddl
//	δλ1       is ddl1
//	δλ2       is ddl2
//	Δλ        is dl, dlamda
//	Δl        is dll
//	Δl1       is dll1
//	Δl2       is dll2
//	𝛼1        is c1
//	𝛼2        is c2
//	𝛼3        is c3

func SolveLinear(K [][]float64, f []float64) (d []float64) {
	ndof := len(K)
	dfinv := npzerosm(ndof)
	dfinv = nplinalginv(K) // TODO : move to separate function for FEM calcs
	d = npdotm(dfinv, f)
	return
}

func ExampleArc2() {
	// # Input of user defined parameters
	th0 := math.Pi / 3
	last_w := 0.25 // w = β/k , see page 24

	// # Define the dimensions of 'load' and 'displacement' vectors
	var ndof int = 2

	// # Iq is the force distribution vector (needs to be defined explicitly)
	𝐪 := npzeros(ndof)
	𝐪[0] = 0.
	𝐪[1] = 1.

	// # Define the b function needed for calculations
	b := func(a float64) float64 {
		// see formula (3.11)
		// B(a1, th0) = 1.0 - 2*a1*sin(th0)+a1*a1
		return 1.0 + a*a - 2.0*a*math.Sin(th0) // use formula
	}

	var dfcn func(a []float64) (df [][]float64)

	// # Define the system of equations in the form F(u)=0
	fcn := func(a []float64, λ float64) (f []float64) {
		// # f is the system of equations (They need to be defined explicitly)
		// [See Function fcn]
		f = npzeros(ndof)
		w := last_w
		bb := b(a[0])

		// use correct name of variables
		f[0] = (1./math.Sqrt(bb)-1.0)*(math.Sin(th0)-a[0]) - w*(a[1]-a[0])
		// formula (3.12)
		f[1] = w*(a[1]-a[0]) - λ
		// formula (3.13)

		return
	}

	// # Define the tangent matrix (stiffness matrix)
	// # It contains the derivatives of the equations w.r.t the variables
	// # The function returns both the matrix as df
	// as well as it's inverse as dfinv
	dfcn = func(a []float64) (df [][]float64) {
		df = npzerosm(ndof)

		bb := b(a[0])
		w := last_w
		// Tangent Matrix
		// look like https://en.wikipedia.org/wiki/
		// Jacobian_matrix_and_determinant#Example_1
		df[0][0] = (1 + w) -
			(1.-math.Pow(math.Sin(th0), 2.0))/(math.Pow(bb, 1.5))
		df[0][1] = -w
		df[1][0] = -w
		df[1][1] = w
		return df
	}

	// TODO : break a substep
	// TODO : break a calculation
	stopStep := func(step int, λ float64, u []float64) bool {
		maxiter := 20000
		return maxiter < step || 2 < λ || 3.5 <= u[1]
	}
	stopSubstep := func(substep int, fcheck float64) bool {
		maxiter := 100
		return maxiter < substep || fcheck < tol
	}

	c := DefaultConfig()
	// c.Radius = 2.0e-3
	data := arcm(dfcn, 𝐪, stopStep, stopSubstep, c)
	printData(data, "data.txt", 𝐪, dfcn, fcn)

	var errorValue float64
	for _, r := range data {
		// print error
		f := fcn(r.u, r.lambda)
		for _, v := range f {
			errorValue = math.Max(errorValue, math.Abs(v))
		}
	}
	fmt.Fprintf(os.Stdout, "error value = %.1e\n", errorValue)

	// Output:
	// error value = 4.8e-04
}

func printData(data []row, filename string,
	q []float64,
	Kt func(a []float64) (df [][]float64),
	F func(x []float64, lambda float64) []float64) {
	// gnuplot graph
	// plot "data.txt" using 2:1 title "rotation", \
	//      "data.txt" using 3:1 title "vertical disp"
	var buf bytes.Buffer
	yintegral := make([]float64, len(q))
	for index, r := range data {
		fmt.Fprintf(&buf, "%.12f", r.lambda)
		for i := range r.u {
			fmt.Fprintf(&buf, " %.12f", r.u[i])
		}
		if F != nil {
			ys := F(r.u, 0) // r.lambda)
			if 0 < index {
				yintegral = summa(
					yintegral,
					npdotm(Kt(r.u), summa(data[index].u, scale(-1, data[index-1].u))))
			}
			for i := range ys {
				fmt.Fprintf(&buf, " %.12f %.12f", yintegral[i], ys[i])
			}
		}
		fmt.Fprintf(&buf, "\n")
	}
	if err := os.WriteFile(filename, buf.Bytes(), 0644); err != nil {
		panic(err)
	}
}

func ExampleArc3() {
	stopStep := func(step int, λ float64, u []float64) bool {
		maxiter := 20000
		return maxiter < step || 2 < λ || 20 <= u[0]
	}
	stopSubstep := func(substep int, fcheck float64) bool {
		maxiter := 100
		return maxiter < substep || fcheck < tol
	}

	q := []float64{120}
	F := func(x []float64, lambda float64) []float64 {
		return []float64{
			-0.06*math.Pow(x[0], 3) + 1.2*math.Pow(x[0], 2) + 3*x[0] - lambda*q[0],
		}
	}
	K := func(u []float64) [][]float64 {
		k := [][]float64{
			{-0.06*3*math.Pow(u[0], 2) + 1.2*2*u[0] + 3},
		}
		return k
	}

	c := DefaultConfig()
	c.Radius = 2.0e-0
	data := arcm(K, q, stopStep, stopSubstep, c)

	var errorValue float64
	for _, r := range data {
		// print error
		f := F(r.u, r.lambda)
		for _, v := range f {
			errorValue = math.Max(errorValue, math.Abs(v))
		}
	}
	fmt.Fprintf(os.Stdout, "error value = %.1e\n", errorValue)

	for i := range data {
		data[i].lambda *= q[0]
	}

	printData(data, "arc3.txt", q, K, F)
	// Output:
	// error value = 5.1e+00
}

type row struct {
	lambda float64
	u      []float64
}

type Config struct {
	// TODO: hyperellipsoid ratio - input data
	Ksi float64

	// TODO : radius
	Radius float64
}

func DefaultConfig() *Config {
	c := Config{
		Ksi:    1.0,
		Radius: 1e-3,
	}
	return &c
}

// TODO : dfcn, 𝐪  dependens of u
// TODO : Uo - initialization deformation
func arcm(Kstiff func([]float64) [][]float64, 𝐪 []float64,
	stopStep func(step int, λ float64, u []float64) bool,
	stopSubstep func(substep int, fcheck float64) bool,
	c *Config,
) (data []row) {
	// TODO : add error handling

	ndof := len(𝐪)
	u := npzeros(ndof)

	// Arc Length Parameters

	// Lambda - load proportionality factor (LPF)
	// the dimensionless ``load'' vector
	var λ float64
	var 𝜓, Δl float64
	if c == nil {
		c = DefaultConfig()
	}
	𝜓, Δl = c.Ksi, c.Radius
	if 𝜓 <= 0 {
		panic("1")
	}
	if Δl <= 0 {
		panic("2")
	}

	data = append(data, row{
		lambda: λ,
		u:      u,
	})

	for step := 0; ; step++ {
		if stopStep(step, λ, u) {
			break
		}

		var (
			Δu = npzeros(ndof)

			δu       []float64
			δu1, δu2 []float64

			det    float64
			Δλ     float64
			fcheck float64

			δλ       float64
			δλ1, δλ2 float64
		)

		begin := func(isFirst bool) {
			Kt := Kstiff(summa(u, Δu))
			// For formula (2.14):
			// δū = -invert[KT](uo+Δu) * (Fint*(uo+Δu)-(λo+Δλ)*𝐪)
			// value of Fint is precision value and for we cannot find them
			// only by Jacobi matrix.
			// Fint*(uo+Δu)-(λo+Δλ)*𝐪 is equal R(uo+Δu), but
			// theoretically R(uo+Δu) = 0, then vector is zero always:
			// TODO : find the solution
			δū := npzeros(ndof)
			δut := SolveLinear(Kt, 𝐪)

			// Formula (2.12):
			// (∆u + δu)T*(∆u + δu) + ψ^2*(∆λ + δλ)^2*(𝐪T * 𝐪) = ∆l^2
			//
			// Formula (2.14)
			// δu = δū + δλ*δut
			//
			// Formula (2.15)
			// 𝛼1*δλ^2 + 𝛼2*δλ + 𝛼3 = 0
			//
			// symbolic math:
			// pow(deltau+(δu_+δλ*δut),2)+ψ2*pow(deltaλ+δλ,2)*(q2)-l2
			//
			// deltau*deltau + 2*deltau*δu_ + δu_*δu_ + 2*deltau*δut*δλ + \
			// ::::::::::::::::::::::::::::::::::::::   ---------------   \
			// 2*δu_*δut*δλ + δut*δut*δλ*δλ + deltaλ*deltaλ*q2*ψ2 +       \
			// -------------  =============   :::::::::::::::::::         \
			// 2*deltaλ*q2*δλ*ψ2 + q2*δλ*δλ*ψ2 - l2                       \
			// -----------------   ============::::                       \
			//
			// 𝛼1 = δutT*δut + ψ^2*(𝐪T * 𝐪)
			// 𝛼2 = 2*(∆u+δū)*δut+2*ψ^2*∆λ*(𝐪T * 𝐪)
			// 𝛼3 = (∆u + δū)T*(∆u + δū)+ψ^2*∆λ^2*(𝐪T * 𝐪)-∆l^2
			//
			var (
				// calculate the coefficients of the polynomial
				𝛼1 = npdot(δut, δut) +
					math.Pow(𝜓, 2.0)*npdot(𝐪, 𝐪)
				𝛼2 = 2.0*npdot(summa(Δu, δū), δut) +
					2.0*Δλ*math.Pow(𝜓, 2.0)*npdot(𝐪, 𝐪)
				𝛼3 = npdot(summa(Δu, δū), summa(Δu, δū)) +
					math.Pow(𝜓, 2.0)*math.Pow(Δλ, 2.0)*npdot(𝐪, 𝐪) -
					math.Pow(Δl, 2.0)

				// determinant
				D = 𝛼2*𝛼2 - 4.*𝛼1*𝛼3
			)
			if 0.0 < D {
				// acceptable 2 solutions
				δλ1 = (-𝛼2 - math.Sqrt(D)) / (2.0 * 𝛼1)
				δλ2 = (-𝛼2 + math.Sqrt(D)) / (2.0 * 𝛼1)
			} else {
				δλ1 = -𝛼2 / (2.0 * 𝛼1)
				δλ2 = -𝛼2 / (2.0 * 𝛼1)
				panic((fmt.Errorf("not implemented: (%e,%e,%e) - %e",
					𝛼1, 𝛼2, 𝛼3, D)))
				// TODO : check coverage for that part of code : D < 0.0
			}
			// After checking - acceptable swap the results, but no need
			// δλ2, δλ1 = δλ1, δλ2

			// Formula (2.14):
			// δu = δū + δλ*δut
			//
			δu1 = summa(δū, scale(δλ1, δut))
			δu2 = summa(δū, scale(δλ2, δut))

			// calculate determinant matrix of stiffiners
			//
			det = nplinalgdet(Kt)
		}
		begin(true)

		if npsign(det) == npsign(δλ1) {
			δu, δλ = δu1, δλ1
		} else {
			δu, δλ = δu2, δλ2
		}

		finish := func() {
			Δu = summa(Δu, δu)
			Δλ = Δλ + δλ
			fcheck = math.Max(nplinalgnorm(δu), math.Abs(δλ))
		}
		finish()

		// Run substeps
		for substep := 1; ; substep++ {
			if stopSubstep(substep, fcheck) {
				break
			}

			begin(false)

			daomag := npdot(Δu, Δu)
			if daomag == 0.0 {
				if npsign(Δλ+δλ1) == npsign(det) {
					δu, δλ = δu1, δλ1
				} else {
					δu, δλ = δu2, δλ2
				}
			} else {
				// Formula (2.16):
				//
				DOT1 := npdot(summa(Δu, δu1), Δu) +
					math.Pow(𝜓, 2)*Δλ*(Δλ+δλ1)*npdot(𝐪, 𝐪)
				DOT2 := npdot(summa(Δu, δu2), Δu) +
					math.Pow(𝜓, 2)*Δλ*(Δλ+δλ2)*npdot(𝐪, 𝐪)

				if DOT1 > DOT2 {
					δu, δλ = δu1, δλ1
				} else {
					δu, δλ = δu2, δλ2
				}
			}
			if δλ1 == δλ2 {
				δu, δλ = δu1, δλ1
			}

			finish()
		}

		u = summa(u, Δu)
		λ += Δλ

		// TODO: add visualization for steps and substeps
		// TODO: add recorder for each step

		data = append(data, row{
			lambda: λ,
			u:      u,
		})
	}

	return
}
