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
	// TODO: one simple initialization for all input templorary allocations for reusing

	// # Input of user defined parameters
	th0 := math.Pi / 3
	last_w := 0.25 // w = β/k , see page 24
	//dll := 2.5e-4

	// # Define the dimensions of 'load' and 'displacement' vectors
	var ndof int = 2

	// # Iq is the force distribution vector (needs to be defined explicitly)
	𝐪 := npzeros(ndof)
	// TODO : KI : FORCE on each direction
	// TODO: strange vector. What happen if ndof more 2
	𝐪[0] = 0.
	𝐪[1] = 1.
	// [0] - rotation
	// [1] - vertical displacement

	// # a is the dimensionless ``displacement'' vector (no need to define)
	// TODO: change to some dimention, not dimensionless
	// u := npzeros(ndof)

	// # df is the tangent matrix to the system of equations
	// (Contains derivatives)
	// # dfinv is the inverse of df
	// # df elements need to be defined explicitly in function dfcn
	// df := npzerosm(ndof)
	// dfinv := npzerosm(ndof)

	// # dao is an araay that stores the last converged ``displacement
	//  correction''
	// Δu := npzeros(ndof)

	// # Define the b function needed for calculations
	b := func(a float64) float64 {
		// TODO see formula (3.11)
		// B(a1, th0) = 1.0 - 2*a1*sin(th0)+a1*a1
		return 1.0 + a*a - 2.0*a*math.Sin(th0) // TODO: use formula
	}

	var dfcn func(a []float64) (df [][]float64)

	// # Define the system of equations in the form F(u)=0
	fcn := func(a []float64, λ float64) (f []float64) {

		// усилия в глобальной системе координат

		f = npzeros(ndof)

		// by formula (3.4):
		// a = u/Lo
		// λ = P/(2*k*Lo)
		// k = E` * Ao` / Lo //
		// β = E  * Ao  / Lo // vertical
		// w = β/k

		// # f is the system of equations (They need to be defined explicitly)
		// [See Function fcn]
		w := last_w
		bb := b(a[0])

		// TODO: use correct name of variables
		f[0] = (1./math.Sqrt(bb)-1.0)*(math.Sin(th0)-a[0]) - w*(a[1]-a[0])
		// TODO formula (3.12)
		f[1] = w*(a[1]-a[0]) - λ
		// TODO formula (3.13)

		//
		//
		// R := Fint + K * d - f = Fint + K * d - lambda * q ---> 0
		//
		// var Fint []float64
		// if len(data) != 0 {
		// 	Fint = data[len(data)-1].Fint
		// 	// last := data[len(data)-1]
		// 	// K := dfcn(last.a)
		// 	// Fint = w3
		// }
		// K := dfcn(a)
		// g := summa(npdotm(K, a), scale(-λ, 𝐪))
		// // for i := range data {
		// // 	K := dfcn(data[i].u)
		// // 	r := summa(npdotm(K, data[i].u), scale(-data[i].lambda, 𝐪))
		// // 	g = summa(g, scale(-1, r))
		// // }
		// g = summa(Fint, scale(-1, g))
		// _ = Fint
		// fmt.Println("> K = ", K)
		// fmt.Println("> λ = ", λ)
		// fmt.Println("> Fint = ", Fint)
		// fmt.Println(">>>", g, f)
		// f = g

		return
	}

	// # Define the tangent matrix (stiffness matrix)
	// # It contains the derivatives of the equations w.r.t the variables
	// # The function returns both the matrix as df
	// as well as it's inverse as dfinv
	dfcn = func(a []float64) (df [][]float64) {

		// Typical truss matrix stiffiners
		//
		//  A * E    | +c*c +c*s -c*c -c*s |
		//  -----  * | +c*s +s*s -c*s -s*s |
		//    L      | -c*c -c*s +c*c +c*s |
		//           | -c*s -s*s +c*s +s*s |
		//
		// angle = 90
		// Ao' * E'   | +c*c +c*s |   Ao  * E    | 0 0 |       | 0 0 |
		// -------- * | +c*s +s*s | = -------- * | 0 1 | = β * | 0 1 |
		//    Lo                         Lo
		//
		// angle = th0 = acos( Lo
		// Ao  * E    | +c*c +c*s |       | +c*c +c*s |
		// -------- * | +c*s +s*s | = k * | +c*s +s*s |
		//    Lo

		df = npzerosm(ndof)

		//dfinv = npzerosm(ndof)
		bb := b(a[0])
		w := last_w
		//      # Tangent Matrix
		// TODO: look like https://en.wikipedia.org/wiki/
		// Jacobian_matrix_and_determinant#Example_1
		df[0][0] = (1 + w) -
			(1.-math.Pow(math.Sin(th0), 2.0))/(math.Pow(bb, 1.5))
		df[0][1] = -w
		df[1][0] = -w
		df[1][1] = w
		// # Inverse of Tangent Matrix
		//dfinv = nplinalginv(df)
		// TODO : move to separate function for FEM calcs
		return df // , dfinv
	}

	data := arcm(dfcn, 𝐪)

	fmt.Printf("ok\n")

	// gnuplot graph
	// plot "data.txt" using 2:1 title "rotation", \
	//      "data.txt" using 3:1 title "vertical disp"
	var buf bytes.Buffer
	var errorValue float64
	for _, r := range data {
		fmt.Fprintf(&buf, "%.12f", r.lambda)
		for i := 0; i < ndof; i++ {
			fmt.Fprintf(&buf, " %.12f", r.u[i])
		}
		// print error
		f := fcn(r.u, r.lambda)
		for _, v := range f {
			fmt.Fprintf(&buf, " %.12f", v)
			errorValue = math.Max(errorValue, math.Abs(v))
		}

		fmt.Fprintf(&buf, "\n")
	}
	if err := os.WriteFile("data.txt", buf.Bytes(), 0644); err != nil {
		panic(err)
	}
	fmt.Printf("error value = %.1e\n", errorValue)

	// TODO : remove output data to specific file

	// Output:
	// ok
	// error value = 4.8e-04
}

type row struct {
	lambda float64
	u      []float64
}

// TODO : dfcn, 𝐪  dependens of u
// TODO : Uo - initialization deformation
func arcm(Kstiff func([]float64) [][]float64, 𝐪 []float64) (data []row) {

	ndof := len(𝐪)

	// # al is the dimensionless ``load'' vector
	λ := 0.0
	u := npzeros(ndof)
	𝜓 := 1.0
	// 	#Arc Length Parameters
	// var (
	//𝜓  = 1.0   // TODO: hyperellipsoid ratio - input data
	Δl := 1.e-3 // TODO : radius
	// )

	// # Define the maximum number of Riks increments
	// var ( // TODO: input data
	// 	// riks    = 20000
	// )
	// TODO : KI : names:
	// Lambda - load proportionality factor (LPF)

	stopStep := func(step int, λ float64, u []float64) bool {
		maxiter := 20000
		return maxiter < step || 2 < λ || 3.5 <= u[1]
	}

	// TODO : break a substep
	// TODO : break a calculation
	stopSubstep := func(substep int, fcheck float64) bool {
		maxiter := 100
		return maxiter < substep || fcheck < tol
	}

	for step := 0; ; step++ {
		if stopStep(step, λ, u) {
			break
		}

		// 	# Increment starts; Set all variations=0
		var (
			// TODO : minimaze allocations
			Δu = npzeros(ndof)

			// δů  []float64
			δu []float64
			// δut       []float64
			δu1, δu2 []float64
			// f   []float64 // TODO: remove

			// df  [][]float64
			det    float64
			Δλ     float64
			fcheck float64

			δλ       float64
			δλ1, δλ2 float64
		)

		stepa := func(isFirst bool) {
			Kt := Kstiff(summa(u, Δu))
			δut := SolveLinear(Kt, 𝐪)
			var δū []float64
			if isFirst {
				δū = npzeros(ndof)
			} else {
				// f = fcn(summa(u, Δu), (λ + Δλ)) //
				// temp := SolveLinear(Kt, f)      //
				// δů = scale(-1, temp)            //
				// fmt.Println(">", δů,
				// 	summa(SolveLinear(Kt, scale(Δλ,  𝐪)),scale(-1,Δu)))
				δū = summa(SolveLinear(Kt, scale(Δλ, 𝐪)), scale(-1, Δu))
			}

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
			// TODO : why if change ddl1 and ddl2 algorithm are fail??
			if D > 0. {
				// acceptable 2 solutions
				δλ1 = (-𝛼2 - math.Sqrt(D)) / (2 * 𝛼1)
				δλ2 = (-𝛼2 + math.Sqrt(D)) / (2 * 𝛼1)
			} else {
				panic(fmt.Errorf("not implemented: (%e,%e,%e) - %e",
					𝛼1, 𝛼2, 𝛼3, D))
				// δλ1 = -𝛼2 / 2 * 𝛼1
				// δλ2 = -𝛼2 / 2 * 𝛼1
				// // TODO : check coverage for that part of code
				// fmt.Println("Possible issue in Arc Length equation")
			}

			// Formula (2.14):
			// δu = δū + δλ*δut
			δu1 = summa(δū, scale(δλ1, δut))
			δu2 = summa(δū, scale(δλ2, δut))

			det = nplinalgdet(Kt)
		}
		stepa(true)

		// df = dfcn(summa(a, Δu))
		// δu = SolveLinear(df, 𝐪)
		// δλ1, δλ2 = square_root(Δu, δů, δu, Δλ, 𝐪)
		// TODO: why?? generally values are zeros
		// δu1 = summa(δů, scale(δλ1, δu))
		// δu2 = summa(δů, scale(δλ2, δu))
		// det = nplinalgdet(df)

		if npsign(det) == npsign(δλ1) {
			δu, δλ = δu1, δλ1
		} else {
			δu, δλ = δu2, δλ2
		}

		finish := func() {
			Δu = summa(Δu, δu)
			Δλ = Δλ + δλ
			// f = fcn(summa(u, Δu), (λ + Δλ))
			// fcheck = nplinalgnorm(f)
			fcheck = math.Max(nplinalgnorm(δu), math.Abs(δλ))
		}
		finish()

		// var iters int = 1 // TODO: in my point of view - it is 1

		// Run substeps
		for substep := 1; ; substep++ {
			if stopSubstep(substep, fcheck) {
				break
			}

			// ; fcheck > tol && iters <= maxiter; iters++ {

			stepa(false)

			// df = dfcn(summa(a, Δu))
			// δu = SolveLinear(df, 𝐪)
			// f = fcn(summa(a, Δu), (al + Δλ))
			// temp := SolveLinear(df, f)
			// δů = scale(-1, temp)
			// δλ1, δλ2 = square_root(Δu, δů, δu, Δλ, 𝐪)
			// Formula (2.14)
			// δu1 = summa(δů, scale(δλ1, δu))
			// δu2 = summa(δů, scale(δλ2, δu))
			// det = nplinalgdet(df)

			daomag := npdot(Δu, Δu)
			if daomag == 0. {
				if npsign(Δλ+δλ1) == npsign(det) {
					δu, δλ = δu1, δλ1
				} else {
					δu, δλ = δu2, δλ2
				}
			} else {
				// see page 14
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
			// Δu = summa(Δu, δu)
			// Δλ = Δλ + δλ
			// f = fcn(summa(a, Δu), (al + Δλ))
			// fcheck = nplinalgnorm(f)
		}

		// if iters > maxiter {
		// 	// TODO: create error description
		// 	panic("Max iteration error")
		// }

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
