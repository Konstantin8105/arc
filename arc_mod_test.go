package arc_test

import (
	"fmt"
	"math"
)

// # The arc length function solves the 2nd order equation w.r.t. ddl
// # Returns two values as ddl1 and ddl2
//
// Formula (2.15)
//
//	𝜓         is psi
//	𝐪         is iq
//	Δu        is da, dao
//	δu~ (δů)  is dab
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
func square_root(Δu, δů, δu []float64, Δλ float64, 𝐪 []float64) (
	δλ1, δλ2 float64) {

	// 	#Arc Length Parameters
	var (
		𝜓  = 1.0   // TODO: hyperellipsoid ratio - input data
		Δl = 1.e-3 // TODO : radius
	)

	// TODO: add comments for each variable
	// TODO: rename in according to arc documentation

	// 	# Calculate the coefficients of the polynomial
	var (
		𝛼1 = npdot(δu, δu) +
			math.Pow(𝜓, 2.0)*npdot(𝐪, 𝐪)
		𝛼2 = 2. * (npdot(summa(Δu, δů), δu) +
			Δλ*math.Pow(𝜓, 2)*npdot(𝐪, 𝐪))
		𝛼3 = npdot(summa(Δu, δů), summa(Δu, δů)) +
			math.Pow(Δλ, 2.0)*math.Pow(𝜓, 2.0)*npdot(𝐪, 𝐪) -
			math.Pow(Δl, 2)
	)

	// TODO : why if change ddl1 and ddl2 algorithm are fail??
	if 𝛼2*𝛼2-4.*𝛼1*𝛼3 > 0. { // TODO: this is determinant
		// # dls will store the 2 solutions from the 2nd order polynomial w.r.t. ddl
		dls := nproots(𝛼1, 𝛼2, 𝛼3)
		δλ1 = dls[0]
		δλ2 = dls[1]
	} else {
		δλ1 = -𝛼2 / 2 * 𝛼1
		δλ2 = -𝛼2 / 2 * 𝛼1
		// TODO : check coverage for that part of code
		fmt.Println("Possible issue in Arc Length equation")
	}

	return //  ddl1,ddl2
}

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
	𝜓 := 1.0
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
	a := npzeros(ndof)

	// # df is the tangent matrix to the system of equations
	// (Contains derivatives)
	// # dfinv is the inverse of df
	// # df elements need to be defined explicitly in function dfcn
	// df := npzerosm(ndof)
	// dfinv := npzerosm(ndof)

	// # dao is an araay that stores the last converged ``displacement
	//  correction''
	// Δu := npzeros(ndof)

	// # al is the dimensionless ``load'' vector
	al := 0.0

	// # Define the b function needed for calculations
	b := func(a1 float64) float64 {
		// TODO see formula (3.11)
		// B(a1, th0) = 1.0 - 2*a1*sin(th0)+a1*a1
		return 1.0 + a1*a1 - 2.0*a1*math.Sin(th0) // TODO: use formula
	}

	// # Define the system of equations in the form F(u)=0
	fcn := func(a []float64, λ float64) (f []float64) {
		// # f is the system of equations (They need to be defined explicitly)
		// [See Function fcn]
		w := last_w
		f = npzeros(ndof)
		bb := b(a[0])

		// by formula (3.4):
		// a = u/Lo
		// λ = P/(2*k*Lo)

		// TODO: use correct name of variables
		f[0] = (1./math.Sqrt(bb)-1.0)*(math.Sin(th0)-a[0]) - w*(a[1]-a[0])
		// TODO formula (3.12)
		f[1] = w*(a[1]-a[0]) - λ
		// TODO formula (3.13)
		return
	}

	// # Define the tangent matrix (stiffness matrix)
	// # It contains the derivatives of the equations w.r.t the variables
	// # The function returns both the matrix as df
	// as well as it's inverse as dfinv
	dfcn := func(x []float64) (df [][]float64) {
		// (df, dfinv [][]float64)
		df = npzerosm(ndof)
		//dfinv = npzerosm(ndof)
		bb := b(x[0])
		y := th0
		w := last_w
		//      # Tangent Matrix
		// TODO: look like https://en.wikipedia.org/wiki/
		// Jacobian_matrix_and_determinant#Example_1
		df[0][0] = (1 + w) -
			(1.-math.Pow(math.Sin(y), 2.0))/(math.Pow(bb, 1.5))
		df[0][1] = -w
		df[1][0] = -w
		df[1][1] = w
		// # Inverse of Tangent Matrix
		//dfinv = nplinalginv(df)
		// TODO : move to separate function for FEM calcs
		return df // , dfinv
	}

	// # Define the maximum number of Riks increments
	var ( // TODO: input data
		// riks    = 20000
		maxiter = 100
	)

	// TODO : KI : names:
	// Lambda - load proportionality factor (LPF)

	for { // i := 0; i < riks; i++ {
		// TODO: add stop factors
		if a[1] >= 3.5 {
			break
		}
		// 	# Increment starts; Set all variations=0
		var (
			// TODO : minimaze allocations
			Δu = npzeros(ndof)

			// δů  []float64
			δu  []float64
			δu1 []float64
			δu2 []float64
			f   []float64

			// df  [][]float64
			det    float64
			Δλ     float64
			fcheck float64

			δλ  float64
			δλ1 float64
			δλ2 float64
		)

		step := func(isFirst bool) {
			Kt := dfcn(summa(a, Δu))
			δu = SolveLinear(Kt, 𝐪)
			var δů []float64
			if isFirst {
				δů = npzeros(ndof)
			} else {
				f = fcn(summa(a, Δu), (al + Δλ)) //
				temp := SolveLinear(Kt, f)       //
				δů = scale(-1, temp)             //
			}
			δλ1, δλ2 = square_root(Δu, δů, δu, Δλ, 𝐪)
			// Formula (2.14)
			δu1 = summa(δů, scale(δλ1, δu))
			δu2 = summa(δů, scale(δλ2, δu))
			det = nplinalgdet(Kt)
		}
		step(true)

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
			f = fcn(summa(a, Δu), (al + Δλ))
			fcheck = nplinalgnorm(f)
		}
		finish()

		var iters int = 1 // TODO: in my point of view - it is 1
		for ; fcheck > tol && iters <= maxiter; iters++ {

			step(false)

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

		if iters > maxiter {
			// TODO: create error description
			panic("Max iteration error")
		}

		a = summa(a, Δu)
		al += Δλ

		// TODO: add visualization for steps and substeps
		// TODO: add recorder for each step
		// fmt.Printf("%.12f", al)
		// for i := 0; i < ndof; i++ {
		// 	fmt.Printf(" %.12f", a[i])
		// }
		// fmt.Printf("\n")
	}
	fmt.Printf("ok\n")

	// TODO : remove output data to specific file

	// Output:
	// ok
}
