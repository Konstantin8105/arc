package arc_test

import (
	"fmt"
	"math"
)

func ExampleNewton2() {
	// # Input of user defined parameters
	const (
		th0 = math.Pi / 3.0
		dl  = 1.5e-2
	)

	var (
		// # a is the dimensionless ``displacement'' vector (no need to define)
		a float64

		// # lambda is the dimensionless ``load'' vector
		lambda float64

		// # Define the b function needed for calculations
		b = func(a float64) float64 {
			return 1.0 - 2.0*a*math.Sin(th0) + a*a
		}

		// # Define the system of equations in the form F(u)=0
		fcn = func(a, lambda float64) float64 {
			bb := b(a) //, θ)
			// # f is the system of equations (They need to be defined explicitly) [See Function fcn]
			f := (1.0/math.Sqrt(bb)-1.0)*(math.Sin(th0)-a) - lambda
			return f
		}

		// # Define the tangent matrix (stiffness matrix)
		// # It contains the derivatives of the equations w.r.t the variables
		// # The function returns both the matrix as df as well as it's inverse as dfinv
		//
		// # df is the tangent matrix to the system of equations (Contains derivatives)
		// # dfinv is the inverse of df
		// # df elements need to be defined explicitly in function dfcn
		dfcn = func(a float64) (_ float64) {
			bb := b(a) // , θ)
			//      # Tangent Matrix
			df := 1 - (1.-math.Pow(math.Sin(th0), 2.0))/(math.Pow(bb, 1.5))
			return df
			// # Inverse of Tangent Matrix
			// dfinv := 1 / df // nplinalginv(df)
			// return dfinv
		}
		norm2 = func(f float64) float64 {
			return math.Sqrt(f * f)
		}
		solve = func(K float64, x *float64, f float64) {
			*x = f / K
		}
	)

	// # Define the maximum number
	for newton := 1000; 0 <= newton; newton-- {
		if a >= 2.5 {
			break
		}
		lambda += dl
		f := fcn(a, lambda)
		da := 0.0
		fcheck := 1.0
		for fcheck > tol {
			if newton < 0 {
				panic(fmt.Errorf("fcheck: %e. lambda = %e", fcheck, lambda))
			}
			newton--
			solve(dfcn(a), &da, f)
			a -= da
			f = fcn(a, lambda)
			// fcheck = norm2(f)
			fcheck = norm2(da)
		}
		fmt.Printf("Newton: %.5f %.12f\n", lambda, a)
	}

	// Output:
	// Newton: 0.01500 0.020180174816
	// Newton: 0.03000 0.040751259053
	// Newton: 0.04500 0.061766145983
	// Newton: 0.06000 0.083289422764
	// Newton: 0.07500 0.105401244875
	// Newton: 0.09000 0.128203021799
	// Newton: 0.10500 0.151826066682
	// Newton: 0.12000 0.176445343109
	// Newton: 0.13500 0.202302531242
	// Newton: 0.15000 0.229747510832
	// Newton: 0.16500 0.259320165489
	// Newton: 0.18000 0.291933974171
	// Newton: 0.19500 0.329378371401
	// Newton: 0.21000 0.376288867399
	// Newton: 0.22500 0.474554314136
	// Newton: 0.24000 2.024111611297
	// Newton: 0.25500 2.041201134099
	// Newton: 0.27000 2.058204961515
	// Newton: 0.28500 2.075128186054
	// Newton: 0.30000 2.091975490037
	// Newton: 0.31500 2.108751187989
	// Newton: 0.33000 2.125459263653
	// Newton: 0.34500 2.142103402451
	// Newton: 0.36000 2.158687020045
	// Newton: 0.37500 2.175213287550
	// Newton: 0.39000 2.191685153860
	// Newton: 0.40500 2.208105365491
	// Newton: 0.42000 2.224476484252
	// Newton: 0.43500 2.240800903028
	// Newton: 0.45000 2.257080859921
	// Newton: 0.46500 2.273318450947
	// Newton: 0.48000 2.289515641454
	// Newton: 0.49500 2.305674276428
	// Newton: 0.51000 2.321796089811
	// Newton: 0.52500 2.337882712927
	// Newton: 0.54000 2.353935682146
	// Newton: 0.55500 2.369956445841
	// Newton: 0.57000 2.385946370727
	// Newton: 0.58500 2.401906747646
	// Newton: 0.60000 2.417838796857
	// Newton: 0.61500 2.433743672872
	// Newton: 0.63000 2.449622468895
	// Newton: 0.64500 2.465476220896
	// Newton: 0.66000 2.481305911358
	// Newton: 0.67500 2.497112472723
	// Newton: 0.69000 2.512896790576
}
