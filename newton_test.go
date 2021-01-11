package arc_test

import (
	"fmt"
	"math"
)

func Example(){	

	// # Define tolerance
	tol := 1.0e-12

	// # Input of user defined parameters
	th0 := math.Pi / 3.0

	// # Define the dimensions of 'load' and 'displacement' vectors
	const ndof int = 1

	// # Iq is the force distribution vector (needs to be defined explicitly)
	iq := [ndof]float64{}
	iq[0] = 1.0

	// # a is the dimensionless ``displacement'' vector (no need to define)
	a := [ndof]float64{}

	// # f is the system of equations (They need to be defined explicitly) [See Function fcn]
	f := [ndof]float64{}

	// # df is the tangent matrix to the system of equations (Contains derivatives)
	// # dfinv is the inverse of df
	// # df elements need to be defined explicitly in function dfcn
	df := [ndof][ndof]float64{}
	dfinv := [ndof][ndof]float64{}

	//
	// # dls will store the 2 solutions from the 2nd order polynomial w.r.t. ddl
	// dls := [2]float64{}

	// # dao is an araay that stores the last converged ``displacement correction''
	// dao := [ndof]float64{}
	// # al is the dimensionless ``load'' vector
	al := 0.0

	// # Define the b function needed for calculations
	b := func(x, y float64) float64 {
		return 1. + x*x - 2.0*x*math.Sin(y)
	}
	// # Define the system of equations in the form F(u)=0
	fcn := func(x [ndof]float64, y, z float64) [ndof]float64 {
		bb := b(x[0], y)
		f[0] = (1./math.Sqrt(bb)-1.0)*(math.Sin(y)-x[0]) - z
		return f
	}

	// # Define the tangent matrix (stiffness matrix)
	// # It contains the derivatives of the equations w.r.t the variables
	// # The function returns both the matrix as df as well as it's inverse as dfinv
	dfcn := func(x [ndof]float64, y, z float64) (_, _ [ndof][ndof]float64) {
		bb := b(x[0], y)
		// 	# Tangent Matrix
		df[0][0] = 1 - (1.-math.Pow(math.Sin(y), 2.0))/(math.Pow(bb, 1.5))
		// # Inverse of Tangent Matrix
		dfinv[0][0] = 1.0 / df[0][0] // TODO: np.linalg.inv(df)
		// TODO: divide by zero
		return df, dfinv
	}

	// # Define the maximum number of Riks increments
	newton := 10000
	maxiter := 100
	//
	summ := func(a, da *[ndof]float64) {
		for i := 0; i < ndof; i++ {
			a[i] += da[i]
		}
	}
	npdot := func(v1, v2 [ndof]float64) (res float64) {
		for i := 0; i < ndof; i++ {
			res += v1[i] * v2[i]
		}
		return
	}
	npdotm := func(m [ndof][ndof]float64, v [ndof]float64) (res [ndof]float64) {
		for i := 0; i < ndof; i++ {
			for k := 0; k < ndof; k++ {
				res[k] += m[k][i] * v[i]
			}
		}
		return
	}
	scale := func(f float64, v [ndof]float64) (res [ndof]float64) {
		for i := 0; i < ndof; i++ {
			res[i] = f * v[i]
		}
		return
	}

	for i := 0; i < newton; i++ {
		if a[0] >= 2.5 {
			break
		}
		// 	# Increment starts; Set all variations=0
		da := [ndof]float64{}
		dl := 1.5e-2
		summ(&a, &da)
		al += dl
		f = fcn(a, th0, al)
		fcheck := math.Sqrt(npdot(f, f))
		if fcheck < tol {
			fmt.Println(a, al)
		} else {
			iters := 0
			for fcheck > tol {
				iters += 1

				_, dfinv = dfcn(a, th0, al)
				da = scale(-1., npdotm(dfinv, f))
				summ(&a, &da)
				f = fcn(a, th0, al)
				fcheck = math.Sqrt(npdot(f, f))

				if iters > maxiter {
					panic(fmt.Errorf("iter : %d. fcheck: %e", iters, fcheck))
				}
			}
			fmt.Printf("%15.8f %15.8f\n", a, al)
			fmt.Println("Newton increment", i, "completed successfully and use" ,iters,"steps")
		}
	}
	fmt.Println("The program completed successfully")

	// Output:
}
