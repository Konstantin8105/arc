package arc_test

import (
	"fmt"
	"math"
)

func summ(a, da []float64) {
	size := len(a)
	for i := 0; i < size; i++ {
		a[i] += da[i]
	}
}

func npdot(v1, v2 []float64) (res float64) {
	ndof := len(v1)
	for i := 0; i < ndof; i++ {
		res += v1[i] * v2[i]
	}
	return
}

func npzeros(ndof int) []float64 {
	return make([]float64, ndof)
}

func npzerosm(ndof int) [][]float64 {
	res := make([][]float64, ndof)
	for i := 0; i < ndof; i++ {
		res[i] = npzeros(ndof)
	}
	return res
}

func npdotm(m [][]float64, v []float64) (res []float64) {
	ndof := len(m)
	res = npzeros(ndof)
	for i := 0; i < ndof; i++ {
		for k := 0; k < ndof; k++ {
			res[k] += m[k][i] * v[i]
		}
	}
	return
}

func scale(f float64, v []float64) (res []float64) {
	ndof := len(v)
	res = npzeros(ndof)
	for i := 0; i < ndof; i++ {
		res[i] = f * v[i]
	}
	return
}

func ExampleNewton() {

	// # Define tolerance
	tol := 1.0e-12

	// # Input of user defined parameters
	th0 := math.Pi / 3.0

	// # Define the dimensions of 'load' and 'displacement' vectors
	var ndof int = 1

	// # Iq is the force distribution vector (needs to be defined explicitly)
	iq := npzeros(ndof)
	iq[0] = 1.0

	// # a is the dimensionless ``displacement'' vector (no need to define)
	a := npzeros(ndof)

	// # f is the system of equations (They need to be defined explicitly) [See Function fcn]
	f := npzeros(ndof)

	// # df is the tangent matrix to the system of equations (Contains derivatives)
	// # dfinv is the inverse of df
	// # df elements need to be defined explicitly in function dfcn
	df := npzerosm(ndof)
	dfinv := npzerosm(ndof)

	// # al is the dimensionless ``load'' vector
	al := 0.0

	// # Define the b function needed for calculations
	b := func(x, y float64) float64 {
		return 1. + x*x - 2.0*x*math.Sin(y)
	}
	// # Define the system of equations in the form F(u)=0
	fcn := func(x []float64, y, z float64) []float64 {
		bb := b(x[0], y)
		f[0] = (1./math.Sqrt(bb)-1.0)*(math.Sin(y)-x[0]) - z
		return f
	}

	// # Define the tangent matrix (stiffness matrix)
	// # It contains the derivatives of the equations w.r.t the variables
	// # The function returns both the matrix as df as well as it's inverse as dfinv

	dfcn := func(x []float64, y, z float64) (_ [][]float64) {
		bb := b(x[0], y)
		//      # Tangent Matrix
		df[0][0] = 1 - (1.-math.Pow(math.Sin(y), 2.0))/(math.Pow(bb, 1.5))
		// # Inverse of Tangent Matrix
		dfinv[0][0] = 1.0 / df[0][0] // TODO: np.linalg.inv(df)
		// TODO: divide by zero
		return dfinv
	}


       // # Define the maximum number of Riks increments
       newton := 10000
       maxiter := 100

	for i := 0; i < newton; i++ {
		if a[0] >= 2.5 {
			break
		}
		//      # Increment starts; Set all variations=0
		da := npzeros(ndof)
		dl := 1.5e-2
		// TODO strange summ(&a, &da)
		al += dl
		f = fcn(a, th0, al)
		fcheck := math.Sqrt(npdot(f, f))

		iters := 0
		for fcheck > tol {
			iters += 1

			dfinv = dfcn(a, th0, al)
			da = scale(-1., npdotm(dfinv, f))
			summ(a, da)
			f = fcn(a, th0, al)
			fcheck = math.Sqrt(npdot(f, f))

			if iters > maxiter {
				panic(fmt.Errorf("iter : %d. fcheck: %e", iters, fcheck))
			}
		}
		fmt.Printf(
			"Newton %3d(steps %3d). Coordinate: %.15f",
			i, iters, al,
		)
		for i := 0; i < ndof; i++ {
			fmt.Printf(" %.12f", a[i])
		}
		fmt.Printf("\n")
	}

	// Output:
	// Newton   0(steps   3). Coordinate: 0.015000000000000 0.020180174816
	// Newton   1(steps   3). Coordinate: 0.030000000000000 0.040751259053
	// Newton   2(steps   3). Coordinate: 0.045000000000000 0.061766145983
	// Newton   3(steps   3). Coordinate: 0.060000000000000 0.083289422764
	// Newton   4(steps   3). Coordinate: 0.075000000000000 0.105401244875
	// Newton   5(steps   3). Coordinate: 0.090000000000000 0.128203021799
	// Newton   6(steps   3). Coordinate: 0.105000000000000 0.151826066682
	// Newton   7(steps   3). Coordinate: 0.120000000000000 0.176445343109
	// Newton   8(steps   3). Coordinate: 0.135000000000000 0.202302531242
	// Newton   9(steps   3). Coordinate: 0.150000000000000 0.229747510832
	// Newton  10(steps   4). Coordinate: 0.165000000000000 0.259320165489
	// Newton  11(steps   4). Coordinate: 0.180000000000000 0.291933974171
	// Newton  12(steps   4). Coordinate: 0.195000000000000 0.329378371401
	// Newton  13(steps   4). Coordinate: 0.210000000000000 0.376288867399
	// Newton  14(steps   7). Coordinate: 0.225000000000000 0.474554314136
	// Newton  15(steps  29). Coordinate: 0.240000000000000 2.024111611297
	// Newton  16(steps   3). Coordinate: 0.255000000000000 2.041201134099
	// Newton  17(steps   3). Coordinate: 0.270000000000000 2.058204961515
	// Newton  18(steps   3). Coordinate: 0.285000000000000 2.075128186054
	// Newton  19(steps   3). Coordinate: 0.300000000000000 2.091975490037
	// Newton  20(steps   3). Coordinate: 0.315000000000000 2.108751187989
	// Newton  21(steps   3). Coordinate: 0.330000000000000 2.125459263653
	// Newton  22(steps   3). Coordinate: 0.345000000000000 2.142103402451
	// Newton  23(steps   3). Coordinate: 0.360000000000000 2.158687020045
	// Newton  24(steps   3). Coordinate: 0.375000000000000 2.175213287550
	// Newton  25(steps   3). Coordinate: 0.390000000000000 2.191685153860
	// Newton  26(steps   3). Coordinate: 0.405000000000000 2.208105365491
	// Newton  27(steps   3). Coordinate: 0.420000000000000 2.224476484252
	// Newton  28(steps   3). Coordinate: 0.435000000000000 2.240800903028
	// Newton  29(steps   3). Coordinate: 0.450000000000000 2.257080859921
	// Newton  30(steps   3). Coordinate: 0.465000000000000 2.273318450947
	// Newton  31(steps   3). Coordinate: 0.480000000000000 2.289515641454
	// Newton  32(steps   3). Coordinate: 0.495000000000000 2.305674276428
	// Newton  33(steps   3). Coordinate: 0.510000000000000 2.321796089811
	// Newton  34(steps   3). Coordinate: 0.525000000000000 2.337882712927
	// Newton  35(steps   3). Coordinate: 0.540000000000000 2.353935682146
	// Newton  36(steps   3). Coordinate: 0.555000000000000 2.369956445841
	// Newton  37(steps   3). Coordinate: 0.570000000000000 2.385946370727
	// Newton  38(steps   3). Coordinate: 0.585000000000000 2.401906747646
	// Newton  39(steps   3). Coordinate: 0.600000000000000 2.417838796857
	// Newton  40(steps   3). Coordinate: 0.615000000000000 2.433743672872
	// Newton  41(steps   3). Coordinate: 0.630000000000000 2.449622468895
	// Newton  42(steps   3). Coordinate: 0.645000000000000 2.465476220896
	// Newton  43(steps   3). Coordinate: 0.660000000000000 2.481305911358
	// Newton  44(steps   3). Coordinate: 0.675000000000000 2.497112472723
	// Newton  45(steps   3). Coordinate: 0.690000000000001 2.512896790576
}
