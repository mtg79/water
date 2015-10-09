#ifndef SHALLOW2D_H
#define SHALLOW2D_H

#include <vector>
#include <cmath>
#include <array>
//added vector also
//ldoc on
/**
 * # Shallow water equations
 * 
 * ## Physics picture
 * 
 * The shallow water equations treat water as incompressible and
 * inviscid, and assume that the horizontal velocity remains constant
 * in any vertical column of water.  The unknowns at each point are
 * the water height and the total horizontal momentum in a water
 * column; the equations describe conservation of mass (fluid is
 * neither created nor destroyed) and conservation of linear momentum.
 * We will solve these equations with a numerical method that also
 * exactly conserves mass and momentum (up to rounding error), though
 * it only approximately conserves energy.
 * 
 * The basic variables are water height ($h$), and the velocity components
 * ($u, v$).  We write the governing equations in the form
 * $$
 *   U_t = F(U)_x + G(U)_y
 * $$
 * where
 * $$
 *   U = \begin{bmatrix} h \\ hu \\ hv \end{bmatrix},
 *   F = \begin{bmatrix} hu \\ h^2 u + gh^2/2 \\ huv \end{bmatrix}
 *   G = \begin{bmatrix} hv \\ huv \\ h^2 v + gh^2/2 \end{bmatrix}
 * $$
 * The functions $F$ and $G$ are called *fluxes*, and describe how the
 * conserved quantities (volume and momentum) enter and exit a region
 * of space.
 * 
 * Note that we also need a bound on the characteristic wave speeds
 * for the problem in order to ensure that our method doesn't explode;
 * we use this to control the Courant-Friedrichs-Levy (CFL) number
 * relating wave speeds, time steps, and space steps.  For the shallow
 * water equations, the characteristic wave speed is $\sqrt{g h}$
 * where $g$ is the gravitational constant and $h$ is the height of the
 * water; in addition, we have to take into account the velocity of
 * the underlying flow.
 * 
 * ## Implementation
 * 
 * Our solver takes advantage of C++ templates to get (potentially)
 * good performance while keeping a clean abstraction between the
 * solver code and the details of the physics.  The `Shallow2D`
 * class specifies the precision of the comptutation (single precision),
 * the data type used to represent vectors of unknowns and fluxes
 * (the C++ `std::array`).  We are really only using the class as 
 * name space; we never create an instance of type `Shallow2D`,
 * and the `flux` and `wave_speed` functions needed by the solver are
 * declared as static (and inline, in the hopes of getting the compiler
 * to optimize for us).
 */

struct Shallow2D {

    // Type parameters for solver
    typedef float real;
	// a single vector,  which contains values for all points of the grid
	typedef std::vector<real> vec;	

	// We initialize a 3 x nxall x ny all vector (might have to chance to vector<vector,3>
	typedef std::array< vec ,3> tvec; //called tvec for three vector

    // Gravitational force (compile time constant)
    static constexpr real g = 9.8;
	
	//!!!!!!!!!!!!!!!!!!!!!!!! MAG have to somehow pass nx_all and ny_all to this  class
	
    // Compute shallow water fluxes F(U), G(U) 
	//MAG: now updates whole tvec vector
    static void flux(tvec& FU, tvec& GU, const tvec& U) {
		for (int index= 0 ; index < u.size(); ++index){  /// u.size() might not give what I want it to ? which is 3
			for (int iy = 1; iy < ny_all-1; ++iy){
				for (int ix = 1; ix < nx_all-1; ++ix) {
					real h = U[0][iy*nx_all+ix], hu = U[1][iy*nx_all+ix], hv = U[2][iy*nx_all+ix];
					FU[0][iy*nx_all+ix] = hu;
					FU[1][iy*nx_all+ix] = hu*hu/h + (0.5*g)*h*h;
					FU[2][iy*nx_all+ix] = hu*hv/h;

					GU[0][iy*nx_all+ix] = hv;
					GU[1][iy*nx_all+ix] = hu*hv/h;
					GU[2][iy*nx_all+ix] = hv*hv/h + (0.5*g)*h*h;

					}
				}
			}	
		}
    }

    // Compute shallow water wave speed
    static void wave_speed(real& cx, real& cy, const vtec& U) {
		
        using namespace std;
		real cell_cx, cell_cy;
		for (int index= 0 ; index < u.size(); ++index){  /// u.size() might not give what I want it to ? which is 3
			for (int iy = 0; iy < ny_all; ++iy){ // can precompute iy*nxall here (later)
				for (int ix = 0; ix < nx_all; ++ix) { 
					real h = U[0][iy*nx_all+ix], hu = U[1][iy*nx_all+ix], hv = U[2][iy*nx_all+ix];
					real root_gh = sqrt(g * h);  // NB: Don't let h go negative!
					cell_cx = abs(hu/h) + root_gh;
					cell_cy = abs(hv/h) + root_gh;										
					cx = max(cx, cell_cx);
					cy = max(cy, cell_cy);
				}
			}	
		}	
		
		
    }
};

//ldoc off
#endif /* SHALLOW2D_H */
