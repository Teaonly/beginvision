// Limited-memory BFGS (L-BFGS) algorithm implementation as described by Nocedal.
// L-BFGS is an unconstrained quasi-Newton optimization method that uses a limited memory variation
// of the Broyden–Fletcher–Goldfarb–Shanno (BFGS) update to approximate the inverse Hessian matrix.
// The implementation is robust as it uses a simple line-search technique (backtracking in one
// direction only) and still works even if the L-BFGS algorithm returns a non descent direction (as 
// it will then restart the optimization starting from the current solution).
// Its robustness enables it to minimize non-smooth functions, such as the hinge loss.
//
// Copyright (c) 2012 Idiap Research Institute, http://www.idiap.ch/
// Written by Charles Dubout <charles.dubout@idiap.ch>
//
// This file is part of LBFGS.
//
// LBFGS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 3 as
// published by the Free Software Foundation.
//
// LBFGS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with LBFGS. If not, see <http://www.gnu.org/licenses/>.

#ifndef _BV_LBFGS_H_
#define _BV_LBFGS_H_

#include <Eigen/Core>

#include <algorithm>
#include <cassert>
#include <iostream>

namespace bv{

/// Limited-memory BFGS (L-BFGS) algorithm implementation as described by Nocedal.
/// L-BFGS is an unconstrained quasi-Newton optimization method that uses a limited memory variation
/// of the Broyden–Fletcher–Goldfarb–Shanno (BFGS) update to approximate the inverse Hessian matrix.
/// The implementation is robust as it uses a simple line-search technique (backtracking in one
/// direction only) and still works even if the L-BFGS algorithm returns a non descent direction (as 
/// it will then restart the optimization starting from the current solution).
/// Its robustness enables it to minimize non-smooth functions, such as the hinge loss.
class MF_LBFGS
{
public:
	/// Callback interface to provide objective function and gradient evaluations.
	class IFunction
	{
	public:
		/// Destructor.
		virtual ~IFunction() {

		}
		
		/// Returns The number of variables.
		virtual int dim() const = 0;
		
		/// Provides objective function and gradient evaluations.
		/// @param[in] x Current solution.
		/// @param[out] g The gradient vector which must be computed for the current solution.
		/// @returns The value of the objective function for the current solution.
		virtual double operator()(const double * x, double * g = 0) const = 0;
		
		/// Provides information about the current iteration.
		/// @param[in] x The current solution.
		/// @param[in] g The gradient vector of the current solution.
		/// @param[in] n The number of variables.
		/// @param[in] fx The current value of the objective function.
		/// @param[in] xnorm The Euclidean norm of the current solution.
		/// @param[in] gnorm The Euclidean norm of the gradient vector.
		/// @param[in] step The line-search step used for this iteration.
		/// @param[in] t The iteration count.
		/// @param[in] ls The number of evaluations called for this iteration.
		virtual void progress(const double * x, const double * g, int n, double fx, double xnorm,
							  double gnorm, double step, int t, int ls) const {

		}
	};
	
public:
	/// Constructor.
	/// @param[in] function Callback function to provide objective function and gradient
	/// evaluations.
	/// @param[in] epsilon Accuracy to which the solution is to be found.
	/// @param[in] maxIterations Maximum number of iterations allowed.
	/// @param[in] int maxLineSearches Maximum number of line-searches per iteration allowed.
	/// @param[in] maxHistory Maximum history length of previous solutions and gradients.
	MF_LBFGS(const IFunction * function = 0, double epsilon = 1e-6, int maxIterations = 400,
		  int maxLineSearches = 20, int maxHistory = 10) {

		assert(!function || (function->dim() > 0));
		assert(epsilon > 0.0);
		assert(maxIterations > 0);
		assert(maxLineSearches > 0);
		assert(maxHistory >= 0);
	}
	
	/// Starts the L-BFGS optimization process.
	/// @param[in,out] x Initial solution on entry. Receives the optimization result on exit.
	/// @returns The final value of the objective function.
	double run(double * argx) const {
		// Define the types ourselves to make sure that the matrices are col-major
		typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor> VectorXd;
		typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixXd;
		
		assert(function_);
		assert(argx);
	
		// Convert the current solution to an Eigen::Map
		Eigen::Map<VectorXd> x(argx, function_->dim());
		
		// Initial value of the objective function and gradient
		VectorXd g(x.rows());
		double fx = (*function_)(x.data(), g.data());
		function_->progress(argx, g.data(), x.rows(), fx, x.norm(), g.norm(), 0.0, 0, 1);
		
		// Histories of the previous solution (required by L-BFGS)
		VectorXd px; // Previous solution x_{t-1}
		VectorXd pg; // Previous gradient g_{t-1}
		MatrixXd dxs(x.rows(), maxHistory_); // History of the previous dx's = x_{t-1} - x_{t-2}, ...
		MatrixXd dgs(x.rows(), maxHistory_); // History of the previous dg's = g_{t-1} - g_{t-2}, ...
		
		// Number of iterations remaining
		int nbIterations = maxIterations_;
		
		for (int i = 0; i < nbIterations; ++i) {
			// Relative tolerance
			const double relativeEpsilon = epsilon_ * std::max(1.0, x.norm());
			
			// Check the norm of the gradient against convergence threshold
			if (g.norm() < relativeEpsilon)
				return fx;
			
			// Get a new descent direction using the L-BFGS algorithm
			VectorXd z = g;
			
			if (i && maxHistory_) {
				// Update the histories
				const int h = std::min(i, maxHistory_); // Current length of the histories
				const int end = (i - 1) % h;
				
				dxs.col(end) = x - px;
				dgs.col(end) = g - pg;
				
				// Initialize the variables
				VectorXd p(h);
				VectorXd a(h);
				
				for (int j = 0; j < h; ++j) {
					const int k = (end - j + h) % h;
					p(k) = 1.0 / dxs.col(k).dot(dgs.col(k));
					a(k) = p(k) * dxs.col(k).dot(z);
					z -= a(k) * dgs.col(k);
				}
				
				// Scaling of initial Hessian (identity matrix)
				z *= dxs.col(end).dot(dgs.col(end)) / dgs.col(end).dot(dgs.col(end));
				
				for (int j = 0; j < h; ++j) {
					const int k = (end + j + 1) % h;
					const double b = p(k) * dgs.col(k).dot(z);
					z += dxs.col(k) * (a(k) - b);
				}
			}
			
			// Save the previous state
			px = x;
			pg = g;
			
			// If z is not a valid descent direction (because of a bad Hessian estimation), restart the
			// optimization starting from the current solution
			double descent = -z.dot(g);
			
			if (descent > -0.0001 * relativeEpsilon) {
				std::cout << "Reset LBFGS after " << i << "iterations." << std::endl;
				z = g;
				nbIterations -= i;
				i = 0;
				descent = -z.dot(g);
			}
			
			// Backtracking using Wolfe's first condition (Armijo condition)
			double step = i ? 1.0 : (1.0 / g.norm());
			
			bool down = false;
			
			int ls;
			
			for (ls = 0; ls < maxLineSearches_; ++ls) {
				// Tentative solution, gradient and loss
				const VectorXd nx = x - step * z;
				VectorXd ng(x.rows());
				const double nfx = (*function_)(nx.data(), ng.data());
				
				if (nfx <= fx + 0.0001 * step * descent) { // First Wolfe condition
					if ((-z.dot(ng) >= 0.9 * descent) || down) { // Second Wolfe condition
						x = nx;
						g = ng;
						fx = nfx;
						break;
					}
					else {
						step *= 2.0;
					}
				}
				else {
					step *= 0.5;
					down = true;
				}
			}
			
			if (ls == maxLineSearches_)
				return fx;
			
			function_->progress(argx, g.data(), x.rows(), fx, x.norm(), g.norm(), step, i + 1, ls + 1);
		}
		
		return fx;	
	}
	
private:
	// Constructor parameters
	const IFunction * function_;
	double epsilon_;
	int maxIterations_;
	int maxLineSearches_;
	int maxHistory_;
};

}
#endif // _BV_LBFGS_H_
