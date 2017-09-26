#ifndef __EIGENMULTIVARIATENORMAL_HPP
#define __EIGENMULTIVARIATENORMAL_HPP

#include <Eigen/Dense>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/normal_distribution.hpp>    

/*
We need a functor that can pretend it's const,
but to be a good random number generator
it needs mutable state.  The standard Eigen function
Random() just calls rand(), which changes a global
variable.
*/
namespace Eigen {
	namespace internal {
		template<typename Scalar>
		struct scalar_normal_dist_op
		{
			static boost::random::lagged_fibonacci607 rng;                        // The uniform pseudo-random algorithm
			mutable boost::normal_distribution<Scalar> norm;  // The gaussian combinator

			EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

			template<typename Index>
			inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
			inline void seed(const uint64_t &s) { rng.seed(s); }
		};

		template<typename Scalar>
		boost::random::lagged_fibonacci607 scalar_normal_dist_op<Scalar>::rng;

		template<typename Scalar>
		struct functor_traits<scalar_normal_dist_op<Scalar> >
		{
			enum { Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false };
		};

	} // end namespace internal
	/**
	Find the eigen-decomposition of the covariance matrix
	and then store it for sampling from a multi-variate normal
	*/
	template<typename Scalar>
	class EigenMultivariateNormal
	{
		Matrix<Scalar, Dynamic, Dynamic> _covar;
		Matrix<Scalar, Dynamic, Dynamic> _transform;
		Matrix< Scalar, Dynamic, 1> _mean;
		internal::scalar_normal_dist_op<Scalar> randN; // Gaussian functor


	public:
		EigenMultivariateNormal(const Matrix<Scalar, Dynamic, 1>& mean, const Matrix<Scalar, Dynamic, Dynamic>& covar, const uint64_t &seed = 1251)
		{
			randN.seed(seed);
			setMean(mean);
			setCovar(covar);
		}

		void setMean(const Matrix<Scalar, Dynamic, 1>& mean) { _mean = mean; }
		void setCovar(const Matrix<Scalar, Dynamic, Dynamic>& covar)
		{
			_covar = covar;

			// Assuming that we'll be using this repeatedly,
			// compute the transformation matrix that will
			// be applied to unit-variance independent normals

			
			Eigen::LLT<Eigen::Matrix<Scalar,Dynamic,Dynamic> > cholSolver(_covar);
			// We can only use the cholesky decomposition if
			// the covariance matrix is symmetric, pos-definite.
			// But a covariance matrix might be pos-semi-definite.
			// In that case, we'll go to an EigenSolver
			if (cholSolver.info()==Eigen::Success) {
			// Use cholesky solver
			_transform = cholSolver.matrixL();
			} else {
			SelfAdjointEigenSolver<Matrix<Scalar, Dynamic, Dynamic> > eigenSolver(_covar);
			_transform = eigenSolver.eigenvectors()*eigenSolver.eigenvalues().cwiseMax(0).cwiseSqrt().asDiagonal();
			}

		}

		/// Draw nn samples from the gaussian and return them
		/// as columns in a Dynamic by nn matrix
		Matrix<Scalar, Dynamic, -1> samples(int nn)
		{
			return (_transform * Matrix<Scalar, Dynamic, -1>::NullaryExpr(_covar.rows(), nn, randN)).colwise() + _mean;
		}
	}; // end class EigenMultivariateNormal
} // end namespace Eigen
#endif