////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	 	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////

#ifndef __TENSOR_H__
#define __TENSOR_H__

#include "Utilities.h"
#include "Algebra.h"
#include <Eigen/Eigen>

using namespace Eigen;

#define XX 0
#define YY 1
#define ZZ 2

#define TENSOR_EPSILON  0.0001f

class Tensor	{
	public:
		///Constructor
		Tensor ()	{matrix_form = Matrix3f::Zero();}
		Tensor(Vector3f const &_eval, Vector3f const &_emax, Vector3f const &_emid, Vector3f const &_emin);
		Tensor(Matrix3f const &_matrix_form);
		///Destructor
		~Tensor();

		///The eigenvector types
		///NOTE: The indices are inverted because Eigen sorts the
		///eigenvalues/vectors in increasing order
		typedef enum EIGENVECTOR_TYPE	{
			T_EMAX = 0,
			T_EMID = 1,
			T_EMIN = 2
		} EIGENVECTOR_TYPE;

		typedef enum EIGENVALUE_DIFF_TYPE	{
			L1MINUSL2 = 0,
			L2MINUSL3 = 1,
			L3 = 2
		} EIGENVALUE_DIFF_TYPE;

		typedef enum TENSOR_GEOMETRIC_TYPE	{
			UNKNOWN = -1,
			SURFACE = 0,
			CURVE = 1,
			JOINT = 2
		} TENSOR_GEOMETRIC_TYPE;

		///Returns the requested eigenvector
		Vector3f getEigenVector(EIGENVECTOR_TYPE _type) const;

		///Returns the eigenvalues
		Vector3f getEigenValues() const;

		///Returns the eigenvalue difference
		float getEigenValueDiffs(EIGENVALUE_DIFF_TYPE _type) const ;

		///Returns the geometric type of the tensor
		Tensor::TENSOR_GEOMETRIC_TYPE getGeometricType() const;

		///Returns the tensor in matrix form
		Matrix3f getTensorInMatrixForm() const;

		///Given a tensor it decomposes it into its eigenvalues and eigenvectors
		static bool decomposeTensor(Matrix3f const &tensor, Vector3f &eval, Vector3f &emax, Vector3f &emid, Vector3f &emin);

		///Compare two tensors and return the result
		friend float tensorComparison(Tensor const &tensor1, Tensor const &tensor2);

		///Computes the matrix form of a matrix
        static Matrix3f computeTensorMatrix(Vector3f const &eval, Vector3f const &e1, Vector3f const &e2, Vector3f const &e3);

        ///Normalizes the matrix corresponding to this tensor
        void normalize();

        void update();
        void add(Tensor const &tensor);
        void scale(float scale_factor);
	private:
		///The eigenvalues
		Vector3f eval;
		///The eigenvectors
		Vector3f evecs[3];
		///The eigenvalue differences
		Vector3f eval_diffs;
		///The geometric type of the tensor
		TENSOR_GEOMETRIC_TYPE geometric_type;
		///The tensor in matrix form
		Eigen::Matrix3f matrix_form;
		///The inverse matrix
		Eigen::Matrix3f inv_sqrt_matrix_form;
};


#endif
