////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	 	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////

#ifndef __TENSOR_CPP__
#define __TENSOR_CPP__

#include "Tensor.h"
#include "float.h"
#include <math.h>


bool Tensor::decomposeTensor(Matrix3f const &tensor, Vector3f &eval,Vector3f &emax,Vector3f &emid,Vector3f &emin)	{
    if (tensor.isZero(TC_EPSILON))  {
        return false;
    }

    SelfAdjointEigenSolver<Matrix3f> eigensolver;
    eigensolver.compute(tensor);

    emax = eigensolver.eigenvectors().col(2);
    emid = eigensolver.eigenvectors().col(1);
    emin = eigensolver.eigenvectors().col(0);

#ifdef DEBUG
    if (isnan(eigensolver.eigenvalues()[0]) || isnan(eigensolver.eigenvalues()[1]) || isnan(eigensolver.eigenvalues()[2])) {
        std::cout << " In tensor.cpp decomposeTensor" << std::endl;
        std::cout << eigensolver.eigenvalues() << std::endl;
        std::cout << tensor<< std::endl;
        exit(0);
    }
#endif
    ///The eigen values
    eval(T_EMAX) = eigensolver.eigenvalues()(2);
    eval(T_EMID) = eigensolver.eigenvalues()(1);
    eval(T_EMIN) = eigensolver.eigenvalues()(0);
    //if (eval.norm() > EPSILON) eval.normalize();
   // std::cout << eval.transpose() << std::endl;
	return true;
}

float tensorComparison(Tensor const &label_tensor, Tensor const &sample_tensor)	{

    Matrix3f A = label_tensor.getTensorInMatrixForm() ;
    Matrix3f B = sample_tensor.getTensorInMatrixForm();

    float a_norm = A.norm();
    float b_norm = B.norm();

    float dist = 1.0f - (A*B).trace()/(a_norm*b_norm);

 return std::max(dist,0.0f);
}


Tensor::Tensor(Vector3f const &_eval, Vector3f const &_emax, Vector3f const &_emid, Vector3f const &_emin)	{
	eval = _eval;
	evecs[T_EMAX]= _emax;
	evecs[T_EMID]= _emid;
	evecs[T_EMIN]= _emin;

	///Compute the differences
	eval_diffs[0] = (eval(0)-eval(1));
	eval_diffs[1] = (eval(1)-eval(2));
	eval_diffs[2] = (eval(2));

	///Save the geometric type
	if (eval_diffs[0] > eval_diffs[1])	{
		geometric_type = SURFACE;
	}
	else	{
		if (eval_diffs[1] > eval_diffs[2])	{
			geometric_type = CURVE;
		}
		else	{
			geometric_type = JOINT;
		}
	}

	///Compute the tensor in matrix form
	matrix_form = computeTensorMatrix(eval, evecs[0], evecs[1], evecs[2]);
    ///Compute the sqrt of the matrix
    //computeInverseOfSqrtMatrix();

}

Tensor::Tensor(Matrix3f const &_matrix_form)	{
    ///Initialize the matrix form
    matrix_form = _matrix_form;

	bool result = decomposeTensor(matrix_form,eval, evecs[T_EMAX], evecs[T_EMID], evecs[T_EMIN]);
    if (!result) geometric_type = UNKNOWN;
    else    {
        ///Compute the differences
        eval_diffs(0) = (eval(0)-eval(1));
        eval_diffs(1) = (eval(1)-eval(2));
        eval_diffs(2) = (eval(2));
        //std::cout << eval_diffs << std::endl;

        ///Save the geometric type
        if (eval_diffs(0) > eval_diffs(1))	{
            geometric_type = SURFACE;
        }
        else	{
            if (eval_diffs(1) > eval_diffs(2))	{
                geometric_type = CURVE;
            }
            else	{
                geometric_type = JOINT;
            }
        }
    }
}

Tensor::~Tensor()	{

}

Matrix3f Tensor::getTensorInMatrixForm() const	{
    return matrix_form;
}

Vector3f Tensor::getEigenVector(EIGENVECTOR_TYPE _type) const	{
	return evecs[_type];
}

Vector3f Tensor::getEigenValues() const {
    return eval;
}

float Tensor::getEigenValueDiffs(EIGENVALUE_DIFF_TYPE _type) const	{
	return eval_diffs(_type);
}

Tensor::TENSOR_GEOMETRIC_TYPE Tensor::getGeometricType() const	{
	return geometric_type;
}

Matrix3f Tensor::computeTensorMatrix(Vector3f const &eval,Vector3f const &e1,Vector3f const &e2,Vector3f const &e3)    {
	//compute the left-hand side matrix
	Matrix3f lhs;
	lhs(0,0)=e1(0);  lhs(0,1)=e2(0);  lhs(0,2)=e3(0);
	lhs(1,0)=e1(1);  lhs(1,1)=e2(1);  lhs(1,2)=e3(1);
	lhs(2,0)=e1(2);  lhs(2,1)=e2(2);  lhs(2,2)=e3(2);

	//compute the lambda matrix
	Matrix3f lm;
	lm(0,0)=eval(0);      lm(0,1)=0.0f;         lm(0,2)=0.0f;
	lm(1,0)=0.0f;         lm(1,1)=eval(1);      lm(1,2)=0.0f;
	lm(2,0)=0.0f;         lm(2,1)=0.0f;         lm(2,2)=eval(2);

	//compute the right-hand side matrix
	Matrix3f rhs;
	rhs(0,0)=e1(0);  rhs(0,1)=e1(1);  rhs(0,2)=e1(2);
	rhs(1,0)=e2(0);  rhs(1,1)=e2(1);  rhs(1,2)=e2(2);
	rhs(2,0)=e3(0);  rhs(2,1)=e3(1);  rhs(2,2)=e3(2);

	//compute the tensor product
	Matrix3f tensor_matrix = lhs * (lm * rhs);

	return tensor_matrix;
}

void Tensor::normalize()    {
    if (matrix_form.norm() > TENSOR_EPSILON)   matrix_form.normalize();
    ///Recompute the decomposition
    update();
    return;
}


void Tensor::update()   {
	///Recompute the decomposition
	bool result = decomposeTensor(matrix_form,eval, evecs[T_EMAX], evecs[T_EMID], evecs[T_EMIN]);
    if (!result) geometric_type = UNKNOWN;
    else    {
        ///Compute the differences
        eval_diffs(0) = (eval(0)-eval(1));
        eval_diffs(1) = (eval(1)-eval(2));
        eval_diffs(2) = (eval(2));
        //std::cout << eval_diffs << std::endl;

        ///Save the geometric type
        if (eval_diffs(0) > eval_diffs(1))	{
            geometric_type = SURFACE;
        }
        else	{
            if (eval_diffs(1) > eval_diffs(2))	{
                geometric_type = CURVE;
            }
            else	{
                geometric_type = JOINT;
            }
        }
    }
    return;
}

void Tensor::add(Tensor const &tensor)  {
    matrix_form += tensor.getTensorInMatrixForm();
    return;
}

void Tensor::scale(float scale_factor)  {
    matrix_form *= scale_factor;
    return;
}
#endif
