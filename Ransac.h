////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////

#ifndef __RANSAC_H__
#define __RANSAC_H__

#include <vector>

#include "RNG.h"
#include "Algebra.h"

#include <Eigen/Eigen>
using namespace Eigen;

using namespace std;


// x         - Data sets to which we are seeking to fit a model M

//fitting_function - Handle to a function that fits a model to s
//                   data from x.  It is assumed that the function is of the
//                   form:
//                   M = fitting_function(x)
//                   M is a vector with the plane coefficients

//evaluate_fit_function - Takes the plane coefficients estimated and the set of points. Finds all the inliers
//                        which are less than a specified threshold

//is_degenerate_function - Checks if a set of points are a degenerate case

//minimum_sampled_for_fitting - Minimum number of samples needed to fit a model

//error_threshold - Specifies the threshold distance of the points from the plane

//output: best_plane_coeffs - coefficients with maximum number of inliers
//        inlier_indices - vector with indices of the inliers

static Vector4f plane_fitting_function( std::vector<Vector3f> const &points, Vector3f const &plane_point)  {

		float sum_xx, sum_xy, sum_xz, sum_yy, sum_yz, sum_zz;

		sum_xx = sum_xy = sum_xz = 0.0f;
		sum_yy = sum_yz = 0.0f;
		sum_zz = 0.0f;

        ///sum the squares of the differences
		for ( unsigned int i=0;i<points.size();i++ )    {
            ///compute the differences
            float diff_x = float ( points[i] ( 0 ) - plane_point ( 0 ) );
            float diff_y = float ( points[i] ( 1 ) - plane_point ( 1 ) );
            float diff_z = float ( points[i] ( 2 ) - plane_point ( 2 ) );
            ///sum the squares
            sum_xx += diff_x * diff_x;
            sum_xy += diff_x * diff_y;
            sum_xz += diff_x * diff_z;

            sum_yy += diff_y * diff_y;
            sum_yz += diff_y * diff_z;

            sum_zz += diff_z * diff_z;
		}

        ///build a matrix
		float mat[3][3];
		mat[0][0] = sum_xx;	mat[0][1] = sum_xy;	mat[0][2] = sum_xz;
		mat[1][0] = sum_xy;	mat[1][1] = sum_yy;	mat[1][2] = sum_yz;
		mat[2][0] = sum_xz;	mat[2][1] = sum_yz;	mat[2][2] = sum_zz;

		float lambda[3];

		eig_sys3d ( mat,lambda );

		Vector3f plane_normal = Vector3f (mat[0][2], mat[1][2], mat[2][2] );

		Vector4f plane_coeffs = Vector4f(plane_normal(0), plane_normal(1), plane_normal(2), -plane_normal.dot(plane_point));


/*
	//FITTING USING ALL POINTS
	//compute the new plane coefficients
	MathFunctions mf;
	double **A = 0x00;
	mf.createCMatrix ( A,points.size(),4 );
	double **s = 0x00;
	double **u = 0x00;
	double **vt = 0x00;
	int row_a = 0;
	Vector<T,4> plane_coeffs;

	//go through the points and create the matrices
	for ( unsigned int i=0;i<points.size();i++ )
	{
		Vector3d pt = convert<double> ( points[i] );
		A[row_a][0] = double ( pt ( 0 )-plane_point ( 0 ) );
		A[row_a][1] = double ( pt ( 1 )-plane_point ( 1 ) );
		A[row_a][2] = double ( pt ( 2 )-plane_point ( 2 ) );
		A[row_a][3] = 1.0;
		row_a++;
	}

	//run the SVD to get the result
	int Y_s,X_s,Y_u,X_u,Y_vt,X_vt;
	mf.computeSVD ( A,points.size(),4,s,Y_s,X_s,u,Y_u,X_u,vt,Y_vt,X_vt );
	if ( !vt )
	{
		std::cout << "vt is empty, number of plane points: " << points.size() << std::endl;
		mf.deleteCMatrix ( A,points.size(),4 );
		mf.deleteCMatrix ( s,Y_s,X_s );
		mf.deleteCMatrix ( u,Y_u,X_u );
		mf.deleteCMatrix ( vt,Y_vt,X_vt );
		return plane_coeffs;
	}
	//the plane coefficients are in the vt (Transpose of V last column->row)
	for ( unsigned int i=0;i<X_vt;i++ )
	{
		plane_coeffs ( i ) = T ( vt[Y_vt-1][i] );
	}

	if ( plane_coeffs == Vector<T,4> ( ( T ) 0.0, ( T ) 0.0, ( T ) 0.0, ( T ) 0.0 ) )
	{
		std::cout << "zero plane coefficients or zero indices returned" <<std::endl;
		mf.deleteCMatrix ( A,points.size(),4 );
		mf.deleteCMatrix ( s,Y_s,X_s );
		mf.deleteCMatrix ( u,Y_u,X_u );
		mf.deleteCMatrix ( vt,Y_vt,X_vt );
		return plane_coeffs;
	}

	//normalize the normal
	Vector<T,3> plane_normal = Vector<T,3> ( plane_coeffs ( 0 ), plane_coeffs ( 1 ), plane_coeffs[2] );
	plane_normal.Normalize();

	plane_coeffs ( 0 ) = plane_normal ( 0 );
	plane_coeffs ( 1 ) = plane_normal ( 1 );
	plane_coeffs ( 2 ) = plane_normal ( 2 );
	plane_coeffs ( 3 ) = -DotProduct ( plane_normal,plane_point );

	mf.deleteCMatrix ( A,points.size(),4 );
	mf.deleteCMatrix ( s,Y_s,X_s );
	mf.deleteCMatrix ( u,Y_u,X_u );
	mf.deleteCMatrix ( vt,Y_vt,X_vt );
	//UP TO HERE
*/

	return plane_coeffs;
}

static void evaluate_plane_fit_function( std::vector<Vector3f> const &points, Vector3f const &plane_point, Vector4f const &plane_coeffs, float const &error_threshold, float *total_error, std::vector<int> *inlier_indices)
{
	( *total_error ) = 0.0f;
	( *inlier_indices ).clear();

	Vector3f plane_normal = Vector3f(plane_coeffs(0), plane_coeffs(1), plane_coeffs(2));
	///go through  the plane points and compute the error
	float squared_error = 0.0f;
	for ( int i=0;i<points.size();i++ ) {
		Vector3f point_to_plane = points[i] - plane_point;
		float distance_to_plane = point_to_plane.dot(plane_normal);
		if ( fabs ( distance_to_plane ) <= error_threshold )    {
			( *inlier_indices ).push_back (i);
		}
		squared_error += distance_to_plane*distance_to_plane;
	}

	( *total_error ) = squared_error/float(points.size());

	return;
}

static bool is_degenerate_function ( std::vector<Vector3f > const &points, int const &minimum_number_of_samples )
{
	int duplicate_points = 0;
	///check if the size is correct
	if ( points.size() < minimum_number_of_samples )    return true;
	///check if there are any duplicate points
	for ( int i=0;i<points.size();i++ ) {
		for (int j=0;j<points.size();j++ )  {
			if ( i==j ) continue;
			if ( points[i] == points[j] ) duplicate_points++;
		}
	}
	if ( duplicate_points!=0 )
	{
		//printf("Degenerate case: %d\n",duplicate_points);
		return true;
	}
	return false;
}

static std::vector<int> randsample ( std::vector<Vector3f > const &points, int number_of_samples )
{
	RNG rng;
	std::vector<int> indices;
	///randomly choose indices
    for ( int i=0;i<number_of_samples;i++ ) {
		double index;
		bool done = false;
		while ( !done )
		{
			done = true;
			index = rng.uniformRand ( 0.0,double ( points.size()-1 ) );
			for ( int j=0;j<indices.size();j++ )
			{
				if ( indices[j] == round ( index ) )
				{
					done = false;
					break;
				}
			}
			if ( done )
			{
//                                 switch(i)       {
//                                         case 0: printf("index: %d ",_round(index));
//                                                 break;
//                                         case 2:printf("%d\n",_round(index));
//                                                                  break;
//                                         default: printf("%d ",_round(index));
//                                                  break;
//                                 }
				indices.push_back ( round ( index ) );
			}
		}
	}
	return indices;
}

static void ransac ( std::vector<Vector3f> const &points,
                     Vector4f ( *fitting_function ) ( std::vector<Vector3f > const &,Vector3f const & ),
                     void ( *evaluate_fit_function ) ( std::vector<Vector3f > const &,Vector3f const &, Vector4f const &, float const &,float *,std::vector<int> * ),
                     bool ( *is_degenerate_function ) ( std::vector<Vector3f > const &, int const & ),
                     int const &minimum_samples_for_fitting,
                     float const &error_threshold,
                     Vector3f const &plane_point,
                     Vector4f *best_plane_coeffs,
                     std::vector<int> *inlier_indices )
{

	float p = 0.95f;    ///Desired probability of choosing at least one sample
	///free from outliers
	int ninliers = 0;
	int maxTrials = 1000;    // Maximum number of trials before we give up.
	int maxDataTrials = 1000; // Max number of attempts to select a non-degenerate data set.
	int trialcount = 0;
	int bestscore =  0;
	std::vector<int> bestinliers;
	Vector4f bestmodel;
	int N = 1;            // Dummy initialisation for number of trials.

	//printf("RANSAC:\n Number of points entered: %d\n",points.size());
	while ( N > trialcount || bestscore == 0 )
	{
		// Select at random s datapoints to form a trial model, M.
		// In selecting these points we have to check that they are not in
		// a degenerate configuration.
		bool is_degenerate = true;
		int count = 1;
		Vector4f plane_coefficients;
		while ( is_degenerate )
		{
			// Generate s random indicies in the range 1..npts
			// (If you do not have the statistics toolbox, or are using Octave,
			// use the function RANDOMSAMPLE from my webpage)
			std::vector<int> indices = randsample ( points, minimum_samples_for_fitting );

			//get the 3d points using the indices
			std::vector<Vector3f > sample_points;
			for (int i=0;i<indices.size();i++ )
			{
				sample_points.push_back ( points[indices[i]] );
			}
			// Test that these points are not a degenerate configuration.
			is_degenerate = is_degenerate_function ( sample_points,minimum_samples_for_fitting );

			if ( !is_degenerate )
			{
				// Fit model to this random selection of data points.
				plane_coefficients = fitting_function ( sample_points,plane_point );
			}
			else
			{
				//Safeguard against being stuck in this loop forever
				count = count + 1;
				if ( count > maxDataTrials )
				{
					printf ( "Unable to select a non-degenerate data set within range.\n" );
					//for (unsigned int i=0;i<points.size();i++)      {
					//	printf("%lf %lf %lf\n",points[i](0),points[i](1),points[i](2));
					//}
					//getchar();
					break;
				}
			}
		}

		// Once we are out here we should have some kind of model...
		// Evaluate distances between points and model returning the indices
		float total_error;
		std::vector<int> indices;
		evaluate_fit_function ( points,plane_point, plane_coefficients,error_threshold,&total_error,&indices );

		//Find the number of inliers to this model.
		ninliers = indices.size();
		if ( ninliers > bestscore )
		{
			// Largest set of inliers so far...
			bestscore = ninliers;  // Record data for this model
			bestinliers = indices;
			bestmodel = plane_coefficients;

			// Update estimate of N, the number of trials to ensure we pick,
			// with probability p, a data set with no outliers.
			float fracinliers =  float ( ninliers ) /float( points.size() );
			float pNoOutliers = ( float ) 1.0 -  pow ( fracinliers,minimum_samples_for_fitting );
			pNoOutliers = max ( float( TC_EPSILON ), pNoOutliers );  // Avoid division by -Inf
			pNoOutliers = min ( float ( 1.0-TC_EPSILON ), pNoOutliers );// Avoid division by 0.
			N = round ( log ( 1.0-p ) /log ( pNoOutliers ) );
			//std::cout << "N " << N << ", pNoOutliers: " << pNoOutliers << std::endl;
		}

		trialcount = trialcount+1;

		// Safeguard against being stuck in this loop forever
		if ( trialcount > maxTrials )
		{
			printf ( "Reached the maximum number of %d trials\n",maxTrials );
			//getchar();
			break;
		}
	}

	if ( bestscore!=0 )
	{
		//We got a solution
		( *best_plane_coeffs ) = bestmodel;
		( *inlier_indices ) = bestinliers;
	}
	else
	{
		printf ( "RANSAC was unable to find a useful solution\n" );
// 	for (unsigned int i=0;i<points.size();i++)	{
// 		std::cout << points[i] << std::endl;
// 	}
// 	exit(0);
	}
	return;
}

#endif
