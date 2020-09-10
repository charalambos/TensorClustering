////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////
#ifndef __GAUSSIAN_DISTRIBUTION_H__
#define __GAUSSIAN_DISTRIBUTION_H__

#include <vector>
using namespace std;

#include <Eigen/Eigen>
#include <iostream>

using namespace Eigen;

#define GD_EPSILON  0.000001

template <typename T, int D>
class GaussianDistribution	{
	public:
		///Constructor
		GaussianDistribution<T,D>()	{
			means =  Matrix<T,D,1>::Zero();
            covariance_matrix = Matrix<T,D,D>::Zero();
            unnorm_covariance_matrix = Matrix<T,D,D>::Zero();
			samples.clear();
			sample_diffs.clear();
			sum_of_samples  = Matrix<T,D,1>::Zero();
            determinant = 0.0f;
            recompute = true;
            constant = 0.0f;
		}


        ///Constructor
        GaussianDistribution<T,D>(Matrix<T,D,1> const &_means, Matrix<T,D,D> const &_covariance_matrix)	{
            means = _means;
            covariance_matrix = _covariance_matrix;
            for (int i=0;i<D;i++)	{
                if (covariance_matrix(i,i) < GD_EPSILON)	{
                    covariance_matrix(i,i) = GD_EPSILON;
                }
            }
            unnorm_covariance_matrix = covariance_matrix;
            covariance_matrix.normalize();
            ///Compute the inverse covariance matrix
            bool invertible;
            covariance_matrix.computeInverseAndDetWithCheck(inv_covariance_matrix,determinant,invertible);
            if (!invertible)	{
                std::cout << "Inverse covariance matrix failed" << std::endl;
                std::cout << covariance_matrix << std::endl;
                std::cout << means << std::endl;
            }
            ///set the flag
            recompute = true;
        }

		///Destructor
		~GaussianDistribution<T,D>()	{
		}

        ///Given a set of samples of type T and dimensions D
        ///this function fits a gaussian of the same dimensions
        void fit(std::vector<Matrix<T,D,1> > const &_samples)	{
            samples.clear();
            sample_diffs.clear();
            sum_of_samples = Matrix<T,D,1>::Zero();
            covariance_matrix = Matrix<T,D,D>::Zero();
            means = Matrix<T,D,1>::Zero();

            for (int i=0;i<_samples.size();i++) {
                    ///Add the new sample to it
                    samples.push_back(_samples[i]);
                    ///Add to the sum of the samples
                    sum_of_samples += _samples[i];
            }
            ///Compute the new means
            means = sum_of_samples/float(samples.size());

            ///Compute the difference from the mean
            for (size_t i=0;i<samples.size();i++)   {
                    Matrix<T,D,1> diffs = samples[i] - means;
                    sample_diffs.push_back(diffs);

                    ///Compute the values for each element of the matrix
                    ///Do it symmetrically
                    for (int j=0;j<D;j++)	{
                        for (int k=j;k<D;k++)	{
                            float diff = diffs(j)*diffs(k);
                            covariance_matrix(j,k) += diff;
                            covariance_matrix(k,j) += diff;
                        }
                    }
            }

            if (samples.size() == 1)    {
                covariance_matrix = Matrix<T,D,D>::Identity();
                return;
            }

            unnorm_covariance_matrix = covariance_matrix;
            //covariance_matrix.normalize();

            ///Compute the inverse covariance matrix
            bool invertible;
            covariance_matrix.computeInverseAndDetWithCheck(inv_covariance_matrix,determinant,invertible);
            if (!invertible)	{
                std::cout << "Inverse covariance matrix failed" << std::endl;
                std::cout << covariance_matrix << std::endl;
                std::cout << means << std::endl;
            }

            ///set the flag
            recompute = true;

            return;
		}

		///This is an update function for processing speed-up. All necessary quantities are kept in memory.
		///So when a new sample is added there is no need to re-fit the distribution on ALL samples
		void update(Matrix<T,D,1> const &sample, bool recalculate)	{

			///Add the new sample to it
			samples.push_back(sample);
			///Add to the sum of the samples
			sum_of_samples += sample;
            ///Recompute the new means
			means = sum_of_samples/float(samples.size());
			///Compute the difference from the mean
            Matrix<T,D,1> diffs = sample - means;
			sample_diffs.push_back(diffs);

            ///Compute the values for each element of the matrix
            ///Do it symmetrically
            for (int j=0;j<D;j++)	{
                for (int k=j;k<D;k++)	{
                    float diff = diffs(j)*diffs(k);
                    unnorm_covariance_matrix(j,k) += diff;
                    unnorm_covariance_matrix(k,j) += diff;
                }
            }

            covariance_matrix = unnorm_covariance_matrix;
            covariance_matrix.normalize();
            ///Compute the inverse covariance matrix
            bool invertible;
            covariance_matrix.computeInverseAndDetWithCheck(inv_covariance_matrix,determinant,invertible);
            if (!invertible)	{
                std::cout << "Inverse covariance matrix failed" << std::endl;
                std::cout << covariance_matrix << std::endl;
                std::cout << means << std::endl;
            }

            ///set the flag
            recompute = true;

            return;
		}

		///Evaluates the gaussian function at the given sample
		float evaluate(Matrix<T,D,1> const &sample)	{

            ///if something changed then the flag is true, so recompute the constant
            if (recompute)  {
                constant = 1.0f/(pow(2.0f*float(M_PI),float(D)/2.0f) * std::sqrt(determinant));
                recompute = false;
            }

			///Compute the difference from the mean and its transpose
            Matrix<T,D,1> diffs = (sample-means);

            Matrix<T,1,1> arg_vec = -T(0.5)*diffs.transpose()*(inv_covariance_matrix*diffs);
            float arg = arg_vec(0,0);

            float val = (float) max(GD_EPSILON,double(constant * exp(arg)));

            return val;
		}

		///Returns the probability of the given sample in this distribution
		///This is the negative log probability value
		float prob(Matrix<T,D,1> const &_sample)	{
			///Evaluate the gaussian function at this point
			float eval = evaluate(_sample);
			///Compute the negative log probability
			float neg_log_prob = -log(eval);
			return neg_log_prob;
		}

        Matrix<T,D,1> getMeans()	{
			return means;
		}

		///HACKING function. DO NOT USE unless you know what you are doing
		void setMeans(Matrix<T,D,1> const &_mean)	{
			means = _mean;
		}
		void setCovarianceMatrix(Matrix<T,D,D> const &_covariance_matrix)	{
			covariance_matrix = _covariance_matrix;
		}


		Matrix<T,D,D> getCovarianceMatrix()	{
			return covariance_matrix;
		}

		void getSamples(std::vector<Matrix<T,D,1> > &_samples)	{
			for (int i=0;i<samples.size();i++) {
                _samples.push_back(samples[i]);
			}
		}

	private:
		///The mean
		Matrix<T,D,1> means;
		///The unbiased variance
		Matrix<T,D,D> covariance_matrix;
        ///The unnormalized unbiased variance
        Matrix<T,D,D> unnorm_covariance_matrix;
        ///The inverse covariance matrix (kept for speed optimization)
        Matrix<T,D,D> inv_covariance_matrix;

        ///the determinant of the covariance matrix
        float determinant;
        ///the recompute flag
        bool recompute;
        ///the constant
        float constant;

		///Private information used for speeding up the process
		///The samples
		std::vector<Matrix<T,D,1> > samples;
		///The diffs of the samples from the means
		std::vector<Matrix<T,D,1> > sample_diffs;
		///The sum of the samples
        Matrix<T,D,1> sum_of_samples;
};

typedef	GaussianDistribution<float,1>	GaussianDistribution1f;
typedef	GaussianDistribution<float,2>	GaussianDistribution2f;
typedef	GaussianDistribution<float,3>	GaussianDistribution3f;


#endif
