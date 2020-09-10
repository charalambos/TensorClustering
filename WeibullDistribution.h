////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////
#ifndef __WEIBULL_DISTRIBUTION_H__
#define __WEIBULL_DISTRIBUTION_H__

#include <vector>
using namespace std;

#include <Eigen/Eigen>
using namespace Eigen;

#define USING_MPFR  1

#ifdef USING_MPFR
    #include <mpreal.h>
    using namespace mpfr;
#endif

#define WD_EPSILON 0.0001f

class WeibullDistribution	{
	public:
		///Constructor
		WeibullDistribution()	{
#ifdef USING_MPFR
            using mpfr::mpreal;

            // Required precision of computations in decimal digits
            // Play with it to check different precisions
            const int digits = 5;

            // Setup default precision for all subsequent computations
            // MPFR accepts precision in bits - so we do the conversion

            mpreal::set_default_prec(mpfr::digits2bits(digits));

            // Compute all the vital characteristics of mpreal (in current precision)
            // Analogous to lamch from LAPACK

            const mpreal one         =    1.0;
            const mpreal zero        =    0.0;
            const mpreal eps         =    std::numeric_limits<mpreal>::epsilon();
            const int    base        =    std::numeric_limits<mpreal>::radix;
            const mpreal prec        =    eps * base;
            const int bindigits      =    std::numeric_limits<mpreal>::digits(); // eqv. to mpfr::mpreal::get_default_prec();
            const mpreal rnd         =    std::numeric_limits<mpreal>::round_error();
            const mpreal maxval      =    std::numeric_limits<mpreal>::max();
            const mpreal minval      =    std::numeric_limits<mpreal>::min();
            const mpreal small       =    one / maxval;
            const mpreal sfmin       =    (small > minval) ? small * (one + eps) : minval;
            const mpreal round       =    std::numeric_limits<mpreal>::round_style();
            const int    min_exp     =    std::numeric_limits<mpreal>::min_exponent;
            const mpreal underflow   =    std::numeric_limits<mpreal>::min();
            const int    max_exp     =    std::numeric_limits<mpreal>::max_exponent;
            const mpreal overflow    =    std::numeric_limits<mpreal>::max();
#endif

			mean = 0.0f;
			variance = 0.0f;

			samples.clear();

            ///shape parameter
			alpha = 1.0f;

            ///scale parameter
			beta = 1.0f;

            sum_log_Xi = 0.0f;
            log_Xi.clear();

		}


		///Destructor
		~WeibullDistribution()	{
		}

        /// Returns f'(x) given a pointer to f, and value x.
#ifdef USING_MPFR
        mpreal derivative(mpreal a)   {
            mpreal da = 0.01f * a;
            if (da < WD_EPSILON)    return 0.0;

            mpreal result = (derivative_of_log_likelihood_weibull(a+da) - derivative_of_log_likelihood_weibull(a))/da;

            return result;
        }
#else
        float derivative(float a)   {

            float da = 0.01f * a;
            if (da < WD_EPSILON)    return 0.0f;

            float result = (derivative_of_log_likelihood_weibull(a+da) - derivative_of_log_likelihood_weibull(a))/da;

            return result;
        }
#endif

        /// Solves f(x) = 0 for x,
        /// Provide: function pointer f, a guess (x0), and the required precision.
#ifdef USING_MPFR
        mpreal NewtonRaphson(mpreal x0, mpreal precision, int remaining_iterations = 5)  {
            if (remaining_iterations==0) {
                    //std::cout << "Reached maximum number of iterations." << std::endl;
                    return x0;
            }

            mpreal dx = derivative(x0);
            if (abs(dx) < WD_EPSILON)    {
                //std::cout << "dx is 0" << std::endl;
                return x0;
            }
            mpreal x1 = x0 - derivative_of_log_likelihood_weibull(x0)/dx;

            if (abs(x1-x0) < precision) {
                    return x1;
            }

            return NewtonRaphson( x1 , precision, remaining_iterations-1);
        }
#else
    float NewtonRaphson(float x0, float precision, int remaining_iterations = 5)  {
            if (remaining_iterations==0) {
                //std::cout << "Reached maximum number of iterations." << std::endl;
                return x0;
            }

            float dx = derivative(x0);
            if (abs(dx) < WD_EPSILON)    {
                //std::cout << "dx is 0" << std::endl;
                return x0;
            }
            float x1 = x0 - derivative_of_log_likelihood_weibull(x0)/dx;

            if (abs(x1-x0) < precision) {
                return x1;
            }

            return NewtonRaphson( x1 , precision, remaining_iterations-1);
        }
#endif

#ifdef USING_MPFR
        mpreal derivative_of_log_likelihood_weibull(mpreal a) {

            mpreal n = mpreal(samples.size());

            mpreal sum_Xi_pow_a = 0.0;
            mpreal sum_Xi_pow_a_log_Xi = 0.0;
            mpreal Xi_pow_a = 0.0;
            mpreal sample_log_Xi = 0.0;
            for (int i=0;i<samples.size();i++)  {
                Xi_pow_a = pow(samples[i], a);
                sample_log_Xi = log_Xi[i];

                sum_Xi_pow_a += Xi_pow_a;
                sum_Xi_pow_a_log_Xi += Xi_pow_a*sample_log_Xi;
            }
//std::cout << sum_log_Xi << " " << sum_Xi_pow_a << " " << sum_Xi_pow_a_log_Xi << std::endl;
           mpreal result = (sum_Xi_pow_a_log_Xi/max(WD_EPSILON,sum_Xi_pow_a))-(1.0/max(WD_EPSILON,a))-sum_log_Xi/n;

            return result;
        }
#else
        float derivative_of_log_likelihood_weibull(float a) {

            float n = float(samples.size());

            float sum_Xi_pow_a = 0.0;
            float sum_Xi_pow_a_log_Xi = 0.0;
            float Xi_pow_a = 0.0;
            float sample_log_Xi = 0.0;
            for (int i=0;i<samples.size();i++)  {
                Xi_pow_a = pow(samples[i], a);
                sample_log_Xi = log_Xi[i];

                sum_Xi_pow_a += Xi_pow_a;
                sum_Xi_pow_a_log_Xi += Xi_pow_a*sample_log_Xi;
            }
    //std::cout << sum_log_Xi << " " << sum_Xi_pow_a << " " << sum_Xi_pow_a_log_Xi << std::endl;
            float result = (sum_Xi_pow_a_log_Xi/std::max(WD_EPSILON,sum_Xi_pow_a))-(1.0/std::max(WD_EPSILON,a))-sum_log_Xi/n;

            return result;
        }
#endif

#ifdef USING_MPFR
		void fit(std::vector<float> const &_samples)	{
            samples.clear();
            log_Xi.clear();
            sum_log_Xi = 0.0;

			///Add the new sample to it
			for (int i=0;i<_samples.size();i++) {
                samples.push_back(_samples[i]);
                ///Precompute log_Xi
                mpreal sample_log_Xi = log(samples[i]);
                log_Xi.push_back(sample_log_Xi);
                sum_log_Xi += sample_log_Xi;
            }
            mpreal n = mpreal(samples.size());

            ///Optimize the alpha - shape
            alpha = max(WD_EPSILON, NewtonRaphson(alpha, WD_EPSILON));

            ///Compute the beta - scale
            mpreal sum_Xi_pow_a = 0.0;
            for (int i=0;i<samples.size();i++)  {
                sum_Xi_pow_a += pow(samples[i], alpha);
            }

            mpreal one_over_alpha = 1.0/alpha;

            beta = max(WD_EPSILON, pow(sum_Xi_pow_a/n, one_over_alpha));

            ///Compute the mean
            mean = beta*gamma(1.0+one_over_alpha);

            ///Compute the variance
            mpreal temp = gamma(1.0+one_over_alpha);
            variance = beta*beta*(gamma(1.0+2.0*one_over_alpha) - temp*temp);

			return;
		}
#else
        void fit(std::vector<float> const &_samples)	{
            samples.clear();
            log_Xi.clear();
            sum_log_Xi = 0.0;

            ///Add the new sample to it
            for (int i=0;i<_samples.size();i++) {
                samples.push_back(_samples[i]);
                ///Precompute log_Xi
                float sample_log_Xi = log(samples[i]);
                log_Xi.push_back(sample_log_Xi);
                sum_log_Xi += sample_log_Xi;
            }
            float n = float(samples.size());

            ///Optimize the alpha - shape
            alpha = max(WD_EPSILON, NewtonRaphson(alpha, WD_EPSILON));

            ///Compute the beta - scale
            float sum_Xi_pow_a = 0.0;
            for (int i=0;i<samples.size();i++)  {
                sum_Xi_pow_a += pow(samples[i], alpha);
            }

            float one_over_alpha = 1.0/alpha;

            beta = max(WD_EPSILON, pow(sum_Xi_pow_a/n, one_over_alpha));

            ///Compute the mean
            mean = beta*gamma(1.0+one_over_alpha);

            ///Compute the variance
            float temp = gamma(1.0+one_over_alpha);
            variance = beta*beta*(gamma(1.0+2.0*one_over_alpha) - temp*temp);

            return;
        }
#endif

		///This is an update function for processing speed-up. All necessary quantities are kept in memory.
		///So when a new sample is added there is no need to re-fit the distribution on ALL samples
#ifdef USING_MPFR
        void update(float sample, bool recalculate)	{

			///Add the new sample to it
			samples.push_back(sample);

            ///Precompute log_Xi
            mpreal sample_log_Xi = log(sample);
            log_Xi.push_back(sample_log_Xi);
            sum_log_Xi += sample_log_Xi;

            if (!recalculate)   return;

			if (samples.size() == 1)    {
                alpha = 0.5f;
                beta = 1.0f;
                alpha_over_beta = alpha/beta;
                variance = sample * 0.5f;
                mean = sample;
                return;
			}

            mpreal n = mpreal(samples.size());

            ///Optimize the alpha - shape
            alpha = max(WD_EPSILON, NewtonRaphson(alpha, WD_EPSILON));

            ///Compute the beta - scale
            mpreal sum_Xi_pow_a = 0.0;
            for (int i=0;i<samples.size();i++)  {
                sum_Xi_pow_a += pow(samples[i], alpha);
            }

            mpreal one_over_alpha = 1.0/alpha;

            beta = max(WD_EPSILON, pow(sum_Xi_pow_a/n, one_over_alpha));

            alpha_over_beta = alpha/beta;

            ///Compute the mean
            mean = beta*gamma(1.0+one_over_alpha);

            ///Compute the variance
            mpreal temp = gamma(1.0+one_over_alpha);
            variance = beta*beta*(gamma(1.0+2.0*one_over_alpha) - temp*temp);

            //std::cout << alpha.toFloat() << " " << beta.toFloat() << " " << mean.toFloat() << " " << variance.toFloat() << std::endl;
			return;
		}
#else
        void update(float sample, bool recalculate = true)	{

            ///Add the new sample to it
            samples.push_back(sample);

            ///Precompute log_Xi
            float sample_log_Xi = log(sample);
            log_Xi.push_back(sample_log_Xi);
            sum_log_Xi += sample_log_Xi;

            if (!recalculate)   return;

            if (samples.size() == 1)    {
                alpha = 0.5f;
                beta = 1.0f;

                variance = sample * 0.5f;
                mean = sample;
                return;
            }

            float n = float(samples.size());

            ///Optimize the alpha - shape
            alpha = max(WD_EPSILON, NewtonRaphson(alpha, WD_EPSILON));

            ///Compute the beta - scale
            float sum_Xi_pow_a = 0.0;
            for (int i=0;i<samples.size();i++)  {
                sum_Xi_pow_a += pow(samples[i], alpha);
            }

            float one_over_alpha = 1.0f/alpha;

            beta = max(WD_EPSILON, pow(sum_Xi_pow_a/n, one_over_alpha));

            ///Compute the mean
            mean = beta*gamma(1.0f+one_over_alpha);

            ///Compute the variance
            float temp = gamma(1.0f+one_over_alpha);
            variance = beta*beta*(gamma(1.0+2.0*one_over_alpha) - temp*temp);

            //std::cout << alpha.toFloat() << " " << beta.toFloat() << " " << mean.toFloat() << " " << variance.toFloat() << std::endl;
            return;
        }
#endif

        ///Evaluates the weibull distribution function at the given sample
#ifdef USING_MPFR
        float evaluate(float sample)	{
            mpreal temp = sample/beta;
			mpreal val = (alpha_over_beta) * pow(temp,alpha-1.0f) * exp(-pow(temp, alpha));
if (isnan(val)) {
            std::cout << val.toFloat() << std::endl;
            std::cout << (alpha/beta).toFloat() << std::endl;
            std::cout << pow(sample/beta,alpha-1.0f).toFloat() << " " <<  sample << " " << alpha.toFloat() << " " << beta.toFloat() << std::endl;
            std::cout << exp(-pow(sample/beta, alpha)).toFloat() << std::endl;

        for (int i=0;i<samples.size();i++)  {
            std::cout << i <<": " <<samples[i].toFloat() << std::endl;
        }
        exit(0);
}
			return val.toFloat();
		}
#else
        float evaluate(float sample)	{

            float val = (alpha/beta) * pow(sample/beta,alpha-1.0f) * exp(-pow(sample/beta, alpha));
            if (isnan(val)) {
                std::cout << val << std::endl;
                std::cout << (alpha/beta) << std::endl;
                std::cout << pow(sample/beta,alpha-1.0f) << " " <<  sample << " " << alpha << " " << beta << std::endl;
                std::cout << exp(-pow(sample/beta, alpha)) << std::endl;

                for (int i=0;i<samples.size();i++)  {
                    std::cout << i <<": " <<samples[i] << std::endl;
                }
                exit(0);
            }
            return val;
        }
#endif

		///Returns the probability of the given sample in this distribution
		///This is the negative log probability value
		float prob(float _sample)	{
			///Evaluate the gaussian function at this point
			float eval = evaluate(_sample);
			///Compute the negative log probability
			float neg_log_prob = -log(eval);
			return neg_log_prob;
		}

#ifdef USING_MPFR
        float getMean()	{
		    return mean.toFloat();
		}
#else
        float getMean()	{
            return mean;
        }
#endif
		///HACKING function. DO NOT USE unless π.χ.
		void setMeans(float _mean)	{
			mean = _mean;
		}
		void setVariance(float _variance)	{
			variance = _variance;
		}

#ifdef USING_MPFR
		float getVariance()	{
			return variance.toFloat();
		}
#else
        float getVariance()	{
            return variance;
        }
#endif

#ifdef USING_MPFR
		void getSamples(std::vector<float> &_samples)	{
            _samples.clear();
            for (int i=0;i<samples.size();i++)  {
                _samples.push_back(samples[i].toFloat());
            }
            return;
		}
#else
        void getSamples(std::vector<float> &_samples)	{
            _samples.clear();
            for (int i=0;i<samples.size();i++)  {
                _samples.push_back(samples[i]);
            }
            return;
        }
#endif

#ifdef USING_MPFR
        float getShape()    {
            return alpha.toFloat();
        }
#else
        float getShape()    {
            return alpha;
        }
#endif

#ifdef USING_MPFR
        float getScale()    {
            return beta.toFloat();
        }
#else
        float getScale()    {
            return beta;
        }
#endif

	private:
		///The mean
#ifdef USING_MPFR
		mpreal mean;
#else
        float mean;
#endif
		///The unbiased variance
#ifdef USING_MPFR
		mpreal variance;
#else
        float variance;
#endif
		///The shape parameter
#ifdef USING_MPFR
		mpreal alpha;
#else
        float alpha;
#endif
		///The scale parameter
#ifdef USING_MPFR
        mpreal beta;
#else
        float beta;
#endif
		///The samples
#ifdef USING_MPFR
		std::vector<mpreal> samples;
#else
        std::vector<float> samples;
#endif

        ///For optimization
#ifdef USING_MPFR
        std::vector<mpreal> log_Xi;
#else
        std::vector<float> log_Xi;
#endif

#ifdef USING_MPFR
        mpreal sum_log_Xi;
#else
        float sum_log_Xi;
#endif

#ifdef  USING_MPFR
        mpreal alpha_over_beta;
#else
        float alpha_over_beta;
#endif

};

#endif

