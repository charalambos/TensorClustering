////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////
#ifndef __GABOR_JET_H__
#define __GABOR_JET_H__

//#include "opencv.hpp"
#include <opencv2/opencv.hpp>
#include "Utilities.h"
#include "GaborFilter.h"
#include "ImageProcessing.h"

class GaborJet	{
        public:
		///Constructor
                GaborJet(int _kernel_size,
                         int _number_of_angles, double _min_angle, double _max_angle,
                         int _number_of_frequencies,double _min_frequency, double _max_frequency,
                         double _gaussian_sigma=1.0, double _phase=0.0, double _gamma=1.0, bool _save_out = false)   {

                        ///save the parameters
                        kernel_size = _kernel_size;
                        gaussian_sigma = _gaussian_sigma;
                        phase = _phase;
                        gamma = _gamma;

                        number_of_angles = _number_of_angles;
                        number_of_frequencies = _number_of_frequencies;

                        min_angle = _min_angle;
                        max_angle = _max_angle;

                        min_frequency = _min_frequency;
                        max_frequency = _max_frequency;

                        verbose = _save_out;
                }

		///Destructor
                ~GaborJet()     {;}

		///Returns a normalized set of images after applying all the filters
                std::vector<Image *> applyFilters(Image *input_image, std::vector<Vector2i> const &boundary_points);


        private:
                int number_of_angles;
                int number_of_frequencies;
                int kernel_size;
                double sigma_x,sigma_y;
                double gaussian_sigma;
                double phase;
                double gamma;
		        double min_angle, max_angle;
                double min_frequency, max_frequency;
                bool verbose;
};


#endif
