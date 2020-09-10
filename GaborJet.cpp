////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////
#ifndef __GABOR_JET_CPP__
#define __GABOR_JET_CPP__

#include "GaborJet.h"
#include <sys/stat.h>
#include <sys/types.h>


std::vector<Image *> GaborJet::applyFilters(Image *input_image, std::vector<Vector2i> const &boundary_points)       {
    int nsize = 7;
    ///Create the output directory
    int result = mkdir("gabor_filters",0777);
    if (result != 0)    {
        std::cout << "Directory already exists" << std::endl;
    }

    std::vector<Image *> result_images;
    if (!input_image)  return result_images;
    ///Get the width and height of the input image
    int width = input_image->getWidth();
    int height = input_image->getHeight();

    ///Convert Image to cv::Mat
    cv::Mat mat_image;
    ImageProcessing::Image2cvMat(input_image, mat_image);
    cv::Mat mat_image_f;
    mat_image.convertTo(mat_image_f,CV_32F, 1.0/255.0);

    for (int a=0;a<number_of_angles;a++)   {
            double max_value = -DBL_MAX;
            double min_value = DBL_MAX;
            //compute the angle
            double angle = min_angle + double(a)*(max_angle-min_angle)/double(number_of_angles);

            for (int f=0;f<number_of_frequencies;f++)      {
                    int kernel_factor = (f+1)*(f+1)*(f+1);
                    if (kernel_factor%2 == 0) kernel_factor++;

                    std::cout << "\r\033[1;32m[" << int(100.0 * ((a)*number_of_frequencies + (f+1))/(number_of_angles*number_of_frequencies)) << "%]\033[0m Applying Gabor Jets";
                    std::cout.flush();
                    ///compute the frequency
                    double frequency = min_frequency + double(f)*(max_frequency-min_frequency)/double(number_of_frequencies);
                    ///Get the kernel
                    cv::Mat kernel_f = getGaborKernel(  cv::Size(kernel_size*double(kernel_factor), kernel_size),
                                                        gaussian_sigma,
                                                        angle, frequency,
                                                        gamma, phase, CV_32F);

                    ///Apply the kernel
                    cv::Mat result_f;
                    cv::filter2D(mat_image_f, result_f, CV_32F, kernel_f);
                    cv::Mat result;
                    result_f.convertTo(result,CV_8UC3, 255.0);

#ifdef FALSE
                    cv::Mat kernel;
                    kernel_f.convertTo(kernel,CV_8UC3, 255.0);
                    cv::imwrite(_format("gabor_filters/gf_%.3d_%.3d.png",a,f), kernel);
#endif

                    ///Save the image
                    Image *response_image;
                    ImageProcessing::cvMat2Image(result, response_image);
                    ///Remove the response of the boundary points and neighbours
                    for (int i=0;i<boundary_points.size();i++)  {
                        for (int y=-nsize;y<=nsize;y++) {
                            for (int x=-nsize;x<=nsize;x++) {
                                if (outOfBounds(response_image,boundary_points[i](0)+x, boundary_points[i](1)+y))   continue;
                                response_image->setPixel(boundary_points[i](0)+x, boundary_points[i](1)+y,Color(0.0f));
                            }
                        }
                    }

                    ///Save out the response image
                    response_image->saveImage(_format("gabor_filters/gf_result_%.3d_%.3d.pfm",a,f).c_str());
                    response_image->saveImage(_format("gabor_filters/gf_result_%.3d_%.3d.png",a,f).c_str());
                    delete response_image;
            }
    }
    std::cout << std::endl;

	///For each orientation set, normalize the images and add them together
	///this will result to a set of images of the size of the number of angles
	for (int i=0;i<number_of_angles;i++)	{
		std::vector<Image *> orientation_set;
		for (int j=0;j<number_of_frequencies;j++)	{
            std::cout << "\r\033[1;32m[" << int(100.0 * ((i)*number_of_frequencies + (j+1))/(number_of_angles*number_of_frequencies)) << "%]\033[0m Combining responses per angle";
            std::cout.flush();

			//read in the response image
			Image *response_image = new Image();
			response_image->loadImage(_format("gabor_filters/gf_result_%.3d_%.3d.pfm",i,j).c_str());
			orientation_set.push_back(response_image);
		}

		///normalize the set
		normalizeSet(orientation_set);
		///add all the images of this set together
		Image *result_image = new Image(width,height);
		result_image->clear(0.0f,0.0f,0.0f);
		for (int j=0;j<number_of_frequencies;j++)	{
			///add this image to the resulting image
			result_image->add(orientation_set[j]);
			delete orientation_set[j];
		}
		orientation_set.clear();
        result_image->normalize();

		///add the resulting image to the result images
		result_images.push_back(result_image);

        if (verbose)    {
            result_image->saveImage(_format("gabor_filters/gf_result_orientation_set%.3d.pfm",i).c_str());
            result_image->saveImage(_format("gabor_filters/gf_result_orientation_set%.3d.png",i).c_str());
        }
	}

    if (verbose)    {
        ///create an image to add all the responses for the highest frequency
        Image *addition = new Image(width,height);
        addition->clear(0.0f,0.0f,0.0f);
        for (unsigned int i=0;i<result_images.size();i++)       {
                addition->add(result_images[i]);
        }
        addition->normalize();
        addition->saveImage("gabor_filters/sum_result_images.pfm");
        addition->saveImage("gabor_filters/sum_result_images.png");

        delete addition;
    }
    std::cout << std::endl;

    return result_images;
}

#endif
