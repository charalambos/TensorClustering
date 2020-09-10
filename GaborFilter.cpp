//
// Copyright © 2006 Charalambos Poullis
//
#ifndef __GABOR_FILTER_CPP__
#define __GABOR_FILTER_CPP__

#include "GaborFilter.h"

void GaborFilter::buildGaborFilter()   {
        int half_sizex = kernel_sizex/2;
        int half_sizey = kernel_sizey/2;

        ///declare some constants
        float cosw = cos(spatial_frequency_angle_w0);
        float sinw = sin(spatial_frequency_angle_w0);
        float cos_theta = cos(gaussian_rotation_angle);
        float sin_theta = sin(gaussian_rotation_angle);

        float u0 = spatial_frequency_magnitude_F0*cosw;
        float v0 = spatial_frequency_magnitude_F0*sinw;

        float two_pi = 2.0f*M_PI;
        float K = 1.0f/(two_pi*sigma_x*sigma_y);
        float a = 1.0f/sigma_x;
        float a2 = a*a;
        float b = 1.0f/sigma_y;
        float b2 = b*b;

        float param;
        float real_value,imag_value;
	for (int y=-half_sizey,local_y=0;y<=half_sizey;y++,local_y++)       {
                for (int x=-half_sizex,local_x=0;x<=half_sizex;x++,local_x++)       {

                        ///compute more contants;
                        float XminusX0r = (float(x)-gaussian_peak(0))*cos_theta + (float(y)-gaussian_peak(1))*sin_theta;
                        float YminusY0r = -(float(x)-gaussian_peak(0))*sin_theta + (float(y)-gaussian_peak(1))*cos_theta;

                        ///first
                        float component1 = K*exp(-M_PI*(a2*XminusX0r*XminusX0r + b2*YminusY0r*YminusY0r));

                        ///second
                        float param2 = two_pi*(u0*float(x)+v0*float(y))+P;
                        float real_component2 = cos(param2);
                        float imag_component2 = sin(param2);

                        real_value = component1 * real_component2;
                        imag_value = component1 * imag_component2;

                        ///real part
                        real_kernel->setPixel(local_x, local_y, Color(real_value));
                        ///imaginary part
                        imag_kernel->setPixel(local_x, local_y, Color(imag_value));
                }
        }

	///Remove the DC coefficient from the real part
	float m = 0.0f;
	for (int y=0;y<kernel_sizey;y++)	{
		for (int x=0;x<kernel_sizex;x++)	{
			m += real_kernel->getPixel(x, y).r();
		}
	}
	m /= float(kernel_sizex*kernel_sizey);
	for (int y=0;y<kernel_sizey;y++)	{
		for (int x=0;x<kernel_sizex;x++)	{
			real_kernel->setPixel(x, y, real_kernel->getPixel(x,y)-Color(m));
		}
	}

        return;
}

void GaborFilter::buildGaborFilterInFrequencyDomain(Image **(&spatial_frequency_domain), Image *(&spatial_contributions),float max_frequency)    {
        int size = spatial_frequency_domain[0]->getWidth();
        int half_size = size/2;

        ///declare some constants
        float cosw = cos(spatial_frequency_angle_w0);
        float sinw = sin(spatial_frequency_angle_w0);
        float cos_theta = cos(gaussian_rotation_angle);
        float sin_theta = sin(gaussian_rotation_angle);

        float u0 = spatial_frequency_magnitude_F0*cosw;
        float v0 = spatial_frequency_magnitude_F0*sinw;

        float two_pi = 2.0f*M_PI;
        float K = 1.0f/(two_pi*sigma_x*sigma_y);
        float a = 1.0f/sigma_x;
        float a2 = a*a;
        float b = 1.0f/sigma_y;
        float b2 = b*b;
        float coeff = K/(a*b);

        ///cosine and sine vary from -1 to 1
        float min_f = -max_frequency;
        float max_f = max_frequency;
        float max_minus_min = max_f - min_f;
        if (max_minus_min == 0.0f)      max_minus_min = 1.0f;
	///max_minus_min is the whole width and height. i want half of that
	float half_width = max_minus_min/2.0f;
	float half_height = max_minus_min/2.0f;
        float u0_scaled = ((u0-min_f)/half_width) - 1.0f;
        float v0_scaled = ((v0-min_f)/half_height) - 1.0f;

        float real_value, imag_value, param1, param2;
        for (int v=-half_size,local_v=0;v<=half_size;v++,local_v++)    {
                for (int u=-half_size,local_u=0;u<=half_size;u++,local_u++)    {
			///bring it to -1 to 1
                        float u_scaled = float(u + half_size)/float(half_size)-1.0f;
                        float v_scaled = float(v + half_size)/float(half_size)-1.0f;

                        float UminusU0r = (u_scaled-u0_scaled)*cos_theta + (v_scaled-v0_scaled)*sin_theta;
                        float VminusV0r =-(u_scaled-u0_scaled)*sin_theta + (v_scaled-v0_scaled)*cos_theta;

                        param1 = -two_pi*(gaussian_peak(0)*(u_scaled-u0_scaled) + gaussian_peak(1)*(v_scaled-v0_scaled)) + P;
                        param2 = exp(-M_PI*(UminusU0r*UminusU0r/a2 + VminusV0r*VminusV0r/b2));

                        ///real part
                        real_value = coeff * cos(param1) * param2;
                        spatial_frequency_domain[0]->add(local_u, local_v, Color(real_value));

                        ///imaginary part
                        imag_value = coeff * sin(param1) * param2;
                        spatial_frequency_domain[1]->add(local_u, local_v, Color(imag_value));

                        spatial_contributions->add(local_u, local_v, Color(real_value*imag_value));
                }
        }

        return;
}

Image *GaborFilter::applyFilter(Image *input_image, float *max_value, float *min_value) {
        ///check if input is ok
        if (!input_image || input_image->getWidth() == 0 || input_image->getHeight()== 0)      return 0x00;

        int half_sizex = real_kernel->getWidth()/2;  ///should be an odd number equal to the half size of the_kernel
        int half_sizey = real_kernel->getHeight()/2;

        ///allocate memory
        Image *result_image = new Image(input_image->getWidth(),input_image->getHeight());
        result_image->clear(0.0f,0.0f,0.0f);

        float global_sum_real = 0.0f;
        float global_sum_imag = 0.0f;

        ///for every pixel in the image
        for (int i=0;i<input_image->getHeight();i++) {
                for (int j=0;j<input_image->getWidth();j++) {

                        float local_sum_real = 0.0f;
                        float local_sum_imag = 0.0f;
                        bool not_good = false;

                        ///apply the filter
                        int local_y = 0;
                        for (int y=-half_sizey;y<=half_sizey;y++)       {
                                if (not_good)   break;
                                int local_x = 0;
                                for (int x=-half_sizex;x<=half_sizex;x++)       {
                                        int index_y = i+y;
                                        int index_x = j+x;
                                        ///boundary check
                                        if (index_y < 0 || index_y >= input_image->getHeight() || index_x < 0 || index_x >= input_image->getWidth())   {
                                                not_good = true;
                                                break;
                                        }

                                        //real part
                                        local_sum_real  += (input_image->getPixel(index_x, index_y).r() * real_kernel->getPixel(local_x,local_y).r()) +
                                        								(input_image->getPixel(index_x, index_y).g() * real_kernel->getPixel(local_x,local_y).g()) +
                                        								(input_image->getPixel(index_x, index_y).b() * real_kernel->getPixel(local_x,local_y).b());

                                        //imag part
                                        local_sum_imag  += (input_image->getPixel(index_x, index_y).r() * imag_kernel->getPixel(local_x, local_y).r()) +
                                        								 (input_image->getPixel(index_x, index_y).g() * imag_kernel->getPixel(local_x, local_y).g()) +
                                        								 (input_image->getPixel(index_x, index_y).b() * imag_kernel->getPixel(local_x, local_y).b());
                                        local_x++;
                                }
                                local_y++;
                        }
                        if (not_good)   continue;

                        global_sum_real += local_sum_real;
                        global_sum_imag += local_sum_imag;

                        //compute result
                        float value = sqrtf(local_sum_real*local_sum_real + local_sum_imag*local_sum_imag);
			//if the value is less than epsilon then make it zero because it
			//will cause problems with the normalization later
			if (value < TC_EPSILON) value = 0.0f;
                        result_image->setPixel(j, i, Color(value));

                        //find max or min
                        if (max_value) {
                                if ((*max_value) < value) (*max_value) = value;
                        }
                        if (min_value) {
                                if ((*min_value) > value) (*min_value) = value;
                        }
                }
        }

        return result_image;

}

void GaborFilter::thinning(Image *(&input_image))       {
        int kernel_size = 3;
        int half_kernel_size = kernel_size/2;

        for (int y=0;y<input_image->getHeight();y++) {
                for (int x=0;x<input_image->getWidth();x++) {
                        ///perform the thinning

                        ///first find the maximum value in the neighbourhood
                        float max_value=0.0f;
                        for (int i=y-half_kernel_size;i<=y+half_kernel_size;i++)        {
                                for (int j=x-half_kernel_size;j<=x+half_kernel_size;j++)        {
                                        if (i < 0 || i >= input_image->getHeight() || j < 0 || j >= input_image->getWidth())   continue;
                                        if (max_value <  input_image->getPixel(j,i).r())     max_value = input_image->getPixel(j,i).r();
                                }
                        }

                        ///then check if the pixel in question if less than the max value. If it is, then set it to zero, otherwize move on
                        if (input_image->getPixel(x,y).r() < max_value)      input_image->setPixel(x, y, Color(0.0f));
                }
        }

        return;
}

#endif


