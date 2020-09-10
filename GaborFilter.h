//
// Copyright © 2006 Charalambos Poullis
//
#ifndef __GABOR_FILTER_H__
#define __GABOR_FILTER_H__

#include "Image.h"
#include "Utilities.h"

#define GABOR_FILTER_SIZE_X  32
#define GABOR_FILTER_SIZE_Y  32

///All angles are specified in radians
class GaborFilter     {
        public:
                GaborFilter(float _sigma_x, float _sigma_y,
                            float _gaussian_rotation_angle, Vector2f _gaussian_peak,
                            float _spatial_frequency_magnitude_F0,float _spatial_frequency_angle_w0,
                            float _P,
                            Image *(&spatial_contributions), Image **(&spatial_frequency_domain),float max_frequency,
                            int _kernel_sizex=GABOR_FILTER_SIZE_X, int _kernel_sizey=GABOR_FILTER_SIZE_Y, const char *basename=0x00) {


                        ///make sure the size of the kernel is odd
                        if (_kernel_sizex%2==0) _kernel_sizex++;
                        if (_kernel_sizey%2==0) _kernel_sizey++;

                        ///allocate memory
                        real_kernel = new Image(_kernel_sizex,_kernel_sizey);
                        real_kernel->clear(0.0f,0.0f,0.0f);
                        imag_kernel = new Image(_kernel_sizex,_kernel_sizey);
                        imag_kernel->clear(0.0f,0.0f,0.0f);

                        ///save the parameters
                        sigma_x = _sigma_x;
                        sigma_y = _sigma_y;
                        gaussian_rotation_angle = _gaussian_rotation_angle;
                        gaussian_peak = _gaussian_peak;
                        spatial_frequency_magnitude_F0 = _spatial_frequency_magnitude_F0;
                        spatial_frequency_angle_w0 = _spatial_frequency_angle_w0;
                        P = _P;
                        kernel_sizex = _kernel_sizex;
                        kernel_sizey = _kernel_sizey;

                        ///build the complex gabor fuction
                        buildGaborFilter();


                        ///build the gabor filter in the frequency domain
                        if (spatial_contributions && spatial_frequency_domain && spatial_frequency_domain[0] && spatial_frequency_domain[1])  {
                                ///if there is input images then add the data to those
                                buildGaborFilterInFrequencyDomain(spatial_frequency_domain,spatial_contributions,max_frequency);
                        }
                        else    {
                                ///otherwise create new ones with just this data only
                                Image **local_spatial_frequency_domain = new Image*[2];
                                int size = 501; //size of the output image for the spatial frequency
                                local_spatial_frequency_domain[0] = new Image(size,size);
                                local_spatial_frequency_domain[0]->clear(0.0f,0.0f,0.0f);
                                local_spatial_frequency_domain[1] = new Image(size,size);
                                local_spatial_frequency_domain[1]->clear(0.0f,0.0f,0.0f);

                                buildGaborFilterInFrequencyDomain(local_spatial_frequency_domain,spatial_contributions,max_frequency);

                                if (basename)   {
                                        local_spatial_frequency_domain[0]->saveImage(_format("%s_filter_frequency_real_local.png",basename));
                                        local_spatial_frequency_domain[1]->saveImage(_format("%s_filter_frequency_imag_local.png",basename));
                                }
                                //clean up
                                ///local_spatial_frequency_domain[0]->Release();
                                ///local_spatial_frequency_domain[1]->Release();
                                delete local_spatial_frequency_domain[0];
                                delete local_spatial_frequency_domain[1];
                                delete [] local_spatial_frequency_domain;
                        }

                        if (basename)   {
                                Image *temp_img = new Image(real_kernel->getWidth(),real_kernel->getHeight());
                                temp_img->copy(real_kernel);
                                temp_img->normalize();
                                temp_img->saveImage(_format("%s_filter_real.png",basename));
				temp_img->saveImage(_format("%s_filter_real.pfm",basename));
                                ///temp_img->Release();
                                delete temp_img;
                                temp_img = new Image(imag_kernel->getWidth(),imag_kernel->getHeight());
                                temp_img->copy(imag_kernel);
                                temp_img->normalize();
                                temp_img->saveImage(_format("%s_filter_imag.png",basename));
				temp_img->saveImage(_format("%s_filter_imag.pfm",basename));
                                ///temp_img->Release();
                                delete temp_img;
                        }
                }

                Image *getRealKernel()     {return real_kernel;}

                Image *getImagKernel()     {return imag_kernel;}

                ///Returns a single image unnormalized
                Image *applyFilter(Image *input_image, float *max_value=0x00, float *min_value=0x00);

                ~GaborFilter()        {
                        ///real_kernel->Release();
                        delete real_kernel;
                        ///imag_kernel->Release();
                        delete imag_kernel;
                }


        private:
                int kernel_sizex,kernel_sizey;
                float sigma_x,sigma_y;
                float gaussian_rotation_angle;
                Vector2f gaussian_peak;
                float spatial_frequency_magnitude_F0,spatial_frequency_angle_w0;     ///polar coordinates
                float P;
                Image *real_kernel;
                Image *imag_kernel;

                ///builds the gabor function
                void buildGaborFilter();

                void buildGaborFilterInFrequencyDomain(Image **(&spatial_frequency_domain), Image *(&spatial_contributions), float max_frequency);

                void thinning(Image *(&input_image));

};

#endif


