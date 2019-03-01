////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	      //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////

#include <sys/stat.h>
#include <Eigen/Eigen>

#include <Eigen/Eigen>
using namespace Eigen;

#include "Patch.h"
#include "GaborJet.h"

#include "ImageProcessing.h"
#include "opencv2/opencv.hpp"
#include "GeometryExporter.h"
#include "GeometryProcessing.h"
#include "Tensor.h"
#include "WeibullDistribution.h"

#define TENSOR_DEBUG
//#define EXPORT_DISTRIBUTION_ANIMATION
//#define EXPORT2R
//#undef TENSOR_DEBUG

//#define EXPORT_TENSOR_STEP

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/print_wavefront.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <math.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Exact_predicates_exact_constructions_kernel K;
//typedef CGAL::Simple_cartesian<double> K;

typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3  Point;
typedef K::Plane_3  Plane;

std::string file_name;              ///The file name of the XYZ map
std::string name;
FILE_TYPE type;


///parameters for the gabor filters
int kernel_size = 5;
int number_of_angles = 16;
int number_of_frequencies = 5;
double min_angle = 0.0;
double max_angle = 2.0*M_PI;
double min_frequency = 2.0;
double max_frequency = 20.0;
double phase = M_PI_2;
double gaussian_sigma = 2.0;
double elongation_gamma = 1.0;
double verbose = true;
FILE *fileptr_export2R = 0x00;
int tensor_count = 0;
std::vector<Patch *> patches;

static Image *loadXYZmap()    {
    Image *xyz_map = new Image();
    if (xyz_map->loadImage(file_name))  {
        return xyz_map;
    }
    else    {
        delete xyz_map;
        xyz_map = 0x00;
        return 0x00;
    }
}

static void applyGaborJets(Image *input_image, std::vector<Vector2i> const &boundary_points, std::vector<Image *> &response_images)	{
	///Apply the gabor jets and return the sum of the response images

	///Create the gabor jet object using the parameters above
	GaborJet *gj = new GaborJet(kernel_size,
                                number_of_angles,min_angle,max_angle,
								number_of_frequencies,min_frequency,max_frequency,
								gaussian_sigma,
								phase, elongation_gamma, verbose);

	response_images = gj->applyFilters(input_image, boundary_points);

	///Create a single image as the sum of the response images
	Image *sum_response_images = new Image(input_image->getWidth(), input_image->getHeight());
	sum_response_images->clear(Color(0.0f));
	for (int i=0;i<response_images.size();i++)	{
		sum_response_images->add(response_images[i]);
	}
	//sum_response_images->blur(2.0f);
	//sum_response_images->saveImage("sum_response_image1.pfm");
	///Normalize
	sum_response_images->normalize();
	sum_response_images->saveImage("sum_response_image.pfm");

	return;
}

static void exportTensor2R(Tensor const &tensor,int x,int y,float z, float scale=5.0f)    {

    Matrix3f mat = tensor.getTensorInMatrixForm();
    fprintf(fileptr_export2R,"Sigma%d <- matrix(c(",tensor_count);
    for (int m=0;m<3;m++)   {
        for (int n=0;n<3;n++)   {
            if (m==2 && n ==2)  fprintf(fileptr_export2R,"%f",mat(m,n));
            else                fprintf(fileptr_export2R,"%f,",mat(m,n));
        }
    }
    fprintf(fileptr_export2R,"), 3,3)\n");
    fprintf(fileptr_export2R,"Mean%d<-c(%f, %f, %f)\n",tensor_count,float(x)*scale,float(y)*scale,z);
    //fprintf(fileptr_export2R,"x%d <- MASS::mvrnorm(1000, Mean%d, Sigma%d)\n",tensor_count,tensor_count,tensor_count);
    fprintf(fileptr_export2R,"plot3d(ellipse3d(Sigma%d, centre=Mean%d), col=c1[runif(1, 1, length(c1))], alpha=0.5, add = TRUE)\n\n",tensor_count,tensor_count);
//arrow3d(c(0,0,0), c(1,0,0), barblen=.2, lwd=3, col="red")


    tensor_count++;

    return;
}

static void encodeResponses2Tensor(Image *input_image, Image *height_and_normal_variation_map, std::vector<Image *> &response_images, Image *normal_map, Image *saliency_map,Image *orientation_map, std::vector<Tensor> &pixels_as_tensors, Image *tensor_index_map)  {
	static Matrix3f zero_matrix;
	zero_matrix.setZero();
	///Converts each pixel from the gabor responses into a tensor and sums all the responses
	///it results in a saliency map with (l1-l2,l2-l3,l3) values and an orientation map

	///Check if there are response images
    if (!response_images.size())    return;

    ///the size of the images
	int sizex = response_images[0]->getWidth();
	int sizey = response_images[0]->getHeight();

    ///The number of orientations.
	float n = float(response_images.size());

	///FIRST STEP IS TO CONVERT THE GABOR JET RESPONSES TO PLATE TENSORS CORRESPONDING TO CURVELS
	///AND CONVERT THE NORMAL MAP TO STICK TENSORS CORRESPONDING TO SURFACE PIXELS
	Image *response_intensity_map = new Image(sizex, sizey);
    ///Get the range of the input image
    Color *range = input_image->range();
    float max_minus_min = range[1].b() - range[0].b();
	///Compute the weighted average of the response images
	for (int y=0;y<sizey;y++)       {
		for (int x=0;x<sizex;x++)      {
			Matrix3f tensor_matrix;
			tensor_matrix.setZero();

            ///Get the e1 vector from the normal map
            Vector3f normal = color2vector3(normal_map->getPixel(x,y));
            ///Use the normal as eigenvector e1
            Vector3f e1 = normal;
            if (e1.norm() > TC_EPSILON) e1.normalize();
            if ( (e1 - Vector3f(0.0f,0.0f,0.0f)).norm() < TC_EPSILON)	{
                //std::cout << x << " " << y << std::endl;
                continue;
            }

            ///Get the height of the pixel
            float height = input_image->getPixel(x,y).b();
            ///Normalize height
            height = 255.0f*(height - range[0].b())/max_minus_min;

            ///Get the height variation of this pixel
            float height_variation = height_and_normal_variation_map->getPixel(x,y).r();

            ///Get the normal variation of this pixel
            float normal_variation = height_and_normal_variation_map->getPixel(x,y).g();

			for (int i=0;i<response_images.size();i++)        {
                std::cout << "\r\033[1;32m[" << int(100.0 * (y*sizex*response_images.size() + x*response_images.size() + i+1)/(sizey*sizex*response_images.size())) << "%]\033[0m Encoding responses to tensors";
                std::cout.flush();
				///Get the orientation angle for this one
				float orientation_angle = min_angle + float(i)*(max_angle-min_angle)/float(number_of_angles);

				///compute the direction of this response image
				float cos_theta = float(cos(orientation_angle));
				float sin_theta = float(sin(orientation_angle));

				///The gabors are offset by 90 so this is tangent should be rotated before used i.e. (sin,cos)
				Vector3f tangent_offset = Vector3f(cos_theta,sin_theta,0.0f);
				Vector3f e3 = Vector3f(tangent_offset(1),tangent_offset(0),0.0f);           ///Flipped due to above

				///The eigenvalues: the gabors detect curve pixels so the eigenvalue will have l1=l2=1 and l3=0
				///actually it has to be l2-l3>l3 & l2-l3>l1-l2
				///convert all the pixels to tensors using the information computed above and perform tensor addition

				Vector3f e2 = e1.cross(e3);
                if (e2.norm() > TC_EPSILON) e2.normalize();
                e3 = e2.cross(e1);
                if (e3.norm() > TC_EPSILON) e3.normalize();

				if (fabs(e1.dot(e2)) > 0.01 || fabs(e1.dot(e3)) > 0.01 || fabs(e2.dot(e3)) > 0.01) {
                    std::cout << e1.dot(e2) << " " << e1.dot(e3) << " " << e2.dot(e3) << std::endl;
                    std::cout << e1 << std::endl;
                    std::cout << e3 << std::endl;
                }

				///Get the response of this pixel to the filter
				float curve_response = response_images[i]->getPixel(x,y).r();

				///Get the so-far aggregated response
				Color val = response_intensity_map->getPixel(x,y);
				///Add the response of the gabor for this orientation
				val.g() += curve_response;
				///Set it back to the map
				response_intensity_map->setPixel(x,y,val);
#ifdef TENSOR_DEBUG
                if (curve_response > 1.0f || height_variation > 1.0f || normal_variation > 1.0f)    {
                    std::cout << curve_response << " " << height_variation << " " << normal_variation << std::endl;
                    exit(0);
                }
#endif
    			///Scale the eigenvalue
    			Vector3f curve_values = Vector3f(curve_response, height_variation, normal_variation);
    			if (curve_values(0) < 0.0f || curve_values(1) < 0.0f || curve_values(2) < 0.0f)   {
                    std::cout << "here3: " << curve_values << std::endl;
                    exit(0);
    			}
    			///Get the magnitude of the curve values
    			float c = curve_values.norm();
    			float m = sqrtf(3.0f);
                Vector3f eval = Vector3f((m-c+c*n)/n,
                                                                  c,
                                                                 0.0f);


                if (eval.norm() > TC_EPSILON)  eval.normalize();

				///Compute the tensor product
				///add it to the existing matrix
				tensor_matrix += Tensor::computeTensorMatrix(eval, e1, e2, e3);
			}
			///Average the pixels as tensors
			pixels_as_tensors.push_back(Tensor(tensor_matrix/float(response_images.size())));
			//pixels_as_tensors[pixels_as_tensors.size()-1].normalize();
			#ifdef EXPORT2R
            ///Export to R
			if (x%EXPORT_TENSOR_STEP==0 && y%EXPORT_TENSOR_STEP==0)    {
                exportTensor2R(pixels_as_tensors[pixels_as_tensors.size()-1],x,y,0.0f);
            }
            #endif
			///Keep track of the index
			tensor_index_map->setPixel(x,y,Color(float(pixels_as_tensors.size()-1)));
		}
	}
	std::cout << std::endl;
	///Normalize the response intensity map
	response_intensity_map->normalize();
	response_intensity_map->saveImage("response_intensity_map.pfm");

	///SECOND STEP IS TO DECOMPOSE THE RESULTING TENSORS AND USE THE EIGENVECTORS AND EIGENVALUES TO CLASSIFY EACH POINT
	///For each pixel and each resulting tensor compute the eigenvectors and eigenvalues

	for (int y=0;y<sizey;y++)      {
		for (int x=0;x<sizex;x++)      {
            int index = int(tensor_index_map->getPixel(x,y).r());
            if (index == -1) continue;

            ///GREEN: If it's a curve then store the saliency in the green channel
			if (pixels_as_tensors[index].getGeometricType() == Tensor::CURVE)	{
				Vector3f emin = pixels_as_tensors[index].getEigenVector(Tensor::T_EMIN);
				///Check if the X is positive
				if (emin(0) < 0.0f) emin = -emin;
				orientation_map->setPixel(x,y, Color(emin(0),emin(1),emin(2),1.0f));
				saliency_map->setPixel(x,y, Color(0.0f, pixels_as_tensors[index].getEigenValueDiffs(Tensor::L2MINUSL3), 0.0f, 1.0f));
			}
			///RED:If it's a surface then store the saliency in the red channel
			else	{
				if (pixels_as_tensors[index].getGeometricType() == Tensor::SURFACE)	{
					Vector3f emax = pixels_as_tensors[index].getEigenVector(Tensor::T_EMAX);
					///Check if the Z is pointing upwards
					if (emax(2) < 0.0f) emax = -emax;
					orientation_map->setPixel(x,y, Color(emax(0),emax(1),emax(2),1.0f));
                    saliency_map->setPixel(x,y, Color(pixels_as_tensors[index].getEigenValueDiffs(Tensor::L1MINUSL2),0.0f,0.0f,1.0f));
				}
				///BLUE: Otherwise it's a ball tensor, blue channel
				else	{
					if (pixels_as_tensors[index].getGeometricType() == Tensor::JOINT)	{
						orientation_map->setPixel(x,y, Color(0.0f,0.0f,0.0f,1.0f));
                        saliency_map->setPixel(x,y, Color(0.0f,0.0f,pixels_as_tensors[index].getEigenValueDiffs(Tensor::L3),1.0f));
					}
					else	{
						std::cout << "UNKNOWN type because the tensor matrix is zero." << std::endl;
						std::cout << x << " " << y << std::endl;
					}
				}
			}
		}
	}

	//Save the maps
	orientation_map->saveImage(_format("decompose_orientation_map.pfm").c_str());
	saliency_map->normalize();
	saliency_map->saveImage(_format("decompose_saliency_map.pfm").c_str());

	///Clean up
	delete response_intensity_map;
	return;
}

static Image *computeHeightAndNormalVariation(Image *xyz_map, Image *normal_map)	{
    int neighbourhood_size = 3;
	int neighbourhood_search_size = neighbourhood_size/2;
	///Create a map the same size as the xyz map
	Image *height_and_normal_variation_map = new Image(xyz_map->getWidth(), xyz_map->getHeight(),0.0f,0.0f,0.0f,1.0f);

	///Go through the xyz map and for each point compute the min max values in a nxn window area
	for (int y=0;y<xyz_map->getHeight();y++)	{
		for (int x=0;x<xyz_map->getWidth();x++)	{
			///Check only the valid points
			if ((xyz_map->getPixel(x,y) == Color(0.0f,0.0f,0.0f)) || (normal_map->getPixel(x,y) == Color(0.0f,0.0f,0.0f))) continue;
			///Get the height of the point in question
			float current_height = xyz_map->getPixel(x,y).b();

			///Get the normal of the point in question
			Vector3f normal = color2vector3(normal_map->getPixel(x,y));

			float min_dot = 1.0f;
			float max_dot = 0.0f;
            //Vector3f avg_normal = Vector3f(0.0f,0.0f,0.0f);
			//float avg_dot = 0.0f;
			//float num_dots = 0.0f;

            float min_height = xyz_map->getPixel(x,y).b();
            float max_height = xyz_map->getPixel(x,y).b();

			///Search in a neighbourhood
			for (int i=y-neighbourhood_search_size;i<=y+neighbourhood_search_size;i++){
				for (int j=x-neighbourhood_search_size;j<=x+neighbourhood_search_size;j++)	{
					///If it's the same continue
					if (i==y && j==x)	continue;
					///Check for out of bounds
					if (outOfBounds(xyz_map, j, i))	continue;
					///if its a valid point
					if (xyz_map->getPixel(j,i) == Color(0.0f,0.0f,0.0f) || (normal_map->getPixel(j,i) == Color(0.0f,0.0f,0.0f))) continue;

					///Check if max or min for the dot product
					float dot_product = std::max(0.0f,1.0f - fabs(normal.dot(color2vector3(normal_map->getPixel(j,i)))));
                    //Vector3f neighbour_normal = color2vector3(normal_map->getPixel(j,i));
                    //avg_normal += neighbour_normal;

                    //float dot_product = (normal - color2vector3(normal_map->getPixel(j,i))).norm();


                    if (dot_product < 0.0f)  {
                        ///Check if it's a negligible value
                        if (dot_product < -TC_EPSILON)  {
                            std::cerr << "Dot product has negative value: " << dot_product << std::endl;
                            exit(0);
                        }
                        else    {
                            dot_product = 0.0f;
                        }
                    }


					if (min_dot > dot_product)	min_dot = dot_product;
					if (max_dot < dot_product)	max_dot = dot_product;

                   //avg_dot += dot_product;
					//num_dots++;
                    ///Check if max or min for the depth
                    if (min_height > xyz_map->getPixel(j,i).b())	min_height = xyz_map->getPixel(j,i).b();
                    if (max_height < xyz_map->getPixel(j,i).b())	max_height = xyz_map->getPixel(j,i).b();
				}
			}

			///Compute the height variation
			float max_minus_min = max_height - min_height;
			if (max_minus_min < TC_EPSILON)    max_minus_min = 1.0f;

			float height_var = (current_height-min_height);///max_minus_min;
            if (height_var < 0.0f) {
                std::cout << "Problem: " << max_height << " " << min_height << " " << current_height << " " << height_var << std::endl;
                exit(0);
			}
			///Compute the normal variation
			//if (num_dots < TC_EPSILON)	num_dots = 1.0f;
			//float max_minus_min_dot = max_dot - min_dot;
			//if (max_minus_min_dot < TC_EPSILON)    max_minus_min_dot = 1.0f;
            //avg_normal.normalize();
			float normal_var = (max_dot - min_dot);//(normal - avg_normal).norm();///max_minus_min_dot;//avg_dot/num_dots;
			if (normal_var < 0.0f) {
                normal_var = 0.0f;
                //std::cout << "Problem: " << max_dot << " " << min_dot << " " << normal_var << std::endl;
                //exit(0);
			}
			///Store it in the variation map
			height_and_normal_variation_map->setPixel(x,y,Color(height_var,normal_var,0.0f));
		}
	}
    ///Normalize
    height_and_normal_variation_map->normalize();
    ///Save it out
    height_and_normal_variation_map->saveImage(_format("%s_height_and_dot_variation_map.pfm", name.c_str()));
/*
	///Go through the map again and compute the second derivative of differences
	for (int y=0;y<height_and_normal_variation_map->getHeight();y++)	{
		for (int x=0;x<height_and_normal_variation_map->getWidth();x++)	{
			///Check only the valid points
			if (fabs(height_and_normal_variation_map->getPixel(x,y).r()) < EPSILON &&
				fabs(height_and_normal_variation_map->getPixel(x,y).g()) < EPSILON)	continue;

			///Search in a neighbourhood
			float sum_of_differences_height = 0.0f;
			float sum_of_differences_normal= 0.0f;
			float neighbours = 0.0f;
			for (int i=y-neighbourhood_search_size;i<=y+neighbourhood_search_size;i++){
				for (int j=x-neighbourhood_search_size;j<=x+neighbourhood_search_size;j++)	{
					///If it's the same continue
					if (i==y && j==x)	continue;
					///Check for out of bounds
					if (outOfBounds(height_and_normal_variation_map, j, i))	continue;
					///if its a valid point
					///Check only the valid points
					if (fabs(height_and_normal_variation_map->getPixel(j,i).r()) < EPSILON &&
						fabs(height_and_normal_variation_map->getPixel(j,i).g()) < EPSILON)	continue;
					sum_of_differences_height += fabs(height_and_normal_variation_map->getPixel(x,y).r()-height_and_normal_variation_map->getPixel(j,i).r());
					sum_of_differences_normal+= fabs(height_and_normal_variation_map->getPixel(x,y).g()-height_and_normal_variation_map->getPixel(j,i).g());
					neighbours++;
				}
			}
			sum_of_differences_height /= std::max(1.0f,neighbours);
			sum_of_differences_normal /= std::max(1.0f,neighbours);

			Color val = height_and_normal_variation_map->getPixel(x,y);
			height_and_normal_variation_map->setPixel(x,y,Color(val.r(), val.g(), sum_of_differences_height, sum_of_differences_normal));
		}
	}

	if (TENSOR_DEBUG)	{
		Image *secondary_image = new Image(height_and_normal_variation_map->getWidth(), height_and_normal_variation_map->getHeight(), 0.0f,0.0f,0.0f,1.0f);
		for (int y=0;y<height_and_normal_variation_map->getHeight();y++)	{
			for (int x=0;x<height_and_normal_variation_map->getWidth();x++)	{
				if (fabs(height_and_normal_variation_map->getPixel(x,y).r()) < EPSILON &&
					fabs(height_and_normal_variation_map->getPixel(x,y).g()) < EPSILON)	continue;
				secondary_image->setPixel(x,y, Color(height_and_normal_variation_map->getPixel(x,y).b(), height_and_normal_variation_map->getPixel(x,y).a(),1.0f));
			}
		}
		secondary_image->saveImage(_format("%s_height_and_dot_variation_map_B.pfm", file_name.c_str()));
		delete secondary_image;
	}
*/
	return height_and_normal_variation_map;
}

Matrix3f computeTensorMatrix(Vector3f const &eval,Vector3f const &e1,Vector3f const &e2,Vector3f const &e3)    {
	///compute the left-hand side matrix
	Matrix3f lhs;
	lhs(0,0)=e1(0);  lhs(0,1)=e2(0);  lhs(0,2)=e3(0);
	lhs(1,0)=e1(1);  lhs(1,1)=e2(1);  lhs(1,2)=e3(1);
	lhs(2,0)=e1(2);  lhs(2,1)=e2(2);  lhs(2,2)=e3(2);

	///compute the lambda matrix
	Matrix3f lm;
	lm(0,0)=eval(0);      lm(0,1)=0.0f;         lm(0,2)=0.0f;
	lm(1,0)=0.0f;         lm(1,1)=eval(1);      lm(1,2)=0.0f;
	lm(2,0)=0.0f;         lm(2,1)=0.0f;         lm(2,2)=eval(2);

	///compute the right-hand side matrix
	Matrix3f rhs;
	rhs(0,0)=e1(0);  rhs(0,1)=e1(1);  rhs(0,2)=e1(2);
	rhs(1,0)=e2(0);  rhs(1,1)=e2(1);  rhs(1,2)=e2(2);
	rhs(2,0)=e3(0);  rhs(2,1)=e3(1);  rhs(2,2)=e3(2);

	///compute the tensor product
	Matrix3f tensor_matrix = lhs * (lm * rhs);

	return tensor_matrix;
}

Image *computePartialDerivativeX(Image *xyz_map)    {
    Image *dx = new Image(xyz_map->getWidth(), xyz_map->getHeight(), 0.0f,0.0f,0.0f,1.0f);
    for (int y=0;y<xyz_map->getHeight();y++)    {
        for (int x=0;x<xyz_map->getWidth();x++) {
            ///Check for out of bounds
            if (outOfBounds(xyz_map, x+1,y) || outOfBounds(xyz_map, x,y+1) ) continue;
            ///Compute the difference
            float diff = (fabs(xyz_map->getPixel(x+1,y).r() - xyz_map->getPixel(x,y).r()) +
                         fabs(xyz_map->getPixel(x,y+1).r() - xyz_map->getPixel(x,y).r()))*0.5f;
            dx->setPixel(x,y,Color(diff));
        }
    }
    return dx;
}

Image *computePartialDerivativeY(Image *xyz_map)    {
    Image *dy = new Image(xyz_map->getWidth(), xyz_map->getHeight(), 0.0f,0.0f,0.0f,1.0f);
    for (int y=0;y<xyz_map->getHeight();y++)    {
        for (int x=0;x<xyz_map->getWidth();x++) {
            ///Check for out of bounds
            if (outOfBounds(xyz_map, x+1,y) || outOfBounds(xyz_map, x,y+1) ) continue;
            ///Compute the difference
            float diff = (fabs(xyz_map->getPixel(x+1,y).g() - xyz_map->getPixel(x,y).g()) +
                         fabs(xyz_map->getPixel(x,y+1).g() - xyz_map->getPixel(x,y).g()))*0.5f;
            dy->setPixel(x,y,Color(diff));
        }
    }
    return dy;
}

Image *computePartialDerivativeZ(Image *xyz_map)    {
    Image *dz = new Image(xyz_map->getWidth(), xyz_map->getHeight(), 0.0f,0.0f,0.0f,1.0f);
    for (int y=0;y<xyz_map->getHeight();y++)    {
        for (int x=0;x<xyz_map->getWidth();x++) {
           ///Check for out of bounds
            if (outOfBounds(xyz_map, x+1,y) || outOfBounds(xyz_map, x,y+1) ) continue;
            ///Compute the difference
            float diff = (fabs(xyz_map->getPixel(x+1,y).b() - xyz_map->getPixel(x,y).b()) +
                         fabs(xyz_map->getPixel(x,y+1).b() - xyz_map->getPixel(x,y).b()))*0.5f;
             dz->setPixel(x,y,Color(diff));
        }
    }
    return dz;
}

void triangulatePoints(Image *xyz_map)  {
	///Triangulate the XYZ map to get an object
	GeometricObject *new_object = GeometryProcessing::triangulate(xyz_map);
	///Export to file
	GeometryExporter ge;
	ge.exportToOBJ("original_mesh", new_object);

	delete new_object;
	return;
	/*
    std::vector<Point> points;
    for (int y=0;y<xyz_map->getHeight();y++)    {
        for (int x=0;x<xyz_map->getWidth();x++) {
            points.push_back(Point(xyz_map->getPixel(x,y).r(), xyz_map->getPixel(x,y).g(), xyz_map->getPixel(x,y).b()));
        }
    }

    Delaunay triangulation(points.begin(), points.end());

    FILE *fileptr_mesh = fopen("original_mesh.obj","w");
    int vertex_pos = 1;
    for(Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin();fit != triangulation.finite_faces_end(); ++fit) {
        //Delaunay::Face_handle face = fit;
        fprintf(fileptr_mesh, "v %f %f %f\n", triangulation.triangle(fit)[0].x(), triangulation.triangle(fit)[0].y(), triangulation.triangle(fit)[0].z());
        fprintf(fileptr_mesh, "v %f %f %f\n", triangulation.triangle(fit)[1].x(), triangulation.triangle(fit)[1].y(), triangulation.triangle(fit)[1].z());
        fprintf(fileptr_mesh, "v %f %f %f\n", triangulation.triangle(fit)[2].x(), triangulation.triangle(fit)[2].y(), triangulation.triangle(fit)[2].z());
	  	fprintf(fileptr_mesh, "f %d %d %d\n", vertex_pos, vertex_pos+1, vertex_pos+2);
        vertex_pos += 3;
    }
    fclose(fileptr_mesh);
*/

    return;
}

bool createTensorField(Image *input_image, Image *height_and_normal_variation_map, Image *tensor_index_map, std::vector<Tensor> &pixels_as_tensors, Image *normal_map)    {

	///The size of the images
	int sizex = input_image->getWidth();
	int sizey = input_image->getHeight();


	Image *xyg_image = new Image();
	xyg_image->copy(input_image);
	Color *min_max_val = xyg_image->range();
	//xyg_image->normalize(min_max_val[1].b() - min_max_val[0].b());
	///Set the X and Y components to match the indices of the map
	///Mark the boundary points
    std::vector<Vector2i> boundary_points;
	for (int y=0;y<xyg_image->getHeight();y++)  {
        for (int x=0;x<xyg_image->getWidth();x++)   {
            Color col = xyg_image->getPixel(x,y);
            xyg_image->setPixel(x,y,col);//Color(float(x),float(y),col.b()));
            if (col.r() < TC_EPSILON && col.g() < TC_EPSILON) {
                boundary_points.push_back(Vector2i(x,y));
            }

        }
    }
	///Export the triangulation of the data
    triangulatePoints(xyg_image);
    /*
    GeometryExporter ge;
	GeometricObject *object = GeometryProcessing::triangulate(xyg_image, normal_map);
	ge.exportToOBJ("xyg_map", object);
	delete object;
    */
    ///FIRST PHASE: GABOR JETS
	Image *zzz_image = new Image();
	zzz_image->copy(input_image);
	zzz_image->grayscale(2);
	zzz_image->normalize();
	zzz_image->blur(1.0f);
	zzz_image->saveImage(_format("zzz_map.pfm").c_str());
	std::vector<Image *> response_images;
	applyGaborJets(zzz_image, boundary_points, response_images);
    delete zzz_image;


	///SECOND PHASE: ENCODE THE RESPONSE IMAGES AND NORMAL MAP TO TENSORS
	Image *saliency_map = new Image(sizex, sizey);
	Image *orientation_map = new Image(sizex, sizey);
	encodeResponses2Tensor(xyg_image, height_and_normal_variation_map, response_images, normal_map, saliency_map, orientation_map, pixels_as_tensors, tensor_index_map);
	///Clean up
	for (int i=0;i<response_images.size();i++)	{
		delete response_images[i];
	}
	response_images.clear();
	delete saliency_map;
	delete orientation_map;
	delete xyg_image;

    //fprintf(fileptr2,"%d %d\n", sizex, sizey);
    for (int y=0;y<sizey;y++)   {
        for (int x=0;x<sizex;x++)   {
            int index = int(tensor_index_map->getPixel(x,y).r());
            if (index < 0) {
                std::cout << index << std::endl;
                exit(0);
            }
            Tensor tensor = pixels_as_tensors[index];
            Vector3f eval = tensor.getEigenValues();
            Vector3f e1 = tensor.getEigenVector(Tensor::T_EMAX);
            Vector3f e2 = tensor.getEigenVector(Tensor::T_EMID);
            Vector3f e3 = tensor.getEigenVector(Tensor::T_EMIN);
            //fprintf(fileptr2,"%f %f %f %f %f %f %f %f %f %f %f %f\n", eval[0], eval[1], eval[2], e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], e3[0], e3[1], e3[2]);
        }
    }

    return true;
}

void tensorClustering(Image *xyz_map, Image *tensor_index_map, Image *patch_index_map, std::vector<Tensor> const &pixels_as_tensors)    {
	///The search area used for growing
	int search_area = 1;
	///The region indexA
	int patch_index = 0;
#ifdef EXPORT_DISTRIBUTION_ANIMATION
    int number_of_patch_indices = 4;
    std::vector<int> lookup_patch_indices;
    lookup_patch_indices.push_back(69);//71);
    lookup_patch_indices.push_back(95);//84);
    lookup_patch_indices.push_back(122);//130);
    lookup_patch_indices.push_back(83);//96);
    std::vector<FILE *> fileptrs;
    std::vector<int> lookup_image_positions;
    for (int i=0;i<number_of_patch_indices;i++) {
        fileptrs.push_back(fopen(_format("distribution_%d.r",i).c_str(),"w"));
        fprintf(fileptrs[fileptrs.size()-1],"library(fitdistrplus)\n");
        lookup_image_positions.push_back(0);
    }
#endif
    ///A map to hold the already processed points in the map
	Image *already_processed = new Image( xyz_map->getWidth(),xyz_map->getHeight(), 0.0f, 0.0f, 0.0f, 1.0f);

#ifdef TENSOR_DEBUG
	int patch_image_counter = 0;
	Image *patch_image = new Image(xyz_map->getWidth(),xyz_map->getHeight(), -1.0f, -1.0f, -1.0f, 1.0f);

    ///Color coded patched
    Image *color_patches = new Image( xyz_map->getWidth(),xyz_map->getHeight(), -1.0f, -1.0f, -1.0f, 1.0f);
    srand (time(NULL));
#endif // TENSOR_DEBUG

	///Begin growing regions
	int pos = 1;
	for ( int y=0;y<xyz_map->getHeight();y++ )	{
		for ( int x=0;x<xyz_map->getWidth();x++ )	{
            std::cout << "\r\033[1;32m[" << int(100.0 * (pos+1)/(xyz_map->getWidth()*xyz_map->getHeight())) << "%]\033[0m Clustering tensors";
            std::cout.flush();

			///if this point is valid
			if ( xyz_map->getPixel(x,y) == Color(0.0f,0.0f,0.0f))	continue;
            ///if its not processed
			if (already_processed->getPixel(x,y) == Color(1.0f,1.0f,1.0f))	continue;
            ///if the index to the tensor field is -1
			int tensor_index = int(tensor_index_map->getPixel(x,y).r());
   			if (tensor_index == -1) {
   			    std::cout << "ERR1: Negative tensor index " << tensor_index << std::endl;
                continue;
   			}

#ifdef TENSOR_DEBUG
            float r= ((double) rand() / (RAND_MAX));
            float g= ((double) rand() / (RAND_MAX));
            float b= ((double) rand() / (RAND_MAX));
#endif // TENSOR_DEBUG

			///get the tensor corresponding to this point
			Tensor tensor_matrix = pixels_as_tensors[tensor_index];

            ///get the point
            Vector3f point = color2vector3(xyz_map->getPixel(x,y));

			///get the point index
			Vector2i point_index = Vector2i(x,y);

			///create a region with this point added
			Patch *a_patch = new Patch(patch_index);
			///add the tensor corresponding to the initial point
			a_patch->update(point, point_index, tensor_matrix);
			///add the point to the region map
			patch_index_map->setPixel(x,y,Color(float(patch_index)));

#ifdef TENSOR_DEBUG
			color_patches->setPixel(x,y,Color(r,g,b));
			///add the point to the patch image
			patch_image->setPixel(x,y,Color(r,g,b));
#endif // TENSOR_DEBUG
			///create a processing list with this point on it
			std::vector<Vector2i> processing_list, marked_points;
			marked_points.push_back(point_index);
			processing_list.push_back ( point_index );
			///mark it as processed
			already_processed->setPixel(x,y, Color(1.0f,1.0f,1.0f));
			pos++;

			///start the processing
			while ( processing_list.size() )	{
			    ///get the current point
				Vector2i current_point_index = processing_list[processing_list.size()-1];
				///remove it from the list
				processing_list.pop_back();

				///check the local neighbourhood
				for (int i=current_point_index(1) - search_area;i<=current_point_index(1) + search_area;i++)	{
					for (int j=current_point_index(0) - search_area;j<=current_point_index(0) + search_area;j++)	{
						///if its the same
						if ( i==y && j==x )	continue;
						///if its inbound
						if (outOfBounds(xyz_map, j, i))	continue;
						///if its not processed
						if (already_processed->getPixel(j,i) == Color(1.0f,1.0f,1.0f))	continue;
						///if its a valid point
						if (xyz_map->getPixel(j,i) == Color(0.0f,0.0f,0.0f))	continue;
                        ///if its a valid index
                        int candidate_tensor_index = tensor_index_map->getPixel(j,i).r();
                        if (candidate_tensor_index == -1)    {
                            std::cout << "ERR2: Negative tensor index " << tensor_index << std::endl;
                            continue;
                        }
						///get the saliency and normal of this point
						Vector3f candidate_point = color2vector3(xyz_map->getPixel(j,i));

						///get the tensor corresponding to the candidate point
						Tensor candidate_tensor = pixels_as_tensors[candidate_tensor_index];

						///check is the distance of this position from the plane
						///and the normal difference is similar to the other points
						if (a_patch->isLikelyToBePartOf(candidate_point, candidate_tensor))	{
							Vector2i candidate_point_index = Vector2i(j,i);
							///add it to the list
							processing_list.push_back(candidate_point_index);
							///add it to the region
							a_patch->update(candidate_point, candidate_point_index, candidate_tensor);

							///mark it as processed
							already_processed->setPixel(j,i,Color(1.0f,1.0f,1.0f));
							pos++;
							///add it to the region map
							patch_index_map->setPixel(j,i,Color(float(patch_index)));
#ifdef TENSOR_DEBUG
							color_patches->setPixel(j,i,Color(r,g,b));
							///Save out an image of the region with the new point. Used for animations
							patch_image->setPixel(j,i,Color(r,g,b));
							//patch_image->saveImage(_format("region_patch_%.06d.png",patch_image_counter++));
#endif // TENSOR_DEBUG

#ifdef EXPORT_DISTRIBUTION_ANIMATION
                            for (int k=0;k<number_of_patch_indices;k++) {
                                if (lookup_patch_indices[k] == patch_index) {
                                    //GaussianDistribution *distr = a_patch->getPatchDistribution();
                                    WeibullDistribution *distr = a_patch->getPatchDistribution();
                                    std::vector<float> samples;
                                    distr->getSamples(samples);

                                    fprintf(fileptrs[k],"#shape <- %f\n", distr->getShape());
                                    fprintf(fileptrs[k],"#scale <- %f\n", distr->getScale());
                                    fprintf(fileptrs[k],"#mean <- %f\n", distr->getMean());
                                    fprintf(fileptrs[k],"#variance <- %f\n", distr->getVariance());
                                    fprintf(fileptrs[k], "x <- c(%f ", samples[0]);
                                    for (int l = 1; l < samples.size(); l++) {
                                        fprintf(fileptrs[k], ", %f", samples[l]);
                                    }
                                    fprintf(fileptrs[k], ")\n");
                                    fprintf(fileptrs[k], "weibull <- try(fitdist(x, distr=\"weibull\",method=\"mle\",lower=c(0,0)))\n");
                                    fprintf(fileptrs[k], "if (!inherits(weibull,\'try-error\')){\n");
                                    fprintf(fileptrs[k], _format("png(filename=\"distribution_animation/distr_%05d.png\")\n", lookup_image_positions[k]++).c_str());
                                    fprintf(fileptrs[k], "summary(weibull)\n");
                                    fprintf(fileptrs[k], "plot(weibull)\n");
                                    fprintf(fileptrs[k], "dev.off()\n");
                                    fprintf(fileptrs[k], "}\n");
                                }
                            }
#endif
							marked_points.push_back(candidate_point_index);
						}
					}
				}
			}

			///add the region to the resulting regions if it's stable otherwise just skip it
			if (a_patch->isStable())	{
				patches.push_back(a_patch);
				patch_index++;
                color_patches->saveImage("color_coded_patches.png");

			}
			else	{
				///Unmark any points that were marked
				for (int i=0;i<marked_points.size();i++)	{
					patch_index_map->setPixel(marked_points[i](0), marked_points[i](1), Color(-1.0f));

				#ifdef TENSOR_DEBUG
					color_patches->setPixel(marked_points[i](0), marked_points[i](1), Color(0.0f));
				#endif // DEBUG
				}
			}
		}
	}
    std::cout << std::endl;
    std::cout << "Number of patches: " << patches.size() << std::endl;
#ifdef TENSOR_DEBUG
/*
    for (int i=0;i<patches.size();i++)  {
        std::cout << patches[i]->getPatchDistribution()->getMeans()  << std::endl;
        std::cout << "-----------------------------------------" << std::endl;
    }
*/
    color_patches->saveImage("color_coded_patches.pfm");
    color_patches->saveImage("color_coded_patches.png");
#endif // TENSOR_DEBUG

	///Clean up
	delete already_processed;
#ifdef TENSOR_DEBUG
	delete color_patches;
	delete patch_image;
#endif // TENSOR_DEBUG
#ifdef EXPORT_DISTRIBUTION_ANIMATION
    for (int i=0;i<number_of_patch_indices;i++)	{
		fclose(fileptrs[i]);
	}
#endif // EXPORT_DISTRIBUTION_ANIMATION

	return;
}

void decomposeTensors(Image *xyz_map)	{
	Image *decomposed_tensors_saliency = new Image(xyz_map->getWidth(), xyz_map->getHeight(), 0.0f,0.0f,0.0f,1.0f);
	Image *decomposed_tensors_orientations = new Image(xyz_map->getWidth(), xyz_map->getHeight(), 0.0f,0.0f,0.0f,1.0f);

	///Go through each patch
	for (int i=0;i<patches.size();i++)	{
		///Get the patch point indices
		std::vector<Vector2i> patch_point_indices = patches[i]->getPointIndices();
		///Get the average tensor's type for this patch
		Tensor::TENSOR_GEOMETRIC_TYPE type = patches[i]->getGeometricType();
        Vector3f surface_normal = patches[i]->getAvgTensor().getEigenVector(Tensor::T_EMAX);
        Vector3f curve_tangent = patches[i]->getAvgTensor().getEigenVector(Tensor::T_EMIN);
		for (int j=0;j<patch_point_indices.size();j++) {
			Vector2i index = patch_point_indices[j];
			switch (type) {
				case Tensor::TENSOR_GEOMETRIC_TYPE::SURFACE:
					decomposed_tensors_saliency->setPixel(index(0),index(1), Color(1.0f,0.0f,0.0f));
                    decomposed_tensors_orientations->setPixel(index(0),index(1), vector2color3(surface_normal));
			        break;
                case Tensor::TENSOR_GEOMETRIC_TYPE::CURVE:
                    decomposed_tensors_saliency->setPixel(index(0),index(1), Color(0.0f,1.0f,0.0f));
                    decomposed_tensors_orientations->setPixel(index(0),index(1), vector2color3(curve_tangent));
                    break;
                case Tensor::TENSOR_GEOMETRIC_TYPE::JOINT:
                    decomposed_tensors_saliency->setPixel(index(0),index(1), Color(0.0f,0.0f,1.0f));
                    decomposed_tensors_orientations->setPixel(index(0),index(1), Color(0.0f,0.0f,0.0f));
                    break;
                case Tensor::TENSOR_GEOMETRIC_TYPE::UNKNOWN:
                    decomposed_tensors_saliency->setPixel(index(0),index(1), Color(1.0f,1.0f,1.0f,1.0f));
                    decomposed_tensors_orientations->setPixel(index(0),index(1), Color(-1.0f,-1.0f,-1.0f));
                    break;
            }
		}
	}


	decomposed_tensors_saliency->saveImage("decomposed_tensors_saliency.pfm");
	decomposed_tensors_orientations->saveImage("decomposed_tensors_orientations.pfm");
	delete decomposed_tensors_orientations;
	delete decomposed_tensors_saliency;
}

///NOTE:(0,0,0) HAS A SPECIAL MEANING IN THE INPUT DATA. IT SIGNIFIES THAT THERE IS
///NO POINT INFORMATION AT THAT PIXEL. ACTUAL DATA OF VALUE (0,0,0) SHOULD BE AVOIDED E.G. OFFSETTED.
int main(int argc, char *argv[])    {
    if (argc < 2)  {
        std::cout << "Number of arguments should be at least 2. Quiting..." << std::endl;
        return -1;
    }
    else   {
		FILE *fileptr_runtime_info = fopen("runtime.info","w");
		fprintf(fileptr_runtime_info, "kappa: %f\n", KAPPA);
		fprintf(fileptr_runtime_info, "Started @ %s\n", timestamp().c_str());
        std::cout << "kappa: " << KAPPA << std::endl;
		std::cout << "Started @ " << timestamp() << std::endl;
        ///If you reach this point then it means you have 2 arguments.
        ///The second argument should be the filename of the HDR XYZ map containing the geometry
        file_name = std::string(argv[1]);
        if (!separate(file_name, name, type)) return 0;

        ///Load the XYZ map
        Image *xyz_map = loadXYZmap();
        if (!xyz_map)   return 0;


        ///Triangulate and export the original mesh
        triangulatePoints(xyz_map);

        ///Compute the normal map corresponding to the xyz map
        Image *normal_map = GeometryProcessing::computeNormalMap(xyz_map, true);
        if (!normal_map)    return 0;
        normal_map->blur(1.0f);
        ///Normals are not normalized after the blur. Repeat normalize.
        normal_map->perPixelNormalization();
        normal_map->saveImage(_format("%s_normal_map.pfm", name.c_str()));

        ///Compute the height variation between the points
        Image *height_and_normal_variation_map = computeHeightAndNormalVariation(xyz_map, normal_map);
        if (!height_and_normal_variation_map)   return 0;

	    ///Compute a surface tensor for each point
	    #ifdef EXPORT2R
			fileptr_export2R = fopen("tensorfield.r","w");
            fprintf(fileptr_export2R,"library(rgl)\n\n");
            fprintf(fileptr_export2R, "c1 <-colors(TRUE)\n");
        #endif
	    ///Initialize the structure to keep a tensor for each point in the xyz_map
        Image *tensor_index_map = new Image(xyz_map->getWidth(), xyz_map->getHeight(), -1.0f,-1.0f,-1.0f,1.0f);

        //fileptr2 = fopen("tensorfield.xyz","w");
        std::vector<Tensor> pixels_as_tensors;
	    createTensorField(xyz_map, height_and_normal_variation_map, tensor_index_map, pixels_as_tensors,normal_map);
        //fclose(fileptr2);
        #ifdef EXPORT2R
        fclose(fileptr);
        #endif

        float max_dist = 0.0f;
        float min_dist = FLT_MAX;
        Image *tensor_distance_variance = new Image(xyz_map->getWidth(), xyz_map->getHeight(), 0.0f,0.0f,0.0f,1.0f);
        for (int y=0;y<xyz_map->getHeight();y++)    {
            for (int x=0;x<xyz_map->getWidth();x++) {
                float dist = 0.0f;
                float contributions = 0.0f;
                for (int i=-1;i<=1;i++) {
                    for (int j=-1;j<=1;j++) {
                        if (outOfBounds(tensor_distance_variance,x+j,y+i))  continue;
                        float this_dist = tensorComparison(pixels_as_tensors[int(tensor_index_map->getPixel(x,y).r())],pixels_as_tensors[int(tensor_index_map->getPixel(x+j,y+i).r())]);
                        if (this_dist > max_dist)   max_dist = this_dist;
                        if (this_dist < min_dist)   min_dist = this_dist;
                        dist += this_dist;
                        contributions++;
                    }
                }
                if (contributions > float(TC_EPSILON)) dist/= contributions;
                tensor_distance_variance->setPixel(x,y,Color(dist));
            }
        }
        tensor_distance_variance->normalize();
        tensor_distance_variance->saveImage("tensor_distance_variance.pfm");
        tensor_distance_variance->saveImage("tensor_distance_variance.png");
        delete tensor_distance_variance;
        //std::cout << "Max dist: " << max_dist << std::endl;
        //std::cout << "Min dist: " << min_dist << std::endl;

        //std::cout << "Performing tensor clustering..." << std::endl;
        ///Perform tensor clustering
        Image *patch_index_map = new Image(xyz_map->getWidth(), xyz_map->getHeight(), 0.0f,0.0f,0.0f,1.0f);
        tensorClustering(xyz_map, tensor_index_map, patch_index_map, pixels_as_tensors);
        patch_index_map->saveImage("patch_index_map.pfm");
        //std::cout << "...done" << std::endl;

		///Decompose into surfaces, curves and junctions
		decomposeTensors(xyz_map);

		//Perform hole filling on the patch indices
		//ImageProcessing::holeFill(patch_index_map,Color(-1.0f,-1.0f,-1.0f),6,2,true);
		//patch_index_map->saveImage("patch_index_map_hf.pfm");


		std::cout << "Finished @ " << timestamp() << std::endl;
		fprintf(fileptr_runtime_info, "Finished @ %s\n", timestamp().c_str());
		fclose(fileptr_runtime_info);
        ///Clean up everything
//        delete refined_depth_map;
        //delete refined_object;
        //delete ge;
        delete xyz_map;
        delete normal_map;
        delete height_and_normal_variation_map;
        delete tensor_index_map;
//        delete patch_index_map;
    }



    return 0;
}



