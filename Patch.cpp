////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////
#ifndef __PATCH_CPP__
#define __PATCH_CPP__

#include <Contour.h>
#include <opencv2/ml.hpp>
#include "Patch.h"

#define WD_OFFSET   1.0f
#define WD_SCALE    1.0f

Patch::Patch(int _cluster_id)	{
#ifdef USE_WEIBULL
    patch_descriptors_distr = new WeibullDistribution();
#else
    patch_descriptors_distr = new GaussianDistribution();
#endif // USE_WEIBULL
    invalid_fit = false;
    min_probability = 0.0f;
    max_allowed_probability = 0.0f;
	geo_box = 0x00;
	cluster_id = _cluster_id;
	surface = 0x00;
	bbox = new BoundingBox();
	stable = false;
    sum_height = 0.0f;
	avg_height = 0.0f;
    previous_mean = previous_variance = 0.0f;
    previous_previous_mean = previous_previous_variance = 0.0f;
    since_last_recalculation = 0;

    fit_quality = -1.0;
    marked_for_removal = false;
    substitute_id = -1;
}

Patch::~Patch()	{
	if (patch_descriptors_distr)    {
        delete patch_descriptors_distr;
	}

	geo_box = 0x00;
	if (bbox)	delete bbox;
}

void Patch::belongsTo(GeospatialBoundingBox *_geo_box)	{
	geo_box = _geo_box;
	return;
}

GeometricObject *Patch::getObject()	{
	return geo_box->getObject();
}

BoundingBox *Patch::getBoundingBox()	{
	return bbox;
}

std::vector<Vector2i> Patch::getPointIndices()	{
	return point_indices;
}

int Patch::getClusterId()	{
	return cluster_id;
}

void Patch::belongsTo(Surface *_surface)	{
	surface = _surface;
	return;
}

Surface *Patch::belongsTo()	{
	return surface;
}

void Patch::update(Vector3f new_point, Vector2i new_point_index, Tensor const &new_tensor)   {
    //if (points.size()%500 == 0)     std::cout << points.size() << " " << timestamp() << std::endl;
	///Add the point, point index and tensor
	points.push_back(new_point);
	point_indices.push_back(new_point_index);
	tensors.push_back(new_tensor);
    sum_tensor.add(new_tensor);
    avg_tensor = Tensor(sum_tensor.getTensorInMatrixForm()/float(tensors.size()));

    if (!stable && point_indices.size() >= MINIMUM_STABLE)	{
		stable = true;
	}
    ///This factor controls the variance at the beginning stages so that it doesn't become zero
    float factor = 0.0f;//1.0f/std::max(1.0f,(float) log(points.size()*points.size()));
	sum_height += new_point(2);
	avg_height = sum_height/float(points.size());

    float dist = tensorComparison(avg_tensor, new_tensor);
//std::cout << dist << std::endl;
if (dist < 0.0f)    {
    std::cout << "dist is neg " << dist << std::endl;
    exit(0);
}

    bool recalculate = (points.size() < MINIMUM_SAMPLES_TO_PAUSE_RECALCULATION) ||
                        (previous_mean/previous_previous_mean < RATIO_FACTOR &&
                         previous_variance/previous_previous_variance < RATIO_FACTOR);

    if (!recalculate)   {
        since_last_recalculation++;
        ///Check how long it's been since last recalculation
        ///If it's been a while then force the recalculation
        if (since_last_recalculation >= MAXIMUM_SAMPLES_TO_FORCE_RECALCULATION)  {
            since_last_recalculation = 0;
            recalculate = true;
        }
    }

    ///Fit the samples to the WDs
    float sample = dist;
    patch_descriptors_distr->update((sample + WD_OFFSET)*WD_SCALE, recalculate);

    float mean = patch_descriptors_distr->getMean();
    previous_previous_mean = previous_mean;
    previous_mean = mean;

    ///The minimum probability will be that of the means
    min_probability = patch_descriptors_distr->prob(mean);
    ///The maximum allowed probability will be that of the mean normal + or -  a factor = K*variance
    float variance= patch_descriptors_distr->getVariance();
    previous_previous_variance = previous_variance;
    previous_variance = variance;
//std::cout << mean << " " << variance << std::endl;

if (variance < -TC_EPSILON|| isnan(variance)) {
    std::cout << "problematic variance: " << variance<< std::endl;
    std::cout << "problematic mean: " << mean << std::endl;
    std::cout << "number of samples: " << tensors.size() << std::endl;
    std::cout << "dist: " << dist << std::endl;
    std::cout << "sample: " << sample << std::endl;
    //std::cout << "alpha: " << patch_descriptors_distr->getShape() << std::endl;
    //std::cout << "beta: " << patch_descriptors_distr->getScale() << std::endl;
    std::vector<float> samples;
    patch_descriptors_distr->getSamples(samples);
    for (int i=0;i<samples.size();i++)  {
        std::cout << samples[i] << std::endl;
    }
    exit(0);
}

    max_allowed_probability = patch_descriptors_distr->prob(mean + KAPPA*sqrtf(variance));
/*
std::cout << mean << std::endl;
std::cout << sqrtf(variance)<< std::endl;
std::cout << patch_descriptors_distr->prob(mean) << " " << patch_descriptors_distr->prob(mean + KAPPA*sqrtf(variance))<< " " << patch_descriptors_distr->evaluate(mean + KAPPA*sqrtf(variance)) << std::endl;
*/
    Vector3f min_pt = bbox->getMinPt();
    Vector3f max_pt = bbox->getMaxPt();
    ///Check if max or min
    if (max_pt(0) < new_point(0))	max_pt(0) = new_point(0);
    if (max_pt(1) < new_point(1))	max_pt(1) = new_point(1);
    if (max_pt(2) < new_point(2))	max_pt(2) = new_point(2);

    if (min_pt(0) > new_point(0))	min_pt(0) = new_point(0);
    if (min_pt(1) > new_point(1))	min_pt(1) = new_point(1);
    if (min_pt(2) > new_point(2))	min_pt(2) = new_point(2);

    bbox->setMinPt(min_pt);
    bbox->setMaxPt(max_pt);


	return;
}

bool Patch::isLikelyToBePartOf(Vector3f const &point, Tensor const &candidate_tensor)	{
    ///If this is the first point of the region then add it.
	if ( points.size() == 0)	return true;

	///Otherwise we need to compute the probability of the new sample (tensor) of being part of this patch
	///Compute the probability of the candidate point's information
	float dist = tensorComparison(avg_tensor, candidate_tensor);
	//std::cout << "dist: " << dist << std::endl;
	float sample = dist;
	float prob = patch_descriptors_distr->prob((sample + WD_OFFSET)*WD_SCALE);
/*
std::cout << "-------------------" << std::endl;
std::cout << min_probability[0] << " " << "p0: " << prob_e1 << " " << max_allowed_probability[0] << std::endl;
std::cout << min_probability[1] << " " << "p1: " << prob_e2 << " " << max_allowed_probability[1] << std::endl;
std::cout << min_probability[2] << " " << "p2: " << prob_e3 << " " << max_allowed_probability[2] << std::endl;
std::cout << min_probability[3] << " " << "p3: " << prob_l << " " << max_allowed_probability[3] << std::endl;
*/
	///If the probability is within the acceptable range then add the point to the patch i.e. return true
	if (prob <= max_allowed_probability)	{
		return true;
	}

	///Otherwise its not part of the patch
	return false;
}

bool Patch::isStable()	{
	return stable;
}

std::vector<Vector3f> Patch::getPoints()	{
	return points;
}

#ifdef USE_WEIBULL
WeibullDistribution *Patch::getPatchDistribution()	{
	return patch_descriptors_distr;
}
#else
GaussianDistribution *Patch::getPatchDistribution()	{
	return patch_descriptors_distr;
}
#endif // USE_WEIBULL

/*
///Set the boundary indices
void Patch::setBoundaryIndices(std::vector<Vector2i> const &_boundary_indices)   {
    boundary_indices.clear();
    for (int i=0;i<_boundary_indices.size();i++)    {
        boundary_indices.push_back(_boundary_indices[i]);
    }
    return;
}

///Set the boundary points
void Patch::setBoundaryPoints(std::vector<Vector3f> const &_boundary_points) {
    boundary_points.clear();
    for (int i =0;i<_boundary_points.size();i++)    {
        boundary_points.push_back(_boundary_points[i]);
    }
    return;
}


void Patch::getBoundaryIndices(std::vector<Vector3f> &_boundary_points) {
    for (size_t i=0;i<boundary_points.size();i++)   {
        _boundary_points.push_back(boundary_points[i]);
    }
    return;
}

void Patch::getBoundaryPoints(std::vector<Vector2i> &_boundary_indices) {
    for (size_t i=0;i<boundary_indices.size();i++)  {
        _boundary_indices.push_back(boundary_indices[i]);
    }
return;
}
*/

Tensor::TENSOR_GEOMETRIC_TYPE Patch::getGeometricType() {
    Tensor norm_avg_tensor = avg_tensor;
    norm_avg_tensor.normalize();
    return norm_avg_tensor.getGeometricType();
}

Tensor Patch::getAvgTensor()    {
    return avg_tensor;
}

Contour &Patch::getBoundaryAt(int index)  {
    if (index < boundaries.size())  {
        return boundaries[index];
    }
    else    {
        std::cout << "Error: Patch::getBoundaryAt" << std::endl;
    }
}

std::vector<Contour> Patch::getBoundaries() {
    return boundaries;
}

void Patch::setPoints(std::vector<Vector3f> const &_points, bool erase_previous) {
    if (erase_previous) points.clear();
    for (int i=0;i<_points.size();i++)  {
        points.push_back(_points[i]);
    }
    pcaData();
    return;
}

void Patch::setPointIndices(std::vector<Vector2i> const &_point_indices, bool erase_previous)    {
    if (erase_previous) point_indices.clear();
    for (int i=0;i<_point_indices.size();i++)   {
        point_indices.push_back(_point_indices[i]);
    }
    return;
}

void Patch::setPlane(Plane3d const &_plane3d, double _fit_quality) {
    plane = _plane3d;
    fit_quality = _fit_quality;
    invalid_fit = std::isnan(fit_quality) || std::isinf(fit_quality);// || !(fit_quality >= 0.5f && fit_quality <= 1.0f);
    return;
}

Plane3d Patch::getPlane()  {
    if (invalid_fit)  {
        std::cerr << "Warning: No plane has been set for this patch." << std::endl;
    }
    return plane;
}

bool Patch::hasInvalidFit()    {
    return invalid_fit;
}

double Patch::getPlaneFitQuality()  {
    if (fit_quality < 0.0)  {
        std::cerr << "Warning: No plane has been set for this patch." << std::endl;
    }
    return fit_quality;
}

void Patch::setMarkedForRemoval(bool flag, int _substitute_id)  {
    marked_for_removal = flag;
    substitute_id = _substitute_id;
    return;
}

bool Patch::getMarkedForRemoval()   {
    return marked_for_removal;
}

int Patch::getSubstituteId()    {
    return substitute_id;
}

void Patch::setBoundaryAt(int index, Contour const &_contour)    {
    if (index < boundaries.size())  {
        boundaries[index] = Contour(_contour);
    }
    else    {
        std::cout << "Error: Patch::setBoundaryAt()" << std::endl;
    }
    return;
}

void Patch::setBoundaries(std::vector<Contour> const &_contours)    {
    boundaries.clear();
    for (int i=0;i<_contours.size();i++)    {
        boundaries.push_back(Contour(_contours[i]));
    }
    return;
}

float Patch::getOrientation() {
    return atan2(eigen_vectors[0](1), eigen_vectors[0](0));
}

void Patch::getEigenValues(float &val0, float &val1) {
    val1 = eigen_values[1];
    val0 = eigen_values[0];
}

void Patch::pcaData(){
    Vector3f centroid_temp = Vector3f(0.0f,0.0f,0.0f);
    
    cv::Mat data = cv::Mat(points.size(), 2, CV_64FC1);
    for (int i = 0; i < data.rows; i++) {
        centroid_temp += points[i];
        
        data.at<double>(i, 0) = points[i](0);
        data.at<double>(i, 1) = points[i](1);
    }
    cv::PCA pca(data, cv::Mat(), cv::PCA::DATA_AS_ROW);
    for(int i = 0; i < 2; i++){
        eigen_vectors.push_back(Vector2f(pca.eigenvectors.at<double>(i, 0), pca.eigenvectors.at<double>(i, 1)));
        eigen_values.push_back(pca.eigenvalues.at<double>(i));
    }
    centroid = centroid_temp/ std::max(1.0f,float(points.size()));
    
}
Vector3f Patch::getCentroid(){
    return centroid;
}

int Patch::patchSize(){
    return points.size();
}

bool Patch::contains(Vector2f p){
    for(int i = 0; i < points.size();i++){
    Vector2f point = Vector2f(points[i](0),points[i](1));
        if (fabs((p-point).norm()) < PATCH_EPSILON) {
//		if(p(0) == points[i](0) && p(1)==points[i](1)){
			return true;
		}
	}
	return false;
}

/*
Vector3f Patch::getCentroid()   {
    Vector3f centroid = Vector3f(0.0f,0.0f,0.0f);
    for (int i=0;i<points.size();i++)   {
        centroid += points[i];
    }
    centroid /= std::max(1.0f,float(points.size()));
    return centroid;
}*/
#endif

