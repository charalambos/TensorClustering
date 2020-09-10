////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////
#ifndef __PATCH_H__
#define __PATCH_H__

#define USE_WEIBULL

#define	KAPPA	    							    1.0f
#define MINIMUM_SAMPLES_TO_PAUSE_RECALCULATION      150
#define MINIMUM_STABLE                              5
#define MAXIMUM_SAMPLES_TO_FORCE_RECALCULATION      100
#define RATIO_FACTOR                                0.95f

#include "Tensor.h"
#include "Utilities.h"
#include "BoundingBox.h"
//#include "GeometricObject.h"
//#include "GeometryExporter.h"
#include "GeometryProcessing.h"
#ifdef USE_WEIBULL
    #include "WeibullDistribution.h"
#else
    #include "GaussianDistribution.h"
#endif // USE_WEIBULL
//#include "DelaunayTriangulation.h"
#include "GeospatialBoundingBox.h"
#include "Contour.h"
//#include "GaussianMixtureDistribution.h"

#define PATCH_EPSILON                                     0.0001f
///Forward declaration
class Surface;

class Patch	{
	public:
		///Constructor
		Patch(int _cluster_id);
		///Destructor
		~Patch();

		///Registers the geospatial bounding box the surface belongs to
		void belongsTo(GeospatialBoundingBox *_geo_box);

		///Get the object the patch belongs to
		GeometricObject *getObject();

		///Returns the point indices
		std::vector<Vector2i> getPointIndices();

		///Returns the points
		std::vector<Vector3f> getPoints();

		///Get the patch's cluster id
		int getClusterId();

		///Registers the surface the patch belongs to
		void belongsTo(Surface *_surface);

		///Returns the belonging surface
		Surface *belongsTo();

		///Returns the centroid of the inlier points
		//Vector3f getInlierCentroid();

        ///Returns the centroid
        Vector3f getCentroid();

		///Returns the bounding box of the patch
		BoundingBox *getBoundingBox();

		///These function only affects the first distribution function to update the region with a new point and update the distributions
		void update(Vector3f new_point, Vector2i new_point_index, Tensor const &new_tensor);

		///Function to check if a point is likely to be part of this region
		bool isLikelyToBePartOf(Vector3f const &point, Tensor const &candidate_tensor);

		///Returns if the patch is stable
		bool isStable();

	/*
		///Set the boundary indices
		void setBoundaryIndices(std::vector<Vector2i> const &_boundary_indices);

		///Set the boundary points
		void setBoundaryPoints(std::vector<Vector3f> const &_boundary_points);

		void getBoundaryIndices(std::vector<Vector3f> &_boundary_points);

		void getBoundaryPoints(std::vector<Vector2i> &_boundary_indices);
*/
		Tensor::TENSOR_GEOMETRIC_TYPE getGeometricType();

		Tensor getAvgTensor();

		///Set the boundaries
		void setBoundaryAt(int index, Contour const &contour);
		void setBoundaries(std::vector<Contour> const &contour);

		///Get the boundaries
		Contour &getBoundaryAt(int index);
		std::vector<Contour> getBoundaries();


		///Used for reading in the patch information from a file
		void setPoints(std::vector<Vector3f> const &_points, bool erase_previous);
		void setPointIndices(std::vector<Vector2i> const &_point_indices, bool erase_previous);

		///Set/get the fitted plane
		void setPlane(Plane3d const &_plane3d, double _fit_quality);
		Plane3d getPlane();
		double getPlaneFitQuality();
		bool hasInvalidFit();

		///Marked for removal
		void setMarkedForRemoval(bool flag, int _substitute_id);
		bool getMarkedForRemoval();
		int getSubstituteId();

		///Returns the orientation of the patch
		float getOrientation();
        ///Returns the ratio of eigen_values
        void getEigenValues(float &, float &);
        
        ///Returns size or area of patch 
        int patchSize();

		///Contains point
		bool contains(Vector2f);

	///Returns the Weibull Distribution
#ifdef USE_WEIBULL
		WeibullDistribution *getPatchDistribution();
#else
		GaussianDistribution1f *getPatchDistribution();
#endif // USE_WEIBULL

	Color class_label;

	private:
		///The distributions of the patch descriptors
#ifdef USE_WEIBULL
		WeibullDistribution *patch_descriptors_distr;
#else
		GaussianDistribution1f *patch_descriptors_distr;
#endif // USE_WEIBULL

		///The fitting errors
		std::vector<float> fitting_errors;
		///The point indices of the patch
		std::vector<Vector2i> point_indices;
		///The boundary
		std::vector<Contour> boundaries;
		///The 3D points
		std::vector<Vector3f> points;
		///The average height
		float avg_height, sum_height;
		///The 3D normals
		std::vector<Vector3f> normals;
		///All the tensors
		std::vector<Tensor> tensors;
		///The previous variance,mean
		float previous_variance, previous_mean;
		///The previous-previous variance, mean
		float previous_previous_variance, previous_previous_mean;
        ///Counter
        int since_last_recalculation;


        ///The running sum of tensors and average tensor
        Tensor avg_tensor, sum_tensor;

		///The geospatial bounding box the linear surface belongs to
		GeospatialBoundingBox *geo_box;
		///The id of the corresponding cluster
		int cluster_id;
		///The surface this patch belongs to
		Surface *surface;
		///The bounding box of the patch
		BoundingBox *bbox;

		///The minimum and maximmum allowed probability which determines if a sample should be part of the distribution or not
		float min_probability, max_allowed_probability;

		///Defines if the patch is stable
		bool stable;

		///The fitted plane
		Plane3d plane;
		double fit_quality;
		bool invalid_fit;

		///Marked for removal
		bool marked_for_removal;
		int substitute_id;
	
		///PCA analysis
		void pcaData();
	
		///Centroid
		Vector3f centroid;
		
		///Eigen Values
		std::vector<Vector2f> eigen_vectors;
		std::vector<double> eigen_values;

};



#endif

