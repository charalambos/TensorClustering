////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////


#ifndef __GEOMETRIC_OBJECT_H__
#define __GEOMETRIC_OBJECT_H__


/**
*	Geometric Object
* This class defines a geometric object with all it's properties
**/
#include <math.h>
#include <iostream>

#include "BoundingBox.h"
#include "Edge.h"
#include "Face.h"
#include <Eigen/Eigen>
using namespace Eigen;

class GeometricObject {
	public:
		///Constructors
		GeometricObject();
		GeometricObject(std::vector<Vector3f> const &_vertices,
				std::vector<Vector3f> const &_normals,
				std::vector<Vector2f> const &_tex_coords,
				std::vector<Face *> const &_faces,
				std::vector<Edge *> const &_edges);

		///Copy constructor
		GeometricObject(GeometricObject const &_object);
		///Destructor
		~GeometricObject();

		///Get methods
		std::vector<Vector3f> getVertices() const ;
		Vector3f getVertexAt(int index) const ;

		std::vector<Vector3f> getNormals() const ;
		Vector3f getNormalAt(int index) const ;

		std::vector<Vector2f> getTextureCoords() const ;
		Vector2f getTextureCoordAt(int index) const ;

		std::vector<Face *> getFaces() const ;
		Face *getFaceAt(int index) const ;

		std::vector<Edge *> getEdges() const ;
		Edge *getEdgeAt(int index) const;

		void setVertices(std::vector<Vector3f> const _vertices, bool erase = true);
		void setNormals(std::vector<Vector3f> const _normals, bool erase = true);
		void setTexCoords(std::vector<Vector2f> const _tex_coords, bool erase = true);
		void addFace(Face *_face);

		int addVertex(Vector3f const &vertex);

		///Get functions for the extrinsic parameters
		Vector3f getTranslation();
		Vector3f getScale();
		Vector3f getRotation();

        ///Returns the name of the object
        std::string getName();

		///Rotate
		bool rotate();
		bool rotate(Vector3f const &rotation);
		bool rotate(Quaternionf const &rotation);
		bool rotate(Matrix3f const &rotation_matrix);

		///Scale
		bool scale();
		bool scale(Vector3f const &scale);

		///Translate
		bool translate();
		bool translate(Vector3f const &translation);

		///Apply all transformations
		bool applyTransformations();

		///Static functions
		static GeometricObject *grid(const float length, const float width, const int nbSub);
		static GeometricObject *axes(const float length);
        static GeometricObject *camera();

		///Return the object name
		std::string getObjectName();

		///Returns the bounding box
		BoundingBox getBBox();
		///Returns the bbox centroid
		Vector3f getBBoxCentroid();
		///Initialize bounding box
		bool setupBBox();
		///Returns the diagonal of bbox
		Vector3f getBBoxDiagonal();

private:
		///A list of vertices
		std::vector<Vector3f> vertices;
		///A list of normals
		std::vector<Vector3f> normals;
		///A list of texture coordinates
		std::vector<Vector2f> texture_coords;
		///A list of edges
		std::vector<Edge *> edges;
		///A list of faces
		std::vector<Face *> faces;

		///Bounding box
		BoundingBox bbox;

        void setGeometry();

		///Flag to indicate if VBO contents changes
		bool update;
        ///Count of indices
        int number_of_indices;

		///The rotation and translation of the object
		Vector3f translation_vector;
		Quaternionf rotation;
		Vector3f scale_vector;

        ///The name of the object
        std::string name;

        ///A global object id counter
        static long obj_id;

        bool set_geometry;
};

#endif
