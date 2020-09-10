////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////


#ifndef __GEOMETRIC_OBJECT_CPP__
#define __GEOMETRIC_OBJECT_CPP__


#include "GeometricObject.h"
#include "GeometryProcessing.h"

long GeometricObject::obj_id = 0;

GeometricObject::GeometricObject()	{
    name = _format("object_%.2d",obj_id++);

    set_geometry = false;
}

GeometricObject::GeometricObject(std::vector<Vector3f> const &_vertices,
				 std::vector<Vector3f> const &_normals,
				 std::vector<Vector2f> const &_tex_coords,
				 std::vector<Face *> const &_faces,
				 std::vector<Edge *> const &_edges):
				 vertices(_vertices),normals(_normals),texture_coords(_tex_coords),faces(_faces),edges(_edges)	{
    name = _format("object_%.2d",obj_id++);
	update = true;
	number_of_indices = 0;

	setupBBox();

    set_geometry = false;
}

GeometricObject::~GeometricObject()	{

	vertices.clear();
	normals.clear();
	texture_coords.clear();

	for (int i=0;i<edges.size();i++)	{
		delete edges[i];
	}
	edges.clear();

	for (int i=0;i<faces.size();i++)	{
		delete faces[i];
	}
	faces.clear();

	if (set_geometry) {
        // glDeleteVertexArrays(1, &VAO);
        // glDeleteBuffers(1, &vertex_VBO);
        // glDeleteBuffers(1, &normal_VBO);
        // glDeleteBuffers(1, &tex_coord_VBO);
        // glDeleteBuffers(1, &EBO);
    }
}

std::vector<Vector3f> GeometricObject::getVertices() const	{
	return vertices;
}

Vector3f GeometricObject::getVertexAt(int index) const	{

	if (index>=0 && index<vertices.size())	{
		return vertices[index];
	}
	else	{
		std::cout << "here 4b: " << index << " out of " << vertices.size() << std::endl;
	}
	return Vector3f();
}

std::vector<Vector3f> GeometricObject::getNormals() const	{

	return normals;
}

Vector3f GeometricObject::getNormalAt(int index) const	{

	if (index>=0 && index<normals.size())	{
		return normals[index];
	}
	return Vector3f();
}

std::vector<Vector2f> GeometricObject::getTextureCoords() const	{

	return texture_coords;
}

Vector2f GeometricObject::getTextureCoordAt(int index) const	{

	if (index>=0 && index<texture_coords.size())	{
		return texture_coords[index];
	}
	return Vector2f();
}

std::vector<Face *> GeometricObject::getFaces() const	{

	return faces;
}

Face *GeometricObject::getFaceAt(int index) const	{

	if (index>=0 && index<faces.size())	{
		return faces[index];
	}
	return 0x00;
}

std::vector<Edge *> GeometricObject::getEdges() const	{

	return edges;
}

Edge *GeometricObject::getEdgeAt(int index) const	{

	if (index>=0 && index<edges.size())	{
		return edges[index];
	}
	return 0x00;
}

void GeometricObject::setVertices(std::vector<Vector3f> const _vertices, bool erase)	{
	if (erase)	vertices.clear();
	for (int i=0;i<_vertices.size();i++)	{
		vertices.push_back(_vertices[i]);
	}
	update = true;

	setupBBox();
	return;
}

void GeometricObject::setNormals(std::vector<Vector3f> const _normals, bool erase)	{
	if (erase)	normals.clear();
	for (int i=0;i<_normals.size();i++)	{
		normals.push_back(_normals[i]);
	}
	update = true;
	return;
}

void GeometricObject::setTexCoords(std::vector<Vector2f> const _tex_coords, bool erase)	{
	if (erase)	texture_coords.clear();
	for (int i=0;i<_tex_coords.size();i++)	{
        texture_coords.push_back(_tex_coords[i]);
	}
	update = true;
	return;
}

int GeometricObject::addVertex(Vector3f const &vertex)  {
    ///search for it
    int index = -1;
    for (int i=0;i<vertices.size();i++) {
        if (vertices[i] == vertex)  {
            ///If it exists then return the index
            return i;
        }
    }

    ///If it's not there then add it an return the new index
    vertices.push_back(vertex);
	update = true;
	setupBBox();
    return vertices.size()-1;
}


void GeometricObject::addFace(Face *_face)	{
	faces.push_back(_face);
	update = true;
	return;
}

void GeometricObject::setGeometry()	{

    if (!set_geometry)  {
        // glGenVertexArrays(1, &VAO);
        // glGenBuffers(1, &vertex_VBO);
        // glGenBuffers(1, &normal_VBO);
        // glGenBuffers(1, &tex_coord_VBO);
        // glGenBuffers(1, &EBO);

        set_geometry = true;
    }

	// glBindVertexArray(VAO);

 //    ///Vertex attribute
	// glBindBuffer(GL_ARRAY_BUFFER, vertex_VBO);
	// glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(Vector3f), &vertices.front(), GL_STATIC_DRAW);
	// GLuint vertex_id = glGetAttribLocation(shader_program->getProgramObject(), "in_vertex");
	// glVertexAttribPointer(vertex_id, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
	// glEnableVertexAttribArray(vertex_id);

 //    ///Normal attribute
 //    glBindBuffer(GL_ARRAY_BUFFER, normal_VBO);
 //    glBufferData(GL_ARRAY_BUFFER, normals.size()*sizeof(Vector3f), &normals.front(), GL_STATIC_DRAW);
 //    GLuint normal_id = glGetAttribLocation(shader_program->getProgramObject(), "in_normal");
 //    glVertexAttribPointer(normal_id, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
 //    glEnableVertexAttribArray(normal_id);

 //    ///Texture coordinate attribute
 //    glBindBuffer(GL_ARRAY_BUFFER, tex_coord_VBO);
 //    glBufferData(GL_ARRAY_BUFFER, texture_coords.size()*sizeof(Vector2f), &texture_coords.front(), GL_STATIC_DRAW);
 //    GLuint tex_coord_id = glGetAttribLocation(shader_program->getProgramObject(), "in_texture_coord");
 //    glVertexAttribPointer(tex_coord_id, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), (GLvoid*)0);
 //    glEnableVertexAttribArray(tex_coord_id);

	std::vector<unsigned int> indices;
	for (int i=0;i<faces.size();i++)	{
        std::vector<int> face_indices = faces[i]->getVertices();
        indices.insert(indices.end(),face_indices.begin(),face_indices.end());
	}
    number_of_indices = indices.size();

	// glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	// glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices.front(), GL_STATIC_DRAW);

 //    glBindVertexArray(0); // Unbind VAO

    update = false;
    
	return;
}

///Rotate
bool GeometricObject::rotate(Vector3f const &_rotation)	{
	//change the rotation angles to the quaternion
	rotation = Quaternionf(AngleAxisf(_rotation(0), Vector3f::UnitX()) * AngleAxisf(_rotation(1), Vector3f::UnitY()) * AngleAxisf(_rotation(2), Vector3f::UnitZ()));
	return true;
}

bool GeometricObject::rotate(Quaternionf const &_rotation)    {
    rotation = _rotation;
    return true;
}

bool GeometricObject::rotate(Matrix3f const &rotation_matrix)	{
	rotation = Quaternionf(rotation_matrix);
	return true;
}

bool GeometricObject::rotate()	{
	for (int i=0;i<vertices.size();i++)	{
		vertices[i] = rotation * vertices[i];
	}
	rotation= Quaternionf(1.0f,0.0f,0.0f,0.0f);
	setupBBox();
	return true;
}

///Scale
bool GeometricObject::scale(Vector3f const &_scale)	{
	scale_vector = _scale;
	return true;
}

bool GeometricObject::scale()	{
	for (int i=0;i<vertices.size();i++)	{
		vertices[i](0) = scale_vector(0) * vertices[i](0);
		vertices[i](1) = scale_vector(1) * vertices[i](1);
		vertices[i](2) = scale_vector(2) * vertices[i](2);
	}
	scale_vector = Vector3f(1.0f,1.0f,1.0f);
	setupBBox();
	return true;
}

///Translate
bool GeometricObject::translate()	{
	for (int i=0;i<vertices.size();i++)	{
		vertices[i] = vertices[i] + translation_vector;
	}
	translation_vector = Vector3f(0.0f,0.0f,0.0f);
	setupBBox();
	return true;
}
bool GeometricObject::translate(Vector3f const &translation)	{
	translation_vector = translation;
	return true;
}

Vector3f GeometricObject::getTranslation()	{
	return translation_vector;
}

Vector3f GeometricObject::getScale()	{
	return scale_vector;
}

Vector3f GeometricObject::getRotation()	{
	Vector3f angles = rotation.toRotationMatrix().eulerAngles(0,1,2);
	return angles;
}

bool GeometricObject::applyTransformations()	{
	rotate();
	scale();
	translate();
	scale_vector = Vector3f(1.0f,1.0f,1.0f);
	translation_vector = Vector3f(0.0f,0.0f,0.0f);
	rotation = Quaternionf(1.0f,0.0f,0.0f,0.0f);
	return true;
}

std::string GeometricObject::getName()  {
    return name;
}

GeometricObject *GeometricObject::grid(const float length, const float width, const int nbSub) {
    std::vector<Vector3f> vertices, normals;
    std::vector<Vector2f> tex_coords;
    std::vector<Face *> faces;
    std::vector<Edge *> edges;

	float step_x = length/nbSub;
	float step_z = width/nbSub;
	float x = -length/2.0f;
	for (int i=0;i<nbSub;i++)	{
		float z = -width/2.0f;
		for (int j=0;j<nbSub;j++)	{
			vertices.push_back(Vector3f(x,0.0f, z));
			normals.push_back(Vector3f(0.0f,1.0f,0.0f));

		    z += step_z;
    	}
		x += step_x;
    }

	for (int i=0,pos=0;pos < vertices.size()-nbSub;i++)	{
		std::vector<int> indices;
		indices.push_back(pos);
		indices.push_back(pos+nbSub);
		indices.push_back(pos+1);

		Face *new_face = new Face();
		new_face->setVertices(indices);
		faces.push_back(new_face);

		Edge *new_edge1 = new Edge(pos, pos+nbSub);
		Edge *new_edge2 = new Edge(pos+nbSub, pos+1);
		Edge *new_edge3 = new Edge(pos+1, pos);
		edges.push_back(new_edge1);
		edges.push_back(new_edge2);
		edges.push_back(new_edge3);

		indices.clear();
		indices.push_back(pos+nbSub);
		indices.push_back(pos+nbSub+1);
		indices.push_back(pos+1);

		new_face = new Face();
		new_face->setVertices(indices);
		faces.push_back(new_face);

		new_edge1 = new Edge(pos+nbSub, pos+nbSub+1);
		new_edge2 = new Edge(pos+nbSub+1, pos+1);
		new_edge3 = new Edge(pos+1, pos+nbSub);
		edges.push_back(new_edge1);
		edges.push_back(new_edge2);
		edges.push_back(new_edge3);

		pos++;
		if ((pos+1)%nbSub == 0) pos++;
	}

    return new GeometricObject(vertices,normals,tex_coords,faces,edges);
}

GeometricObject *GeometricObject::axes(float length) {
    std::vector<Vector3f> vertices, normals;
    std::vector<Vector2f> tex_coords;
    std::vector<Face *> faces;
    std::vector<Edge *> edges;
    std::vector<int> indices;

    vertices.push_back(Vector3f(0.0f,0.0f,0.0f));
    vertices.push_back(Vector3f(length,0.0f,0.0f));
    Edge *new_edge1 = new Edge(0,1);

    vertices.push_back(Vector3f(0.0f,length,0.0f));
    Edge *new_edge2 = new Edge(0,2);

    vertices.push_back(Vector3f(0.0f,0.0f,length));
    Edge *new_edge3 = new Edge(0,3);

    edges.push_back(new_edge1);
    edges.push_back(new_edge2);
    edges.push_back(new_edge3);

    indices.clear();
    indices.push_back(0);
    indices.push_back(1);
    indices.push_back(0);
//    indices.push_back(2);
    Face *face1 = new Face();
    face1->setVertices(indices);

    indices.clear();
    indices.push_back(0);
    indices.push_back(2);
    indices.push_back(0);
//    indices.push_back(3);
    Face *face2 = new Face();
    face2->setVertices(indices);

    indices.clear();
    indices.push_back(0);
    indices.push_back(3);
    indices.push_back(0);
//    indices.push_back(1);
    Face *face3 = new Face();
    face3->setVertices(indices);

    faces.push_back(face1);
    faces.push_back(face2);
    faces.push_back(face3);

    return new GeometricObject(vertices,normals,tex_coords,faces,edges);
}

std::string GeometricObject::getObjectName()	{
	return name;
}

BoundingBox GeometricObject::getBBox()	{
	return bbox;
}

Vector3f GeometricObject::getBBoxCentroid()	{
	return 0.5f*(bbox.getMaxPt()+bbox.getMinPt());
}

bool GeometricObject::setupBBox()	{
	Vector3f min_pt = Vector3f(FLT_MAX, FLT_MAX, FLT_MAX);
	Vector3f max_pt = Vector3f(FLT_MIN, FLT_MIN, FLT_MIN);
	for (int i=0;i<vertices.size();i++)	{
		if (vertices[i](0) > max_pt(0)) max_pt(0) = vertices[i](0);
		if (vertices[i](1) > max_pt(1)) max_pt(1) = vertices[i](1);
		if (vertices[i](2) > max_pt(2)) max_pt(2) = vertices[i](2);

		if (vertices[i](0) < min_pt(0)) min_pt(0) = vertices[i](0);
		if (vertices[i](1) < min_pt(1)) min_pt(1) = vertices[i](1);
		if (vertices[i](2) < min_pt(2)) min_pt(2) = vertices[i](2);
	}
	return true;
}

Vector3f GeometricObject::getBBoxDiagonal()	{
	return bbox.getMaxPt() - bbox.getMinPt();
}
#endif

