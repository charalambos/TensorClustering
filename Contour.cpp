//-----------------------------------------------------------------------------
// Polygon Triangulation and Extrusion Sample
// Copyright � 2001 Michael F�tsch.
//
// This file contains the implementation of the Contour class,
// which encapsulates a contour as passed to the GLU tesselator.
//
// Visit my site: http://www.geocities.com/foetsch
// Contact me: foetsch@crosswinds.net
//-----------------------------------------------------------------------------

#include "Contour.h"

Contour::Contour()  {

}

Contour::Contour(std::vector<Vector2i> const &_indices, std::vector<Vector3f> const &_positions)    {
    setContour(_indices, _positions);
}

Contour::~Contour()
{
}

int Contour::getNumberOfPoints()   {
    return positions.size();
}

int Contour::getNumberOfPoints() const  {
    return positions.size();
}

Vector3f &Contour::getPositionAt(int index) {
    return positions[index];
}

Vector2i &Contour::getIndexAt(int index) {
    return indices[index];
}

Vector3f Contour::getPositionAt(int index) const {
    return positions[index];
}

Vector2i Contour::getIndexAt(int index) const {
    return indices[index];
}

bool Contour::addContourPosIndex(Vector2i const &_index, Vector3f const &_position) {
    for (int i=positions.size()-1;i>=0;i--)	{
        if ((positions[i] - _position).norm() < CP_EPSILON && (indices[i] - _index).norm() < CP_EPSILON)	{
            return false;
        }
    }
    positions.push_back(_position);
    indices.push_back(_index);
    return true;
}

Vector3f Contour::getCentroidPos()   {
    Vector3f centroid = Vector3f(0.0f,0.0f,0.0f);
    for (int i=0;i<positions.size();i++)    {
        centroid += positions[i];
    }
    centroid /= float(positions.size());
    return centroid;
}

Vector2f Contour::getCentroidIndex()    {
    Vector2f centroid = Vector2f(0.0f,0.0f);
    for (int i=0;i<indices.size();i++)  {
        centroid += Vector2f(float(indices[i](0)), float(indices[i](1)));
    }
    centroid /= float(indices.size());
    return centroid;
}


Vector3f Contour::EdgeNormal(int vertexNum) const
{
    Vector3f temp = Vector3f(GetVertexCE(vertexNum + 1) - GetVertexCE(vertexNum));
    Vector3f temp_normal = Vector3f(temp(1), -temp(0), 0.0f);
    temp_normal.normalize();
    return temp_normal;
}


Vector3f Contour::VertexNormal(int vertexNum) const     {
    Vector3f v1 = Vector3f(GetVertexCE(vertexNum - 1)-GetVertexCE(vertexNum));
    Vector3f v1_temp = Vector3f(v1(1),-v1(0),0.0f);
    v1_temp.normalize();
    Vector3f v2 = Vector3f(GetVertexCE(vertexNum + 1)-GetVertexCE(vertexNum));
    Vector3f v2_temp = Vector3f(v2(1),-v2(0),0.0f);
    v2_temp.normalize();

    Vector3f out = (v1_temp + v2_temp*-1.0f);
    out.normalize();
    out *= -1.0f;
    return out;
}


double Contour::VertexCosAngle(int vertexNum) const
{
    Vector3f v1 = GetVertexCE(vertexNum-1) - GetVertexCE(vertexNum);
    Vector3f v2 = GetVertexCE(vertexNum+1) - GetVertexCE(vertexNum);
    return CosAngle(v1,v2);
}


Vector3f Contour::GetVertexCE(int vertexNum) const  {
    if (vertexNum < 0)
    {
        vertexNum = positions.size() - (-vertexNum % positions.size());
    }
    return positions.at(vertexNum % positions.size());
}


