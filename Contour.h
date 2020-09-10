//-----------------------------------------------------------------------------
// Polygon Triangulation and Extrusion Sample
// Copyright � 2001 Michael F�tsch.
//
// This file contains the declaration of the Contour class,
// which encapsulates a contour as passed to the GLU tesselator.
//
// Visit my site: http://www.geocities.com/foetsch
// Contact me: foetsch@crosswinds.net
//-----------------------------------------------------------------------------

#if !defined(ContourH)
#define ContourH

#include <Eigen/Eigen>
#include <vector>
using namespace Eigen;

#define CP_EPSILON                  0.0001f

class Contour    {
public:
    Contour();
    Contour(std::vector<Vector2i> const &_indices, std::vector<Vector3f> const &_positions);
    virtual ~Contour();

    int getNumberOfPoints();
    int getNumberOfPoints() const;

    Vector3f getCentroidPos();
    Vector2f getCentroidIndex();

    Vector3f &getPositionAt(int index);
    Vector2i &getIndexAt(int index);
    Vector3f getPositionAt(int index) const;
    Vector2i getIndexAt(int index) const;

    void addPosition(Vector3f const &_position);
    void addIndex(Vector2i const &_index);

    bool addContourPosIndex(Vector2i const &_index, Vector3f const &_position);

        /**
         *  retrieves the normal vector of the edge starting at vertex #vertexNum
         */
    Vector3f EdgeNormal(int vertexNum) const;

    /**
     *  VertexNormal
     *  Retrieves the normal vector of the vertex #vertexNum.
     */
    Vector3f VertexNormal(int vertexNum) const;

    /**
     *  VertexCosAngle
     *  retrieves the cosine of the angle at vertex #vertexNum
     */
    double VertexCosAngle(int vertexNum) const;

    /**
     *  GetVertexCE
     *  Retrieves the vertex #vertexNum. If vertexNum i s out of the range
     *  of vertices, the method performs a circular extension, i.e.,
     *  the index is wrapped around to the other end of the vector.
     */
    Vector3f GetVertexCE(int vertexNum) const;

    //std::vector<Point> contour;

    ///set the positions and indices
    void setContour(std::vector<Vector2i> const &_indices, std::vector<Vector3f> const &_positions) {
        assert(_indices.size() == _positions.size());
        positions.clear();
        positions.insert(positions.begin(), _positions.begin(), _positions.end());

        indices.clear();
        indices.insert(indices.begin(), _indices.begin(), _indices.end());
    }

    double CosAngle(Vector3f const &v1, Vector3f const &v2) const
    {
        double divisor = v1.norm() * v2.norm();
        if (divisor == 0.0)
        {
            return 0.0;
        }
        else
        {
            return v1.dot(v2) / divisor;
        }
    }



private:

    std::vector<Vector3f> positions;

    std::vector<Vector2i> indices;
};

#endif // !defined(ContourH)
