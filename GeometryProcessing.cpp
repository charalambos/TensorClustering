////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
////////////////////////////////////////////////////////////////////////////////////
#ifndef __GEOMETRY_PROCESSING_CPP__
#define __GEOMETRY_PROCESSING_CPP__

#include "GeometryProcessing.h"
#include "Color.h"
#include "Ransac.h"

GeometricObject *GeometryProcessing::triangulate(Image *xyz_map, Image *_normal_map)	{

	if (!xyz_map)	{
		return 0x00;
	}

	///Compute the normal map from the xyz map
	Image *normal_map = 0x00;
	if (_normal_map)	{
		normal_map = _normal_map;
	}
	else	{
		normal_map = GeometryProcessing::computeNormalMap(xyz_map, true);
	}

	///create a connected mesh of the same size as the xyz map
	std::vector<Vector3f> new_vertices,normals;
	std::vector<Vector2f> tex_coords;
	std::vector<Face *> new_faces;
	int vertex_count = 0 ;
	Image *already_added = new Image(xyz_map->getWidth(),xyz_map->getHeight(), -1.0f,-1.0f,-1.0f,1.0f);
	for (int y=0;y<xyz_map->getHeight();y++)	{
		for (int x=0;x<xyz_map->getWidth();x++)	{
			bool stop_it = false;
			int offset_x=1;
			int offset_y=1;
			if (color2vector3(xyz_map->getPixel(x,y)).norm() < GP_EPSILON) continue;
			if (y+offset_y >= xyz_map->getHeight() || x+offset_x >= xyz_map->getWidth()) continue;
			while (color2vector3(xyz_map->getPixel(x+offset_x,y)).norm() < GP_EPSILON)	{
				offset_x++;
				if (x+offset_x>=xyz_map->getWidth()) {
					stop_it = true;
					break;
				}
			}
			if (stop_it) continue;
			while (color2vector3(xyz_map->getPixel(x,y+offset_y)).norm() < GP_EPSILON)	{
				offset_y++;
				if (y+offset_y>=xyz_map->getHeight()) {
					stop_it = true;
					break;
				}
			}
			if (stop_it) continue;
			if (y+offset_y >= xyz_map->getHeight() || x+offset_x >= xyz_map->getWidth()) continue;
			if (color2vector3(xyz_map->getPixel(x+offset_x, y+offset_y)).norm() < GP_EPSILON) continue;
			std::vector<int> vertex_indices;
			int index1,index2,index3,index4;
			//add the points
			if (already_added->getPixel(x,y) == Color(-1.0f,-1.0f,-1.0f))	{
				Vector3f point1 = color2vector3(xyz_map->getPixel(x,y));
				new_vertices.push_back(point1);
				Vector3f normal1 = color2vector3(normal_map->getPixel(x,y));
				normals.push_back(normal1);
				Vector2f tex_coord1 = Vector2f(float(x)/float(xyz_map->getWidth()), float(y)/float(xyz_map->getHeight()));
				tex_coords.push_back(tex_coord1);
				index1 = vertex_count;
				already_added->setPixel(x,y,Color(float(vertex_count)));
				vertex_count++;
			}
			else	{
				index1 = int(already_added->getPixel(x,y).r());
			}
			if (already_added->getPixel(x+offset_x,y) == Color(-1.0f,-1.0f,-1.0f))	{
				Vector3f point2 = color2vector3(xyz_map->getPixel(x+offset_x,y));
				new_vertices.push_back(point2);
				Vector3f normal2 = color2vector3(normal_map->getPixel(x+offset_x,y));
				normals.push_back(normal2);
				Vector2f tex_coord2 = Vector2f(float(x+offset_x)/float(xyz_map->getWidth()), float(y)/float(xyz_map->getHeight()));
				tex_coords.push_back(tex_coord2);
				index2 = vertex_count;
				already_added->setPixel(x+offset_x,y,Color(float(vertex_count)));
				vertex_count++;
			}
			else	{
				index2 = int(already_added->getPixel(x+offset_x,y).r());
			}
			if (already_added->getPixel(x,y+offset_y)==Color(-1.0f,-1.0f,-1.0f))	{
				Vector3f point3 = color2vector3(xyz_map->getPixel(x,y+offset_y));
				new_vertices.push_back(point3);
				Vector3f normal3 = color2vector3(normal_map->getPixel(x,y+offset_y));
				normals.push_back(normal3);
				Vector2f tex_coord3 = Vector2f(float(x)/float(xyz_map->getWidth()), float(y+offset_y)/float(xyz_map->getHeight()));
				tex_coords.push_back(tex_coord3);
				index3 = vertex_count;
				already_added->setPixel(x,y+offset_y, Color(float(vertex_count)));
				vertex_count++;
			}
			else	{
				index3 = int(already_added->getPixel(x,y+offset_y).r());
			}
			if (already_added->getPixel(x+offset_x, y+offset_y) == Color(-1.0f,-1.0f,-1.0f))	{
				Vector3f point4 = color2vector3(xyz_map->getPixel(x+offset_x, y+offset_y));
				new_vertices.push_back(point4);
				Vector3f normal4 = color2vector3(normal_map->getPixel(x+offset_x,y+offset_y));
				normals.push_back(normal4);
				Vector2f tex_coord4 = Vector2f(float(x+offset_x)/float(xyz_map->getWidth()), float(y+offset_y)/float(xyz_map->getHeight()));
				tex_coords.push_back(tex_coord4);
				index4 = vertex_count;
				already_added->setPixel(x+offset_x,y+offset_y,Color(float(vertex_count)));
				vertex_count++;
			}
			else	{
				index4 = int(already_added->getPixel(x+offset_x, y+offset_y).r());
			}
			//Create 2 faces for each of the triangles form
			//printf("%d %d %d %d\n",index1,index2,index3,index4);
			vertex_indices.push_back(index1);
			vertex_indices.push_back(index3);
			vertex_indices.push_back(index2);
			//create a new face
			Face *new_face_a = new Face();
			new_face_a->setVertices(vertex_indices);
			new_faces.push_back(new_face_a);

			vertex_indices.clear();
			vertex_indices.push_back(index2);
			vertex_indices.push_back(index3);
			vertex_indices.push_back(index4);
			//create a new face
			Face *new_face_b = new Face();
			new_face_b->setVertices(vertex_indices);
			new_faces.push_back(new_face_b);
		}
	}

	///add a new object
	std::vector<Edge *> edges;

	GeometricObject *new_object = new GeometricObject(new_vertices,normals,tex_coords,new_faces,edges);

	delete already_added;
	if (_normal_map == 0x00)	{
		delete normal_map;
	}

 	return new_object;
}

Image *GeometryProcessing::computeNormalMap(Image *xyz_map, bool z_up)      {
         Image *normal_map = new Image(xyz_map->getWidth(),xyz_map->getHeight(),0.0f,0.0f,0.0f,1.0f);
        for (int y=0;y<xyz_map->getHeight();y++) {
                for (int x=0;x<xyz_map->getWidth();x++) {
                		if (xyz_map->getPixel(x,y) == Color(0.0f,0.0f,0.0f))	continue;
               			normal_map->setPixel(x,y,vector2color3(computeLocalNormal(xyz_map,Vector2i(x,y),z_up)));
                }
        }
        return normal_map;
}

Vector3f GeometryProcessing::computeLocalNormal(Image *xyz_map, Vector2i const &index, bool z_up)        {

        //get an average normal for this point
        //check the  8 neighbours of the pixel in counter clockwise order
		std::vector<Vector2i> indices;
		indices.push_back(Vector2i(index(0)-1,index(1)+1));
		indices.push_back(Vector2i(index(0)-1,index(1)));
		indices.push_back(Vector2i(index(0)-1,index(1)-1));
		indices.push_back(Vector2i(index(0),index(1)-1));
		//index5 is the pixel in question
		indices.push_back(Vector2i(index(0)+1,index(1)-1));
		indices.push_back(Vector2i(index(0)+1,index(1)));
		indices.push_back(Vector2i(index(0)+1,index(1)+1));
		indices.push_back(Vector2i(index(0),index(1)+1));

		std::vector<Vector3f> good_points;
		for (int i=0;i<indices.size();i++)	{
			if (outOfBounds(xyz_map, indices[i](0),indices[i](1)))	continue;
			if (xyz_map->getPixel(indices[i](0), indices[i](1)) == Color(0.0f,0.0f,0.0f))	continue;
			good_points.push_back(color2vector3(xyz_map->getPixel(indices[i](0),indices[i](1))));
		}

		std::vector<Vector3f> diff_vectors;
		for (int i=0;i<good_points.size();i++)	{
			diff_vectors.push_back(good_points[i] - color2vector3(xyz_map->getPixel(index(0),index(1))));
		}

		Vector3f normal = Vector3f(0.0f,0.0f,0.0f);
		for (int i=0;i<diff_vectors.size();i++)	{
			int next = (i+1)%diff_vectors.size();
			Vector3f current_normal = diff_vectors[i].cross(diff_vectors[next]);
			if (current_normal.norm() > TC_EPSILON)    current_normal.normalize();
	        ///THIS IS ADDED FOR AERIAL LIDAR DATA ONLY. THE Z SHOULD POINT UP
			if (z_up && current_normal(2) < 0.0f) {
				current_normal = -current_normal;
			}
			normal += current_normal;
		}

        normal /= float(std::max(1,int(diff_vectors.size())));
        if (normal.norm() > TC_EPSILON)    {
            normal.normalize();
        }

	return normal;
}

///Performs linear plane fitting on a set of points
Vector4f GeometryProcessing::linearPlaneFit ( std::vector<Vector3f > const &points, bool use_ransac,  float *fitting_error,Vector3f *plane_point)   {

    Vector4f plane_coeffs = Vector4f(0.0f,0.0f,0.0f,0.0f);

    Vector3f plane_normal = Vector3f(0.0f,0.0f,0.0f);

    ///compute the centroid
    Vector3f centroid = Vector3f(0.0f,0.0f,0.0f);
    for ( int i=0;i<points.size();i++ ) {
        centroid+= points[i];
    }
    if ( points.size() )	centroid/=float ( points.size() );

    if ( plane_point ) ( *plane_point ) = centroid;

    ///if you have too many points then use ransac
    if ( use_ransac || points.size() > 100 )  {
        ///USING RANSAC
        std::vector<int> indices;
        ///ransac
        ransac ( points, plane_fitting_function, evaluate_plane_fit_function, is_degenerate_function, 3, 0.1f, centroid, &plane_coeffs,&indices );

        plane_normal = Vector3f(plane_coeffs(0),plane_coeffs(1),plane_coeffs(2));
    }
    else    {

        float sum_xx, sum_xy, sum_xz, sum_yy, sum_yz, sum_zz;

        sum_xx = sum_xy = sum_xz = 0.0f;
        sum_yy = sum_yz = 0.0f;
        sum_zz = 0.0f;

        ///sum the squares of the differences
        for ( int i=0;i<points.size();i++ ) {
            ///compute the differences
            float diff_x = float ( points[i] ( 0 ) - centroid ( 0 ) );
            float diff_y = float ( points[i] ( 1 ) - centroid ( 1 ) );
            float diff_z = float ( points[i] ( 2 ) - centroid ( 2 ) );
            ///sum the squares
            sum_xx += diff_x * diff_x;
            sum_xy += diff_x * diff_y;
            sum_xz += diff_x * diff_z;

            sum_yy += diff_y * diff_y;
            sum_yz += diff_y * diff_z;

            sum_zz += diff_z * diff_z;
        }

        ///build a matrix
        float mat[3][3];

        mat[0][0] = sum_xx;	mat[0][1] = sum_xy;	mat[0][2] = sum_xz;
        mat[1][0] = sum_xy;	mat[1][1] = sum_yy;	mat[1][2] = sum_yz;
        mat[2][0] = sum_xz;	mat[2][1] = sum_yz;	mat[2][2] = sum_zz;

        float lambda[3];

        eig_sys3d ( mat,lambda);

        plane_normal = Vector3f(mat[0][2], mat[1][2], mat[2][2] );
        plane_coeffs = Vector4f(plane_normal(0), plane_normal(1), plane_normal(2), -plane_normal.dot(centroid));
    }

    plane_coeffs(3) /= plane_normal.norm();
    plane_normal/= plane_normal.norm();

    ///if there was an error fitting variable then compute the error
    if ( fitting_error )    {
        ///so you compute the distance of the candidate point to the plane
        float squared_error = 0.0f;
        for (int i=0;i<points.size();i++ )    {
            Vector3f point_to_plane = points[i] - centroid;
            float distance_to_plane = point_to_plane.dot(plane_normal);
            squared_error += distance_to_plane*distance_to_plane;
        }

        (*fitting_error) =  squared_error/float(points.size());
    }

    return plane_coeffs;
}

#endif
