////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
//////////////////////////////////////////////////////////////////////////////////// 
#ifndef __ONB_H__
#define __ONB_H__

#include <Eigen/Eigen>
using namespace Eigen;

class ONB
{
public:
   ONB() {};

   ONB(const Vector3f& a, const Vector3f& b, const Vector3f& c)
   { U = a; V = b; W = c; }

   void initFromU( const Vector3f& u );
   void initFromV( const Vector3f& v );
   void initFromW( const Vector3f& w );

   void set(const Vector3f& a, const Vector3f& b, const Vector3f& c)
   { U = a; V = b; W = c; }

   // Calculate the ONB from two vectors
   // The first one is the Fixed vector (it is just normalized)
   // The second is normalized and its direction can be adjusted
   void  initFromUV( const Vector3f& u, const Vector3f& v );
   void  initFromVU( const Vector3f& v, const Vector3f& u );

   void  initFromUW( const Vector3f& u, const Vector3f& w );
   void  initFromWU( const Vector3f& w, const Vector3f& u );

   void  initFromVW( const Vector3f& v, const Vector3f& w );
   void  initFromWV( const Vector3f& w, const Vector3f& v );
   
   friend std::istream &operator>>(std::istream &is, ONB &t);
   friend std::ostream &operator<<(std::ostream &os, const ONB &t);
   friend bool  operator==(const ONB& o1, const ONB &o2);

   // Get a component from the ONB basis
   Vector3f u() const { return U; }
   Vector3f v() const { return V; }
   Vector3f w() const { return W; }

   Vector3f U,V,W;
};

#endif // _ONB_H_
