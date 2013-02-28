// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
///
/// \file		cySpatial.h 
/// \author		Cem Yuksel
/// \version	1.0
/// \date		June 7, 2007
///
/// \brief Spatial vector algebra classes
/// 
/// This file includes spatial vector algebra classes intended for the
/// implementation of Featherstone's articulated rigid body dynamics method.
/// cySpatialVector6f class is both for spatial motion vectors and spatial
/// force vectors, cySpatialTrans6f is a spatial matrix class for coordinate
/// transformations only, and cySpatialMatrix6f is the general spatial
/// matrix class.
///
//-------------------------------------------------------------------------------

#ifndef _CY_SPATIAL_H_INCLUDED_
#define _CY_SPATIAL_H_INCLUDED_

//-------------------------------------------------------------------------------

#include "cyPoint.h"
#include "cyMatrix3.h"

//-------------------------------------------------------------------------------

/// 6D spatial vector (for 3D).
/// This class is both for spatial motion vectors and spatial force vectors.

class cySpatialVector6f
{
public:
	cyPoint3f a, b;

	///@name Constructors
	cySpatialVector6f() {}
	cySpatialVector6f( const cyPoint3f &p1, const cyPoint3f &p2 ) { Set(p1,p2); }
	cySpatialVector6f( float a1, float a2, float a3, float b1, float b2, float b3 ) { Set(a1,a2,a3,b1,b2,b3); }
	cySpatialVector6f( const cySpatialVector6f &v ) { a=v.a; b=v.b; }

	///@name Initialization methods
	void Set( const cyPoint3f &p1, const cyPoint3f &p2 ) { a = p1; b = p2; }
	void Set( float a1, float a2, float a3, float b1, float b2, float b3 ) { a.Set(a1,a2,a3); b.Set(b1,b2,b3); }
	void Zero() { a.Zero(); b.Zero(); }

	///@name Transpose methods
	void SetTranspose() { cyPoint3f p=a; a=b; b=p; }
	cySpatialVector6f Transpose() const { return cySpatialVector6f( b, a ); }

	///@name Unary operators
	cySpatialVector6f operator-() const { return cySpatialVector6f(-a,-b); } 

	///@name Binary operators
	cySpatialVector6f	operator + ( const cySpatialVector6f &s ) const { return cySpatialVector6f(a+s.a, b+s.b); }
	cySpatialVector6f	operator - ( const cySpatialVector6f &s ) const { return cySpatialVector6f(a-s.a, b-s.b); }
	cySpatialVector6f	operator * ( float t ) const { return cySpatialVector6f( a*t, b*t ); }

	/// Scalar product of two vectors.
	/// Note that one of the vectors should be motion vector ant the other should be a force vector.
	/// Otherwise, scalar product is not defined in spatial vector algebra.
	float operator * ( const cySpatialVector6f &s ) const { return a.Dot(s.a) + b.Dot(s.b); }

	///@name Assignment operators
	void operator =  ( const cySpatialVector6f &v ) { a=v.a; b=v.b; }
	void operator += ( const cySpatialVector6f &s ) { a+=s.a; b+=s.b; }
	void operator -= ( const cySpatialVector6f &s ) { a-=s.a; b-=s.b; }
	void operator *= ( float t ) { a*=t; b*=t; }

};


//-------------------------------------------------------------------------------

/// 6D spatial matrix for coordinate transforms.
/// This is a special case for cySpatialMatrix6f class,
/// where the matrix represents a coordinate transformation.
/// In this case, instead of keeping a full 6x6 matrix values,
/// we can keep a 3x3 matrix for rotation, and a 3D point
/// for translation. This compact representation simplifies
/// some computations, therefore you should use this class
/// instead of cySpatialMatrix6f whenever you represent
/// a coordinate transformation. However, for general matrix
/// operations, you have to use cySpatialMatrix6f.
///

class cySpatialTrans6f
{
public:

	// |    R     0 |
	// | -r x R   R |

	cyMatrix3f R;	///< Rotation matrix
	cyPoint3f r;	///< Transformation

	///@name Constructors
	cySpatialTrans6f() {}
	cySpatialTrans6f( const cySpatialTrans6f &mat ) { R=mat.R; r=mat.r; }
	cySpatialTrans6f( const cyMatrix3f &_R, const cyPoint3f &_r ) { Set(_R,_r); }

	///@name Initialization methods
	void Set( const cyMatrix3f &_R, const cyPoint3f &_r ) { R=_R; r=_r; }
	void SetIdentity() { R.SetIdentity(); r.Zero(); }

	///@name Unary operators
	cySpatialTrans6f operator - () const { return cySpatialTrans6f( -R, -r ); }


	///@name Binary operators
	cySpatialVector6f operator * ( const cySpatialVector6f &p ) const { cyPoint3f Ra = R*p.a; return cySpatialVector6f( Ra, cyMatrix3f(-r)*Ra + R*p.b ); }

	cySpatialTrans6f operator * ( const cySpatialTrans6f &mat ) const { return cySpatialTrans6f( R*mat.R, r + R*mat.r ); }
	cySpatialTrans6f operator + ( const cySpatialTrans6f &mat ) const { return cySpatialTrans6f( R + mat.R, r + mat.r ); }
	cySpatialTrans6f operator - ( const cySpatialTrans6f &mat ) const { return cySpatialTrans6f( R - mat.R, r - mat.r ); }

	cySpatialTrans6f operator * ( float t ) const { return cySpatialTrans6f(R*t,r*t); }
	cySpatialTrans6f operator / ( float t ) const { float d=1.0f/t; return *this * d; }


	///@name Assignment operators
	void operator *= ( const cySpatialTrans6f &mat ) { *this = *this * mat; }
	void operator += ( const cySpatialTrans6f &mat ) { R+=mat.R; r+=mat.r; }
	void operator -= ( const cySpatialTrans6f &mat ) { R-=mat.R; r-=mat.r; }
	void operator *= ( float t ) { *this = *this * t; }

};

//-------------------------------------------------------------------------------

/// 6D spatial matrix.
/// This is the general class for 6D spatial matrices.
/// For representing coordinate transformation matrices
/// use cySpatialTrans6f instead, since it is more efficient.
/// However, cySpatialTrans6f cannot be used for general
/// matrix operations that do not correspond to a
/// coordinate transformation.

class cySpatialMatrix6f
{
public:

	// | m[0]  m[2] |
	// | m[1]  m[3] |

	cyMatrix3f m[4];	///< Matrix data in column major order


	///@name Constructors
	cySpatialMatrix6f() {}
	cySpatialMatrix6f( const cySpatialMatrix6f &mat ) { m[0]=mat.m[0]; m[1]=mat.m[1]; m[2]=mat.m[2]; m[3]=mat.m[3]; }
	cySpatialMatrix6f( const cyMatrix3f &_R, const cyPoint3f &_r ) { Set(_R,_r); }
	cySpatialMatrix6f( const cyMatrix3f &m11, const cyMatrix3f &m21, const cyMatrix3f &m12, const cyMatrix3f &m22 ) { Set(m11,m21,m12,m22); }
	cySpatialMatrix6f( const cySpatialVector6f &p1, const cySpatialVector6f &p2 ) { Set(p1,p2); }	// Outer product of two spatial vectors
	cySpatialMatrix6f( const cySpatialTrans6f &tm ) { Set(tm); }


	///@name Initialization methods
	void Set( const cyMatrix3f &_R, const cyPoint3f &_r ) { m[0]=m[3]=_R; m[1]=cyMatrix3f(-_r)*_R; m[2].Zero(); }
	void Set( const cyMatrix3f &m11, const cyMatrix3f &m21, const cyMatrix3f &m12, const cyMatrix3f &m22 ) { m[0]=m11; m[1]=m21; m[2]=m12; m[3]=m22; }
	void Set( const cySpatialTrans6f &tm ) { m[0]=m[3]=tm.R; m[1]=cyMatrix3f(-tm.r)*tm.R; m[2].Zero(); }

	/// Sets the matrix as the outer product of two vectors.
	void Set( const cySpatialVector6f &p1, const cySpatialVector6f &p2 )
	{
			SetMatrix( m[0], p1.a, p2.a );
			SetMatrix( m[1], p1.b, p2.a );
			SetMatrix( m[2], p1.a, p2.b );
			SetMatrix( m[3], p1.b, p2.b );
	}

	void SetIdentity() { m[0].SetIdentity(); m[1].Zero(); m[2].Zero(); m[3].SetIdentity(); }
	void Zero() { m[0].Zero(); m[1].Zero(); m[2].Zero(); m[3].Zero(); }


	///@name Unary operators
	cySpatialMatrix6f operator - () const { return cySpatialMatrix6f( -m[0], -m[1], -m[2], -m[3] ); }


	///@name Unary operators

	cySpatialVector6f operator * ( const cySpatialVector6f &p ) const { return cySpatialVector6f( m[0]*p.a + m[2]*p.b, m[1]*p.a + m[3]*p.b ); }

	cySpatialMatrix6f operator * ( const cySpatialMatrix6f &mat ) const { return cySpatialMatrix6f( m[0]*mat.m[0]+m[2]*mat.m[1], m[1]*mat.m[0]+m[3]*mat.m[1], m[0]*mat.m[2]+m[2]*mat.m[3], m[1]*mat.m[2]+m[3]*mat.m[3] ); }
	cySpatialMatrix6f operator + ( const cySpatialMatrix6f &mat ) const { return cySpatialMatrix6f( m[0]+mat.m[0], m[1]+mat.m[1], m[2]+mat.m[2], m[3]+mat.m[3] ); }
	cySpatialMatrix6f operator - ( const cySpatialMatrix6f &mat ) const { return cySpatialMatrix6f( m[0]-mat.m[0], m[1]-mat.m[1], m[2]-mat.m[2], m[3]-mat.m[3] ); }

	cySpatialMatrix6f operator * ( float t ) const { return cySpatialMatrix6f(m[0]*t,m[1]*t,m[2]*t,m[3]*t); }
	cySpatialMatrix6f operator / ( float t ) const { float d=1.0f/t; return *this * d; }


	///@name Assignment operators
	void	operator *= ( const cySpatialMatrix6f &mat ) { *this = *this * mat; }
	void	operator += ( const cySpatialMatrix6f &mat ) { m[0]+=mat.m[0]; m[1]+=mat.m[1]; m[2]+=mat.m[2]; m[3]+=mat.m[3]; }
	void	operator -= ( const cySpatialMatrix6f &mat ) { m[0]-=mat.m[0]; m[1]-=mat.m[1]; m[2]-=mat.m[2]; m[3]-=mat.m[3]; }
	void	operator *= ( float t ) { *this = *this * t; }


protected:

	/// \internal

	// Sets the given matrix as the outer product of the given two vectors.
	void SetMatrix( cyMatrix3f &m, const cyPoint3f &p1, const cyPoint3f &p2 )
	{ 
		float val[] = {p1.x * p2.x, p1.y * p2.x, p1.z * p2.x,   p1.x * p2.y, p1.y * p2.y, p1.z * p2.y,   p1.x * p2.z, p1.y * p2.z, p1.z * p2.z};
		m.Set( val ); 
	}

};

//-------------------------------------------------------------------------------


namespace cy {
	typedef cySpatialVector6f SpatialVector6f;
	typedef cySpatialTrans6f  SpatialTrans6f;
	typedef cySpatialMatrix6f SpatialMatrix6f;
}


//-------------------------------------------------------------------------------

#endif
