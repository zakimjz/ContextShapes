#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

//index bbox of surface triangles
#include "CGAL/Cartesian.h"
#include "CGAL/Segment_tree_k.h"
#include "CGAL/Range_segment_tree_traits.h"

//used to index triangle, to find closest pairs
#include "CGAL/basic.h"
#include "CGAL/Point_3.h"
#include "CGAL/Range_tree_k.h"

typedef CGAL::Cartesian<double> UINT_Representation;
typedef CGAL::Range_tree_map_traits_3<UINT_Representation, unsigned int> UINT_Traits;
typedef CGAL::Range_tree_3<UINT_Traits> UINT_Range_tree_3_type;
typedef UINT_Traits::Key UINT_Key;
typedef UINT_Traits::Pure_key UINT_Pure_key;
typedef UINT_Traits::Interval UINT_Interval;

typedef unsigned int    TID_VERT;
typedef unsigned int    TID_VERT_TYPE;
typedef unsigned int    TID_ATOM;
typedef unsigned int  CONTEXT_RAY_INT_TYPE;

#endif
