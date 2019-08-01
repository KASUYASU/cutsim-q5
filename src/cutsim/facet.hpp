#ifndef FACET_H
#define FACET_H

#include <cassert>

#include "glvertex.hpp"

namespace cutsim {

class Facet {

   public:
	   Facet(GLVertex n, GLVertex p1, GLVertex p2, GLVertex p3) { normal = n; v1 = p1; v2 = p2; v3 = p3; }
	   GLVertex	normal;
	   GLVertex	v1, v2, v3;

       static GLVertex facetCenter(const Facet *facet) { GLVertex center;
                                                  center.x = (facet->v1.x + facet->v2.x + facet->v3.x) * 1.0/3.0;
                                                  center.y = (facet->v1.y + facet->v2.y + facet->v3.y) * 1.0/3.0;
                                                  center.z = (facet->v1.z + facet->v2.z + facet->v3.z) * 1.0/3.0;
                                                  return center;
                                                }

};

}

#endif // FACET_H

