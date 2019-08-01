/*
 *  Copyright 2010 Anders Wallin (anders.e.e.wallin "at" gmail.com)
 *  Copyright 2015      Kazuyasu Hamada (k-hamada@gifu-u.ac.jp)
 *
 *  This file is part of OpenCAMlib.
 *
 *  OpenCAMlib is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  OpenCAMlib is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with OpenCAMlib.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <cassert>
#include <cmath>

#include <src/app/cutsim_def.hpp>

#include "volume.hpp"

namespace cutsim {

//************* Sphere **************/

/// sphere at center
SphereVolume::SphereVolume() {
	type = SPHERE_VOLUME;
    center = GLVertex(0,0,0);
    radius = 1.0;
    calcBB();
}

double SphereVolume::dist(const GLVertex& p ) const {
    double d = (center-p).norm();
    return radius-d; // positive inside. negative outside.
}

/// set the bounding box values
void SphereVolume::calcBB() {
    bb.clear();
    GLVertex maxpt = GLVertex(center.x + radius, center.y + radius, center.z + radius);
    GLVertex minpt = GLVertex(center.x - radius, center.y - radius, center.z - radius);
    bb.addPoint( maxpt );
    bb.addPoint( minpt );
}

//************* Rectangle **************/

RectVolume::RectVolume() {
	type = RECTANGLE_VOLUME;
    corner = GLVertex(0,0,0);
    v1 = GLVertex(1,0,0);
    v2 = GLVertex(0,1,0);
    v3 = GLVertex(0,0,1);
}

// FIXME??
void RectVolume::calcBB() {
    bb.clear();
    GLVertex maxp;
    GLVertex minp;
    double max_x, max_y, max_z;
    double min_x, min_y, min_z;

    max_x = fmax(fmax(fmax(corner.x, v1.x), v2.x), v3.x);
    max_y = fmax(fmax(fmax(corner.y, v1.y), v2.y), v3.y);
    max_z = fmax(fmax(fmax(corner.z, v1.z), v2.z), v3.z);
    maxp = GLVertex(max_x, max_y, max_z);
    min_x = fmin(fmin(fmin(corner.x, v1.x), v2.x), v3.x);
    min_y = fmin(fmin(fmin(corner.y, v1.y), v2.y), v3.y);
    min_z = fmin(fmin(fmin(corner.z, v1.z), v2.z), v3.z);
    minp = GLVertex(min_x, min_y, min_z);
    bb.addPoint( maxp );
    bb.addPoint( minp );
}

double RectVolume::dist(const GLVertex& p) const {
    // translate to origo
    double max_x = corner.x + v1.x;
    double min_x = corner.x;
    double max_y = corner.y + v2.y;
    double min_y = corner.y;
    double max_z = corner.z + v3.z;
    double min_z = corner.z;
    double dOut = 0.0;

    if ( (min_x <= p.x) && (p.x <= max_x) && (min_y <= p.y) && (p.y <= max_y) && (min_z <= p.z) && (p.z <= max_z) )   {
        double xdist,ydist,zdist;
        if ( (p.x-min_x) > (max_x-p.x) )
            xdist = max_x-p.x;
        else
            xdist = p.x-min_x;

        if ( (p.y-min_y) > (max_y-p.y) )
            ydist = max_y-p.y;
        else
            ydist = p.y-min_y;

        if ( (p.z-min_z) > (max_z-p.z) )
            zdist = max_z-p.z;
        else
            zdist = p.z-min_z;

        if ( xdist <= ydist && xdist <= zdist )
            dOut = -xdist;
        else if ( ydist < xdist && ydist < zdist )
            dOut = -ydist;
        else if ( zdist < xdist && zdist< xdist )
            dOut = -zdist;
        else {
            assert(0);
            return -1;
        }
    } else if ( (min_y <= p.y) && (p.y <= max_y) && (min_z <= p.z) && (p.z <= max_z) )   {
        if (p.x < min_x) {
            dOut = min_x - p.x;
        } else if ( p.x > max_x ) {
            dOut = p.x-max_x;
        }
    } else if ( (min_x <= p.x) && (p.x <= max_x) && (min_z <= p.z) && (p.z <= max_z) )   {
        if (p.y < min_y) {
            dOut = min_y - p.y;
        } else if ( p.y > max_y ) {
            dOut = p.y-max_y;
        }
    } else if ( (min_x <= p.x) && (p.x <= max_x) && (min_y <= p.y) && (p.y <= max_y) )   {
        if (p.z < min_z) {
            dOut = min_z - p.z;
        } else if ( p.z > max_z ) {
            dOut = p.z-max_z;
        }
    } else if ( (p.x > max_x) && (p.y > max_y))
        dOut = sqrt((p.x - max_x)*(p.x - max_x)+(p.y - max_y)*(p.y - max_y));
    else if ( (p.x > max_x) && (p.z < min_z))
        dOut = sqrt((p.x - max_x)*(p.x - max_x)+(min_z - p.z)*(min_z - p.z));
    else if ( (p.x < min_x) && (p.y > max_y))
        dOut = sqrt((min_x - p.x)*(min_x - p.x)+(p.y - max_y)*(p.y - max_y));
    else if ( (p.y > max_y) && (p.z > max_z))
        dOut = sqrt((p.y - max_y)*(p.y - max_y)+(p.z - max_z)*(p.z - max_z));
    else if ( (p.x > max_x) && (p.z > max_z))
        dOut = sqrt((p.x - max_x)*(p.x - max_x)+(p.z - max_z)*(p.z - max_z));
    else if ( (p.x > max_x) && (p.y < min_y))
        dOut = sqrt((p.x - max_x)*(p.x - max_x)+(min_y - p.y)*(min_y - p.y));
    else if ( (p.x < min_x) && (p.y < min_y))
        dOut = sqrt((min_x - p.x)*(min_x - p.x)+(p.y - max_y)*(p.y - max_y));
    else if ( (p.y < min_y) && (p.z > max_z))
        dOut = sqrt((min_y - p.y)*(min_y - p.y)+(p.z - max_z)*(p.z - max_z));
    else if ( (p.x < min_x) && (p.z < min_z) )
        dOut = sqrt((p.x - max_x)*(p.x - max_x)+(min_z - p.z)*(min_z - p.z));
    else if ( (p.y > max_y) && (p.z < min_z))
        dOut = sqrt((p.y - max_y)*(p.y - max_y)+(min_y - p.y)*(min_y - p.y));
    else if ( (p.x < min_x) && (p.z > max_z))
        dOut = sqrt((min_x - p.x)*(min_x - p.x)+(p.z - max_z)*(p.z - max_z));
    else if ( (p.y < min_y) && (p.z < min_z))
        dOut = sqrt((min_y - p.y)*(min_y - p.y)+(min_z - p.z)*(min_z - p.z));

    return -dOut;
}

RectVolume2::RectVolume2() {
	type = RECTANGLE_VOLUME;
    corner = GLVertex(0,0,0);
    width = 1.0;
    length = 1.0;
    hight = 1.0;
    center = GLVertex(corner.x+width*0.5, corner.y+length*0.5, corner.z+hight*0.5);	// center is located at the center of box
    rotationCenter = GLVertex(0, 0, 0);
    angle = GLVertex(0, 0, 0);
}

void RectVolume2::calcBB() {
    GLVertex p[8] = { center + GLVertex( width * 0.5,  length * 0.5,    0.0 ),
					  center + GLVertex(-width * 0.5,  length * 0.5,    0.0 ),
					  center + GLVertex(-width * 0.5, -length * 0.5,    0.0 ),
					  center + GLVertex( width * 0.5, -length * 0.5,    0.0 ),
					  center + GLVertex( width * 0.5,  length * 0.5,  hight ),
					  center + GLVertex(-width * 0.5,  length * 0.5,  hight ),
					  center + GLVertex(-width * 0.5, -length * 0.5,  hight ),
					  center + GLVertex( width * 0.5, -length * 0.5,  hight ) };

    for (int n = 0; n < 8; n++) {
		GLVertex rotated_p = p[n] - rotationCenter;
		p[n] = rotated_p.rotateABC(angle.x, angle.y, angle.z) + rotationCenter;
	}

    GLVertex maxpt;
    GLVertex minpt;
    maxpt.x = fmax(fmax(fmax(fmax(fmax(fmax(fmax(p[0].x, p[1].x), p[2].x), p[3].x), p[4].x), p[5].x), p[6].x), p[7].x) + TOLERANCE;
    maxpt.y = fmax(fmax(fmax(fmax(fmax(fmax(fmax(p[0].y, p[1].y), p[2].y), p[3].y), p[4].y), p[5].y), p[6].y), p[7].y) + TOLERANCE;
    maxpt.z = fmax(fmax(fmax(fmax(fmax(fmax(fmax(p[0].z, p[1].z), p[2].z), p[3].z), p[4].z), p[5].z), p[6].z), p[7].z) + TOLERANCE;
    minpt.x = fmin(fmin(fmin(fmin(fmin(fmin(fmin(p[0].x, p[1].x), p[2].x), p[3].x), p[4].x), p[5].x), p[6].x), p[7].x) - TOLERANCE;
    minpt.y = fmin(fmin(fmin(fmin(fmin(fmin(fmin(p[0].y, p[1].y), p[2].y), p[3].y), p[4].y), p[5].y), p[6].y), p[7].y) - TOLERANCE;
    minpt.z = fmin(fmin(fmin(fmin(fmin(fmin(fmin(p[0].z, p[1].z), p[2].z), p[3].z), p[4].z), p[5].z), p[6].z), p[7].z) - TOLERANCE;
    bb.clear();
    bb.addPoint( maxpt );
    bb.addPoint( minpt );
}

double RectVolume2::dist(const GLVertex& p) const {
    GLVertex rotated_p = p - rotationCenter;
    rotated_p = rotated_p.rotateCBA(-angle.x, -angle.y, -angle.z) + rotationCenter;
    double max_x = corner.x + width;
    double min_x = corner.x;
    double max_y = corner.y + length;
    double min_y = corner.y;
    double max_z = corner.z + hight;
    double min_z = corner.z;
    double dOut = 0.0;

    if ( (min_x <= rotated_p.x) && (rotated_p.x <= max_x) && (min_y <= rotated_p.y) && (rotated_p.y <= max_y) && (min_z <= rotated_p.z) && (rotated_p.z <= max_z) )   {
        double xdist,ydist,zdist;
        if ( (rotated_p.x-min_x) > (max_x-rotated_p.x) )
            xdist = max_x-rotated_p.x;
        else
            xdist = rotated_p.x-min_x;

        if ( (rotated_p.y-min_y) > (max_y-rotated_p.y) )
            ydist = max_y-rotated_p.y;
        else
            ydist = rotated_p.y-min_y;

        if ( (rotated_p.z-min_z) > (max_z-rotated_p.z) )
            zdist = max_z-rotated_p.z;
        else
            zdist = rotated_p.z-min_z;

        if ( xdist <= ydist && xdist <= zdist )
            dOut = -xdist;
        else if ( ydist < xdist && ydist < zdist )
            dOut = -ydist;
        else if ( zdist < xdist && zdist< xdist )
            dOut = -zdist;
        else {
            assert(0);
            return -1;
        }
    } else if ( (min_y <= rotated_p.y) && (rotated_p.y <= max_y) && (min_z <= rotated_p.z) && (rotated_p.z <= max_z) )   {
        if (rotated_p.x < min_x) {
            dOut = min_x - rotated_p.x;
        } else if ( rotated_p.x > max_x ) {
            dOut = rotated_p.x-max_x;
        }
    } else if ( (min_x <= rotated_p.x) && (rotated_p.x <= max_x) && (min_z <= rotated_p.z) && (rotated_p.z <= max_z) )   {
        if (rotated_p.y < min_y) {
            dOut = min_y - rotated_p.y;
        } else if ( rotated_p.y > max_y ) {
            dOut = rotated_p.y-max_y;
        }
    } else if ( (min_x <= rotated_p.x) && (rotated_p.x <= max_x) && (min_y <= rotated_p.y) && (rotated_p.y <= max_y) )   {
        if (rotated_p.z < min_z) {
            dOut = min_z - rotated_p.z;
        } else if ( rotated_p.z > max_z ) {
            dOut = rotated_p.z-max_z;
        }
    } else if ( (rotated_p.x > max_x) && (rotated_p.y > max_y))
        dOut = sqrt((rotated_p.x - max_x)*(rotated_p.x - max_x)+(rotated_p.y - max_y)*(rotated_p.y - max_y));
    else if ( (rotated_p.x > max_x) && (rotated_p.z < min_z))
        dOut = sqrt((rotated_p.x - max_x)*(rotated_p.x - max_x)+(min_z - rotated_p.z)*(min_z - rotated_p.z));
    else if ( (rotated_p.x < min_x) && (rotated_p.y > max_y))
        dOut = sqrt((min_x - rotated_p.x)*(min_x - rotated_p.x)+(rotated_p.y - max_y)*(rotated_p.y - max_y));
    else if ( (rotated_p.y > max_y) && (rotated_p.z > max_z))
        dOut = sqrt((rotated_p.y - max_y)*(rotated_p.y - max_y)+(rotated_p.z - max_z)*(rotated_p.z - max_z));
    else if ( (rotated_p.x > max_x) && (rotated_p.z > max_z))
        dOut = sqrt((rotated_p.x - max_x)*(rotated_p.x - max_x)+(rotated_p.z - max_z)*(rotated_p.z - max_z));
    else if ( (rotated_p.x > max_x) && (rotated_p.y < min_y))
        dOut = sqrt((rotated_p.x - max_x)*(rotated_p.x - max_x)+(min_y - rotated_p.y)*(min_y - rotated_p.y));
    else if ( (rotated_p.x < min_x) && (rotated_p.y < min_y))
        dOut = sqrt((min_x - rotated_p.x)*(min_x - rotated_p.x)+(rotated_p.y - max_y)*(rotated_p.y - max_y));
    else if ( (rotated_p.y < min_y) && (rotated_p.z > max_z))
        dOut = sqrt((min_y - rotated_p.y)*(min_y - rotated_p.y)+(rotated_p.z - max_z)*(rotated_p.z - max_z));
    else if ( (rotated_p.x < min_x) && (rotated_p.z < min_z) )
        dOut = sqrt((rotated_p.x - max_x)*(rotated_p.x - max_x)+(min_z - rotated_p.z)*(min_z - rotated_p.z));
    else if ( (rotated_p.y > max_y) && (rotated_p.z < min_z))
        dOut = sqrt((rotated_p.y - max_y)*(rotated_p.y - max_y)+(min_y - rotated_p.y)*(min_y - rotated_p.y));
    else if ( (rotated_p.x < min_x) && (rotated_p.z > max_z))
        dOut = sqrt((min_x - rotated_p.x)*(min_x - rotated_p.x)+(rotated_p.z - max_z)*(rotated_p.z - max_z));
    else if ( (rotated_p.y < min_y) && (rotated_p.z < min_z))
        dOut = sqrt((min_y - rotated_p.y)*(min_y - rotated_p.y)+(min_z - rotated_p.z)*(min_z - rotated_p.z));

    return -dOut;
}

//************* Cylinder **************/

CylinderVolume::CylinderVolume() {
	type = CYLINDER_VOLUME;
    radius = 1.0;
    length = 1.0;
    center = GLVertex(0, 0, 0);	// center is located at the bottom of cylinder
    rotationCenter = GLVertex(0, 0, 0);
    angle  = GLVertex(0, 0, 0);
}

void CylinderVolume::calcBB() {
	GLVertex p[8] = { center + GLVertex( radius,  radius,    0.0 ),
					  center + GLVertex(-radius,  radius,    0.0 ),
					  center + GLVertex(-radius, -radius,    0.0 ),
					  center + GLVertex( radius, -radius,    0.0 ),
					  center + GLVertex( radius,  radius, length ),
					  center + GLVertex(-radius,  radius, length ),
					  center + GLVertex(-radius, -radius, length ),
					  center + GLVertex( radius, -radius, length ) };

    for (int n = 0; n < 8; n++) {
		GLVertex rotated_p = p[n] - rotationCenter;
		p[n] = rotated_p.rotateABC(angle.x, angle.y, angle.z) + rotationCenter;
	}

	GLVertex maxpt;
	GLVertex minpt;
	maxpt.x = fmax(fmax(fmax(fmax(fmax(fmax(fmax(p[0].x, p[1].x), p[2].x), p[3].x), p[4].x), p[5].x), p[6].x), p[7].x) + TOLERANCE;
	maxpt.y = fmax(fmax(fmax(fmax(fmax(fmax(fmax(p[0].y, p[1].y), p[2].y), p[3].y), p[4].y), p[5].y), p[6].y), p[7].y) + TOLERANCE;
	maxpt.z = fmax(fmax(fmax(fmax(fmax(fmax(fmax(p[0].z, p[1].z), p[2].z), p[3].z), p[4].z), p[5].z), p[6].z), p[7].z) + TOLERANCE;
	minpt.x = fmin(fmin(fmin(fmin(fmin(fmin(fmin(p[0].x, p[1].x), p[2].x), p[3].x), p[4].x), p[5].x), p[6].x), p[7].x) - TOLERANCE;
	minpt.y = fmin(fmin(fmin(fmin(fmin(fmin(fmin(p[0].y, p[1].y), p[2].y), p[3].y), p[4].y), p[5].y), p[6].y), p[7].y) - TOLERANCE;
	minpt.z = fmin(fmin(fmin(fmin(fmin(fmin(fmin(p[0].z, p[1].z), p[2].z), p[3].z), p[4].z), p[5].z), p[6].z), p[7].z) - TOLERANCE;
	bb.clear();
    bb.addPoint( maxpt );
    bb.addPoint( minpt );
}

double CylinderVolume::dist(const GLVertex& p) const {
    GLVertex rotated_p = p - rotationCenter;
    rotated_p = rotated_p.rotateCBA(-angle.x, -angle.y, -angle.z) + rotationCenter;
    GLVertex tb = rotated_p - center;
    GLVertex tt = rotated_p - (center + GLVertex(0.0, 0.0, length));
    double d = (rotated_p - GLVertex(center.x, center.y, rotated_p.z)).norm();

    if (tb.z >= 0.0 && tt.z <= 0.0)
    	return ((radius - d < tb.z) && (radius - d < -tt.z)) ? radius - d : (tb.z < -tt.z) ? tb.z : -tt.z;  // positive inside. negative outside.
    else if (tb.z < 0.0) {
		// if we are under the cylinder, then return distance to flat cylinder bottom
    	if (d < radius)
    		return tb.z;
    	else {
    		// outside the cylinder, return a distance to the outer lower "ring" of the cylinder
    		GLVertex n = GLVertex(tb.x, tb.y, 0.0);
    		n = n * (radius / d);  // 1/d means normalization
			return -((tb - n).norm());
    	}
    } else {
		 // if we are above the cylinder, then return distance to flat cylinder top
    	if (d < radius)
			    return -tt.z;
    	else {
			  // outside the cylinder, return a distance to the outer upper "ring" of the cylinder
  		    GLVertex n = GLVertex(tt.x, tt.y, 0.0);
  		    n = n * (radius / d);  // 1/d means normalization
  		    return -((tt - n).norm());
  	   }
   }
}

//************* STL **************/

StlVolume::StlVolume() {
	type = STL_VOLUME;
    center = GLVertex(0, 0, 0);	// center is treated as the origin's offset of STL
    rotationCenter = GLVertex(0, 0, 0);
    angle  = GLVertex(0, 0, 0);
    cube_resolution = 0.0;
}

void StlVolume::calcBB() {
    for (int i = 0; i < (int)facets.size(); i++) {
    	facets[i]->v1 += center; facets[i]->v2 += center; facets[i]->v3 += center;
    	facets[i]->normal = facets[i]->normal.rotateABC(angle.x, angle.y, angle.z);
std::cout << "normal x: " << facets[i]->normal.x << " y: " << facets[i]->normal.y << " z: " << facets[i]->normal.z << "\n";
    	GLVertex v1p = facets[i]->v1 - rotationCenter;
    	facets[i]->v1 = v1p.rotateABC(angle.x, angle.y, angle.z) + rotationCenter;
std::cout << "vertex v1: " << facets[i]->v1.x << " y: " << facets[i]->v1.y << " z: " << facets[i]->v1.z << "\n";
    	GLVertex v2p = facets[i]->v2 - rotationCenter;
    	facets[i]->v2 = v2p.rotateABC(angle.x, angle.y, angle.z) + rotationCenter;
std::cout << "vertex v2: " << facets[i]->v2.x << " y: " << facets[i]->v2.y << " z: " << facets[i]->v2.z << "\n";
    	GLVertex v3p = facets[i]->v3 - rotationCenter;
    	facets[i]->v3 = v3p.rotateABC(angle.x, angle.y, angle.z) + rotationCenter;
std::cout << "vertex v3: " << facets[i]->v3.x << " y: " << facets[i]->v3.y << " z: " << facets[i]->v3.z << "\n";

        GLVertex normal = (facets[i]->v2 - facets[i]->v1).cross(facets[i]->v3 - facets[i]->v1);
        normal.normalize();
//assert ((facets[i]->normal - normal).norm() < CALC_TOLERANCE);
        facets[i]->normal = normal;
    }
    if (facets.size()) {
        maxpt.x = fmax(fmax(facets[0]->v1.x, facets[0]->v2.x),facets[0]->v3.x);
        maxpt.y = fmax(fmax(facets[0]->v1.y, facets[0]->v2.y),facets[0]->v3.y);
        maxpt.z = fmax(fmax(facets[0]->v1.z, facets[0]->v2.z),facets[0]->v3.z);
        minpt.x = fmin(fmin(facets[0]->v1.x, facets[0]->v2.x),facets[0]->v3.x);
        minpt.y = fmin(fmin(facets[0]->v1.y, facets[0]->v2.y),facets[0]->v3.y);
        minpt.z = fmin(fmin(facets[0]->v1.z, facets[0]->v2.z),facets[0]->v3.z);
    }
    for (int i = 0; i < (int)facets.size(); i++) {
        maxpt.x = fmax(fmax(fmax(facets[i]->v1.x, facets[i]->v2.x),facets[i]->v3.x), maxpt.x);
        maxpt.y = fmax(fmax(fmax(facets[i]->v1.y, facets[i]->v2.y),facets[i]->v3.y), maxpt.y);
        maxpt.z = fmax(fmax(fmax(facets[i]->v1.z, facets[i]->v2.z),facets[i]->v3.z), maxpt.z);
        minpt.x = fmin(fmin(fmin(facets[i]->v1.x, facets[i]->v2.x),facets[i]->v3.x), minpt.x);
        minpt.y = fmin(fmin(fmin(facets[i]->v1.y, facets[i]->v2.y),facets[i]->v3.y), minpt.y);
        minpt.z = fmin(fmin(fmin(facets[i]->v1.z, facets[i]->v2.z),facets[i]->v3.z), minpt.z);
        V21.push_back(facets[i]->v2 - facets[i]->v1);
        V21invV21dotV21.push_back((facets[i]->v2 - facets[i]->v1) * (1.0/(facets[i]->v2 - facets[i]->v1).dot(facets[i]->v2 - facets[i]->v1)));
        V32.push_back(facets[i]->v3 - facets[i]->v2);
        V32invV32dotV32.push_back((facets[i]->v3 - facets[i]->v2) * (1.0/(facets[i]->v3 - facets[i]->v2).dot(facets[i]->v3 - facets[i]->v2)));
        V13.push_back(facets[i]->v1 - facets[i]->v3);
        V13invV13dotV13.push_back((facets[i]->v1 - facets[i]->v3) * (1.0/(facets[i]->v1 - facets[i]->v3).dot(facets[i]->v1 - facets[i]->v3)));
    }
    bb.clear();
    maxpt += GLVertex(cube_resolution, cube_resolution, cube_resolution) * 2.0;
    minpt -= GLVertex(cube_resolution, cube_resolution, cube_resolution) * 2.0;
std::cout << "STL maxpt x:" << maxpt.x << " y: " << maxpt.y << " z:" << maxpt.z  << "\n";
std::cout << "STL minpt x:" << minpt.x << " y: " << minpt.y << " z:" << minpt.z  << "\n";
    bb.addPoint( maxpt );
    bb.addPoint( minpt );

    maxlength = fmax(fmax(maxpt.x - minpt.x, maxpt.y - minpt.y), maxpt.z - minpt.z);
std::cout << "maxlength:" << maxlength << "\n";

std::cout << "cube_resolution:" << cube_resolution << "\n";

	typedef struct cube {
		double upper_limit;
		double indexcubesize;
		int max_x, max_y, max_z;
		int ratio = 2;
		bool from_facets = false;
	} CUBE;

	std::vector<CUBE> cubes;
	CUBE cube;

	int start_ratio = ROUGH_RATIO;
	assert(start_ratio >= 0); assert((start_ratio % 2) == 0);
	int end_ratio = 0;
	assert(end_ratio >= 0); assert((end_ratio == 1) || (end_ratio % 2) == 0);
	assert(start_ratio > end_ratio);
	int devide_ratio = 2;
	assert((devide_ratio % 2) == 0);
	for (int i = start_ratio; i > end_ratio; i /= devide_ratio) {
		cube.upper_limit = 0.0;
		for (int j = i; j > end_ratio; j /= devide_ratio) {
			cube.upper_limit += maxlength / (MAX_INDEX/j);
		}
		cube.upper_limit += cube_resolution * 2.0;
		cube.upper_limit *= 2.0 + TOLERANCE;
		cube.indexcubesize = maxlength / (MAX_INDEX/i);
		cube.ratio = devide_ratio;
		cube.max_x = MAX_X/i; cube.max_y = MAX_Y/i; cube.max_z = MAX_Z/i;
		assert(cube.max_x > 0); assert(cube.max_y > 0); assert(cube.max_z > 0);
		cube.from_facets = false;
		cubes.push_back(cube);
	}
	cubes[0].from_facets = true;

#ifdef MULTI_THREAD_STL_NEIGHBOR
	int threadNum = QThreadPool::globalInstance()->maxThreadCount();
	std::cout << "Thread Num: " << threadNum << "\n";
	QFuture<void> future[threadNum];
#endif

    for (int i = 0; i < (int)cubes.size(); i++) {
		double upper_limit = cubes[i].upper_limit;
		indexcubesize = cubes[i].indexcubesize;
		ratio = cubes[i].ratio;
		int max_x = cubes[i].max_x, max_y = cubes[i].max_y, max_z = cubes[i].max_z;
		bool from_facets = cubes[i].from_facets;

        std::cout << "STL Facets Searching Level " << i+1 << "/" << (int)cubes.size();
		swapIndex();
emit signalProgressFeature(QString("STL Facets Searching Level ") + QString::number(i+1) + QString("/") + QString::number((int)cubes.size()), 0, max_x);
setProgress(0);
#ifdef MULTI_THREAD_STL_NEIGHBOR
        for (int index_x = 0; index_x < max_x; index_x++) {
            for (int index_y = 0; index_y < max_y; index_y++) {
                for (int index_z = 0; index_z < max_z;) {
                    for (int i = 0; i < threadNum; i++)
					{
						future[i] = QtConcurrent::run(this, &StlVolume::calcNeighborhoodIndex, index_x, index_y, index_z++, upper_limit, from_facets);
					}
                    for (int i = 0; i < threadNum; i++)
					{
						future[i].waitForFinished();
					}
#else
        for (int index_x = 0; index_x < max_x; index_x++) {
            for (int index_y = 0; index_y < max_y; index_y++) {
                for (int index_z = 0; index_z < max_z; index_z++) {
					calcNeighborhoodIndex(index_x, index_y, index_z, upper_limit, from_facets);
#endif
				}
			}
			std::cout << "." << std::flush;
accumlateProgress(1);
sendProgress();
		}
		std::cout << "done.\n" << std::flush;
	}

	invcubesize = 1.0 / indexcubesize;

    for (int index_x = 0; index_x < MAX_X; index_x++)
        for (int index_y = 0; index_y < MAX_Y; index_y++)
            for (int index_z = 0; index_z < MAX_Z; index_z++)
				neighborhoodIndex[src_index][index_x][index_y][index_z].resize(0);
}

void StlVolume::calcNeighborhoodIndex(int index_x, int index_y, int index_z, double upper_limit, bool from_facets)
{
    double x = minpt.x + indexcubesize * index_x + indexcubesize / 2.0;
    if (x > maxpt.x + indexcubesize / 2.0) return;
    double y = minpt.y + indexcubesize * index_y + indexcubesize / 2.0;
    if (y > maxpt.y + indexcubesize / 2.0) return;
    double z = minpt.z + indexcubesize * index_z + indexcubesize / 2.0;
    if (z > maxpt.z + indexcubesize / 2.0) return;

	bool find = false;
	double min = 1.0e+6;
	double ret;
	GLVertex p = GLVertex(x, y, z);
	typedef struct { int index; double dist; } CANDIDATE;
	std::vector<CANDIDATE> candidate;

	if (from_facets == true) {
        for (int i = 0; i < (int)facets.size(); i++) {
			ret = distance(p, i);
    		if (ret <= upper_limit) {
    			find = true;
    			neighborhoodIndex[dst_index][index_x][index_y][index_z].push_back(i);
    		} else {
    			if (ret <= min) {
    				min = ret;
                }
                if (ret < (min + indexcubesize)) {
    				CANDIDATE data = { i, ret };
    				candidate.push_back(data);
    			}
    		}
		}
	} else {
		neighborhoodIndex[dst_index][index_x][index_y][index_z].resize(0);

		int index_size = (int)neighborhoodIndex[src_index][index_x/ratio][index_y/ratio][index_z/ratio].size();
		if (index_size == 0) return;
        for (int ic = 0; ic < index_size; ic++) {
            int i = neighborhoodIndex[src_index][index_x/ratio][index_y/ratio][index_z/ratio][ic];
            ret = distance(p, i);
			if (ret <= upper_limit) {
				find = true;
				neighborhoodIndex[dst_index][index_x][index_y][index_z].push_back(i);
			} else {
				if (ret <= min) {
					min = ret;
                }
                if (ret < (min + indexcubesize)) {
					CANDIDATE data = { i, ret };
					candidate.push_back(data);
				}
			}
		}
	}
	if (find == false) {
		double search_distance = sqrt(min * min + indexcubesize * indexcubesize) + TOLERANCE;
        for (int ic = 0; ic < (int)candidate.size(); ic++) {
            int i = candidate[ic].index;
			if (candidate[ic].dist <= search_distance) {
				neighborhoodIndex[dst_index][index_x][index_y][index_z].push_back(i);
			}
		}
	}
}

double StlVolume::distance(const GLVertex& p, int index) {
	GLVertex q, r;
	GLVertex n1, n2, n3;
    double s12, s23, s31;
	double d, u, abs_d;

	u = (p - facets[index]->v1).dot(V21invV21dotV21[index]);
	q = facets[index]->v1 + V21[index] * u;
	d = (q - p).dot(facets[index]->normal);
	r = p + facets[index]->normal * d;
	n1 = (r - facets[index]->v1).cross(V13[index]);
	n2 = (r - facets[index]->v2).cross(V21[index]);
	n3 = (r - facets[index]->v3).cross(V32[index]);
    s12 = n1.dot(n2); s23 = n2.dot(n3); s31 = n3.dot(n1);

    if ((s12 > 0.0) && (s23 > 0.0) && (s31 > 0.0)) {
        return fabs(d);
	}

	double abs_d12, abs_d13, abs_d32;
	GLVertex q12, q13, q32;

	if (u <= 0.0)
		q12 = facets[index]->v1;
	else if (u >= 1.0)
		q12 = facets[index]->v2;
	else
		/* q = facets[index]->v1 + V21[index] * u */
		q12 = q;
	abs_d12 = (q12 - p).norm();

	u = (p - facets[index]->v3).dot(V13invV13dotV13[index]);
	if (u <= 0.0)
		q13 = facets[index]->v3;
	else if (u >= 1.0)
		q13 = facets[index]->v1;
	else
		q13 = facets[index]->v3 + V13[index] * u;
	abs_d13 = (q13 - p).norm();

	u = (p - facets[index]->v2).dot(V32invV32dotV32[index]);
	if (u <= 0.0)
		q32 = facets[index]->v2;
	else if (u >= 1.0)
		q32 = facets[index]->v3;
	else
		q32 = facets[index]->v2 + V32[index] * u;
	abs_d32 = (q32 - p).norm();

	if ((abs_d12 <= abs_d13) && (abs_d12 <= abs_d32)) {
		abs_d = abs_d12;
	} else if ((abs_d13 < abs_d12) && (abs_d13 <= abs_d32)) {
		abs_d = abs_d13;
	} else {
		abs_d = abs_d32;
	}

	return abs_d;
}

double StlVolume::dist(const GLVertex& p) const {
	int index_x, index_y, index_z;
	double ret = -1.0e+6;

	index_x = (int)((p.x - minpt.x) * invcubesize);
	index_y = (int)((p.y - minpt.y) * invcubesize);
	index_z = (int)((p.z - minpt.z) * invcubesize);

	if ((index_x < 0 || index_x > MAX_X-1) || (index_y < 0 || index_y > MAX_Y-1) || (index_z < 0 || index_z > MAX_Z-1))
		return ret;

	int index_size = (int)neighborhoodIndex[dst_index][index_x][index_y][index_z].size();
	if (index_size == 0) return ret;

	GLVertex q, r;
	GLVertex n1, n2, n3;
	double s12, s23, s31;
    double min = 1.0e+6, dir, u, abs_d;

	double abs_d12, abs_d13, abs_d32;
	GLVertex q12, q13, q32;
	bool correction = false;

	enum StlSide  { INSIDE, OUTSIDE, UNDECIDED };
	StlSide side = UNDECIDED;
	typedef struct {int index; StlSide side; double abs_d; GLVertex q;} EDGE;
	EDGE selected = { 0, UNDECIDED, 1.0e+6, };
	EDGE second	  = { 0, UNDECIDED, 1.0e+6, };
	EDGE third	  = { 0, UNDECIDED, 1.0e+6, };
	std::vector<EDGE> more;
    int normal_vec_calc_count = 0;

//struct { double x, y, z; } target = {-11.500,23.000,23.000};
//struct { double x, y, z; } target = {-28.000,22.000,20.000};
//struct { double x, y, z; } target = {-21.500,29.000,26.500};
//struct { double x, y, z; } target = {16.000,22.000,30.000};
//struct { double x, y, z; } target = {-25.000,0.000,25.000};
//struct { double x, y, z; } target = {-29.8828,3.16406,29.1797};
//struct { double x, y, z; } target = {-22.5,-4.92188,23.9062};
struct { double x, y, z; } target = {-15.000, 7.000, 16.500};
GLVertex test_vector;


//#define DETAIL
//#define GRID    (5.0)
#define GRID    (2.0)
#ifdef DETAIL
if ((target.x-TOLERANCE < p.x) && (p.x < target.x+TOLERANCE) && (target.y-TOLERANCE < p.y) && (p.y < target.y+TOLERANCE) && (target.z-TOLERANCE < p.z) && (p.z < target.z+TOLERANCE)) {
#else
if ((target.x-GRID < p.x) && (p.x < target.x+GRID) && (target.y-GRID < p.y) && (p.y < target.y+GRID) && (target.z-GRID < p.z) && (p.z < target.z+GRID)) {
#endif
{
    std::cout << "target p(" << p.x << "," << p.y << "," << p.z << ") index_size: " << index_size << " " << std::flush;
    for (int ic = 0; ic < index_size; ic++) {
        int i = neighborhoodIndex[dst_index][index_x][index_y][index_z][ic];
        std::cout << "(" << i << ")";
    }
    std::cout << " " << std::flush;
}
}

    for (int ic = 0; ic < index_size; ic++) {
        int i = neighborhoodIndex[dst_index][index_x][index_y][index_z][ic];
		u = (p - facets[i]->v1).dot(V21invV21dotV21[i]);
		q = facets[i]->v1 + V21[i] * u;
        dir = (q - p).dot(facets[i]->normal);
        if ((abs_d = fabs(dir)) > min) continue;
        r = p + facets[i]->normal * dir;
		n1 = (r - facets[i]->v1).cross(V13[i]);
		n2 = (r - facets[i]->v2).cross(V21[i]);
		n3 = (r - facets[i]->v3).cross(V32[i]);
		s12 = n1.dot(n2); s23 = n2.dot(n3); s31 = n3.dot(n1);

		if ((s12 > 0.0) && (s23 > 0.0) && (s31 > 0.0)) {
            double candidate_min = abs_d - CALC_TOLERANCE;
            {
				double d1 = (facets[i]->v1 - p).norm(); double d2 = (facets[i]->v2 - p).norm(); double d3 = (facets[i]->v3 - p).norm();
				if ((candidate_min < d1) && (candidate_min < d2) && (candidate_min < d3)) { // sanity check..
					min = candidate_min;
                    ret = dir;
					selected.index = i;
					if (ret > 0.0)
						selected.side = INSIDE;
					else
						selected.side = OUTSIDE;
					selected.q = q;
					correction = false;
					continue;
				}
			}
		}

		if (u <= 0.0)
			q12 = facets[i]->v1;
		else if (u >= 1.0)
			q12 = facets[i]->v2;
		else
			/* q = facets[i]->v1 + V21[i] * u */
			q12 = q;
		abs_d12 = (q12 - p).norm();

		u = (p - facets[i]->v3).dot(V13invV13dotV13[i]);
		if (u <= 0.0)
			q13 = facets[i]->v3;
		else if (u >= 1.0)
			q13 = facets[i]->v1;
		else
			q13 = facets[i]->v3 + V13[i] * u;
		abs_d13 = (q13 - p).norm();

		u = (p - facets[i]->v2).dot(V32invV32dotV32[i]);
		if (u <= 0.0)
			q32 = facets[i]->v2;
		else if (u >= 1.0)
			q32 = facets[i]->v3;
		else
			q32 = facets[i]->v2 + V32[i] * u;
		abs_d32 = (q32 - p).norm();

		if ((abs_d12 <= abs_d13) && (abs_d12 <= abs_d32)) {
			q = q12;
			abs_d = abs_d12;
		} else if ((abs_d13 < abs_d12) && (abs_d13 <= abs_d32)) {
			q = q13;
			abs_d = abs_d13;
		} else {
			q = q32;
			abs_d = abs_d32;
		}

		if (abs_d >= min && abs_d > second.abs_d)
			continue;

        dir = (q - p).dot(facets[i]->normal);
        if (dir > 0.0)
			side = INSIDE;
		else
			side = OUTSIDE;

		if (abs_d < min) {
            min = abs_d;
            if (side == INSIDE) {
                ret = abs_d;
			} else {
                ret = -abs_d;
			}

            if (third.side != UNDECIDED)
				more.push_back(third);
			third.index = second.index;
			third.side = second.side;
			third.abs_d = second.abs_d;
			third.q = second.q;
			second.index = selected.index;
			second.side = selected.side;
			second.abs_d = selected.abs_d;
			second.q = selected.q;

			selected.index = i;
			selected.side = side;
			selected.abs_d = abs_d;
			selected.q = q;
			correction = true;
		} else {
			if (third.side != UNDECIDED)
				more.push_back(third);
			third.index = second.index;
			third.side = second.side;
			third.abs_d = second.abs_d;
			third.q = second.q;
			second.index = i;
			second.side = side;
			second.abs_d = abs_d;
			second.q = q;
		}
	}

    if ((second.side != UNDECIDED) && (third.side != UNDECIDED)) {
        if ((selected.q - third.q).norm() < (selected.q - second.q).norm()) {
            std::swap(second, third);
        }
        if (!more.empty()) {
            sort(more.begin(), more.end(), [selected](const EDGE& e1, const EDGE& e2) { return (selected.q - e1.q).norm() < (selected.q - e2.q).norm();} );
            if ((selected.q - more[0].q).norm() < (selected.q - third.q).norm()) {
                std::swap(third,  more[0]);
            }
        }
    }

	if (correction == true) {
        if (((second.side != UNDECIDED) && (second.side != selected.side)) || ((third.side != UNDECIDED) && (third.side != selected.side)) || !more.empty()) {
            if ((selected.q - second.q).norm() < CALC_TOLERANCE*100.0) {
//				std::cout << "correction ";
                GLVertex outer_vector = ((selected.q - Facet::facetCenter(facets[selected.index])).normalize() + (selected.q - Facet::facetCenter(facets[second.index])).normalize()).normalize();
                GLVertex normal_avg_vector = facets[selected.index]->normal + facets[second.index]->normal;
                normal_vec_calc_count++;
                if ((third.side != UNDECIDED) && ((selected.q - third.q).norm() < CALC_TOLERANCE*100.0)) {
                    outer_vector += (selected.q - Facet::facetCenter(facets[third.index])).normalize();
                    outer_vector.normalize();
                    normal_avg_vector += facets[third.index]->normal;
                    normal_vec_calc_count++;
                    if (more.empty()) {
//						std::cout << "triple\n";
                    } else {
//                        std::cout << "more triple...(" << more.size() << " facets)\n";
                        for (int i = 0; i < (int)more.size(); i++) {
                            if ((selected.q - more[i].q).norm() < CALC_TOLERANCE*100.0) {
                                outer_vector += (selected.q - Facet::facetCenter(facets[more[i].index])).normalize();
                                outer_vector.normalize();
                                normal_avg_vector += facets[more[i].index]->normal;
                                normal_vec_calc_count++;
                            }
                        }
					}
				} else {
//                    std::cout << "double\n" << std::flush;
				}
                normal_avg_vector.normalize();
                if (normal_avg_vector.dot(outer_vector) < 0.0) {
                    outer_vector.x = -outer_vector.x;
                    outer_vector.y = -outer_vector.y;
                    outer_vector.z = -outer_vector.z;
                }
test_vector = outer_vector;
                if (!(normal_vec_calc_count == 1 && (second.side == selected.side))) {
                    if (outer_vector.dot(selected.q - p) < 0.0) {
                        if (selected.side == INSIDE) {
//          				std::cout << "SIDE WRONG!! INSIDE -> OUTSIDE p(" << p.x <<"," << p.y << "," << p.z << ")\n" <<  std::flush;
                            ret = -selected.abs_d;
                        }
                    } else {
                        if (selected.side == OUTSIDE) {
//                          std::cout << "SIDE WRONG!! OUTSIDE -> INSIDE p(" << p.x <<"," << p.y << "," << p.z << ")\n" <<  std::flush;
                            ret = selected.abs_d;
                        }
                    }
                }
            } else {

            }
        }
    }

#ifdef DETAIL
if ((target.x-TOLERANCE < p.x) && (p.x < target.x+TOLERANCE) && (target.y-TOLERANCE < p.y) && (p.y < target.y+TOLERANCE) && (target.z-TOLERANCE < p.z) && (p.z < target.z+TOLERANCE))
#else
if ((target.x-GRID < p.x) && (p.x < target.x+GRID) && (target.y-GRID < p.y) && (p.y < target.y+GRID) && (target.z-GRID < p.z) && (p.z < target.z+GRID))
#endif
{
    std::cout << "normal_vec_calc_count :" << normal_vec_calc_count << "\n" << std::flush;

    if (third.side != UNDECIDED) {
        std::cout << "side(" << ((ret > 0) ? "INSIDE)" : "OUTSIDE):") << " facets:(" << selected.index << ")("<< second.index << ")("<< third.index << ")\n" << std::flush;
        std::cout << "(selected.q - second.q).norm :" << (selected.q - second.q).norm() << " selected.ads_d :" << selected.abs_d << " second.abs_d :" << second.abs_d << "\n" << std::flush;
        std::cout << "(selected.q - third.q).norm :" << (selected.q - third.q).norm() << " selected.ads_d :" << selected.abs_d << " third.abs_d :" << third.abs_d << "\n" << std::flush;
        std::cout << " selected.q:(" << selected.q.x << "," << selected.q.y << "," << selected.q.z << ")\n" << std::flush;
        std::cout << " second.q:  (" << second.q.x << "," << second.q.y << "," << second.q.z << ")\n" << std::flush;
        std::cout << " third.q:   (" << third.q.x << "," << third.q.y << "," << third.q.z << ")\n" << std::flush;
        std::cout << " test_vector:(" << test_vector.x << "," << test_vector.y << "," << test_vector.z << ")\n" << std::flush;
        std::cout << " facets normal:(" << selected.index << ")(" << facets[selected.index]->normal.x << "," << facets[selected.index]->normal.y << "," << facets[selected.index]->normal.z << ")\n" << std::flush;
        std::cout << " facets normal:(" << second.index << ")(" << facets[second.index]->normal.x << "," << facets[second.index]->normal.y << "," << facets[second.index]->normal.z << ")\n" << std::flush;
        std::cout << " facets normal:(" << third.index << ")(" << facets[third.index]->normal.x << "," << facets[third.index]->normal.y << "," << facets[third.index]->normal.z << ")\n" << std::flush;
        std::cout << " facets:(" << selected.index << ")(" << facets[selected.index]->v1.x << "," << facets[selected.index]->v1.y << "," << facets[selected.index]->v1.z << ")(" << facets[selected.index]->v2.x << "," << facets[selected.index]->v2.y << "," << facets[selected.index]->v2.z << ")(" << facets[selected.index]->v3.x << "," << facets[selected.index]->v3.y << "," << facets[selected.index]->v3.z << ")\n" << std::flush;
        std::cout << " facets:(" << second.index << ")(" << facets[second.index]->v1.x << "," << facets[second.index]->v1.y << "," << facets[second.index]->v1.z << ")(" << facets[second.index]->v2.x << "," << facets[second.index]->v2.y << "," << facets[second.index]->v2.z << ")(" << facets[second.index]->v3.x << "," << facets[second.index]->v3.y << "," << facets[second.index]->v3.z << ")\n" << std::flush;
        std::cout << " facets:(" << third.index << ")(" << facets[third.index]->v1.x << "," << facets[third.index]->v1.y << "," << facets[third.index]->v1.z << ")(" << facets[third.index]->v2.x << "," << facets[third.index]->v2.y << "," << facets[third.index]->v2.z << ")(" << facets[third.index]->v3.x << "," << facets[third.index]->v3.y << "," << facets[third.index]->v3.z << ")\n" << std::flush;
        if (!more.empty()) {
            for (unsigned int i = 0; i < more.size(); i++) {
                std::cout << "  (selected.q - more[" << i << "].q).norm :" << (selected.q - more[i].q).norm() << " selected.ads_d :" << selected.abs_d << " more[" << i << "].abd_d :" << more[i].abs_d << "\n" << std::flush;
                std::cout << "   facets normal: more[" << i << "](" << more[i].index << ")(" << facets[more[i].index]->normal.x << "," << facets[more[i].index]->normal.y << "," << facets[more[i].index]->normal.z << ")\n" << std::flush;
                std::cout << "   facets: more[" << i << "](" << more[i].index << ")("<< facets[more[i].index]->v1.x << "," << facets[more[i].index]->v1.y << "," << facets[more[i].index]->v1.z << ")(" << facets[more[i].index]->v2.x << "," << facets[more[i].index]->v2.y << "," << facets[more[i].index]->v2.z << ")(" << facets[more[i].index]->v3.x << "," << facets[more[i].index]->v3.y << "," << facets[more[i].index]->v3.z << ")\n" << std::flush;
            }
        }
    } else {
        std::cout << "side(" << ((ret > 0) ? "INSIDE)" : "OUTSIDE):") << " facets:(" << selected.index << ")("<< second.index << ")\n" << std::flush;
        std::cout << "(selected.q - second.q).norm :" << (selected.q - second.q).norm() << " selected.ads :" << selected.abs_d << " second.abs_d :" << second.abs_d << "\n" << std::flush;
        std::cout << " selected.q:(" << selected.q.x << "," << selected.q.y << "," << selected.q.z << ")\n" << std::flush;
        std::cout << " second.q:  (" << second.q.x << "," << second.q.y << "," << second.q.z << ")\n" << std::flush;
        std::cout << " test_vector:(" << test_vector.x << "," << test_vector.y << "," << test_vector.z << ")\n" << std::flush;
        std::cout << " facets normal:(" << selected.index << ")(" << facets[selected.index]->normal.x << "," << facets[selected.index]->normal.y << "," << facets[selected.index]->normal.z << ")\n" << std::flush;
        std::cout << " facets normal:(" << second.index << ")(" << facets[second.index]->normal.x << "," << facets[second.index]->normal.y << "," << facets[second.index]->normal.z << ")\n" << std::flush;
        std::cout << " facets:(" << selected.index << ")(" << facets[selected.index]->v1.x << "," << facets[selected.index]->v1.y << "," << facets[selected.index]->v1.z << ")(" << facets[selected.index]->v2.x << "," << facets[selected.index]->v2.y << "," << facets[selected.index]->v2.z << ")(" << facets[selected.index]->v3.x << "," << facets[selected.index]->v3.y << "," << facets[selected.index]->v3.z << ")\n" << std::flush;
        std::cout << " facets:(" << second.index << ")(" << facets[second.index]->v1.x << "," << facets[second.index]->v1.y << "," << facets[second.index]->v1.z << ")(" << facets[second.index]->v2.x << "," << facets[second.index]->v2.y << "," << facets[second.index]->v2.z << ")(" << facets[second.index]->v3.x << "," << facets[second.index]->v3.y << "," << facets[second.index]->v3.z << ")\n" << std::flush;
    }
}
	return ret;		// positive inside. negative outside.
}


//************* CutterVolume **************/

CutterVolume::CutterVolume() {
    radius = 0.0;
    length = 0.0;
    enableholder = false;
    holderradius = 0.0;
    holderlength = 0.0;
}

void CutterVolume::calcBBHolder() {
    bbHolder.clear();
    GLVertex maxpt = GLVertex(center.x + holderradius + TOLERANCE, center.y + holderradius + TOLERANCE, center.z + length + holderlength + TOLERANCE);
    GLVertex minpt = GLVertex(center.x - holderradius - TOLERANCE, center.y - holderradius - TOLERANCE, center.z + length - TOLERANCE);
    bbHolder.addPoint( maxpt );
    bbHolder.addPoint( minpt );
}

//************* CylCutterVolume **************/

CylCutterVolume::CylCutterVolume() {
    radius = 0.0;
    length = 0.0;
}

void CylCutterVolume::calcBB() {
    bb.clear();
    GLVertex maxpt = GLVertex(center.x + maxradius + TOLERANCE, center.y + maxradius + TOLERANCE, center.z + length + TOLERANCE);
    GLVertex minpt = GLVertex(center.x - maxradius - TOLERANCE, center.y - maxradius - TOLERANCE, center.z - TOLERANCE);
    bb.addPoint( maxpt );
    bb.addPoint( minpt );
    if (enableholder)
        calcBBHolder();
}

double CylCutterVolume::dist(const GLVertex& p) const {
#ifdef MULTI_AXIS
	GLVertex rotated_p = p;
	rotated_p = rotated_p.rotateAC(angle.x, angle.z);
	GLVertex t = rotated_p - center;
	double d = (rotated_p - GLVertex(center.x, center.y, rotated_p.z)).norm();
#else
    GLVertex t = p - center;
    double d = (p - GLVertex(center.x, center.y, p.z)).norm();
#endif

    if (t.z >= 0.0) {
    	  return t.z > length ? holderradius - d : radius - d;  // positive inside. negative outside.
    } else {
		 // if we are under the cutter, then return distance to flat cutter bottom
		if (d < radius)
			return t.z;
		else {
			// outside the cutter, return a distance to the outer "ring" of the cutter
			GLVertex n = GLVertex(t.x, t.y, 0.0);
			n = n * (radius / d);  // 1/d means normalization
			return -((t - n).norm());
		}
	}
}

Cutting CylCutterVolume::dist_cd(const GLVertex& p) const {
#ifdef MULTI_AXIS
    GLVertex rotated_p = p;
    rotated_p = rotated_p.rotateAC(angle.x, angle.z);
    GLVertex t = rotated_p - center;
    double d = (rotated_p - GLVertex(center.x, center.y, rotated_p.z)).norm();
#else
    GLVertex t = p - center;
    double d = (p - GLVertex(center.x, center.y, p.z)).norm();
#endif
    Cutting result = { t.z, NO_COLLISION, 1 };
    double rdiff = radius - d;

    if (t.z >= 0.0) {
    	  result.collision |= ((t.z > flutelength) && ((rdiff = neckradius  - d)  > COLLISION_TOLERANCE)) ? NECK_COLLISION   : NO_COLLISION;
    	  result.collision |= ((t.z > reachlength) && ((rdiff = shankradius - d)  > COLLISION_TOLERANCE)) ? SHANK_COLLISION  : NO_COLLISION;
    	  result.collision |= ((t.z > length)      && ((rdiff = holderradius - d) > COLLISION_TOLERANCE)) ? HOLDER_COLLISION : NO_COLLISION;
    	  result.f = rdiff < t.z ? rdiff : t.z;  // positive inside. negative outside.
if (rdiff > effective_radius) result.count = 4;
    	  return result;
    } else {
		    // if we are under the cutter, then return distance to flat cutter bottom
    	  if (d < radius)
			   return result;
    	  else {
			  // outside the cutter, return a distance to the outer "ring" of the cutter
    		   GLVertex n = GLVertex(t.x, t.y, 0.0);
    		   n = n * (radius / d);  // 1/d means normalization
    		   result.f = -((t - n).norm());
    		   return result;
    	  }
   }
}

//************* BallCutterVolume **************/

BallCutterVolume::BallCutterVolume() {
    radius = 0.0;
    length = 0.0;
}

void BallCutterVolume::calcBB() {
    bb.clear();
    GLVertex maxpt = GLVertex(center.x + maxradius + TOLERANCE, center.y + maxradius + TOLERANCE, center.z + length + TOLERANCE);
    GLVertex minpt = GLVertex(center.x - maxradius - TOLERANCE, center.y - maxradius - TOLERANCE, center.z - radius - TOLERANCE);
    bb.addPoint( maxpt );
    bb.addPoint( minpt );
    if (enableholder)
        calcBBHolder();
}

double BallCutterVolume::dist(const GLVertex& p) const {
#ifdef MULTI_AXIS
    GLVertex rotated_p = p;
    rotated_p = rotated_p.rotateAC(angle.x, angle.z);
    GLVertex t = rotated_p - center;
#else
    GLVertex t = p - center;
#endif

    if (t.z < 0.0)
#ifdef MULTI_AXIS
      return radius - (rotated_p - GLVertex(center.x, center.y, center.z)).norm(); // positive inside. negative outside.
#else
      return radius - (center - p).norm(); // positive inside. negative outside.
#endif
    else
#ifdef MULTI_AXIS
    return radius - (rotated_p - GLVertex(center.x, center.y, rotated_p.z)).norm();
#else
    return radius - (p - GLVertex(center.x, center.y, p.z)).norm();
#endif
}

Cutting BallCutterVolume::dist_cd(const GLVertex& p) const {
#ifdef MULTI_AXIS
    GLVertex rotated_p = p;
    rotated_p = rotated_p.rotateAC(angle.x, angle.z);
    GLVertex t = rotated_p - center;
#else
    GLVertex t = p - center;
#endif
    Cutting result = { 0.0, NO_COLLISION, 1 };

    if (t.z < 0) {
    	result.f = radius - t.norm();
    	return result;
    } else {
#ifdef MULTI_AXIS
    	double d = (rotated_p - GLVertex(center.x, center.y, rotated_p.z)).norm();
#else
    	double d = (p - GLVertex(center.x, center.y, p.z)).norm();
#endif
    	result.f = radius - d;
    	result.collision |= ((t.z > flutelength) && ((result.f = neckradius   - d) > COLLISION_TOLERANCE)) ? NECK_COLLISION   : NO_COLLISION;
    	result.collision |= ((t.z > reachlength) && ((result.f = shankradius  - d) > COLLISION_TOLERANCE)) ? SHANK_COLLISION  : NO_COLLISION;
    	result.collision |= ((t.z > length)      && ((result.f = holderradius - d) > COLLISION_TOLERANCE)) ? HOLDER_COLLISION : NO_COLLISION;
    	return result;
    }
}

//************* BullCutterVolume **************/

BullCutterVolume::BullCutterVolume() {
    radius = 0.0;
    length = 0.0;
    r1 = r2 = 0.0;
}

void BullCutterVolume::calcBB() {
    bb.clear();
    GLVertex maxpt = GLVertex(center.x + maxradius + TOLERANCE, center.y + maxradius + TOLERANCE, center.z + length + TOLERANCE);
    GLVertex minpt = GLVertex(center.x - maxradius - TOLERANCE, center.y - maxradius - TOLERANCE, center.z - TOLERANCE);
    bb.addPoint( maxpt );
    bb.addPoint( minpt );
    if (enableholder)
        calcBBHolder();
}

double BullCutterVolume::dist(const GLVertex& p) const {
#ifdef MULTI_AXIS
	GLVertex rotated_p = p;
	rotated_p = rotated_p.rotateAC(angle.x, angle.z);
	GLVertex t = rotated_p - center;
	double d = (rotated_p - GLVertex(center.x, center.y, rotated_p.z)).norm();
#else
    GLVertex t = p - center;
    double d = (p - GLVertex(center.x, center.y, p.z)).norm();
#endif

    if (t.z >= r2) {  // cylindrical part, above toroid
    	  return t.z > length ? holderradius - d : radius - d;  // positive inside. negative outside.
    } else {
		 // if we are under the cutter, then return distance to flat cutter bottom
		if (d < r1)  // cylindrical part, inside toroid
			return t.z;
		else {
			// toroid
			GLVertex n = GLVertex(t.x, t.y, 0.0);
			n = n * (r1 / d);  // 1/d means normalization
			n += GLVertex(0.0, 0.0, r2);
			return -((t - n).norm()) + r2;
		}
	}
}

Cutting BullCutterVolume::dist_cd(const GLVertex& p) const {
#ifdef MULTI_AXIS
	GLVertex rotated_p = p;
	rotated_p = rotated_p.rotateAC(angle.x, angle.z);
	GLVertex t = rotated_p - center;
	double d = (rotated_p - GLVertex(center.x, center.y, rotated_p.z)).norm();
#else
    GLVertex t = p - center;
    double d = (p - GLVertex(center.x, center.y, p.z)).norm();
#endif
    Cutting result = { t.z, NO_COLLISION, 1 };
    double rdiff = radius - d;

    if (t.z >= r2) {  // cylindrical part, above toroid
  	  result.collision |= ((t.z > flutelength) && ((rdiff = neckradius   - d) > COLLISION_TOLERANCE)) ? NECK_COLLISION   : NO_COLLISION;
  	  result.collision |= ((t.z > reachlength) && ((rdiff = shankradius  - d) > COLLISION_TOLERANCE)) ? SHANK_COLLISION  : NO_COLLISION;
  	  result.collision |= ((t.z > length)      && ((rdiff = holderradius - d) > COLLISION_TOLERANCE)) ? HOLDER_COLLISION : NO_COLLISION;
  	  result.f = rdiff < t.z ? rdiff : t.z;  // positive inside. negative outside.
  	  return result;
    } else {
    	// if we are under the cutter, then return distance to flat cutter bottom
			if (d < r1)  // cylindrical part, inside toroid
				return result;
			else {
				// toroid
				GLVertex n = GLVertex(t.x, t.y, 0.0);
				n = n * (r1 / d);  // 1/d means normalization
				n += GLVertex(0.0, 0.0, r2);
				result.f = -((t - n).norm()) + r2;
				return result;
			}
		}
}

//************* ConeCutterVolume **************/

ConeCutterVolume::ConeCutterVolume() {
    radius = 0.0;
    length = 0.0;
    flutelength = 0.0;
    r1 = r2 = 0.0;
    incline_coff = 0.0;
}

void ConeCutterVolume::calcBB() {
    bb.clear();
    GLVertex maxpt = GLVertex(center.x + maxradius + TOLERANCE, center.y + maxradius + TOLERANCE, center.z + length + TOLERANCE);
    GLVertex minpt = GLVertex(center.x - maxradius - TOLERANCE, center.y - maxradius - TOLERANCE, center.z - TOLERANCE);
    bb.addPoint( maxpt );
    bb.addPoint( minpt );
    if (enableholder)
        calcBBHolder();
}

double ConeCutterVolume::dist(const GLVertex& p) const {
#ifdef MULTI_AXIS
	GLVertex rotated_p = p;
	rotated_p = rotated_p.rotateAC(angle.x, angle.z);
	GLVertex t = rotated_p - center;
	double d = (rotated_p - GLVertex(center.x, center.y, rotated_p.z)).norm();
#else
    GLVertex t = p - center;
    double d = (p - GLVertex(center.x, center.y, p.z)).norm();
#endif
    double rdiff;

    if (0.0 <= t.z && t.z <= flutelength)
    	rdiff = (incline_coff * t.z + r1) - d;
    else
    	rdiff = radius - d;

    if (t.z >= 0.0) {
    	  return t.z > length ? holderradius - d : rdiff;  // positive inside. negative outside.
    } else {
		 // if we are under the cutter, then return distance to flat cutter bottom
		if (d < r1)
			return t.z;
		else {
			// toroid
			GLVertex n = GLVertex(t.x, t.y, 0.0);
			n = n * (r1 / d);  // 1/d means normalization
			n += GLVertex(0.0, 0.0, r2);
			return -((t - n).norm());
		}
	}
}

Cutting ConeCutterVolume::dist_cd(const GLVertex& p) const {
#ifdef MULTI_AXIS
	GLVertex rotated_p = p;
	rotated_p = rotated_p.rotateAC(angle.x, angle.z);
	GLVertex t = rotated_p - center;
	double d = (rotated_p - GLVertex(center.x, center.y, rotated_p.z)).norm();
#else
    GLVertex t = p - center;
    double d = (p - GLVertex(center.x, center.y, p.z)).norm();
#endif
    Cutting result = { t.z, NO_COLLISION, 1 };
    double rdiff;

    if (0.0 <= t.z && t.z <= flutelength)
    	rdiff = (incline_coff * t.z + r1) - d;
    else
    	rdiff = radius - d;

    if (t.z >= 0.0) {
    	result.collision |= ((t.z > flutelength) && ((rdiff = neckradius   - d) > COLLISION_TOLERANCE)) ? NECK_COLLISION   : NO_COLLISION;
    	result.collision |= ((t.z > reachlength) && ((rdiff = shankradius  - d) > COLLISION_TOLERANCE)) ? SHANK_COLLISION  : NO_COLLISION;
    	result.collision |= ((t.z > length)      && ((rdiff = holderradius - d) > COLLISION_TOLERANCE)) ? HOLDER_COLLISION : NO_COLLISION;
    	result.f = rdiff < t.z ? rdiff : t.z;  // positive inside. negative outside.
    	return result;
    } else {
		// if we are under the cutter, then return distance to flat cutter bottom
    	if (d < r1)
    		return result;
    	else {
    		// outside the cutter, return a distance to the outer "ring" of the cutter
    		GLVertex n = GLVertex(t.x, t.y, 0.0);
    		n = n * (r1 / d);  // 1/d means normalization
    		result.f = -((t - n).norm());
    		return result;
    	}
    }
}

//************* DrillVolume **************/

DrillVolume::DrillVolume() {
    radius = 0.0;
    length = 0.0;
    flutelength = 0.0;
    tip_hight = 0.0;
    tip_angle = 118.0;
    incline_coff = 0.0;
}

void DrillVolume::calcBB() {
    bb.clear();
    GLVertex maxpt = GLVertex(center.x + maxradius + TOLERANCE, center.y + maxradius + TOLERANCE, center.z + length + TOLERANCE);
    GLVertex minpt = GLVertex(center.x - maxradius - TOLERANCE, center.y - maxradius - TOLERANCE, center.z - TOLERANCE);
    bb.addPoint( maxpt );
    bb.addPoint( minpt );
    if (enableholder)
        calcBBHolder();
}

double DrillVolume::dist(const GLVertex& p) const {
#ifdef MULTI_AXIS
	GLVertex rotated_p = p;
	rotated_p = rotated_p.rotateAC(angle.x, angle.z);
	GLVertex t = rotated_p - center;
	double d = (rotated_p - GLVertex(center.x, center.y, rotated_p.z)).norm();
#else
    GLVertex t = p - center;
    double d = (p - GLVertex(center.x, center.y, p.z)).norm();
#endif
    double rdiff;

    if (0.0 <= t.z && t.z <= tip_hight)
    	rdiff = incline_coff * t.z - d;
    else
    	rdiff = radius - d;

    if (t.z >= 0.0) {
    	return t.z > length ? holderradius - d : rdiff;  // positive inside. negative outside.
    } else {
    	return t.z;
	}
}

Cutting DrillVolume::dist_cd(const GLVertex& p) const {
#ifdef MULTI_AXIS
    GLVertex rotated_p = p;
    rotated_p = rotated_p.rotateAC(angle.x, angle.z);
    GLVertex t = rotated_p - center;
    double d = (rotated_p - GLVertex(center.x, center.y, rotated_p.z)).norm();
#else
    GLVertex t = p - center;
    double d = (p - GLVertex(center.x, center.y, p.z)).norm();
#endif
    Cutting result = { t.z, NO_COLLISION, 1 };
    double rdiff = radius - d;

    if (0.0 <= t.z && t.z <= tip_hight)
    	rdiff = incline_coff * t.z - d;

    if (t.z >= 0.0) {
    	  result.collision |= ((t.z > flutelength) && ((rdiff = neckradius  - d)  > COLLISION_TOLERANCE)) ? NECK_COLLISION   : NO_COLLISION;
    	  result.collision |= ((t.z > reachlength) && ((rdiff = shankradius - d)  > COLLISION_TOLERANCE)) ? SHANK_COLLISION  : NO_COLLISION;
    	  result.collision |= ((t.z > length)      && ((rdiff = holderradius - d) > COLLISION_TOLERANCE)) ? HOLDER_COLLISION : NO_COLLISION;
    	  result.f = rdiff < t.z ? rdiff : t.z;  // positive inside. negative outside.
//if (rdiff > effective_radius) result.count = 4;
    	  return result;
    } else {
		    // if we are under the cutter, then return distance to flat cutter bottom
    	  if (d < tip_radius)
			   return result;
    	  else {
			  // outside the cutter, return a distance to the outer "ring" of the cutter
    		   GLVertex n = GLVertex(t.x, t.y, 0.0);
    		   n = n * (tip_radius / d);  // 1/d means normalization
    		   result.f = -((t - n).norm());
    		   return result;
    	  }
   }
}

} // end namespace
// end of file volume.cpp
