/*
 *  Copyright 2010-2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
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

#ifndef VOLUME_H
#define VOLUME_H

#include <iostream>
#include <list>
#include <cassert>
#include <float.h>

#include "bbox.hpp"
#include "facet.hpp"
#include "glvertex.hpp"
#include "gldata.hpp"
#include "stl.hpp"

#ifdef MULTI_THREAD
#include <QtCore>
#include <QFuture>
#endif

#include <QObject>
#include <QtConcurrent>

#include <algorithm>

namespace cutsim {

/// base-class for defining implicit volumes from which to build octrees
/// an implicit volume is defined as a function dist(Point p)
/// which returns a positive value inside the volume and a negative value outside.
///
/// the "positive inside, negative outside" sign-convetion means that boolean operations can be done with:
///
///  A U B ('union' or 'sum') =  max( d(A),  d(B) )
///  A \ B ('diff' or 'minus') = min( d(A), -d(B) )
///  A int B ('intersection') =  min( d(A),  d(B) )
///
/// reference: Frisken et al. "Designing with Distance Fields"
/// http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.69.9025
///
/// iso-surface extraction using standard marching-cubes requires just the distance
/// field to be stored at each corner vertex of an octree leaf-node.
///
/// advanced iso-surface extraction using extended-marching-cubes or dual-contouring may require
/// more information such as normals of the distance field or exact
/// intersection points and normals. This is similar to a tri-dexel representation.
/// In multi-material simulation a material-index can be stored.
/// Each cutter may also cut the material with a color of its own (new vertices have the color of the cutter).

typedef enum {
       NO_VOLUME				= 0,
       RECTANGLE_VOLUME			= 1,
       CYLINDER_VOLUME			= 2,
       SPHERE_VOLUME			= 3,
       STL_VOLUME				= 4,
} VolumeType;

class Volume : public QObject {
Q_OBJECT
    public:
        /// default constructor
        Volume(){}
        virtual ~Volume(){}

        VolumeType type;

        /// return signed distance from volume surface to Point p
        /// Points p inside the volume should return positive values.
        /// Points p outside the volume should return negative values.
        virtual double dist(const GLVertex& p) const = 0;
/// update the bounding-box
virtual void calcBB() {}
        /// bounding-box. This holds the maximum(minimum) points along the X,Y, and Z-coordinates
        /// of the volume (i.e. the volume where dist(p) returns negative values)
        Bbox bb;
        /// the color of this Volume
        Color color;
        /// set the color
        void setColor(GLfloat r, GLfloat g, GLfloat b) {
            color.r = r; color.g = g; color.b = b;
        }


void setProgress(int value) { progress = value; };
void accumlateProgress(int value) { progress += value; };
void sendProgress() { emit signalProgressValue(progress); }
int Progress() { return progress; }
signals:
void signalProgressValue(int progress);
void signalProgressFeature(QString aMessage, int mim_value, int max_value);

private:
int progress;
};

/// sphere centered at center
class SphereVolume: public Volume {

    public:
        /// default constructor
        SphereVolume();
        virtual ~SphereVolume() {}
        /// set radius of sphere
        void setRadius(double r) {
            radius = r;
            calcBB();
        }
        /// set the centerpoint of the sphere
        void setCenter(GLVertex v) {
            center = v;
            calcBB();
        }
        /// update the Bbox
        void calcBB();
        double dist(const GLVertex& p) const;

        /// center Point of sphere
        GLVertex center;
        /// radius of sphere
        double radius;
};

/// box-volume
/// from corner, go out in the three directions v1,v2,v3
/// interior points = corner + a*v1 + b*v2 + c*v3
/// where a, b, c are in [0,1]
class RectVolume : public Volume {

    public:
        /// default constructor
        RectVolume();
        virtual ~RectVolume() {}
        /// one corner of the box
        GLVertex corner;
        /// first vector from corner, to span the box
        GLVertex v1;
        /// second vector
        GLVertex v2;
        /// third vector
        GLVertex v3;
        /// update the bounding-box
        void calcBB();
        double dist(const GLVertex& p) const;
};

/// box-volume2
/// corner is located at the left bottom one.
/// from corner, width(x) length(y) hight(z)
class RectVolume2 : public Volume {

    public:
        /// default constructor
        RectVolume2();
        virtual ~RectVolume2() {}
        /// set the corner of box
        void setCorner(GLVertex c) {
        	corner = c;
        	center = GLVertex(corner.x + width * 0.5, corner.y + length * 0.5, corner.z);
        	calcBB();
        }
        /// set the length of box
        void setLength(double l) {
            length = l;
            center = GLVertex(corner.x + width * 0.5, corner.y + length * 0.5, corner.z);
            calcBB();
        }
        /// set the width of box
        void setWidth(double w) {
            width = w;
            center = GLVertex(corner.x + width * 0.5, corner.y + length * 0.5, corner.z);
            calcBB();
        }
        /// set the hight of box
        void setHight(double h) {
            hight = h;
            calcBB();
        }
        /// set the center of box
        void setCenter(GLVertex c) {
            center = c;
            corner = GLVertex(center.x - width * 0.5, center.y - length * 0.5, center.z);
            calcBB();
        }
        /// set the rotation center of box
        void setRotationCenter(GLVertex c) {
            rotationCenter = c;
        }
        /// set the angle of box
        void setAngle(GLVertex a) {
            angle = a;
        }
        /// update the bounding-box
        void calcBB();
        double dist(const GLVertex& p) const;

    private:
        /// one corner of the left bottom of box
        GLVertex corner;
        /// width of the box
        double width;
        /// length of the box
        double length;
        /// hight of the box
        double hight;
        /// box center of the bottom
        GLVertex center;
        /// center of rotation
        GLVertex rotationCenter;
        /// box angle
        GLVertex angle;
};

/// cylinder volume
class CylinderVolume: public Volume {

    public:
        CylinderVolume();
        virtual ~CylinderVolume() {}
        /// set radius of cylinder
        void setRadius(double r) {
            radius = r;
            calcBB();
		}
        /// set the length of cylinder
        void setLength(double l) {
            length = l;
            calcBB();
        }
        /// set the center of cylinder
        void setCenter(GLVertex v) {
            center = v;
            calcBB();
        }
        /// set the rotation center of clynder
        void setRotationCenter(GLVertex c) {
            rotationCenter = c;
        }
        /// set the angle of cylinder
        void setAngle(GLVertex a) {
            angle = a;
        }
        /// update the bounding-box of cylinder
        void calcBB();
        double dist(const GLVertex& p) const;

    private:
        /// cylinder radius
        double radius;
        /// cylinder length
        double length;
        /// cylinder center
        GLVertex center;
        /// center of rotation
        GLVertex rotationCenter;
        /// cylinder angle
        GLVertex angle;
};

/// STL volume
class StlVolume: public Stl, public Volume {

	public:
        StlVolume();
        virtual ~StlVolume() {
        	V21.resize(0); V21invV21dotV21.resize(0);
        	V32.resize(0); V32invV32dotV32.resize(0);
        	V13.resize(0); V13invV13dotV13.resize(0);
        }

        /// set the center of STL
        void setCenter(GLVertex v) {
            center = v;
        }
        /// set the rotation center of STL
        void setRotationCenter(GLVertex c) {
            rotationCenter = c;
        }
        /// set the angle of STL
        void setAngle(GLVertex a) {
            angle = a;
        }
        /// set the resplution
        void setCubeResolution(double resolution) {
        	cube_resolution = resolution;
        }
        /// update the bounding-box of STL
        void calcBB();
        double dist(const GLVertex& p) const;

        int readStlFile(QString file) {  int retval = Stl::readStlFile(file); /*calcBB();*/ return retval; }

   private:
        // V21[i] = facets[i]->v2 - facets[i]->v1
        std::vector<GLVertex> V21;
        // V21[i]/<V21[i], V21[i]>
        std::vector<GLVertex> V21invV21dotV21;
        // V32[i] = facets[i]->v3 - facets[i]->v2
        std::vector<GLVertex> V32;
        // V32[i]/<V32[i], V32[i]>
        std::vector<GLVertex> V32invV32dotV32;
        // V13[i] = facets[i]->v1 - facets[i]->v3
        std::vector<GLVertex> V13;
        // V13[i]/<V13[i], V13[i]>
        std::vector<GLVertex> V13invV13dotV13;
        /// STL center
        GLVertex center;
        /// center of rotation
        GLVertex rotationCenter;
        /// STL angle
        GLVertex angle;

enum  { MAX_INDEX = 128, MAX_X = MAX_INDEX, MAX_Y = MAX_INDEX, MAX_Z = MAX_INDEX, ROUGH_RATIO = 16 };

        GLVertex maxpt;
        GLVertex minpt;
        double maxlength;
        int src_index = 0, dst_index = 1;
        void swapIndex() { int tmp = src_index; src_index = dst_index; dst_index = tmp; assert(src_index^dst_index); }
        std::vector<int>  neighborhoodIndex[2][MAX_X][MAX_Y][MAX_Z];
        int  ratio = 2;
        double indexcubesize;
        double invcubesize;
        void calcNeighborhoodIndex(int index_x, int index_y, int index_z, double upper_limit, bool from_facets);
        double distance(const GLVertex& p, int index);
        double cube_resolution;
};


typedef enum {
       NO_TOOL				= 0,
       CYLINDER				= 1,
       BALL					= 2,
       BULL					= 3,
       CONE					= 4,
       DRILL				= 10,
} CutterType;

typedef enum {
		NO_COLLISION 		= 0x0,
		NECK_COLLISION 		= 0x10000,
		SHANK_COLLISION 	= 0x20000,
		HOLDER_COLLISION 	= 0x40000,
		PARTS_COLLISION		= 0x80000,
} CollisionType;

typedef struct cutting {
	double	f;
	int		collision;
	int		count;
} Cutting;

/// cylindrical cutter volume

class CutterVolume: public Volume {

	public:
        CutterVolume();
        /// cutter type
        CutterType cuttertype;
        /// cutter radius
        double radius;
        /// cutter length
        double length;
        /// flute length
        double flutelength;
        /// neck radius
        double neckradius;
        /// reach length
        double reachlength;
        /// shank radius
        double shankradius;
        /// max radius of cutter, neck and shank radius
        double maxradius;
        /// cutter center
        GLVertex center;
        /// cutter angle
        GLVertex angle;

        /// Holder variables
        bool enableholder;
        double holderradius;
        double holderlength;
        void enableHolder(bool flag) {
            enableholder = flag;
            if (enableholder && holderradius <= 0.0)
            	holderradius = DEFAULT_HOLDER_RADIUS;
            if (enableholder && holderlength <= 0.0)
            	holderlength = DEFAULT_HOLDER_LENGTH;
        }
        void setHolderRadius(double r) {
            holderradius = r;
            if (holderradius > 0.0)
            	enableholder = true;
            else
            	enableholder = false;
            calcBBHolder();
        }
        void setHolderLength(double l) {
            holderlength = l;
            calcBBHolder();
        }

        Bbox bbHolder;
        /// update the Bbox
        void calcBBHolder();

        virtual void setRadius(double r) {}
        virtual void setAngle(GLVertex a) {}
        virtual void setCenter(GLVertex v) {}
        virtual void setTool(CutterVolume* cutter) {}
        virtual GLVertex getCenter() { return center; }
        virtual GLVertex getAngle() { return angle; }
        virtual double dist(const GLVertex& p) const { return 0.0; }
        virtual Cutting dist_cd(const GLVertex& p) const { Cutting r = { 0.0, NO_COLLISION, 0 }; return r; }
};

/// cylindrical cutter volume
class CylCutterVolume: public CutterVolume {

	public:
        CylCutterVolume();
        /// set the radius of Cylindrical Cutter
        void setRadius(double r) {
            radius = r;
            neckradius = r;
            shankradius = r;
            maxradius = r;
            effective_radius = 0.25 * radius;
            calcBB();
        }
        /// set the length of Cylindrical Cutter
        void setLength(double l) {
            length = l;
            flutelength = l;
            reachlength = l;
            calcBB();
        }
        /// set the centerpoint of Cylindrical Cutter
        void setCenter(GLVertex v) {
            center = v;
            calcBB();
        }        /// cutter position
        /// set the angle of Cylindrical Cutter
        void setAngle(GLVertex a) {
            angle = a;
            bb.angle = a;
            if (enableholder)
            	bbHolder.angle = a;
         }
        /// set the flute length of Cylindrical Cutter
        void setFluteLength(double fl) {
        	flutelength = fl;
        }
        /// set the neck radius of Cylindrical Cutter
        void setNeckRadius(double nr) {
            neckradius = nr;
        }
        /// set the reach length of Cylindrical Cutter
        void setReachLength(double rl) {
        	reachlength = rl;
        }
        /// set the shank radius of Cylindrical Cutter
        void setShankRadius(double sr) {
            shankradius = sr;
            if (shankradius > radius) {
            	maxradius = sr;
            	calcBB();
            }
        }
        /// get the centerpoint of Cylindrical Cutter
        GLVertex getCenter() { return center; }
        /// get the angle of Cylindrical Cutter
        GLVertex getAngle()  { return angle; }
        /// update the Bbox
        void calcBB();
        double dist(const GLVertex& p) const;
        Cutting dist_cd(const GLVertex& p) const;

    private:
        double effective_radius;
};

/// ball-nose cutter volume

class BallCutterVolume: public CutterVolume {

    public:
        BallCutterVolume();
        /// set the radius of Ball Cutter
        void setRadius(double r) {
            radius = r;
            neckradius = r;
            shankradius = r;
            maxradius = r;
            calcBB();
        }
        /// set the length of Ball Cutter
        void setLength(double l) {
            length = l - radius;
            flutelength = length;
            reachlength = length;
            calcBB();
        }
        /// set the centerpoint of Ball Cutter
        void setCenter(GLVertex v) {
            center = v + GLVertex(0.0, 0.0, radius);
            calcBB();
        }
        /// set the angle of Ball Cutter
        void setAngle(GLVertex a) {
            angle = a;
            bb.angle = a;
            if (enableholder)
            	bbHolder.angle = a;
        }
        /// set the flute length of Ball Cutter
        void setFluteLength(double fl) {
        	flutelength = fl - radius;
        }
        /// set the neck radius of Ball Cutter
        void setNeckRadius(double nr) {
            neckradius = nr;
        }
        /// set the reach length of Ball Cutter
        void setReachLength(double rl) {
        	reachlength = rl - radius;
        }
        /// set the shank radius of Ball Cutter
        void setShankRadius(double sr) {
            shankradius = sr;
            if (shankradius > radius) {
            	maxradius = sr;
            	calcBB();
            }
        }
        /// get the centerpoint of Ball Cutter
        GLVertex getCenter() { return center; }
        /// get the angle of Ball Cutter
        GLVertex getAngle()  { return angle; }
        /// update bounding box
        void calcBB();
        double dist(const GLVertex& p) const;
        Cutting dist_cd(const GLVertex& p) const;

    private:
        double effective_hight;
};

/// bull-nose cutter volume

class BullCutterVolume: public CutterVolume {

    public:
        BullCutterVolume();

        /// radius of cylinder-part
        double r1;
        /// radius of torus
        double r2;

        /// set the radius of Bull Cutter
        void setRadius(double r) {
            radius = r;
            neckradius = r;
            shankradius = r;
            maxradius = r;
            calcBB();
        }
        /// set the nose radius of Bull Cutter
        void setBullRadius(double r) {
        	if (radius - r > 0.0) {
        		r1 = radius - r;
        		r2 = r;
            }
        }
        /// set the length of Bull Cutter
        void setLength(double l) {
            length = l;
            flutelength = l;
            reachlength = l;
            calcBB();
        }
        /// set the centerpoint of Bull Cutter
        void setCenter(GLVertex v) {
            center = v;
            calcBB();
        }        /// cutter position
        /// set the angle of Bull Cutter
        void setAngle(GLVertex a) {
            angle = a;
            bb.angle = a;
            if (enableholder)
            	bbHolder.angle = a;
         }
        /// set the flute length of Bull Cutter
        void setFluteLength(double fl) {
        	flutelength = fl;
        }
        /// set the neck radius of Bull Cutter
        void setNeckRadius(double nr) {
            neckradius = nr;
        }
        /// set the reach length of Bull Cutter
        void setReachLength(double rl) {
        	reachlength = rl;
        }
        /// set the shank radius of Bull Cutter
        void setShankRadius(double sr) {
            shankradius = sr;
            if (shankradius > radius) {
            	maxradius = sr;
            	calcBB();
            }
        }
        /// get the centerpoint of Bull Cutter
        GLVertex getCenter() { return center; }
        /// get the angle of Bull Cutter
        GLVertex getAngle()  { return angle; }
        /// update bounding box
        void calcBB();
        double dist(const GLVertex& p) const;
        Cutting dist_cd(const GLVertex& p) const;

    private:
        double effective_radius;
};

/// cone-nose cutter volume

class ConeCutterVolume: public CutterVolume {

    public:
		ConeCutterVolume();

        /// radius of tip-part
        double r1;
        /// radius of base-part
        double r2;

        /// set the radius of Cone Cutter
        void setRadius(double r) {
            radius = r;
            neckradius = r;
            shankradius = r;
            maxradius = r;
            r1 = r2 = r;
            calcBB();
        }
        /// set the nose tip radius of Cone Cutter
        void setTipRadius(double r) {
        	r1 = r;
        	if (maxradius < r1) {
                maxradius = r1;
                calcBB();
        	}
        	if (flutelength > 0.0)
        		incline_coff = (r2 - r1) / flutelength;
        }
        /// set the nose base radius of Cone Cutter
        void setBaseRadius(double r) {
        	r2 = r;
        	if (maxradius < r2) {
                maxradius = r2;
                calcBB();
        	}
        	if (flutelength > 0.0)
        		incline_coff = (r2 - r1) / flutelength;
        }
        /// set the length of Cone Cutter
        void setLength(double l) {
            length = l;
            flutelength = l;
            reachlength = l;
            calcBB();
        	if (flutelength > 0.0)
        		incline_coff = (r2 - r1) / flutelength;
        }
        /// set the centerpoint of Cone Cutter
        void setCenter(GLVertex v) {
            center = v;
            calcBB();
        }        /// cutter position
        /// set the angle of Cone Cutter
        void setAngle(GLVertex a) {
            angle = a;
            bb.angle = a;
            if (enableholder)
            	bbHolder.angle = a;
         }
        /// set the flute length of Cone Cutter
        void setFluteLength(double fl) {
        	flutelength = fl;
        	if (flutelength > 0.0)
        		incline_coff = (r2 - r1) / flutelength;
        }
        /// set the neck radius of Cone Cutter
        void setNeckRadius(double nr) {
            neckradius = nr;
        }
        /// set the reach length of Cone Cutter
        void setReachLength(double rl) {
        	reachlength = rl;
        }
        /// set the shank radius of Cone Cutter
        void setShankRadius(double sr) {
            shankradius = sr;
            if (shankradius > radius) {
            	maxradius = sr;
            	calcBB();
            }
        }
        /// get the centerpoint of Cone Cutter
        GLVertex getCenter() { return center; }
        /// get the angle of Cone Cutter
        GLVertex getAngle()  { return angle; }
        /// update bounding box
        void calcBB();
        double dist(const GLVertex& p) const;
        Cutting dist_cd(const GLVertex& p) const;

    private:
        double incline_coff;
        double effective_radius;
};

/// drill volume

class DrillVolume: public CutterVolume {

    public:
		DrillVolume();

        /// hight of tip-part
		double tip_hight;

		/// set the radius of Drill
        void setRadius(double r) {
            radius = r;
            neckradius = r;
            shankradius = r;
            maxradius = r;
            assert(tip_angle > 0.0);
            tip_hight = radius / tan(tip_angle/(2.0 * 180.0) * PI);
            if (tip_hight > 0.0)
            	incline_coff = radius / tip_hight;
            tip_radius = 0.1 * radius;
            calcBB();
        }
        /// set the nose tip angle of Drill
        void setTipAngle(double a) {
        	if (a > 0.0) {
        		tip_angle = a;
        		tip_hight = radius / tan(tip_angle/(2.0 * 180.0) * PI);
                if (tip_hight > 0.0)
        			incline_coff = radius / tip_hight;
        	}
        }
        /// set the length of Drill
        void setLength(double l) {
            length = l;
            flutelength = l;
            reachlength = l;
            calcBB();
        }
        /// set the centerpoint of Drill
        void setCenter(GLVertex v) {
            center = v;
            calcBB();
        }        /// cutter position
        /// set the angle of Drill
        void setAngle(GLVertex a) {
            angle = a;
            bb.angle = a;
            if (enableholder)
            	bbHolder.angle = a;
         }
        /// set the flute length of Drill
        void setFluteLength(double fl) {
        	flutelength = fl;
        }
        /// set the neck radius of Drill
        void setNeckRadius(double nr) {
            neckradius = nr;
        }
        /// set the reach length of Drill
        void setReachLength(double rl) {
        	reachlength = rl;
        }
        /// set the shank radius of Drill
        void setShankRadius(double sr) {
            shankradius = sr;
            if (shankradius > radius) {
            	maxradius = sr;
            	calcBB();
            }
        }
        /// get the centerpoint of Drill
        GLVertex getCenter() { return center; }
        /// get the angle of Drill
        GLVertex getAngle()  { return angle; }
        /// update bounding box
        void calcBB();
        double dist(const GLVertex& p) const;
        Cutting dist_cd(const GLVertex& p) const;

    private:
        double tip_angle;
        double tip_radius;
        double incline_coff;
        double effective_hight;
};

} // end namespace
#endif
// end file volume.hpp
