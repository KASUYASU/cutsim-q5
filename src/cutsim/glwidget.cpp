/*
 *  Copyright 2010-2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
 *  Copyright 2015      Kazuyasu Hamada (k-hamada@gifu-u.ac.jp)
 *
 *  This file is part of Cutsim / OpenCAMlib.
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

#include <iostream>

#include <QObject>
#include <QTimer>

#include "glwidget.hpp"

namespace cutsim {

GLWidget::GLWidget( unsigned int sceneRadius, QWidget *parent, char *name ) {
    setSceneRadius(sceneRadius);
    showEntireScene();
    file_number = 0;
    corner_axis = true;
    draw_tool = true;
    enable_animate = true;
    tool.x = tool.y = 0.0;
    tool.z = DEFAULT_MAX_Z_LIMIT;
    spindleradius = DEFAULT_SPINDLE_RADIUS;
    spindlelength = DEFAULT_SPINDLE_LENGTH;
}

/// add new GLData object and return pointer to it.
GLData* GLWidget::addGLData() {
    GLData* g = new GLData();
    glObjects.push_back(g);
    return g;
}

/// loop through glObjects and for each GLData draw it using VBO
void GLWidget::draw()  {

    glPushMatrix();
#ifdef MULTI_AXIS
    glRotatef(tool.a*180.0/PI, 1.0f, 0.0f, 0.0f);
    glRotatef(tool.c*180.0/PI, 0.0f, 0.0f, 1.0f);
#endif

    BOOST_FOREACH( GLData* g, glObjects ) { // draw each object

        // apply a transformation-matrix here !?
        QMutexLocker locker( &(g->renderMutex) );
        glPolygonMode( g->polygonFaceMode(), g->polygonFillMode() );
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_COLOR_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        // http://www.opengl.org/sdk/docs/man/xhtml/glNormalPointer.xml
        glNormalPointer( GLData::coordinate_type, sizeof( GLData::vertex_type ),    ((GLbyte*)g->getVertexArray()  + GLData::normal_offset ) );
        // http://www.opengl.org/sdk/docs/man/xhtml/glColorPointer.xml
        glColorPointer(  3, GLData::color_type     , sizeof( GLData::vertex_type ), ((GLbyte*)g->getVertexArray()  + GLData::color_offset  ) );
        glVertexPointer( 3, GLData::coordinate_type, sizeof( GLData::vertex_type ), ((GLbyte*)g->getVertexArray()  + GLData::vertex_offset  ) );
        // http://www.opengl.org/sdk/docs/man/xhtml/glDrawElements.xml
        glDrawElements( g->GLType() , g->indexCount() , GLData::index_type, g->getIndexArray());
        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
    }
    glPopMatrix();

if (selectedName() >= 0)
{
	// Draw the intersection line only.
	glColor3f(0.0, 0.0, 1.0); // blue
	glBegin(GL_LINES);
	glVertex3fv(orig);
	glVertex3fv(selectedPoint - dir);
	glEnd();
	// Draw (approximated) intersection point on selected object
	glColor3f(0.9f, 0.2f, 0.1f);
	glPointSize(5.0);
	glBegin(GL_POINTS);
	glVertex3fv(selectedPoint);
	glEnd();
} else {
	// Draw the intersection line only.
	glColor3f(0.0, 0.0, 1.0); // blue
	glBegin(GL_LINES);
	glVertex3fv(orig);
	glVertex3fv(orig + 100.0*dir);
	glEnd();
}

    lastFrameTime = QTime::currentTime();
}

void GLWidget::postDraw() {
    QGLViewer::postDraw();
    if (draw_tool)
#ifdef MULTI_AXIS
    drawTool(tool.x, tool.y, tool.z, tool.a, tool.b, tool.c, tool.cutter);
#else
    drawTool(tool.x, tool.y, tool.z, tool.cutter);
#endif
    if (corner_axis)
        drawCornerAxis();
}

void GLWidget::slotWriteScreenshot() {
//    QImage img = grabFrameBuffer();
//    QString file_name = "frame_" + QString::number(file_number) + ".png";
//    img.save( file_name );
//    file_number++;
}

void GLWidget::keyPressEvent(QKeyEvent *e) {
    const Qt::KeyboardModifiers modifiers = e->modifiers();
    // A simple switch on e->key() is not sufficient if we want to take state key into account.
    // With a switch, it would have been impossible to separate 'F' from 'CTRL+F'.
    // That's why we use imbricated if...else and a "handled" boolean.
    bool handled = false;

    if ((e->key() == Qt::Key_A) && (modifiers == Qt::NoButton)) {
    	 enable_animate = !enable_animate;
    	 if (enable_animate == true)
    		 statusBarMessage( tr("Animation Enable") );
    	 else
    		 statusBarMessage( tr("Animation Disable") );

         update();
         handled=true;
    }

    if ((e->key() == Qt::Key_C) && (modifiers == Qt::NoButton)) {
    	corner_axis = !corner_axis;
        update();
    	handled = true;
    }

    if ((e->key() == Qt::Key_T) && (modifiers == Qt::NoButton)) {
    	 draw_tool = !draw_tool;
    	 if (draw_tool == true)
    		 statusBarMessage( tr("Draw Tool Enable") );
    	 else
    		 statusBarMessage( tr("Draw Tool Disable") );

         update();
    	 handled=true;
    }

    if ((e->key() == Qt::Key_S) && (modifiers == Qt::NoButton)) {
    	std::cout << tree->str();
    	handled = true;
    }

    if ((e->key() == Qt::Key_D) && (modifiers == Qt::NoButton)) {
    	tree->g->print();
    	handled = true;
    }

    if ((e->key() == Qt::Key_R) && (modifiers == Qt::NoButton)) {
    	tree->setInvalid();
    	requestRedraw(tr("Refresh Screen"));
    	handled = true;
    }

if ((e->key() == Qt::Key_F) && (modifiers == Qt::NoButton)) {
//tree->treeTransfer(GLVertex(0.0, 0.0, -6.0), X_AXIS);
tree->treeTransfer(GLVertex(0.0, 0.0, -6.0), X_AXIS, true);
requestRedraw(tr("Flip Stock!!"));
handled = true;
}

    if (!handled)
        QGLViewer::keyPressEvent(e);
}

void GLWidget::drawWithNames()
{
    glPushMatrix();
#ifdef MULTI_AXIS
    glRotatef(tool.a*180.0/PI, 1.0f, 0.0f, 0.0f);
    glRotatef(tool.c*180.0/PI, 0.0f, 0.0f, 1.0f);
#endif

    int i = 0;

    BOOST_FOREACH( GLData* g, glObjects ) { // draw each object
        glPushName(i++);
        // apply a transformation-matrix here !?
        QMutexLocker locker( &(g->renderMutex) );
        glPolygonMode( g->polygonFaceMode(), g->polygonFillMode() );
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_COLOR_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        // http://www.opengl.org/sdk/docs/man/xhtml/glNormalPointer.xml
        glNormalPointer( GLData::coordinate_type, sizeof( GLData::vertex_type ),    ((GLbyte*)g->getVertexArray()  + GLData::normal_offset ) );
        // http://www.opengl.org/sdk/docs/man/xhtml/glColorPointer.xml
        glColorPointer(  3, GLData::color_type     , sizeof( GLData::vertex_type ), ((GLbyte*)g->getVertexArray()  + GLData::color_offset  ) );
        glVertexPointer( 3, GLData::coordinate_type, sizeof( GLData::vertex_type ), ((GLbyte*)g->getVertexArray()  + GLData::vertex_offset  ) );
        // http://www.opengl.org/sdk/docs/man/xhtml/glDrawElements.xml
        glDrawElements( g->GLType() , g->indexCount() , GLData::index_type, g->getIndexArray());
        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
        glPopName();
    }
    glPopMatrix();
}

#include <qmessagebox.h>
void GLWidget::postSelection(const QPoint& point)
{
  // Compute orig and dir, used to draw a representation of the intersecting line
  camera()->convertClickToLine(point, orig, dir);

  // Find the selectedPoint coordinates, using camera()->pointUnderPixel().
  bool found;
  selectedPoint = camera()->pointUnderPixel(point, found);
  selectedPoint -= 0.01f*dir; // Small offset to make point clearly visible.
  // Note that "found" is different from (selectedObjectId()>=0) because of the size of the select region.

  if (selectedName() == -1)
    QMessageBox::information(this, "No selection",
			     "No object selected under pixel " + QString::number(point.x()) + "," + QString::number(point.y()));
  else
    QMessageBox::information(this, "Selection",
			     "Spiral number " + QString::number(selectedName()) + " selected under pixel " +
//			     QString::number(point.x()) + "," + QString::number(point.y()));
			     QString::number(selectedPoint.x) + "," + QString::number(selectedPoint.y)+ "," + QString::number(selectedPoint.z));
}

void GLWidget::drawCornerAxis() {
    int viewport[4];
    int scissor[4];
    // The viewport and the scissor are changed to fit the lower left
    // corner. Original values are saved.
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetIntegerv(GL_SCISSOR_BOX, scissor);
    // Axis viewport size, in pixels
    const int size = 150;
    glViewport(0,0,size,size);
    glScissor(0,0,size,size);
    // The Z-buffer is cleared to make the axis appear over the
    // original image.
    glClear(GL_DEPTH_BUFFER_BIT);
    // Tune for best line rendering
    glDisable(GL_LIGHTING);
    glLineWidth(3.0);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(-1, 1, -1, 1, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMultMatrixd( camera()->orientation().inverse().matrix() );
#ifdef GLVIEWER_AXIS
    QGLViewer::drawAxis(1.0);
#else
    glBegin(GL_LINES);
        glColor3f(1.0, 0.0, 0.0); // red X-axis
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(1.0, 0.0, 0.0);
        glColor3f(0.0, 1.0, 0.0); // green Y-axis
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 1.0, 0.0);
        glColor3f(0.0, 0.0, 1.0); // blue Z-axis
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, 1.0);
    glEnd();
#endif
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glEnable(GL_LIGHTING);
    // The viewport and the scissor are restored.
    glScissor(scissor[0],scissor[1],scissor[2],scissor[3]);
    glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
}

#ifdef MULTI_AXIS
void GLWidget::drawTool(double x, double y, double z, double a, double b, double c, CutterVolume* cutter)
#else
void GLWidget::drawTool(double x, double y, double z, CutterVolume* cutter)
#endif
{
	if (cutter == NULL) return;

	GLUquadricObj *quad = gluNewQuadric();

	glColor3f(TOOL_BODY_COLOR);

	glPushMatrix();
	switch (cutter->cuttertype) {
	case CYLINDER:
		// draw side wall
		glTranslated(x, y, z);
		if (cutter->flutelength < cutter->length) {
			glPushMatrix();
		    glColor3f(TOOL_FLUTE_COLOR);
			gluCylinder(quad, cutter->radius, cutter->radius, cutter->flutelength, 30, 1);
			glTranslated(0.0, 0.0, cutter->flutelength);
			gluDisk(quad, 0.0, cutter->radius, 30, 1);
		    glColor3f(TOOL_BODY_COLOR);
			if (cutter->reachlength < cutter->length) {
				gluCylinder(quad, cutter->neckradius, cutter->neckradius, cutter->reachlength - cutter->flutelength, 30, 1);
				glTranslated(0.0, 0.0, cutter->reachlength - cutter->flutelength);
				if (cutter->length - cutter->reachlength > 4.0) {
					gluCylinder(quad, cutter->neckradius, cutter->shankradius, 4.0, 30, 1);
					glTranslated(0.0, 0.0, 4.0);
					gluCylinder(quad, cutter->shankradius, cutter->shankradius, cutter->length - cutter->reachlength - 4.0, 30, 1);
				} else
					gluCylinder(quad, cutter->neckradius, cutter->shankradius, cutter->length - cutter->reachlength, 30, 1);
			} else
				gluCylinder(quad, cutter->neckradius, cutter->neckradius, cutter->length - cutter->flutelength, 30, 1);
			glPopMatrix();
		} else
			gluCylinder(quad, cutter->radius, cutter->radius, cutter->length, 30, 1);

		// draw bottom wall
		gluQuadricOrientation(quad, GLU_INSIDE);
		gluDisk(quad, 0.0, cutter->radius, 30, 1);

		// draw upper wall
		gluQuadricOrientation(quad, GLU_OUTSIDE);
		glTranslated(0.0, 0.0, cutter->length);
		gluDisk(quad, 0.0, cutter->shankradius, 30, 1);

		break;

	case BALL:
	// draw side wall
		glTranslated(x, y, z + cutter->radius);
		if (cutter->flutelength < cutter->length) {
			glPushMatrix();
		    glColor3f(TOOL_FLUTE_COLOR);
			gluCylinder(quad, cutter->radius, cutter->radius, cutter->flutelength, 30, 1);
			glTranslated(0.0, 0.0, cutter->flutelength);
			gluDisk(quad, 0.0, cutter->radius, 30, 1);
		    glColor3f(TOOL_BODY_COLOR);
			if (cutter->reachlength < cutter->length) {
				gluCylinder(quad, cutter->neckradius, cutter->neckradius, cutter->reachlength - cutter->flutelength, 30, 1);
				glTranslated(0.0, 0.0, cutter->reachlength - cutter->flutelength);
				if (cutter->length - cutter->reachlength > 4.0) {
					gluCylinder(quad, cutter->neckradius, cutter->shankradius, 4.0, 30, 1);
					glTranslated(0.0, 0.0, 4.0);
					gluCylinder(quad, cutter->shankradius, cutter->shankradius, cutter->length - cutter->reachlength - 4.0, 30, 1);
				} else
					gluCylinder(quad, cutter->neckradius, cutter->shankradius, cutter->length - cutter->reachlength, 30, 1);
			} else
				gluCylinder(quad, cutter->neckradius, cutter->neckradius, cutter->length - cutter->flutelength, 30, 1);
			glPopMatrix();
		}else
			gluCylinder(quad, cutter->radius, cutter->radius, cutter->length, 30, 1);

		// draw bottom wall
		glColor3f(TOOL_FLUTE_COLOR);
		gluSphere(quad, cutter->radius, 30, 30);

		// draw upper wall
		glColor3f(TOOL_BODY_COLOR);
		glTranslated(0.0, 0.0, cutter->length);
		gluDisk(quad, 0.0, cutter->shankradius, 30, 1);

		break;

	case BULL:
		// draw side wall
		glTranslated(x, y, z);
		if (cutter->flutelength < cutter->length) {
			glPushMatrix();
		    glColor3f(TOOL_FLUTE_COLOR);
			glTranslated(0.0, 0.0, ((BullCutterVolume*)cutter)->r2);
		    gluCylinder(quad, cutter->radius, cutter->radius, cutter->flutelength - ((BullCutterVolume*)cutter)->r2, 30, 1);
			glTranslated(0.0, 0.0, cutter->flutelength - ((BullCutterVolume*)cutter)->r2);
			gluDisk(quad, 0.0, cutter->radius, 30, 1);
		    glColor3f(TOOL_BODY_COLOR);
			if (cutter->reachlength < cutter->length) {
				gluCylinder(quad, cutter->neckradius, cutter->neckradius, cutter->reachlength - cutter->flutelength, 30, 1);
				glTranslated(0.0, 0.0, cutter->reachlength - cutter->flutelength);
				if (cutter->length - cutter->reachlength > 4.0) {
					gluCylinder(quad, cutter->neckradius, cutter->shankradius, 4.0, 30, 1);
					glTranslated(0.0, 0.0, 4.0);
					gluCylinder(quad, cutter->shankradius, cutter->shankradius, cutter->length - cutter->reachlength - 4.0, 30, 1);
				} else
					gluCylinder(quad, cutter->neckradius, cutter->shankradius, cutter->length - cutter->reachlength, 30, 1);
			} else
				gluCylinder(quad, cutter->neckradius, cutter->neckradius, cutter->length - cutter->flutelength, 30, 1);
			glPopMatrix();
		} else
			gluCylinder(quad, cutter->radius, cutter->radius, cutter->length, 30, 1);

		// draw bottom wall
		glColor3f(TOOL_FLUTE_COLOR);
		gluQuadricOrientation(quad, GLU_INSIDE);
		gluDisk(quad, 0.0, ((BullCutterVolume*)cutter)->r1, 30, 1);

		// draw quarter torus
		{
			int numc = 5;
			int numt = 20;
			double s, t, x, y, z, twopi = 2 * PI, halfdpi = PI / 2 + 0.5;
			double r1 = ((BullCutterVolume*)cutter)->r1;
			double r2 = ((BullCutterVolume*)cutter)->r2;

			for (int i = 0; i < numc; i++) {
				glBegin(GL_QUAD_STRIP);
				for (int j = 0; j <= numt; j++) {
					for (int k = 1; k >= 0; k--) {
						s = (i + k) % numc;
						t = j % numt;

						x = (r1 + r2 * cos(s * halfdpi/numc)) * cos(t * twopi/numt);
						y = (r1 + r2 * cos(s * halfdpi/numc)) * sin(t * twopi/numt);
						z = r2 - r2 * sin(s * halfdpi/numc);
						glVertex3f(x, y, z);
					}
				}
				glEnd();
			}
		}

		// draw upper wall
		glColor3f(TOOL_BODY_COLOR);
		gluQuadricOrientation(quad, GLU_OUTSIDE);
		glTranslated(0.0, 0.0, cutter->length);
		gluDisk(quad, 0.0, cutter->shankradius, 30, 1);

		break;

	case CONE:
		// draw side wall
		glTranslated(x, y, z);
		if (cutter->flutelength < cutter->length) {
			glPushMatrix();
		    glColor3f(TOOL_FLUTE_COLOR);
			gluCylinder(quad, ((ConeCutterVolume*)cutter)->r1, ((ConeCutterVolume*)cutter)->r2, cutter->flutelength, 30, 1);
			glTranslated(0.0, 0.0, cutter->flutelength);
			gluDisk(quad, 0.0, ((ConeCutterVolume*)cutter)->r2, 30, 1);
		    glColor3f(TOOL_BODY_COLOR);
			if (cutter->reachlength < cutter->length) {
				gluCylinder(quad, cutter->neckradius, cutter->neckradius, cutter->reachlength - cutter->flutelength, 30, 1);
				glTranslated(0.0, 0.0, cutter->reachlength - cutter->flutelength);
				if (cutter->length - cutter->reachlength > 4.0) {
					gluCylinder(quad, cutter->neckradius, cutter->shankradius, 4.0, 30, 1);
					glTranslated(0.0, 0.0, 4.0);
					gluCylinder(quad, cutter->shankradius, cutter->shankradius, cutter->length - cutter->reachlength - 4.0, 30, 1);
				} else
					gluCylinder(quad, cutter->neckradius, cutter->shankradius, cutter->length - cutter->reachlength, 30, 1);
			} else
				gluCylinder(quad, cutter->neckradius, cutter->neckradius, cutter->length - cutter->flutelength, 30, 1);
			glPopMatrix();
		} else
			gluCylinder(quad, cutter->radius, cutter->radius, cutter->length, 30, 1);

		// draw bottom wall
		gluQuadricOrientation(quad, GLU_INSIDE);
		gluDisk(quad, 0.0, ((ConeCutterVolume*)cutter)->r1, 30, 1);

		// draw upper wall
		gluQuadricOrientation(quad, GLU_OUTSIDE);
		glTranslated(0.0, 0.0, cutter->length);
		gluDisk(quad, 0.0, cutter->shankradius, 30, 1);

		break;

	case DRILL:
		// draw side wall
		glTranslated(x, y, z);
		if (cutter->flutelength < cutter->length) {
			glPushMatrix();
		    glColor3f(TOOL_FLUTE_COLOR);
			gluCylinder(quad, 0.0, cutter->radius, ((DrillVolume*)cutter)->tip_hight, 30, 1);
			glTranslated(0.0, 0.0, ((DrillVolume*)cutter)->tip_hight);
			gluCylinder(quad, cutter->radius, cutter->radius, cutter->flutelength - ((DrillVolume*)cutter)->tip_hight, 30, 1);
			glTranslated(0.0, 0.0, cutter->flutelength - ((DrillVolume*)cutter)->tip_hight);
		    glColor3f(TOOL_BODY_COLOR);
			if (cutter->reachlength < cutter->length) {
				gluCylinder(quad, cutter->neckradius, cutter->neckradius, cutter->reachlength - cutter->flutelength, 30, 1);
				glTranslated(0.0, 0.0, cutter->reachlength - cutter->flutelength);
				if (cutter->length - cutter->reachlength > 4.0) {
					gluCylinder(quad, cutter->neckradius, cutter->shankradius, 4.0, 30, 1);
					glTranslated(0.0, 0.0, 4.0);
					gluCylinder(quad, cutter->shankradius, cutter->shankradius, cutter->length - cutter->reachlength - 4.0, 30, 1);
				} else
					gluCylinder(quad, cutter->neckradius, cutter->shankradius, cutter->length - cutter->reachlength, 30, 1);
			} else
				gluCylinder(quad, cutter->neckradius, cutter->neckradius, cutter->length - cutter->flutelength, 30, 1);
			glPopMatrix();
		} else
			gluCylinder(quad, cutter->radius, cutter->radius, cutter->length, 30, 1);

		// draw upper wall
		gluQuadricOrientation(quad, GLU_OUTSIDE);
		glTranslated(0.0, 0.0, cutter->length);
		gluDisk(quad, 0.0, cutter->shankradius, 30, 1);

		break;

	default:
		glTranslated(x, y, z);
	}

	if (cutter->enableholder) {
		glColor3f(DEFAULT_HOLDER_COLOR);
		gluQuadricOrientation(quad, GLU_INSIDE);
		gluDisk(quad, 0.0, cutter->holderradius, 30, 1);
		gluQuadricOrientation(quad, GLU_OUTSIDE);
		gluCylinder(quad, cutter->holderradius, cutter->holderradius, cutter->holderlength, 30, 1);
		glTranslated(0.0, 0.0, cutter->holderlength);
		gluDisk(quad, 0.0, cutter->holderradius, 30, 1);

		glColor3f(DEFAULT_SPINDLE_COLOR);
		gluCylinder(quad, spindleradius, spindleradius, spindlelength, 30, 1);
	}

	glPopMatrix();

	gluDeleteQuadric(quad);
}

} // end cutsim namespace
