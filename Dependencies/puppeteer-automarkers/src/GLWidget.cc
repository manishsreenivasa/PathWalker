/* 
 * Puppeteer - A Motion Capture Mapping Tool
 * Copyright (c) 2013-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>.
 * All rights reserved.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE* 
 */

#include "GL/glew.h"
#include "GLWidget.h"
#include <GL/glu.h>

#include <QtGui>
#include <QDebug>

#include <algorithm>
#include <iostream>
#include <cmath>

#include <assert.h>
#include "timer.h"
#include "MeshVBO.h"

using namespace std;

TimerInfo timer_info;
double draw_time = 0.;
int draw_count = 0;

Vector4f light_ka (0.2f, 0.2f, 0.2f, 1.0f);
Vector4f light_kd (0.7f, 0.7f, 0.7f, 1.0f);
Vector4f light_ks (0.7f, 0.7f, 0.7f, 1.0f);
//Vector4f light_ks (0.7f, 0.7f, 0.7f, 1.0f);
Vector4f light_position (0.f, 0.f, 0.f, 1.f);

MeshVBO grid_mesh = CreateGrid (4, 4, Vector3f (0.f, 0.f, 1.f), Vector3f (0.1f, 0.1f, 0.1), Vector3f (0.8f, 0.8f, 0.8f));

GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(parent),
		draw_base_axes (false),
		draw_grid (true),
		camera (Camera()),
		scene (NULL),
		colorPickingFrameBuffer(NULL),
		opengl_initialized (false)
{
	setFocusPolicy(Qt::StrongFocus);
	setMouseTracking(true);
}

GLWidget::~GLWidget() {
	cerr << "DESTRUCTOR: drawing time: " << draw_time << "(s) count: " << draw_count << " ~" << draw_time / draw_count << "(s) per draw" << endl;
	delete colorPickingFrameBuffer;

	makeCurrent();
}

/****************
 * Slots
 ****************/
void GLWidget::toggle_draw_grid (bool status) {
	draw_grid = status;
}

void GLWidget::toggle_draw_base_axes (bool status) {
	draw_base_axes = status;
}

void GLWidget::toggle_draw_orthographic (bool status) {
	camera.orthographic = status;

	resizeGL (static_cast<int>(windowWidth), static_cast<int>(windowHeight));
}

void GLWidget::set_front_view () {
	camera.setFrontView();
	emit camera_changed();
}

void GLWidget::set_side_view () {
	camera.setSideView();
	emit camera_changed();
}

void GLWidget::set_top_view () {
	camera.setTopView();
	emit camera_changed();
}

QSize GLWidget::minimumSizeHint() const
{
    return QSize(50, 50);
}

QSize GLWidget::sizeHint() const
{
    return QSize(400, 400);
}

void GLWidget::initializeGL()
{
	GLenum err = glewInit();
	if (GLEW_OK != err) {
		cerr << "Error initializing GLEW: " << glewGetErrorString(err) << endl;
		exit(1);
	}

	qDebug() << "Using GLEW     : " << (const char*) glewGetString (GLEW_VERSION);
	qDebug() << "OpenGL Version : " << (const char*) glGetString (GL_VERSION);
	qDebug() << "GLSL Version   : " << (const char*) glGetString (GL_SHADING_LANGUAGE_VERSION);

	if (!GLEW_ARB_shadow) {
		qDebug() << "Error: ARB_shadow not supported!";
		exit (1);
	}
  if (!GLEW_ARB_depth_texture) {
		qDebug() << "Error: ARB_depth_texture not supported!";
		exit(1);
	}

	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();

	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity();

	glShadeModel (GL_SMOOTH);
	glClearColor (0.f, 0.f, 0.f, 0.f);
	glColor4f (1.f, 1.f, 1.f, 1.f);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	glClearDepth(1.0f);
  glDepthFunc(GL_LEQUAL);
	glEnable (GL_DEPTH_TEST);
	
	glEnable (GL_NORMALIZE);

	glColorMaterial (GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glEnable (GL_COLOR_MATERIAL);
  glMaterialfv(GL_FRONT, GL_SPECULAR, Vector4f (1.f, 1.f, 1.f, 1.f).data());
  glMaterialf(GL_FRONT, GL_SHININESS, 16.0f);

	// initialize lights
	glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ka.data());
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_kd.data());
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_ks.data());

	light_position.set (3.f, 6.f, 3.f, 1.f);
	glLightfv (GL_LIGHT0, GL_POSITION, light_position.data());

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	if (scene)
		scene->initShaders();

	opengl_initialized = true;
}

void GLWidget::drawScene() {
	if (draw_grid)
		grid_mesh.draw(GL_TRIANGLES);

	timer_start (&timer_info);

	if (scene)
		scene->draw();

	draw_time += timer_stop(&timer_info);
	draw_count++;
}

void GLWidget::paintGL() {
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity();

	camera.update(width(), height());

	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLightfv(GL_LIGHT0, GL_POSITION, light_position.data());
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_kd.data());
	glEnable(GL_LIGHT0);	
	glEnable(GL_LIGHTING);

	glBegin (GL_LINES);
	glColor3f (1.f, 0.f, 0.f);
	glVertex3f (0.f, 0.f, 0.f);
	glVertex3f (1.f, 0.f, 0.f);
	glColor3f (0.f, 1.f, 0.f);
	glVertex3f (0.f, 0.f, 0.f);
	glVertex3f (0.f, 1.f, 0.f);
	glColor3f (0.f, 0.f, 1.f);
	glVertex3f (0.f, 0.f, 0.f);
	glVertex3f (0.f, 0.f, 1.f);
	glEnd();

	drawScene();

	glDisable(GL_LIGHTING);

	GLenum gl_error = glGetError();
	if (gl_error != GL_NO_ERROR) {
		cout << "OpenGL Error: " << gluErrorString(gl_error) << endl;
		abort();
	}
}

void GLWidget::paintColorPickingFrameBuffer() {
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity();

	camera.update(width(), height());

	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glDisable(GL_LIGHTING);

	if (scene)
		scene->drawForColorPicking();

	GLenum gl_error = glGetError();
	if (gl_error != GL_NO_ERROR) {
		cout << "OpenGL Error: " << gluErrorString(gl_error) << endl;
		abort();
	}
}

void GLWidget::resizeGL(int width, int height)
{
//	qDebug() << "resizing to" << width << "x" << height;

	if (height == 0)
		height = 1;

	if (width == 0)
		width = 1;

	glViewport (0, 0, width, height);

	windowWidth = width;
	windowHeight = height;

	if (colorPickingFrameBuffer)
		delete colorPickingFrameBuffer;

	QGLFramebufferObjectFormat buffer_format;
	buffer_format.setInternalTextureFormat (GL_RGBA);
	buffer_format.setAttachment(QGLFramebufferObject::Depth);
	colorPickingFrameBuffer = new QGLFramebufferObject(width, height, buffer_format);
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event) {
	if ((mousePressPos - event->pos()).manhattanLength() < 5) {
		if (event->button() == Qt::LeftButton) {
			if (!event->modifiers().testFlag(Qt::ControlModifier)) {
				scene->clearSelection();
			}

			if (scene->objectIsSelected (scene->mouseOverObjectId)) {
				scene->unselectObject (scene->mouseOverObjectId);
				emit object_unselected (scene->mouseOverObjectId);
			} else {
				scene->selectObject (scene->mouseOverObjectId);
				emit object_selected (scene->mouseOverObjectId);
			}
		}
	}
}
	
void GLWidget::mousePressEvent(QMouseEvent *event)
{
	mousePressPos = event->pos();
	lastMousePos = event->pos();
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
	float dx = static_cast<float>(event->x() - lastMousePos.x());
	float dy = static_cast<float>(event->y() - lastMousePos.y());

	if (event->buttons().testFlag(Qt::MiddleButton)
			|| ( event->buttons().testFlag(Qt::LeftButton) && event->buttons().testFlag(Qt::RightButton))) {
		camera.move (dx, dy);
		emit camera_changed();
	} else if (event->buttons().testFlag(Qt::LeftButton)) {
		// rotate
		camera.rotate (dx, dy);
		emit camera_changed();
	} else if (event->buttons().testFlag(Qt::RightButton)) {
		// zoom
		camera.zoom (dy);
		emit camera_changed();
	}

	lastMousePos = event->pos();
	updateGL();

	colorPickingFrameBuffer->bind();
	paintColorPickingFrameBuffer();
	Vector4f color;
	glReadPixels (lastMousePos.x(), static_cast<int>(windowHeight) - lastMousePos.y(), 1, 1, GL_RGBA, GL_FLOAT, color.data());
	if (scene) {
		scene->mouseOverObjectId = vector4_to_object_id(color);
	}

	colorPickingFrameBuffer->release();
}

QImage GLWidget::renderContentOffscreen (int image_width, int image_height, bool use_alpha) {
	makeCurrent();

	// set up the actual format (alpha channel, depth buffer)
	QGLFramebufferObjectFormat buffer_format;
	if (use_alpha)
		buffer_format.setInternalTextureFormat(GL_RGBA);
	else
		buffer_format.setInternalTextureFormat(GL_RGB);
	buffer_format.setAttachment (QGLFramebufferObject::Depth );

	// create the buffer object
	QGLFramebufferObject *fb = new QGLFramebufferObject(image_width, image_height, buffer_format);

	// future drawing shall be performed into this buffer
	fb->bind();

	// resize to the desired size, draw, release, and resize again to the
	// previous size

	int old_width = width();
	int old_height = height(); 

	camera.setSize (image_width, image_height);
	resizeGL(image_width, image_height);
	paintGL();

	fb->release();

	resizeGL (old_width, old_height);

	// now grab the buffer
	QImage result = fb->toImage();

	delete fb;

	return result;
}

