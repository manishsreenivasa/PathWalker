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

#include <iostream>
#include <assert.h>

#include "GL/glew.h"
#include <GL/glu.h>

#include "Scene.h"
#include "Shader.h"

using namespace std;

Vector4f object_id_to_vector4 (int id) {
	Vector4f result (0.f, 0.f, 0.f, 1.f);

	unsigned char byte = ((id + 1) & 0x0000ff);
	result[2] = static_cast<float>(byte) / 255.f;
	byte = ((id + 1) >> 8 & 0x0000ff);
	result[1] = static_cast<float>(byte) / 255.f;
	byte = ((id + 1) >> 16 & 0x0000ff);
	result[0] = static_cast<float>(byte) / 255.f;
	
	return result;
}

int vector4_to_object_id (const Vector4f &color) {
	int result = 0;

	if (color[0] != 0.f) {
		result += (static_cast<int> (color[0] * 255.f) << 16);
	}
	
	if (color[1] != 0.f) {
		result += (static_cast<int> (color[1] * 255.f) << 8);
	}
	
	if (color[2] != 0.f) {
		result += static_cast<int> (color[2] * 255.f);
	}

	return result - 1;
}

void Scene::initShaders() {
//	defaultShader = ShaderProgram::createFromFiles ("shaders/vertex_shader.glsl", "shaders/fragment_shader.glsl");
}

void Scene::drawSceneObjectStyled (const SceneObject *object, DrawStyle style) {
	if (style == DrawStyleHidden || object->noDraw)
		return;

//	defaultShader.printLog();
//	glUseProgram(defaultShader.program_id);
//	defaultShader.setUniformFloat("alpha", 1.0);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glPushMatrix();
	glMultMatrixf (object->transformation.toGLMatrix().data());

	glEnable (GL_DEPTH_TEST);

	if (style == DrawStyleSelected) {
		glDisable(GL_LIGHTING);
		glEnable(GL_CULL_FACE);
		glCullFace(GL_FRONT);
		glPushMatrix();
		glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth (3.f);
		glColor3f (1.f, 0.f, 0.f);
		const_cast<MeshVBO*>(&(object->mesh))->draw(GL_TRIANGLES);
		glPopMatrix();
		glCullFace(GL_BACK);
		glDisable(GL_CULL_FACE);
		glLineWidth (1.f);
		glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

		if (!object->noLighting)
			glEnable(GL_LIGHTING);
		glColor4fv (object->color.data());
		const_cast<MeshVBO*>(&(object->mesh))->draw(GL_TRIANGLES);
	} else if (style == DrawStyleHighlighted) {
		glDisable(GL_LIGHTING);
		glEnable(GL_CULL_FACE);
		glCullFace(GL_FRONT);
		glPushMatrix();
		glScalef (1.03f, 1.03f, 1.03f);
		glColor4f (0.9, 0.9, 0.3, object->color[3]);
		const_cast<MeshVBO*>(&(object->mesh))->draw(GL_TRIANGLES);
		glPopMatrix();
		glCullFace(GL_BACK);
		glDisable(GL_CULL_FACE);
		
		if (!object->noLighting)
			glEnable(GL_LIGHTING);

		glColor4f (0.8, 0.8, 0.2, object->color[3]);
		const_cast<MeshVBO*>(&(object->mesh))->draw(GL_TRIANGLES);
	} else {
		if (object->noLighting)
			glDisable(GL_LIGHTING);
		else
			glEnable(GL_LIGHTING);

		glColor4fv (object->color.data());
		const_cast<MeshVBO*>(&(object->mesh))->draw(GL_TRIANGLES);
	}
	glDisable(GL_BLEND);

	glPopMatrix();
}

void Scene::draw() {
	std::vector<SceneObject*> depth_ignoring_objects;

	for (size_t i = 0; i < objects.size(); i++) {
		if (objects[i]->noDepthTest) {
			depth_ignoring_objects.push_back (objects[i]);
			continue;
		}
		if (objectIsSelected(objects[i]->id)) {
			drawSceneObjectStyled (objects[i], DrawStyleSelected);
		} else if (objects[i]->id == mouseOverObjectId) {
			drawSceneObjectStyled (objects[i], DrawStyleHighlighted);
		} else {
			drawSceneObjectStyled (objects[i], DrawStyleNormal);
		}
	}

	// Draw the outline of the selected objects using stencil buffers
	glClear (GL_DEPTH_BUFFER_BIT);

	glEnable (GL_STENCIL_TEST);
	glClearStencil(0);
	glClear (GL_STENCIL_BUFFER_BIT);

	std::list<int>::iterator selected_iter = selectedObjectIds.begin();

	while (selected_iter != selectedObjectIds.end()) {
		if (*selected_iter == -1)
			break;
		SceneObject* object = getObject<SceneObject>(*selected_iter);

		glPushMatrix();
		glMultMatrixf (object->transformation.toGLMatrix().data());

		// 1st draw: draw the wireframe model of the back and set the stencil
		// buffer to 1
		glDisable (GL_LIGHTING);
		glDisable (GL_BLEND);

		glPolygonMode (GL_BACK, GL_LINE);
		glLineWidth (5.f);
		glColor3f (1.f, 0.f, 0.f);

		glStencilFuncSeparate(GL_BACK, GL_ALWAYS, 1, 1);
		glStencilOpSeparate (GL_BACK, GL_KEEP, GL_KEEP, GL_KEEP);
		const_cast<MeshVBO*>(&(object->mesh))->draw(GL_TRIANGLES);

		// 2nd draw: draw the regular model and set the stencil buffer to 0
		glEnable (GL_LIGHTING);
		glLineWidth (1.f);
		glPolygonMode (GL_FRONT, GL_FILL);
		glColor4fv (object->color.data());
		glStencilFuncSeparate(GL_FRONT, GL_ALWAYS, 0, 1);
		glStencilOpSeparate (GL_FRONT, GL_REPLACE, GL_REPLACE, GL_REPLACE);
		const_cast<MeshVBO*>(&(object->mesh))->draw(GL_TRIANGLES);

		glPopMatrix();

		selected_iter++;
	}

	glClear (GL_DEPTH_BUFFER_BIT);
	glPolygonMode (GL_BACK, GL_FILL);
	glDisable (GL_STENCIL_TEST);

	for (size_t i = 0; i < depth_ignoring_objects.size(); i++) {
		drawSceneObjectStyled (depth_ignoring_objects[i], DrawStyleNormal);
		if (objectIsSelected(depth_ignoring_objects[i]->id)) {
			drawSceneObjectStyled (depth_ignoring_objects[i], DrawStyleSelected);
		} else if (depth_ignoring_objects[i]->id == mouseOverObjectId) {
			drawSceneObjectStyled (depth_ignoring_objects[i], DrawStyleHighlighted);
		} else {
			drawSceneObjectStyled (depth_ignoring_objects[i], DrawStyleNormal);
		}
	}

}

void Scene::drawForColorPicking() {
	glDisable(GL_LIGHTING);

	std::vector<SceneObject*> depth_ignoring_objects;

	for (size_t i = 0; i < objects.size(); i++) {
		if (objects[i]->noDraw)
			continue;

		if (objects[i]->noDepthTest) {
			depth_ignoring_objects.push_back (objects[i]);
			continue;
		}

		glPushMatrix();
		glMultMatrixf (objects[i]->transformation.toGLMatrix().data());

		glColor4fv (object_id_to_vector4 (objects[i]->id).data());
		const_cast<MeshVBO*>(&(objects[i]->mesh))->draw(GL_TRIANGLES);

		glPopMatrix();
	}

	glClear (GL_DEPTH_BUFFER_BIT);

	for (size_t i = 0; i < depth_ignoring_objects.size(); i++) {
		glPushMatrix();
		glMultMatrixf (depth_ignoring_objects[i]->transformation.toGLMatrix().data());

		glColor4fv (object_id_to_vector4 (depth_ignoring_objects[i]->id).data());
		const_cast<MeshVBO*>(&(depth_ignoring_objects[i]->mesh))->draw(GL_TRIANGLES);

		glPopMatrix();
	}
}

void Scene::selectObject (const int id) {
	if (objectIsSelected(id))
		return;

	selectedObjectIds.push_back(id);
}

void Scene::unselectObject (const int id) {
	list<int>::iterator iter = selectedObjectIds.begin();
	do {
		if (*iter == id) {
			selectedObjectIds.erase(iter);
			return;
		}
		iter++;
	} while (iter != selectedObjectIds.end());
}

bool Scene::objectIsSelected (const int id) const {
	for (list<int>::const_iterator iter = selectedObjectIds.begin(); iter != selectedObjectIds.end(); iter++) {
		if (*iter == id) {
			return true;
		}
	}

	return false;
}

void Scene::unregisterSceneObject (const int id) {
	std::vector<SceneObject*>::iterator obj_iter = objects.begin();

	unselectObject (id);

	do {
		if ((*obj_iter)->id == id)  {
			objects.erase (obj_iter); 
			return;
		}

		obj_iter++;
	} while (obj_iter != objects.end());

	cerr << "Error deleting object with id " << id << ": object not found." << endl;
	abort();
}
