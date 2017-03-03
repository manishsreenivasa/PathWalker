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

#include "Scene.h"
#include "MarkerData.h"
#include "c3dfile.h"

#include <limits>

using namespace std;

MarkerData::~MarkerData() {
	if (c3dfile) {
		delete c3dfile;
		c3dfile = NULL;
	}
	for (size_t i = 0; i < markers.size(); i++) {
		scene->destroyObject<MarkerObject>(markers[i]);
	}
}

bool MarkerData::loadFromFile(const char *filename) {
	if (c3dfile) {
		delete c3dfile;
		markers.clear();
	}

	c3dfile = new C3DFile;
	if (!c3dfile->load (filename)) {
		cerr << "Error loading marker data from file '" << filename << "'!" << endl;
		abort();
		return false;
	}

	currentFrame = getFirstFrame();

	if (!scene)
		return true;

	for (size_t i = 0; i < markerNames.size(); i++) {
	  std::string tmp = markerNames.at(i);
	  enableMarker(tmp.c_str(), Vector3f(0.f, 0.f, 1.f));
	}
	return true;
}

void MarkerData::clearMarkers () {
	for (unsigned int i = 0; i < markers.size(); i++) {
		scene->destroyObject<MarkerObject>(markers[i]);
	}

	markers.clear();
}

void MarkerData::enableMarker (const char* marker_name, const Vector3f &color) {
	assert (c3dfile);

	if (markerExists(marker_name)) {
	MarkerObject* scene_marker = scene->createObject<MarkerObject>();
	scene_marker->color.block<3,1>(0,0) = color;

	Vector3f position = getMarkerCurrentPosition(marker_name);
	scene_marker->transformation.translation = position;
	scene_marker->mesh = CreateUVSphere (4, 8);
	scene_marker->transformation.scaling = Vector3f (0.02f, 0.02f, 0.02f);
	scene_marker->noDepthTest = true;
	scene_marker->markerName = marker_name;

	markers.push_back (scene_marker);
	} else {
	  std::cout << "!! WARNING::Marker " << marker_name << " does not exist" << std::endl;
	}
}

bool MarkerData::markerExists(const char* marker_name) {
	std::string point_label(marker_name);
	
	point_label = point_label.substr(0, point_label.find_last_not_of(" ") + 1);

	if (c3dfile->label_point_map.find(point_label) == c3dfile->label_point_map.end()) {
		return false;
	}

	return true;
}

Vector3f MarkerData::getMarkerCurrentPosition(const char * marker_name) {
	FloatMarkerData marker_traj = c3dfile->getMarkerTrajectories (marker_name);

	int index = currentFrame - getFirstFrame();

	if (rotateZ) 
		return Vector3f (-marker_traj.x[index], -marker_traj.y[index], marker_traj.z[index]) * 1.0e-3;
	
	return Vector3f (marker_traj.x[index], marker_traj.y[index], marker_traj.z[index]) * 1.0e-3;
}

std::string MarkerData::getMarkerName (int object_id) {
	for (size_t i = 0; i < markers.size(); i++) {
		if (markers[i]->id == object_id) 
			return markers[i]->markerName;
	}

	cerr << "Error: could not find marker with object id " << object_id << "!" << endl;
	abort();

	return "Error";
}

int MarkerData::getFirstFrame () {
	assert (c3dfile);

	return static_cast<int>(c3dfile->header.first_frame);
}

int MarkerData::getLastFrame () {
	assert (c3dfile);

	return static_cast<int>(c3dfile->header.last_frame);
}

float MarkerData::getFrameRate () {
	assert (c3dfile);

	return c3dfile->header.video_sampling_rate;
}

void MarkerData::setCurrentFrameNumber (int frame_number) {
	assert (frame_number >= getFirstFrame());
	assert (frame_number <= getLastFrame());

	currentFrame = frame_number;

	updateMarkerSceneObjects();
}

void MarkerData::updateMarkerSceneObjects() {
	for (size_t i = 0; i < markers.size(); i++) {
		Vector3f position = getMarkerCurrentPosition(markers[i]->markerName.c_str());
		markers[i]->transformation.translation = position;
	}
}

void MarkerData::calcDataBoundingBox(Vector3f &min, Vector3f &max) {
	min = Vector3f (std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	max = -min;

	int current_frame_temp = currentFrame;

	for (int frame = getFirstFrame(); frame <= getLastFrame(); frame++) {
		setCurrentFrameNumber (frame);
		
		for (size_t mi = 0; mi < markers.size(); mi++) {
			Vector3f pos = getMarkerCurrentPosition (markers[mi]->markerName.c_str());

			for (size_t i = 0; i < 2; i++) {
				min[i] = std::min(pos[i], min[i]);
				max[i] = std::max(pos[i], max[i]);
			}
		}
	}

	setCurrentFrameNumber(current_frame_temp);
}
