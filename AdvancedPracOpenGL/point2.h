#pragma once
#include "GL/glut.h"

// Point2 subclass
class Point2 {	// A simple 2-D point
public:
	// Default constructor
	Point2() { x = 0.0; y = 0.0; }

	// Constructor
	Point2(GLfloat x, GLfloat y)
	{
		this->x = x;
		this->y = y;
	}

	void set(GLfloat x, GLfloat y)
	{
		this->x = x;
		this->y = y;
	}

	GLfloat setx(GLfloat x) { this->x = x; }
	GLfloat sety(GLfloat y) { this->y = y; }

	GLfloat getx() { return this->x; }
	GLfloat gety() { return this->y; }

private:
	GLfloat x;
	GLfloat y;
};