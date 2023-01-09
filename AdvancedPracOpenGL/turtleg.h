#pragma once
#include "point2.h"
#include "GL/glut.h"

// Turtle Graphics drawing tool
class TurtleG {
public:
	//Default constructor
	TurtleG()
	{
		this->CP = Point2(0.0, 0.0);
		this->CD = 0.0;
	}

	// Constructor from Point2
	TurtleG(GLfloat x, GLfloat y, GLfloat angle)
	{
		this->CP = Point2(x, y);
		this->CD = angle;
	}

	// Draw a line from CP to the new vertex
	void lineTo(float x, float y)
	{
		glBegin(GL_LINES);
		glVertex2f((GLfloat)CP.getx(), (GLfloat)CP.gety());
		glVertex2f((GLfloat)x, (GLfloat)y);
		glEnd();

		CP.set(x, y);	// Update the CP
		glFlush();		// Flush the buffer
	}

	void lineTo(Point2 p)
	{
		glBegin(GL_LINES);
		glVertex2f((GLfloat)CP.getx(), (GLfloat)CP.gety());
		glVertex2f((GLfloat)p.getx(), (GLfloat)p.gety());
		glEnd();

		CP.set(p.getx(), p.gety());	// Update the CP
		glFlush();		// Flush the buffer
	}

	// Update the CP
	void moveTo(float x, float y) { CP.set(x, y); }

	void moveTo(Point2 p) { CP.set(p.getx(), p.gety()); }

	// turn the turtle to given angle
	void turnTo(float angle)
	{
		// Convert radians to degrees
		// REPLACE with PI constant
		this->CD = - (angle * (180.0 / 3.141592653589793238463));
	}
	
	// turn the turtle given number of degrees
	void turn(float angle)
	{
		// Convert radians to degrees
		CD += 20 + (angle * (180.0 / 3.141592653589793238463));
	}

	// move turtle forward in a straight line from CP
	void forward(float dist, bool isVisible)
	{
		// Determine endpoint based on radial distance
		const float radPerDeg = 0.017453393;
		float x = CP.getx() + dist * cos(radPerDeg * CD);
		float y = CP.gety() + dist * sin(radPerDeg * CD);

		if (isVisible) lineTo(x, y);	// Move CP
		else moveTo(x, y);
	}

	void setCP(GLfloat x, GLfloat y)
	{
		this->CP = Point2(x, y);
	}

	Point2 getCP()
	{
		return this->CP;
	}

	void setCD(GLfloat angle)
	{
		this->CD = angle;
	}

	GLfloat getCD()
	{
		return this->CD;
	}

private:
	// Current position
	Point2 CP = Point2(0.0, 0.0);

	// Current direction
	GLfloat CD = 0.0;
};

