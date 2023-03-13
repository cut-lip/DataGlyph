#pragma once
#include "point2.h"
#include "GL/glut.h"

// Constant value for converting radians to degrees
const float RADIAN_PER_DEGREE = 0.017453393;

/* TurtleG: Turtle Graphics drawing tool
*  This object is used to perform line drawing
*  operations using the "Turtle Graphics"
*  concept, as described in the OpenGl Red Book.
*  Chapter 3.
* 
*  This object draws lines based on initial and final
*  position of the "turtle", and the currect angle
*  (in degrees) the turtle is pointed torward.
*
*/
class TurtleG {
public:
	/*
	* TurtleG
	* Default constructor. Construct new
	* TurtleG object with current position (CP)
	* set to (0.0, 0.0).
	*/
	TurtleG()
	{
		this->CP = Point2(0.0, 0.0);
		this->CD = 0.0;
	}

	/*
	* TurtleG
	* Constructor. Construct new TurtleG object
	* with specified current position (CP) and
	* current direction (CD).
	* TurtleG object with current position (CP)
	* set to (0.0, 0.0)
	*
	* @param	x		x coordinate of CP
	*			y		y coordinate of CP
	*			angle	initial direction of CD
	*/
	TurtleG(GLfloat x, GLfloat y, GLfloat angle)
	{
		this->CP = Point2(x, y);
		this->CD = angle;
	}

	/*
	* turnTo
	* Turn the turtle to given angle from horizontal.
	* In SPC-SF context, can be considered static angles.
	* 
	* @param angle Value in degrees to  be turned to
	*/
	void turnTo(GLfloat angle)
	{
		this->CD = angle;
	}
	
	/*
	* turn the turtle to given angle relative to current direction
	* In SPC-SF context, can be considered dynamic angles
	* 
	* @param angle Value in degrees to be turned to
	*/
	void turn(GLfloat angle)
	{
		this->CD += angle;
	}

	/*
	* forward
	* Move the turtle forward a given distance,
	* drawing a line from start to finish.
	*
	* @param	distance	the distance to be traveled
	*			isVisible	draw the line or not
	*/
	void forward(float distance, bool isVisible)
	{
		// Determine endpoint based on radial distance
		float x = CP.getx() + distance * cos(RADIAN_PER_DEGREE * CD);
		float y = CP.gety() + distance * sin(RADIAN_PER_DEGREE * CD);

		if (isVisible) lineTo(x, y);	// Move CP
		else moveTo(x, y);
	}

	/*
	* setCP
	* Set the Turtle's Current Position (CP).
	*
	* @param	x	x coordinate of new CP
	*			y	y coordinate of new CP
	*/
	void setCP(GLfloat x, GLfloat y)
	{
		this->CP = Point2(x, y);
	}

	/*
	* getCP
	* Get the Turtle's Current Position (CP).
	*
	* @return	Point2	turtle's current position
	*/
	Point2 getCP()
	{
		return this->CP;
	}

	/*
	* setCD
	* Set the Turtle's Current Direction (CD).
	*
	* @param	angle	direction of the turtle
	*/
	void setCD(GLfloat angle)
	{
		this->CD = angle;
	}

	/*
	* getCD
	* Get the Turtle's Current Direction (CD).
	*
	* @return	GLfloat	current direction of turtle
	*/
	GLfloat getCD()
	{
		return this->CD;
	}

private:

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

	// Current position
	// (defines current location of turtle in viewport coordinates)
	Point2 CP = Point2(0.0, 0.0);

	// Current direction
	// (defines current angle pointed in by turtle, from 0 to 360 degrees)
	GLfloat CD = 0.0;
};

