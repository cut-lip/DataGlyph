#pragma once
#include <vector>
#include "turtleg.h"
#include "myglfuncs.h"

const float SF_SEGMENT_CONSTANT = 0.1;
const float SF_ANGLE_SCALE = 2;
const float GLYPH_SCALE_FACTOR = 0.25;

// An individual SPC-SF glyph located in 2-D space
class SpcSfGlyph {
public:
	// Constructor
	SpcSfGlyph(GLfloat x, GLfloat y)
	{
		this->x = x;
		this->y = y;
	}

	// Draw the three SPC coordinate pairs
	void drawAxesSPC(Point2 pos2, Point2 pos3, std::vector<float>::iterator axes, GLfloat axisLength)
	{
		// Isolate axes attributes
		float x1 = *axes * GLYPH_SCALE_FACTOR;
		float x2 = *++axes * GLYPH_SCALE_FACTOR;
		float x3 = *++axes * GLYPH_SCALE_FACTOR;
		float x4 = *++axes * GLYPH_SCALE_FACTOR;
		float x5 = *++axes * GLYPH_SCALE_FACTOR;
		float x6 = *++axes * GLYPH_SCALE_FACTOR;

		// Enable line stipple
		glPushAttrib(GL_ENABLE_BIT);
		// glPushAttrib is done to return everything to normal after drawing

		glLineStipple(3, 0xAAAA);  // [1]
		glEnable(GL_LINE_STIPPLE);

		// 1st pair: origin (0 - X1, 0 - X2)
		glColor3f(0.0, 1.0, 0.0);

		glBegin(GL_LINES);
		glVertex2f(-x1, -x2);
		glVertex2f(-x1, -x2 + (axisLength * GLYPH_SCALE_FACTOR));
		glEnd();
		glBegin(GL_LINES);
		glVertex2f(-x1, -x2);
		glVertex2f(-x1 + (axisLength * GLYPH_SCALE_FACTOR), -x2);
		glEnd();

		// 2nd pair: origin ( )
		glColor3f(0.0, 0.7, 0.0);

		glBegin(GL_LINES);
		glVertex2f(pos2.getx() - x3, pos2.gety() - x4);
		glVertex2f(pos2.getx() - x3, pos2.gety() - x4 + (axisLength * GLYPH_SCALE_FACTOR));
		glEnd();
		glBegin(GL_LINES);
		glVertex2f(pos2.getx() - x3, pos2.gety() - x4);
		glVertex2f(pos2.getx() - x3 + (axisLength * GLYPH_SCALE_FACTOR), pos2.gety() - x4);
		glEnd();

		// 3nd pair: origin ( )
		glColor3f(0.0, 0.3, 0.0);

		glBegin(GL_LINES);
		glVertex2f(pos3.getx() - x5, pos3.gety() - x6);
		glVertex2f(pos3.getx() - x5, pos3.gety() - x6 + (axisLength * GLYPH_SCALE_FACTOR));
		glEnd();
		glBegin(GL_LINES);
		glVertex2f(pos3.getx() - x5, pos3.gety() - x6);
		glVertex2f(pos3.getx() - x5 + (axisLength * GLYPH_SCALE_FACTOR), pos3.gety() - x6);
		glEnd();

		glPopAttrib();
	}

	// Draw 2-segment SF (Stick Figure) glyph
	void drawGlyphSF(Point2* pos2, Point2* pos3,
		std::vector<float>::iterator stick, bool benign, TurtleG turt)
	{
		// First segment
		if (benign) glColor3f(0.0, 0.0, 1.0);
		else glColor3f(1.0, 0.0, 0.0);
		turt.turnTo(*stick * GLYPH_SCALE_FACTOR);
		turt.forward((*++stick * GLYPH_SCALE_FACTOR) + SF_SEGMENT_CONSTANT, true);
		*pos2 = turt.getCP();

		//Second segment
		if (benign) glColor3f(0.0, 0.0, 0.8);
		else glColor3f(0.8, 0.0, 0.0);
		turt.turn(SF_ANGLE_SCALE * (*++stick * GLYPH_SCALE_FACTOR));
		turt.forward((*++stick * GLYPH_SCALE_FACTOR) + SF_SEGMENT_CONSTANT, true);
		*pos3 = turt.getCP();
	}

private:
	// 2-D position
	GLfloat x;
	GLfloat y;

};