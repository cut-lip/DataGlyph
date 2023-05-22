#pragma once
#include <vector>
#include "turtleg.h"
#include "myglfuncs.h"

float COLOR_0 = 0.0;
float COLOR_1 = 1.0;
float ALPHA_CROSS = 0.7;

float CROSS_SCALE = 0.2;
int CROSS_REDUCE = 2;

float SF_SEGMENT_CONSTANT = 0.5;
float SF_ANGLE_SCALE = 2.0;
float GLYPH_SCALE_FACTOR = 0.25;
float CROSS_SCALE_FACTOR = 0.07;

float R_MIN = 0.1;
float R_MAX = 1.0;
float T_MIN = 10.0;
float T_MAX = 85.0;

// An individual SPC-SF glyph located in 2-D space
class SpcSfGlyph {
public:
	// Constructor
	SpcSfGlyph() {}

	// Draw the three SPC coordinate pairs
	void drawAxesSPC(Point2 pos2, Point2 pos3, std::vector<float>::iterator axes, GLfloat axisLength,
		GLfloat glyphScaleFactor, bool dottedAxes)
	{
		// Initialize scaling params
		GLYPH_SCALE_FACTOR = glyphScaleFactor;

		// Isolate axes attributes
		float x1 = *axes;// * GLYPH_SCALE_FACTOR;
		float x2 = *++axes;// * GLYPH_SCALE_FACTOR;
		float x3 = *++axes;// * GLYPH_SCALE_FACTOR;
		float x4 = *++axes;// * GLYPH_SCALE_FACTOR;
		float x5 = *++axes;// * GLYPH_SCALE_FACTOR;
		float x6 = *++axes;// * GLYPH_SCALE_FACTOR;

		if (dottedAxes)
		{
			// Enable line stipple
			glPushAttrib(GL_ENABLE_BIT);
			// glPushAttrib is done to return everything to normal after drawing

			glLineStipple(3, 0xAAAA);  // [1]
			glEnable(GL_LINE_STIPPLE);
		}

		// 1st pair: origin (0 - X1, 0 - X2)
		glColor4f(COLOR_0, COLOR_1, COLOR_0, ALPHA_CROSS);	// Set color to light green
		glBegin(GL_LINES);	// Draw first cross vertical
		glVertex2f(-x1,
			-x2 - (((axisLength * GLYPH_SCALE_FACTOR) + CROSS_SCALE) / CROSS_REDUCE) + CROSS_SCALE_FACTOR);
		glVertex2f(-x1,
			-x2 + (((axisLength * GLYPH_SCALE_FACTOR) + CROSS_SCALE) / CROSS_REDUCE) + CROSS_SCALE_FACTOR);
		glEnd();	// end horizontal line

		glBegin(GL_LINES);	// Draw first cross horizontal
		glVertex2f(-x1 - (((axisLength * GLYPH_SCALE_FACTOR) + CROSS_SCALE) / CROSS_REDUCE) + CROSS_SCALE_FACTOR,
			-x2);
		glVertex2f(-x1 + (((axisLength * GLYPH_SCALE_FACTOR) + CROSS_SCALE) / CROSS_REDUCE) + CROSS_SCALE_FACTOR,
			-x2);
		glEnd();	// end vertical line

		// 2nd pair: origin ( )
		glColor4f(0.0, 0.7, 0.0, 0.7);	// Set color to medium green
		glBegin(GL_LINES);	// Draw second cross vertical
		glVertex2f(pos2.getx() - x3, pos2.gety() - x4 - (((axisLength * GLYPH_SCALE_FACTOR) + 0.2) / 2)+ CROSS_SCALE_FACTOR);
		glVertex2f(pos2.getx() - x3, pos2.gety() - x4 + (((axisLength * GLYPH_SCALE_FACTOR) + 0.2) / 2)+ CROSS_SCALE_FACTOR);
		glEnd();
		glBegin(GL_LINES);	// Draw second cross horizontal
		glVertex2f(pos2.getx() - x3 - (((axisLength * GLYPH_SCALE_FACTOR) + 0.2) / 2)+ CROSS_SCALE_FACTOR, pos2.gety() - x4);
		glVertex2f(pos2.getx() - x3 + (((axisLength * GLYPH_SCALE_FACTOR) + 0.2) / 2)+ CROSS_SCALE_FACTOR, pos2.gety() - x4);
		glEnd();

		// 3nd pair: origin ( )
		glColor4f(0.0, 0.3, 0.0, 0.7);	// Set color to dark green
		glBegin(GL_LINES);	// Draw third cross vertical
		glVertex2f(pos3.getx() - x5, pos3.gety() - x6 - (((axisLength * GLYPH_SCALE_FACTOR) + 0.2) / 2) + CROSS_SCALE_FACTOR);
		glVertex2f(pos3.getx() - x5, pos3.gety() - x6 + (((axisLength * GLYPH_SCALE_FACTOR) + 0.2) / 2) + CROSS_SCALE_FACTOR);
		glEnd();
		glBegin(GL_LINES);	// Draw third cross horizontal
		glVertex2f(pos3.getx() - x5 - (((axisLength * GLYPH_SCALE_FACTOR) + 0.2) / 2) + CROSS_SCALE_FACTOR, pos3.gety() - x6);
		glVertex2f(pos3.getx() - x5 + (((axisLength * GLYPH_SCALE_FACTOR) + 0.2) / 2) + CROSS_SCALE_FACTOR, pos3.gety() - x6);
		glEnd();

		if (dottedAxes)
		{
			glPopAttrib();
		}
	}

	// Draw 2-segment SF (Stick Figure) glyph
	void drawGlyphSF(Point2* pos2, Point2* pos3,
		std::vector<GLfloat>::iterator stick, bool benign, TurtleG turt, bool dynamicAngles, bool posAngle,
		GLfloat glyphScaleFactor, GLfloat sfSegmentConstant, GLfloat sfAngleScale, bool ANGLE_FOCUS, bool birdFocus,
		float colors[6] )
	{
		// Initialize scaling params
		GLYPH_SCALE_FACTOR = glyphScaleFactor;
		SF_SEGMENT_CONSTANT = sfSegmentConstant;
		SF_ANGLE_SCALE = sfAngleScale;

		// Determine color by class
		if (benign)
		{	// If class is benign
			if (birdFocus)
			{	// If bird focus is on, color
				glColor4f(colors[0], colors[1], colors[2], 1.0);
			}
			else
			{	// If bird focus is off, grey out
				glColor4f(0.0, 0.0, 0.0, 0.3);
			}
		}
		else
		{	// If class is malignant
			if (birdFocus)
			{	// If bird focus is on, color
				glColor4f(colors[0], colors[1], colors[2], 1.0);
			}
			else
			{	// If bird focus is off, grey out
				glColor4f(0.0, 0.0, 0.0, 0.3);
			}
		}

		GLfloat ang1 = 0;
		GLfloat dist1 = 0;
		GLfloat ang2 = 0;
		GLfloat dist2 = 0;
		if (ANGLE_FOCUS)
		{
			ang1 = *stick;
			dist1 = *(stick + 1);
			ang2 = *(stick + 2);
			dist2 = *(stick + 3);
		}
		else if (!ANGLE_FOCUS)
		{
			dist1 = *(stick + 2);
			ang1 = *(stick + 3);
			dist2 = *stick;
			ang2 = *(stick + 1);

		}

		// First segment
		if (posAngle)
		{
			turt.turnTo((normalizeAngle(ang1)));	// Positive angle
		}
		else
		{
			turt.turnTo((-normalizeAngle(ang1)));	// Negative angle
		}

		turt.forward(dist1, true);
		*pos2 = turt.getCP();

		//Determine color by class
		//if (benign) glColor4f(0.0, 0.0, 0.8, 0.7);
		//else glColor4f(0.8, 0.0, 0.0, 0.7);

		// Determine color by class
		if (benign)
		{
			if (birdFocus)
			{
				glColor4f(colors[3], colors[4], colors[5], 0.7);
			}
			else
			{
				glColor4f(0.0, 0.0, 0.0, 0.3);
			}
		}
		else
		{
			if (birdFocus)
			{
				glColor4f(colors[3], colors[4], colors[5], 0.7);
			}
			else
			{
				glColor4f(0.0, 0.0, 0.0, 0.3);
			}
		}

		//Second segment
		if (dynamicAngles)
		{
			turt.turn(normalizeAngle(ang2));	// Dynamic angles
		}
		else 
		{
			turt.turnTo(normalizeAngle(ang2));	// Static angles
		}

		turt.forward(dist2, true);
		*pos3 = turt.getCP();
	}

private:

	// Normalize data value to angle value
	GLfloat normalizeAngle(GLfloat m)
	{
		// Range 0-85 degrees
		GLfloat rmin = R_MIN;		// Min value data range
		GLfloat rmax = R_MAX;		// Max value data range
		GLfloat tmin = T_MIN;		// Min value target range
		GLfloat tmax = T_MAX;	// Max value target range

		// Normalization formula
		return (((m - rmin) / (rmax - rmin)) * (tmax - tmin)) + tmin;
	}
};