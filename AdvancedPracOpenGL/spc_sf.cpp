// spc_sf.cpp : This file implements the hybrid SPC-SF data visualization method
//				described in 'Visual Knowledge Discovery and Machine Learning'
//				by Dr. Boris Kovalerchuk, Figure 13.2
// AUTHOR:	Nicholas Cutlip (Central Washington University Computer Science)
// DATE:	v 0.0: 29 December, 2022
//			v 0.5: 27 February, 2023

	/* ********************************* WBC Attribute Indices ************************************
		#  Attribute                     Domain
		-- ---------------------------------------- -
		0. Clump Thickness               1 - 10
		1. Uniformity of Cell Size       1 - 10
		2. Uniformity of Cell Shape      1 - 10
		3. Marginal Adhesion             1 - 10
		4. Single Epithelial Cell Size   1 - 10
		5. Bare Nuclei                   1 - 10
		6. Bland Chromatin               1 - 10
		7. Normal Nucleoli               1 - 10
		8. Mitoses                       1 - 10
		Class:                          (2 for benign, 4 for malignant)
	***********************************************************************************************/


/**************** INCLUDES ********************/
#include <iostream>		/*  */
#include <sstream>		/*  */
#include <fstream>		/*  */
#include <string>		/*  */
#include <vector>		/*  */
#include "GL/glut.h"	/*  */
#include "turtleg.h"	/*  */
#include "spc_sf.h"		/*  */
#include <cmath>
#include <iomanip>
#include <sstream>
#include <algorithm>

// NEW PARAMS
const int HYPERBLOCK_SIZE = 345;
const int HYPERBLOCK_DATA_SIZE = 8;

const int SEED_DATASET_SIZE = 210;
const int SEED_DATA_SIZE = 6;


/************************* DATA AND DIMENSION CONSTANTS *******************************/
int SCREEN_WIDTH;						/* Screen Width */
int SCREEN_HEIGHT;						/* Screen Height */
int VIEWPORT_SCALE = 8;					/* viewport scaling constant */
float THRESHOLD_VALUE = 4.0;			/* threshold calculation value */
float MIN_THRESHOLD = 2.0;
float ALLOWED_DIFFERENCES = 2;
float AXIS_LENGTH = 1.0;				/* SPC axis length constant 8 */
unsigned int DATA_SIZE = 683;			/* cardinality of data set */
unsigned int DATA_INDEX = 0;			/* current index of data set */
unsigned int MAX_SIG_INDEX = 9;			/* maximum significant data index */
const unsigned int HEIGHT_SCALE = 10;	/* height scaling constant */
const unsigned int WIDTH_SCALE = 4;		/* width scaling constant */
const float GRID_MARGIN = 4.0;			/* glyph grid margin*/
const float MAR = 20.0;					/* general use margin */

/***************************** DISPLAY FLAGS *******************************************/
bool DRAW_EDGES = true;				/* toggle drawing edges between glyphs in PC-SPC-SF */
bool DRAW_AXES = true;				/* toggle glyph SPC axes on / off*/
bool DISPLAY_ALL = false;			/* toggle entire dataset / single neighborhood views */
bool DISPLAY_HYPERCUBES = true;		/* toggle displaying hypercubes overlaying PC-SPC-SF */
bool DISPLAY_SELECTOR = true;		/* toggle PC-SPC-SF and glyph grid views */
bool DYNAMIC_ANGLES = false;		/* toggle computing second SF angle dynamically */
bool POS_ANGLE = false;				/* toggle (+)/(-) first SF angle */
bool PC_OFF = true;					/* toggle Paired Coords view */
bool CLASS_SEPERATION_MODE = true;	/* toggle which class to display first in PC_SPC_SF */
bool DOTTED_AXES = false;			/* toggle dotted / solid SPC axes in glyphs */
bool REPS_OFF = true;				/* toggle grid of representative glyphs */
bool ANGLE_FOCUS = true;			/* toggle focus on divergent glyph lengths or angles */
bool BIRD_FOCUS = true;				/* toggle grey / colored SF birds */

/************************* GLYPH GRID CONSTANTS  *******************************/
int NUM_ROWS = 4;						/* row in glyph grid */
int NUM_COLUMNS = 9;					/* columns in glyph grid */
float GLYPH_SCALE = 1.6;				/* glyph scale factor */
float GLYPH_TRANSLATE_FACTOR = -0.1;	/* glyph translate factor */

/********************** REPRESENTATIVE GLYPH CONSTANTS  ***************************/
std::vector<std::string> labels{};			/* vector of labels of representative glyphs */
std::vector<std::vector<GLfloat>> reps{};	/* vector of representative glyphs */
std::vector<bool> repsClass{};				/* vector of classes of representative glyphs */
bool REPS_COLLECTED = false;				/* flag if rep glyphs have already been collected */
bool IDEAL_COLLECTED = false;				/* flag if ideal glyphs h   ave already been collected */

/************** DEBUG DATA STRUCTS (remove later) *******************/
std::vector<std::vector<GLfloat>> mixedHood{};
std::vector<bool> mixedClass{};
std::vector<std::string> mixedLabels{};

// IDEAL BENIGN POINT
std::vector<GLfloat> idealBenign = { 1, 3, 2, 3, 1, 1, 1, 2, 1, 1 };
// IDEAL MALIGNANT POINT
std::vector<GLfloat> idealMalig = { 7, 7, 6, 7, 6, 8, 7, 5, 6, 3 };

unsigned int CLUSTER_CENTER_FOCUS = 2;
const int CLUSTER = 0;

// collect cluster distance labels
std::vector<std::string> distanceLabels{};

/*
displayHypercubes
This function overlays transparent rectangles over
the SPC-SF visualization. These rectangles are
representative of hypercubes representing branches
of the decision tree implemented for WBC data.
@param			none
@return			void
*/
void displayHypercubes()
{
	glPushMatrix();
	glEnable(GL_BLEND); //Enable blending.
	glDepthMask(GL_FALSE);

	glViewport(0.0, 0.0, SCREEN_WIDTH, SCREEN_HEIGHT);
	gluOrtho2D(0.0, SCREEN_WIDTH, 0.0, SCREEN_HEIGHT);

	// Draw hypercubes based on DT paper
	// Cube a
	glColor4f(0.0, 0.0, 0.0, 0.3);
	glRectf(MAR, MAR, ((SCREEN_WIDTH / 4) + (SCREEN_WIDTH / 24)) * 0.25, (SCREEN_HEIGHT - 10));

	// Cube b
	glColor4f(0.0, 0.0, 1.0, 0.3);
	glRectf(((SCREEN_WIDTH / 4) + (SCREEN_WIDTH / 24)) * 0.25, (SCREEN_HEIGHT - 10) * .45, (SCREEN_WIDTH / 4) + (SCREEN_WIDTH / 24), SCREEN_HEIGHT - 10);

	// Cube c
	glColor4f(1.0, 0.0, 0.0, 0.3);
	glRectf(((SCREEN_WIDTH / 4) + (SCREEN_WIDTH / 24)) * 0.25, MAR, (SCREEN_WIDTH / 4) + (SCREEN_WIDTH / 24), (SCREEN_HEIGHT - 10) * .45);

	// Cube d
	glColor4f(0.0, 0.0, 1.0, 0.3);
	glRectf(((SCREEN_WIDTH) / 3) + MAR, MAR,
		((((SCREEN_WIDTH) / 3) + MAR) + ((((7 * SCREEN_WIDTH) / 12) + (SCREEN_WIDTH / 24)) - (((SCREEN_WIDTH) / 3) + MAR)) * .15),
		SCREEN_HEIGHT - 10);

	// Cube e
	glColor4f(1.0, 0.0, 0.0, 0.3);
	glRectf(((((SCREEN_WIDTH) / 3) + MAR) + ((((7 * SCREEN_WIDTH) / 12) + (SCREEN_WIDTH / 24)) - (((SCREEN_WIDTH) / 3) + MAR)) * .15),
		(SCREEN_HEIGHT - 10) * .45, ((7 * SCREEN_WIDTH) / 12) + (SCREEN_WIDTH / 24), SCREEN_HEIGHT - 10);


	// Cube f
	glColor4f(0.0, 0.0, 0.0, 0.3);
	glRectf(((((SCREEN_WIDTH) / 3) + MAR) + ((((7 * SCREEN_WIDTH) / 12) + (SCREEN_WIDTH / 24)) - (((SCREEN_WIDTH) / 3) + MAR)) * .15),
		MAR, ((7 * SCREEN_WIDTH) / 12) + (SCREEN_WIDTH / 24), (SCREEN_HEIGHT - 10) * .45);

	// Cube g
	glColor4f(1.0, 0.0, 0.0, 0.3);
	glRectf(((2 * SCREEN_WIDTH) / 3) + MAR, SCREEN_HEIGHT * 0.6, ((11 * SCREEN_WIDTH) / 12) + (SCREEN_WIDTH / 24), SCREEN_HEIGHT - 10);

	// Cube h
	glColor4f(0.0, 0.0, 1.0, 0.3);
	glRectf(((2 * SCREEN_WIDTH) / 3) + MAR, MAR,
		((((2 * SCREEN_WIDTH) / 3) + MAR) + ((((11 * SCREEN_WIDTH) / 12) + (SCREEN_WIDTH / 24)) - (((2 * SCREEN_WIDTH) / 3) + MAR)) * .35),
		SCREEN_HEIGHT * 0.6);

	// Cube i
	glColor4f(1.0, 0.0, 0.0, 0.3);
	glRectf(((((2 * SCREEN_WIDTH) / 3) + MAR) + ((((11 * SCREEN_WIDTH) / 12) + (SCREEN_WIDTH / 24)) - (((2 * SCREEN_WIDTH) / 3) + MAR)) * .35),
		MAR, ((11 * SCREEN_WIDTH) / 12) + (SCREEN_WIDTH / 24), (SCREEN_HEIGHT) * 0.6);

	// Disable blending and resume depth mask
	glDepthMask(GL_TRUE);
	glDisable(GL_BLEND);
}

/*
drawLocatedGlyphs
This is the driver function for displaying a set of
three located SPC-SF glyphs in the three paired
coordinate axes of the visualization (PC-SPC-SF)
@param		normalData	Data point to be visualized
			classify	Class of the data point
@return					void
*/
void drawLocatedGlyphs(std::vector<GLfloat>* normalData, bool classify)
{
	glPushMatrix();
	// Bare nucleus
	GLfloat bn = *(normalData->begin() + 6);

	// Uniformity of cell size
	GLfloat uc = *(normalData->begin() + 2);

	// Clump thickness
	GLfloat cl = *(normalData->begin() + 1);

	// Bland chromatin
	GLfloat bc = *(normalData->begin() + 7);

	// Marginal adhesion
	GLfloat mg = *(normalData->begin() + 4);

	std::vector<GLfloat> axesSPC{};
	
		axesSPC.push_back(*((normalData)->begin() + 3));	// X1
		axesSPC.push_back(*((normalData)->begin() + 2));	// Y1
		axesSPC.push_back(*((normalData)->begin() + 3));	// X2
		axesSPC.push_back(*((normalData)->begin() + 4));	// Y2
		axesSPC.push_back(*((normalData)->begin() + 8));	// X3
		axesSPC.push_back(*((normalData)->begin() + 7));	// Y3

		// Encode angles with most meaningful attributes
		std::vector<GLfloat> stickFig{};
		// Populate: CL (angle), UC (length), BN (angle), BC (length)
		stickFig.push_back(*((normalData)->begin()));			// Angle 1
		stickFig.push_back(*((normalData)->begin() + 1));	// Length 1
		stickFig.push_back(*((normalData)->begin() + 5));	// Angle 2
		stickFig.push_back(*((normalData)->begin() + 6));	// Length 2

	// Arrange greater SPC position attributes from optimal positioning described in Worland, Wagle, and Kovalerchuk
	std::vector<GLfloat> position{ uc, bn, bc, cl, bn, mg };
	std::vector<GLfloat>::iterator positionIt = position.begin();

	// Glyph 1
	std::vector<GLfloat>::iterator axes1 = normalData->begin();			// First six attributes used for SPC (6/10)
	std::vector<GLfloat>::iterator stick1 = normalData->begin() + 6;	// Last four attributes used for SF  (4/10)

	// Glyph 2
	std::vector<GLfloat>::iterator axes2 = normalData->begin();			// First six attributes used for SPC (6/10)
	std::vector<GLfloat>::iterator stick2 = normalData->begin() + 6;	// Last four attributes used for SF  (4/10)

	// Glyph 3
	std::vector<GLfloat>::iterator axes3 = normalData->begin();			// First six attributes used for SPC (6/10)
	std::vector<GLfloat>::iterator stick3 = normalData->begin() + 6;	// Last four attributes used for SF  (4/10)

	// *********************** DRAW SF GLYPHS ***********************
	// Note that we are currently converting 'radians to degrees'
	// Draw Stick Figure (save vertex positions)

	// Construct glyph tool
	SpcSfGlyph glyph = SpcSfGlyph();
	// Construct turtle tool
	TurtleG* turt = new TurtleG();

	Point2* pos2 = new Point2();
	Point2* pos3 = new Point2();

	// Locate the lower-left corner of viewport for glyph drawing
	GLfloat x1 = ((SCREEN_WIDTH / WIDTH_SCALE) * *positionIt) - ((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2));
	GLfloat y1 = ((SCREEN_HEIGHT - (SCREEN_HEIGHT / HEIGHT_SCALE)) * *++positionIt) - ((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2));

	glPushMatrix();

	// Place the viewport
	glViewport(
		x1,
		y1,
		SCREEN_WIDTH / VIEWPORT_SCALE,
		SCREEN_WIDTH / VIEWPORT_SCALE
	);

	//Glyph1
	glLineWidth(4.0);
	glyph.drawGlyphSF(pos2, pos3, stickFig.begin(), classify, *turt, DYNAMIC_ANGLES, POS_ANGLE, GLYPH_SCALE_FACTOR, 0.1, 2.0, ANGLE_FOCUS, BIRD_FOCUS);
	// RESET THE CP AND CD
	turt->setCP(0.0, 0.0);
	// *********************** DRAW SPC AXES ***********************
	if (DRAW_AXES)
	{
		glLineWidth(2.0);
		glyph.drawAxesSPC(*pos2, *pos3, axesSPC.begin(), AXIS_LENGTH, GLYPH_SCALE_FACTOR, DOTTED_AXES);
	}


	// Locate the viewport for glyph drawing
	GLfloat x2 = ((SCREEN_WIDTH / WIDTH_SCALE) * *++positionIt) + (SCREEN_WIDTH / 3) - ((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2));
	GLfloat y2 = ((SCREEN_HEIGHT - (SCREEN_HEIGHT / HEIGHT_SCALE)) * *++positionIt) - ((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2));

	// Locate the viewport for glyph drawing
	glViewport(
		x2,
		y2,
		SCREEN_WIDTH / VIEWPORT_SCALE,
		SCREEN_WIDTH / VIEWPORT_SCALE
	);

	//Glyph2
	glLineWidth(4.0);
	glyph.drawGlyphSF(pos2, pos3, stickFig.begin(), classify, *turt, DYNAMIC_ANGLES, POS_ANGLE, 0.25, 0.1, 2.0, ANGLE_FOCUS, BIRD_FOCUS);
	// RESET THE CP AND CD
	turt->setCP(0.0, 0.0);

	// *********************** DRAW SPC AXES ***********************
	if (DRAW_AXES)
	{
		glLineWidth(2.0);
		glyph.drawAxesSPC(*pos2, *pos3, axesSPC.begin(), AXIS_LENGTH, 0.25, DOTTED_AXES);
	}

	// Locate the viewport for glyph drawing
	GLfloat x3 = ((SCREEN_WIDTH / WIDTH_SCALE) * *++positionIt) + (((2 * SCREEN_WIDTH) / 3)) - ((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2));
	GLfloat y3 = ((SCREEN_HEIGHT - (SCREEN_HEIGHT / HEIGHT_SCALE)) * *++positionIt) - ((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2));

	// Locate the viewport for glyph drawing
	glViewport(
		x3,
		y3,
		SCREEN_WIDTH / VIEWPORT_SCALE,
		SCREEN_WIDTH / VIEWPORT_SCALE
	);

	//Glyph3
	glLineWidth(4.0);
	glyph.drawGlyphSF(pos2, pos3, stickFig.begin(), classify, *turt, DYNAMIC_ANGLES, POS_ANGLE, 0.25, 0.1, 2.0, ANGLE_FOCUS, BIRD_FOCUS);
	// RESET THE CP AND CD
	turt->setCP(0.0, 0.0);
	// *********************** DRAW SPC AXES ***********************
	if (DRAW_AXES)
	{
		glLineWidth(2.0);
		glyph.drawAxesSPC(*pos2, *pos3, axesSPC.begin(), AXIS_LENGTH, 0.25, DOTTED_AXES);
	}

	glPopMatrix();
	glViewport(0.0, 0.0, SCREEN_WIDTH, SCREEN_HEIGHT);
	gluOrtho2D(0.0, SCREEN_WIDTH, 0.0, SCREEN_HEIGHT);
	glColor4f(0.0, 0.0, 0.0, 7.0);

	// Draw SPC axis frame * 3

	drawArrow(Point2(MAR, MAR), Point2(MAR, SCREEN_HEIGHT - 10), SCREEN_HEIGHT);
	drawArrow(Point2(MAR, MAR), Point2((SCREEN_WIDTH / 4) + (SCREEN_WIDTH / 24), MAR), SCREEN_HEIGHT);

	drawArrow(Point2(((SCREEN_WIDTH) / 3) + MAR, MAR), Point2((SCREEN_WIDTH / 3) + MAR, SCREEN_HEIGHT - 10), SCREEN_HEIGHT);
	drawArrow(Point2((SCREEN_WIDTH / 3) + MAR, MAR), Point2(((7 * SCREEN_WIDTH) / 12) + (SCREEN_WIDTH / 24), MAR), SCREEN_HEIGHT);

	drawArrow(Point2(((2 * SCREEN_WIDTH) / 3) + MAR, MAR), Point2(((2 * SCREEN_WIDTH) / 3) + MAR, SCREEN_HEIGHT - 10), SCREEN_HEIGHT);
	drawArrow(Point2(((2 * SCREEN_WIDTH) / 3) + MAR, MAR), Point2(((11 * SCREEN_WIDTH) / 12) + (SCREEN_WIDTH / 24), MAR), SCREEN_HEIGHT);

	// Draw edges between glyphs
	if (DRAW_EDGES)
	{	// Determine which class/color the edge belongs to
		if (classify)	glColor4f(0.0, 0.0, 1.0, 7.0);
		else			glColor4f(1.0, 0.0, 0.0, 7.0);
		glLineWidth(1.0);
		glPushMatrix();	// Push new modelview matrix for translation
		glTranslatef(-((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2)), -((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2)), 0);

		// Draw first edge
		drawArrow(Point2(x1 + (SCREEN_WIDTH / VIEWPORT_SCALE), y1 + (SCREEN_WIDTH / VIEWPORT_SCALE)),
			Point2(x2 + (SCREEN_WIDTH / VIEWPORT_SCALE), y2 + (SCREEN_WIDTH / VIEWPORT_SCALE)), SCREEN_HEIGHT);

		// Draw second edge
		drawArrow(Point2(x2 + (SCREEN_WIDTH / VIEWPORT_SCALE), y2 + (SCREEN_WIDTH / VIEWPORT_SCALE)),
			Point2(x3 + (SCREEN_WIDTH / VIEWPORT_SCALE), y3 + (SCREEN_WIDTH / VIEWPORT_SCALE)), SCREEN_HEIGHT);

		glPopMatrix(); // Pop modelview matrix used for translation
	}

	// Reset modelview matrix and flush buffer
	glPopMatrix();
	glFlush();
}

/*	hammingDistance
Return the Hamming distance between two vectors,
by taking the sqrt of the sum of the squares of the differences
using matching attributes of each vector.
@param		vec1	first vector to compare
			vec2	second vector to compare

@return		the Hamming distance between the two vectors
*/
int hammingDistance(std::vector<float>* vec1, std::vector<float>* vec2)
{
	int hammingSum = 0;	// initialize sum variable
	// initialize iterator to vec2 begin
	std::vector<float>::iterator it2 = vec2->begin();

	// for each attributes in vec1
	for (auto& attr : *vec1)
	{
		// if no match, increment Hamming sum
		if (((int)(attr * 10)) != ((int)(*it2 * 10)))
		{
			++hammingSum;
		}
		++it2;	// increment iterator
	}

	return hammingSum;
}
/*	euclideanDistance
Return the Euclidean distance between two vectors,
by taking the sqrt of the sum of the squares of the differences
using matching attributes of each vector.
@param		vec1	first vector to compare
			vec2	second vector to compare

@return		the Euclidean distance between the two vectors
*/
float euclideanDistance(std::vector<float>* vec1, std::vector<float>* vec2)
{
	float sum = 0.0;	// initialize sum variable
	// initialize iterator to vec2 begin
	std::vector<float>::iterator it2 = vec2->begin();

	// for each attributes in vec1
	for (auto& attr : *vec1)
	{
		// Sum the squares of the differences of the matching attributes between vectors
		sum += pow(attr - *it2, 2);
		++it2;	// increment iterator
	}

	// return sqrt of sum
	return sqrt(sum);
}

/*
isClose
This function is a helper function for the isClose function.
It is used to determine whether two given data points are within
the given threshold of each other, and thus are in the same
neighborhood.
@param		vec1		the first data point to be compared
			vec2		the second data point to be compared
@return					boolean value, is vec2 within threshold
						of vec1
*/
bool isClose(std::vector<float>* vec1, std::vector<float>* vec2)
{
	std::vector<float>::iterator it1 = vec1->begin();
	std::vector<float>::iterator it2 = vec2->begin();

	// Check each pair of elements for difference within threshold value
	for (it1; it1 != vec1->end(); ++it1)
	{
		if (abs(*it1 - *it2) > THRESHOLD_VALUE)
			return false;
		it2++;
	}

	// Return true if all attribtes are within threshold
	return true;
}

/*
computeAllDistances
This function is used to determine which data points in the
data set are within the threshold of the data point identified
by DATA_INDEX
@param		curr		data point at DATA_INDEX
			data		pointer to data set
@return					boolean vector identifying threshold points
*/
std::vector<bool> computeAllDistances(std::vector<GLfloat>* curr, std::vector<std::vector<GLfloat>>* data)
{
	// This variable is for temp extracted data to compare
	//std::vector<float> data0();
	std::vector<bool> close{};

	// Populate boolean vector to determine
	// which data points are within threshold of current point
	for (std::vector<std::vector<GLfloat>>::iterator iter = data->begin(); iter != data->end(); ++iter)
	{
		close.push_back(isClose(curr, &(*iter)));
	}

	// Close file streams
	//output.close();

	return close;
}

/*
loadConfig
This function loads global parameters from the file
config.config contained in the current directory.
@param			none
@return			void
*/
void loadConfig()
{
	std::string line = "";
	std::ifstream myFile("config.config");

	if (!myFile)	// Check that file was opened successfully
	{
		std::cout << "Error: Config File did not open.";
	}

	// Check for incorrect input in config file
	//

	// Set Parameters
	getline(myFile, line);
	// AXIS_LENGTH
	std::istringstream sin1(line.substr(line.find("=") + 1));
	sin1 >> AXIS_LENGTH;

	getline(myFile, line);
	// SCREEN_WIDTH
	std::istringstream sin2(line.substr(line.find("=") + 1));
	sin2 >> SCREEN_WIDTH;

	getline(myFile, line);
	// SCREEN_HEIGHT
	std::istringstream sin3(line.substr(line.find("=") + 1));
	sin3 >> SCREEN_HEIGHT;

	/*
	getline(myFile, line);
	// GLYPH_SIZE
	std::istringstream sin4(line.substr(line.find("=") + 1));
	sin4 >> GLYPH_SCALE_FACTOR;

	getline(myFile, line);
	// NUM_ROWS
	std::istringstream sin5(line.substr(line.find("=") + 1));
	sin5 >> NUM_ROWS;

	getline(myFile, line);
	// NUM_COLUMNS
	std::istringstream sin6(line.substr(line.find("=") + 1));
	sin6 >> NUM_COLUMNS;
	*/
}

/*
myIdle
This is the OpenGL idle callback function.
@param			none
@return			void
*/
void myIdle() {}

/*
mouse_button_callback
This is the OpenGL mouse event callback function,
used to register responses to mouse events.
@param		button	mouse button pressed
			state	state of pressed button
			x		x coord of mouse cursor
			y		y coord of mouse cursor
@return				void
*/
void mouse_button_callback(int button, int state, int x, int y)
{
	// Left mouse button decrements the current data index
	if (DATA_INDEX > 0 && state == GLUT_DOWN && button == GLUT_LEFT_BUTTON)
	{
		DATA_INDEX -= 1;
		glutPostRedisplay();
	}
	// Right mouse button increments the current data index
	else if (DATA_INDEX < DATA_SIZE && state == GLUT_DOWN && button == GLUT_RIGHT_BUTTON)
	{
		DATA_INDEX += 1;
		glutPostRedisplay();
	}
}

/*
keyboard_special
This is the OpenGL special keyboard callback function,
used to register callbacks to special key events.
(i.e. arrow keys, function keys)
@param		key		keyboard key pressed to initiate callback
			x		x coord of mouse cursor
			y		y coord of mouse cursor
@return				void
*/
void keyboard_special(int key, int x, int y)
{
	// Move backwards in the list of data points
	if (DATA_INDEX > 0 && key == GLUT_KEY_LEFT)
	{
		glClearColor(1.0, 1.0, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT);

		DATA_INDEX -= 1;
	}
	// Right arrow key increments the current data index
	else if (DATA_INDEX < DATA_SIZE && key == GLUT_KEY_RIGHT)
	{
		glClearColor(1.0, 1.0, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT);

		DATA_INDEX += 1;
	}

	// F1 key response:
	// toggle display all / cycle individual
	else if (key == GLUT_KEY_F1)
	{
		DISPLAY_ALL = !DISPLAY_ALL;
	}

	// F2 key response:
	// toggle which class is drawn first
	else if (key == GLUT_KEY_F2)
	{
		CLASS_SEPERATION_MODE = !CLASS_SEPERATION_MODE;
	}

	// F3 key response:
	// toggle hypercube overlay
	else if (key == GLUT_KEY_F3)
	{
		DISPLAY_HYPERCUBES = !DISPLAY_HYPERCUBES;
	}

	// F4 key response:
	// Switch from SPC-SF to tiled glyphs
	else if (key == GLUT_KEY_F4)
	{
		DISPLAY_SELECTOR = !DISPLAY_SELECTOR;
	}

	// F5 key response:
	// Switch between dynamic / static SF angles
	else if (key == GLUT_KEY_F5)
	{
		DYNAMIC_ANGLES = !DYNAMIC_ANGLES;
	}

	// F6 key response:
	// Switch between positive / negative first SF angle
	else if (key == GLUT_KEY_F6)
	{
		POS_ANGLE = !POS_ANGLE;
	}

	// F7 key response:
	// Switch between grid of glyphs / PC view
	else if (key == GLUT_KEY_F7)
	{
		PC_OFF = !PC_OFF;
	}

	// F8 key response:
	// Switch to grid of representative glyphs
	else if (key == GLUT_KEY_F8)
	{
		REPS_OFF = !REPS_OFF;
	}

	// F9 key response:
	// Toggle dashed / dotted axes
	else if (key == GLUT_KEY_F9)
	{
		DOTTED_AXES = !DOTTED_AXES;
	}

	// F10 key response:
	// Toggle angle / length focus
	else if (key == GLUT_KEY_F10)
	{
		ANGLE_FOCUS = !ANGLE_FOCUS;
	}

	// F11 key response:
	// Toggle grey / color birds
	else if (key == GLUT_KEY_F11)
	{
		BIRD_FOCUS = !BIRD_FOCUS;
	}

	// Up arrow response:
	// toggle drawing edges between glyphs
	else if (key == GLUT_KEY_UP)
	{
		DRAW_EDGES = !DRAW_EDGES;
	}

	// Down arrow response:
	// toggle drawing SPC axes for glyphs
	else if (key == GLUT_KEY_DOWN)
	{
		DRAW_AXES = !DRAW_AXES;
	}

	// Redisplay with updated parameters
	glutPostRedisplay();
}

/*
myKeyboard
This is the OpenGL keyboard callback function.
@param			none
@return			void
*/
void myKeyboard(unsigned char key, int x, int y)
{


	// Redisplay with updated DATA_INDEX
	//glutPostRedisplay();
	//glutIdleFunc(myIdle);
}

// Draw tile grid, set to size determined by
// NUM_ROWS and NUM_COLUMNS
void drawGrid()
{
	glClearColor(1.0, 1.0, 1.0, 0.0);	// Clear window to white
	glClearDepth(1.0f);					// Clear depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);	// Load GL_PROJECTION matrix
	glLoadIdentity();				// Load identity matrix
	glPushMatrix();					// Push new matrix

	// Set orthographic projection to SCREEN_WIDTH x SCREEN_HEIGHT
	gluPerspective(0, float(SCREEN_WIDTH) / float(SCREEN_HEIGHT), 0.1, 100.0);
	glViewport(0.0, 0.0, SCREEN_WIDTH, SCREEN_HEIGHT);
	gluOrtho2D(0.0, SCREEN_WIDTH, 0.0, SCREEN_HEIGHT);

	// ************** DRAW TILE FRAME (NUM_ROWS x NUM_COLUMNS) ****************
	glColor4f(0.0, 0.0, 0.0, 1.0);				// Draw in black
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);	// Draw outlines

	// Draw the tile frame to contain the glyphs
	glLineWidth(2.0);		// Set line width to 2
	glPushMatrix();			// Push new matrix
	for (unsigned int rowNum = 0; rowNum < NUM_ROWS; ++rowNum)			// Rows
	{
		for (unsigned int colNum = 0; colNum < NUM_COLUMNS; ++colNum)	// Columns
		{
			glRectf(		// Draw rectangles for grid of glyphs
				0.0 + GRID_MARGIN,
				SCREEN_HEIGHT - (SCREEN_HEIGHT / NUM_ROWS) * rowNum,
				((SCREEN_WIDTH / NUM_COLUMNS) * (colNum + 1)) + GRID_MARGIN,
				SCREEN_HEIGHT - (SCREEN_HEIGHT / NUM_ROWS) * (rowNum + 1)
			);
		}
	}
	glPopMatrix();		// Pop a matrix from gl matrix stack
	glPopMatrix();		// Pop a matrix from gl matrix stack
}

/*
getIdealGLyphs
This function computes an ideal glyph for each class,
using the class average for each attribute,
and adds the glyphs to the rep glyph vector.
@param		allData		copy of data set
@return		classify	data set class labels
*/
void getIdealGlyphs(std::vector<std::vector<GLfloat>>* all_Data,
	std::vector<bool>* classify)
{
	// Make deep copy of all data
	std::vector<std::vector<GLfloat>> allData(*all_Data);
	std::vector<bool> classifyCopy(*classify);

	// Set iterators to copies
	std::vector<std::vector<GLfloat>>::iterator dataIt = allData.begin();
	std::vector<bool>::iterator classIt = classifyCopy.begin();

	// Initialize temp class vectors
	std::vector<std::vector<GLfloat>> malSet{};
	std::vector<std::vector<GLfloat>> benSet{};

	// Iterate through all data
	for (dataIt; dataIt != allData.end(); ++dataIt)
	{
		if (*classIt)
		{	// add point to class 1
			benSet.push_back(*dataIt);
		}
		else
		{	// add point to class 2
			malSet.push_back(*dataIt);
		}
		// increment class iterator
		++classIt;
	}

	// Initialize temporary benign vec
	std::vector<GLfloat> tempBen{};

	// Compute average points for benign class
	for (unsigned int index = 0; index < MAX_SIG_INDEX; ++index)
	{	// For the current index of each vector
		GLfloat tempAttr = 0.0;
		for (auto& vec : benSet)
		{	// For each vector
			tempAttr += vec[index];
		}
		// Compute average value from sum
		tempAttr /= benSet.size();
		// Add to representative point
		tempBen.push_back(tempAttr);
	}

	// Initialize temporary malig vec
	std::vector<GLfloat> tempMal{};
	// Compute average points for malignant class
	for (unsigned int index = 0; index < MAX_SIG_INDEX; ++index)
	{	// For the current index of each vector
		GLfloat tempAttr = 0.0;
		for (auto& vec : malSet)
		{	// For each vector
			tempAttr += vec[index];
		}
		// Compute average value from sum
		tempAttr /= malSet.size();
		// Add to representative point
		tempMal.push_back(tempAttr);
	}

	// Push ideal glyphs to rep vector
	reps.push_back(tempBen);
	reps.push_back(tempMal);

	// Push ideal glyph labels to label vector
	labels.push_back("ideal B");
	labels.push_back("ideal M");

	// Push ideal glyph class to class vector
	repsClass.push_back(true);
	repsClass.push_back(false);

	// Set condition for ideal glyphs collected
	IDEAL_COLLECTED = true;
}

// Retrieve vector of representative glyphs for each neighborhood
void getRepresentativeGlyphs(std::vector<std::vector<GLfloat>>* all_Data,
	std::vector<bool>* classify)
{
	unsigned int dataIndex = 0;

	// Make deep copy of all data
	std::vector<std::vector<GLfloat>> allData(*all_Data);
	std::vector<bool> classifyCopy(*classify);
	std::vector<bool> thresholds{};

	bool addThis = false;

	int benCount = 0;
	int malCount = 0;
	int hoodCount = 0;
	// Loop through data, adding one point from each neighborhood to the REPS array
	// Until all points have been processed
	while (allData.size() != 0)
	{
		// Compute points within threshold
		thresholds = computeAllDistances(&allData[dataIndex], &allData);

		// Initialize temporary vector for calculating average glyph
		std::vector< std::vector<GLfloat>> tempData{};

		// Initialize vector iterators
		std::vector<std::vector<GLfloat>>::iterator dataIt = allData.begin();
		std::vector<bool>::iterator classIt = classifyCopy.begin();
		std::vector<bool>::iterator threshIt = thresholds.begin();

		// Save first point in cluster for comparison
		std::vector<float> currPoint = (*dataIt);
		//std::vector<float>::iterator it1End = (*dataIt).end();
		/*
		if (hoodCount == 1)
		{	// Add rep glyph to mixed view
			mixedHood.push_back(*dataIt);
			mixedClass.push_back(false);
		}
		*/
		addThis = false;
		// Use threshold values to remove points
		while (dataIt != allData.end())
		{
			addThis = false;
			// If current point is in current neighborhood
			if (*threshIt)
			{
				// Add class to class count
				if (*classIt)
				{
					++benCount;
				}
				else if (!*classIt)
				{
					++malCount;
				}

				// If currently analyzing chosen focus cluster,
				// and current point is off the same class as the center of the cluster.
				if (hoodCount == CLUSTER && (*classIt || !*classIt))
				{	// Save the point and its class for analysis
					// extract outlying benign points hood[0]
					std::vector<float>::iterator it2 = (*dataIt).begin();

					int difCount = 0;	// track number of divergent attributes in current point
					// Check each pair of elements from each vector for difference within threshold value
					for (std::vector<float>::iterator it1 = currPoint.begin(); it1 != currPoint.end(); ++it1)
					{	
						// If the absolute value of the difference of the pair of attributes
						// is below the allowed minimum threshold value
						if (abs(*it1 - *it2) >= MIN_THRESHOLD)
						{
							++difCount;	// Increment the difference count
							if (difCount > ALLOWED_DIFFERENCES)
							{	
								// If the number of allowed differences has been exceeded
								addThis = true;	// flag this data point as divergent from the center
							}
						}
						// increment comparison vector iterator
						it2++;
					}
					
					// addThis = true;
					// save for debugging
					
					// if flag was set to add this point
					if (addThis)
					{
						mixedHood.push_back(*dataIt);	// add data point to analysis data vector
						mixedClass.push_back(*classIt);	// add class to analysis class vector
					}
				}

				// Add point to temp vector
				tempData.push_back(*dataIt);

				// Remove current point from vectors
				dataIt = allData.erase(dataIt);
				classIt = classifyCopy.erase(classIt);
				threshIt = thresholds.erase(threshIt);
			}
			else
			{	// Increment iterators
				++dataIt;
				++classIt;
				++threshIt;
			}

			if (dataIt == allData.end())
			{	// Exit loop
				break;
			}
		}

		// identify count of dominant class of cluster
		float domCount = std::max(benCount, malCount);
		// itentify total count of datapoints in cluster
		float totalCount = benCount + malCount;

		// calculate percent purity of cluster
		float percentPure = (domCount / totalCount) * 100.0;

		// round label to three decimal places
		std::string pureLabel1 = std::to_string(percentPure);
		std::string pureLabel = pureLabel1.substr(0, 5);

		// Add label to vector
		labels.push_back(pureLabel + " pure");
		 

		//labels.push_back(std::to_string(benCount) + " ben., " + std::to_string(malCount) + " mal.");
		repsClass.push_back(benCount > malCount);

		// Initialize new rep vector
		std::vector<GLfloat> repVec{};

		/* Compute average glyph of neighborhood by taking tempData,
		*  the set of vectors saved during clustering, and computing
		*  the average value for each attribute from the data points
		*  in the set.
		*/
		for (unsigned int index = 0; index < MAX_SIG_INDEX; ++index)
		{	// For the current index of each vector
			GLfloat tempAttr = 0.0;
			for (auto& vec : tempData)
			{	// For each vector
				tempAttr += vec[index];
			}
			// Compute average value from sum
			tempAttr /= tempData.size();
			// Add to representative point
			repVec.push_back(tempAttr);
		}

		// Add representative vector to set
		reps.push_back(repVec);

		// Reset counters
		benCount = 0;
		malCount = 0;
		++hoodCount;
	}

	// compute distance labels within cluster
	std::vector<GLfloat>	vecEDist{};
	std::vector<int>		vecHDist{};
	// index here determines new center of cluster
	std::vector<GLfloat>	newCenterPoint = mixedHood[CLUSTER_CENTER_FOCUS];

	// Compute Hamming distance between all points and chosen focus point
	for (std::vector<std::vector<GLfloat>>::iterator iter = mixedHood.begin(); iter != mixedHood.end(); ++iter)
	{
		// compute Hamming distance
		vecHDist.push_back(hammingDistance(&newCenterPoint, &(*iter)));
	}

	// Normalize data in focus cluster vector
	// Normalize data in mixed neighborhood vector vector
	int i = 0;	// loop counter
	for (auto& vec : mixedHood)
	{	// For each data point in the mixed neighborhood
		const unsigned int scaleFactor = 10;
		std::vector<float> normalData;		// Normalize data to [0, 1]
		for (std::vector<float>::iterator iter = vec.begin(); iter < vec.end(); iter++)
		{	// Add normalized data point to copy vector
			normalData.push_back(*iter / scaleFactor);
		}
		// Copy vector to reprentative glyph array
		mixedHood[i] = normalData;
		++i;
	}

	// reset center point with normalized values
	newCenterPoint = mixedHood[CLUSTER_CENTER_FOCUS];

	// Compute Euclidean distance between all points and chosen focus point
	for (std::vector<std::vector<GLfloat>>::iterator iter = mixedHood.begin(); iter != mixedHood.end(); ++iter)
	{
		// compute Euclidean distance
		vecEDist.push_back(euclideanDistance(&newCenterPoint, &(*iter)));
	}

	// collect cluster distance labels
	std::vector<int>::iterator iterH = vecHDist.begin();	// initialize iterator to distances
	for (std::vector<GLfloat>::iterator iterE = vecEDist.begin(); iterE != vecEDist.end(); ++iterE)
	{	// loop through distance vectors, creating label
		std::string startLabel = "E:";	// Euclidean label
		std::stringstream stream;		// initialize stringstream

		// Round euc distance float to 2 decimal places
		// TAKE A LOOK AT THIS!!!!!! MIGHT BE CHANGING VALUES
		stream << std::fixed << std::setprecision(3) << *iterH;
		std::string eucLabel1 = std::to_string(*iterE);
		std::string eucLabel = eucLabel1.substr(0, 4);

		// concatenate euclidean and hamming distances into single label
		std::string tempLabel = startLabel + eucLabel + " H:" + std::to_string((*iterH));
		distanceLabels.push_back(tempLabel);	// add label to distance label vector

		// increment iterator
		++iterH;
	}

	// Insert cluster representative glyph
	//mixedHood.insert(mixedHood.begin(), reps[CLUSTER + 2]);
	// Insert clas of representative glyph
	//mixedClass.insert(mixedClass.begin(), repsClass[CLUSTER + 2]);

	
	// Normalize data in REPS vector
	i = 0;
	for (auto& vec : reps)
	{
		const unsigned int scaleFactor = 10;
		std::vector<float> normalData;		// Normalize data tDo [0, 1]
		for (std::vector<float>::iterator iter = vec.begin(); iter < vec.end(); iter++)
		{	// Add nomrlized data point to copy vector
			normalData.push_back(*iter / scaleFactor);
		}
		// Copy vector to reprentative glyph array
		reps[i] = normalData;
		++i;
	}

	// Set condition for representative glyphs collected
	REPS_COLLECTED = true;
	// Global representative glyph vector is now initialized
}

// create hyperblocks using MHyper algorithm
void mergerHyperblocks()
{

}


// initialize OpenGL state
void openGLInit()
{
	glEnable(GL_DEPTH_TEST);	// Enable depth testing
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //Set blending function.
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(0, float(SCREEN_WIDTH) / float(SCREEN_HEIGHT), 0.1, 100.0);
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClearDepth(1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

// IMPORT DATA
void importData(std::vector<std::vector<GLfloat>>* allData, std::vector<bool>* classify)
{
	// Read data file
	std::string line = "";
	std::ifstream myFile("breast-cancer-wisconsin.DATA");

	for (unsigned int i = 0; i < DATA_SIZE; i++)
	{
		getline(myFile, line);
		int vecCount = 10;

		// Split string into a vector
		std::vector<int> dataInt;
		std::stringstream ss(line);

		while (ss.good()) {		// Replace "?" in data
			std::string substr = "";
			getline(ss, substr, ',');

			if (substr == "?") substr = std::to_string(*dataInt.begin());
			dataInt.push_back(stoi(substr));
			continue;
		}

		std::vector<GLfloat> data(dataInt.begin(), dataInt.end());

		// Set color determined by class
		if (*--(data.end()) == 4) classify->push_back(false);
		else					  classify->push_back(true);

		// Remove labels from data
		data.erase(data.begin());
		data.erase(--data.end());

		// Duplicate data attributes x 4 to move from 9-D to 36-D
		//data.insert(data.end(), data.begin(), data.end());
		//data.insert(data.end(), data.begin(), data.end());

		(*allData)[i] = data;
	}
	myFile.close();
}

// IMPORT SEED DATA
void importSeedData(std::vector<std::vector<GLfloat>>* allData)
{
	std::string line = "";
	std::ifstream myFile("seeds_dataset.txt");

	// For each line in the hyperblock data file
	for (unsigned int i = 0; i < SEED_DATASET_SIZE; i++)
	{
		getline(myFile, line);
		int vecCount = SEED_DATA_SIZE;

		// Split string into a vector
		std::vector<int> dataFloat;
		std::stringstream ss(line);

		while (ss.good()) {
			std::string substr = "";
			getline(ss, substr, '\t');	// extract line from file
			dataFloat.push_back(stof(substr));	// convert to float and push to vector
			continue;
		}

		std::vector<GLfloat> data(dataFloat.begin(), dataFloat.end());
		(*allData)[i] = data;	// add data point to hyperblock vector
	}
		myFile.close();
}


// IMPORT LINCOLN'S HYPERBLOCK DATA FROM CSV FILES
void importHyperblockData(std::vector<std::vector<GLfloat>>* allData)
{
	// Read data file
	std::string line = "";
	std::ifstream myFile("HB1.csv");

	// For each line in the hyperblock data file
	for (unsigned int i = 0; i < HYPERBLOCK_SIZE; i++)
	{
		getline(myFile, line);
		int vecCount = HYPERBLOCK_DATA_SIZE;

		// Split string into a vector
		std::vector<int> dataFloat;
		std::stringstream ss(line);

		while (ss.good()) {
			std::string substr = "";
			getline(ss, substr, ',');	// extract line from file
			dataFloat.push_back(stof(substr));	// convert to float and push to vector
			continue;
		}

		std::vector<GLfloat> data(dataFloat.begin(), dataFloat.end());
		(*allData)[i] = data;	// add data point to hyperblock vector
	}
	myFile.close();		// close file stream
}

// *********************** Display SPC-SF Hybrid Visualization ***********************
// 
// Data from UCI Machine Learning Repository
// breast-cancer-wisconsin.DATA
void myDisplay()
{
	// hyperblock container
	std::vector<std::vector<std::vector<GLfloat>>> hyperblocks{};

	// seed container
	std::vector<std::vector<GLfloat>> seeds(SEED_DATASET_SIZE);
	importSeedData(&seeds);

	/************************** OpenGl Set Up ********************************/
	openGLInit();

	// Process data into a 2-D vector for ease of use
	std::vector<std::vector<GLfloat>> allData(DATA_SIZE);
	std::vector<bool> classify{};	//  vector containing class of data points

	// import data from csv file
	importData(&allData, &classify);

	//*****************************************************************
	// import HB1 (first hyperblock from Lincoln)
	std::vector<std::vector<GLfloat>> hb1(HYPERBLOCK_SIZE);
	importHyperblockData(&hb1);

	if (!IDEAL_COLLECTED)	// Collect ideal class glyphs, if not done
		getIdealGlyphs(&allData, &classify);

	if (!REPS_COLLECTED)	// Collect representative glyphs, if not already done
		getRepresentativeGlyphs(&allData, &classify);

	// Compute points within threshold
	std::vector<bool> close = computeAllDistances(&allData[DATA_INDEX], &allData);
	// Normalize seed data
	unsigned int  index = 0;
	for (auto& vec : seeds)
	{
		float seedNormalizeFactors[] = { 21.2, 17.3, 6.7, 4.1, 6.6 };

		std::vector<float> normalData;		// Normalize data to [0, 1]
		for (unsigned int i = 0; i < vec.size(); ++i)
		{
			normalData.push_back( (vec[i]) / (seedNormalizeFactors[i]));
		}
		seeds[index] = normalData;
		++index;
	}

	// Normalize data
	index = 0;
	for (auto& vec : allData)
	{
		const unsigned int scaleFactor = 10;
		std::vector<float> normalData;		// Normalize data to [0, 1]
		for (std::vector<float>::iterator iter = vec.begin(); iter < vec.end(); iter++)
		{
			normalData.push_back(*iter / scaleFactor);
		}
		allData[index] = normalData;
		++index;
	}

	// normalize hyperblock data
	index = 0;
	for (auto& vec : hb1)
	{
		const unsigned int scaleFactor = 10;
		std::vector<float> normalData;		// Normalize data to [0, 1]
		for (std::vector<float>::iterator iter = vec.begin(); iter < vec.end(); iter++)
		{
			normalData.push_back(*iter / scaleFactor);
		}
		hb1[index] = normalData;
		++index;
	}



	std::vector<bool>::iterator threshold = close.begin();		// initialize threshold iterator
	std::vector<bool>::iterator classVec = classify.begin();	// initialize class iterator

	// ************************************* DISPLAY SPC-SF GRAPH ***************************************
	if (DISPLAY_SELECTOR)
	{
		// Display points within threshold
		// Pass through twice: Draw classes sequentially depending on
		// CLASS_SEPERATION_MODE flag
		for (std::vector<std::vector<GLfloat>>::iterator it = allData.begin();
			it < allData.end(); ++it)
		{
			// Display the point if DISPLAY_ALL flag is set,
			// or if it is in the threshold of the current point
			// if (DISPLAY_ALL || *threshold == true)
			// Use exclusive or here?
			if ((!*classVec && CLASS_SEPERATION_MODE) && (DISPLAY_ALL || *threshold == true))
			{
				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();
				gluPerspective(0, float(SCREEN_WIDTH) / float(SCREEN_HEIGHT), 0.1, 100.0);

				drawLocatedGlyphs(&(*it), *classVec);
			}
			++classVec;
			++threshold;
		}

		// Reset iterators for second pass through data
		threshold = close.begin();
		classVec = classify.begin();

		// Note that doing it this way is doubling runtime....
		// So, inefficient brute force approach
		for (std::vector<std::vector<GLfloat>>::iterator it = allData.begin();
			it < allData.end(); ++it)
		{
			// Display the point if DISPLAY_ALL flag is set,
			// or if it is in the threshold of the current point
			if ((*classVec && CLASS_SEPERATION_MODE) && (DISPLAY_ALL || *threshold == true))
			{
				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();
				gluPerspective(0, float(SCREEN_WIDTH) / float(SCREEN_HEIGHT), 0.1, 100.0);

				drawLocatedGlyphs(&(*it), *classVec);
			}
			++classVec;
			++threshold;
		}

		// If flag is set, display hypercubes
		if (DISPLAY_HYPERCUBES)
		{
			displayHypercubes();
			glPopMatrix();
		}

		glPopMatrix();
		//glFlush();
	}
	// ************************************* DISPLAY GRID OF GLYPHS *************************************
	else if (!DISPLAY_SELECTOR && PC_OFF && REPS_OFF)    // Draw grid of glyphs (GRID_ROWS x GRID_COLUMNS)
	{
		drawGrid();

		// Organization of Data for Glyph (10-D, CL replicated) =
			//	{ SPC:[ (UC, CL), (BC, CL), (MG, BN) ] SF:[ UCShape, SE, NN, M ] }


		// Construct glyph tool
		SpcSfGlyph glyph = SpcSfGlyph();
		// Construct turtle tool
		TurtleG* turt = new TurtleG();
		Point2* pos2 = new Point2();
		Point2* pos3 = new Point2();

		// Initialize data iterator
		std::vector<std::vector<GLfloat>>::iterator it = allData.begin();

		glPushMatrix();		// Scale and translate glyph
		//glScalef(GLYPH_SCALE, GLYPH_SCALE, 0.0);
		//glTranslatef(GLYPH_TRANSLATE_FACTOR, GLYPH_TRANSLATE_FACTOR, 0.0);

		unsigned int rowNum = 0;	// Current Row index
		unsigned int colNum = 0;	// Current Column index

		while (rowNum < NUM_ROWS)		// Grid Rows loop
		{
			while (colNum < NUM_COLUMNS)	// Grid Columns loop
			{
				// If within threshold of current point,
				// display glyph at current grid index
				if (it != allData.end() && *threshold)
				{
					// Draw polygon outlines
					glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

					// Pointer to current data point
					std::vector<GLfloat> processedData = *it;

					// Encode angles with most meaningful attributes
					std::vector<GLfloat> stickFig{};
					// Populate: CL (angle), UC (length), BN (angle), BC (length)
					stickFig.push_back(*processedData.begin());			// Angle 1
					stickFig.push_back(*(processedData.begin() + 1));	// Length 1
					stickFig.push_back(*(processedData.begin() + 5));	// Angle 2
					stickFig.push_back(*(processedData.begin() + 6));	// Length 2

					std::vector<GLfloat> axesSPC{};
				
					// Populate: UCsh, MA, SIN, MA, NN, MIT
					// (MA, UCsh) (MA, SIN) (MIT, NN)
					axesSPC.push_back(*(processedData.begin() + 3));	// X1
					axesSPC.push_back(*(processedData.begin() + 2));	// Y1
					axesSPC.push_back(*(processedData.begin() + 3));	// X2
					axesSPC.push_back(*(processedData.begin() + 4));	// Y2
					axesSPC.push_back(*(processedData.begin() + 8));	// X3
					axesSPC.push_back(*(processedData.begin() + 7));	// Y3
					
					GLfloat maxAtr = 0.0;	// Initialize max attribute variable
					// Check for max shift
					for (std::vector<GLfloat>::iterator maxValIt = (*it).begin();
						maxValIt < (*it).begin() + 6; ++maxValIt)
					{
						maxAtr = std::max(maxAtr, *maxValIt);
					}

					glViewport(		// (rowNum x colNum)
						// Encode shift based off of first SPC axis horizontal/vertical shift.
						// Shift calculation:
							// the LARGER the value, the SMALLER the shift
							// maybe just subtract attribute value from (SCREEN_WIDTH / NUM_COLUMNS) / 2)
						((SCREEN_WIDTH / NUM_COLUMNS) * (colNum)),// - (maxAtr * ((SCREEN_WIDTH / NUM_COLUMNS) / 2)),
						((SCREEN_HEIGHT / NUM_ROWS) * (rowNum)),// - (maxAtr * ((SCREEN_HEIGHT / NUM_ROWS) / 2)),
						(SCREEN_WIDTH / NUM_COLUMNS),// * (1.0 + (maxAtr / 2)),
						(SCREEN_HEIGHT / NUM_ROWS)// * (1.0 + (maxAtr / 2))
					);

					glPushMatrix();		// Push new matrix
					glMatrixMode(GL_PROJECTION);
					glLineWidth(4.0);	// Line width = 4.0
					// Translate glyph based on value of first SPC x-coordinate
					glTranslatef(-(1.0 - *(processedData.begin() + 3)), 0.0, 0.0);

					// *********************** DRAW STICK FIGURE ***********************
					glyph.drawGlyphSF(pos2, pos3, stickFig.begin(), *classVec, *turt, DYNAMIC_ANGLES, POS_ANGLE,
						GLYPH_SCALE_FACTOR, SF_SEGMENT_CONSTANT, SF_ANGLE_SCALE, ANGLE_FOCUS, BIRD_FOCUS);

					// RESET THE CP AND CD
					turt->setCP(0.0, 0.0);

					// *********************** DRAW SPC AXES ***********************
					if (DRAW_AXES)
					{
						glLineWidth(2.0);	// Line width = 2.0
						glyph.drawAxesSPC(*pos2, *pos3, axesSPC.begin(),
							AXIS_LENGTH, 0.25, DOTTED_AXES);
					}

					glPopMatrix();	// Pop from gl matrix stack
					++colNum;		// Increment column index
				}

				if (it == allData.end())
				{	// If at end of data, set loop exit conditions
					colNum = NUM_COLUMNS;
					rowNum = NUM_ROWS;
				}
				else
				{	// Otherwise increment iterators
					++it;
					++classVec;
					++threshold;
				}
			}

			// If exit conditions have not been met
			if (rowNum != NUM_ROWS)
			{
				colNum = 0;		// Reset column index
				++rowNum;		// Increment row index
			}
		}

		// Draw filled polygons
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		// Return to original modelview matrix
		glPopMatrix();
	}
	// ******************************* DISPLAY GRID OF REPRESENTATIVE GLYPHS ****************************
	else if (!DISPLAY_SELECTOR && PC_OFF && !REPS_OFF)
	{
		// Draw grid of representative glyphs
		drawGrid();

		// Construct glyph tool
		SpcSfGlyph glyph = SpcSfGlyph();
		// Construct turtle tool
		TurtleG* turt = new TurtleG();
		Point2* pos2 = new Point2();
		Point2* pos3 = new Point2();

		// Initialize data iterators
		std::vector<std::vector<GLfloat>>::iterator repsIt = reps.begin();
		std::vector<bool>::iterator classIt = repsClass.begin();
		//std::vector<std::vector<GLfloat>>::iterator repsIt = mixedHood.begin();
		//std::vector<bool>::iterator classIt = mixedClass.begin();
		std::vector<std::string>::iterator labelsIt = labels.begin();
		//std::vector<std::string>::iterator labelsIt = distanceLabels.begin();

		glPushMatrix();		// Scale and translate glyph

		unsigned int rowNum = 0;	// Current Row index
		unsigned int colNum = 0;	// Current Column index
		int repCount = 0;
		while (rowNum < NUM_ROWS)		// Grid Rows loop
		{
			while (colNum < NUM_COLUMNS)	// Grid Columns loop
			{
				// If within threshold of current point,
				// display glyph at current grid index
				if (repsIt != reps.end())
				{
					// Draw polygon outlines
					glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

					// Pointer to current data point
					std::vector<GLfloat> processedData = *repsIt;

					// Vector holding Stick Figure attributes
					std::vector<GLfloat> stickFig{};
					// Populate: BN (angle), UCShape (length), UCSize (angle), NN (length)
					stickFig.push_back(*(processedData.begin()));	// Angle 1
					stickFig.push_back(*(processedData.begin() + 5));	// Length 1
					stickFig.push_back(*(processedData.begin() + 7));	// Angle 2
					stickFig.push_back(*(processedData.begin() + 1));	// Length 2

					// Vector holding SPC attributes
					std::vector<GLfloat> axesSPC{};
					// Populate: (CL, BN) (MA, SIN) (MIT, BC)
					axesSPC.push_back(*(processedData.begin()));		// X1
					axesSPC.push_back(*(processedData.begin() + 5));	// Y1
					axesSPC.push_back(*(processedData.begin() + 3));	// X2
					axesSPC.push_back(*(processedData.begin() + 4));	// Y2
					axesSPC.push_back(*(processedData.begin() + 8));	// X3
					axesSPC.push_back(*(processedData.begin() + 6));	// Y3

					glViewport(		// (rowNum x colNum)
						// Encode shift based off of first SPC axis horizontal/vertical shift.
						// Shift calculation:
							// the LARGER the value, the SMALLER the shift
							// maybe just subtract attribute value from (SCREEN_WIDTH / NUM_COLUMNS) / 2)
						((SCREEN_WIDTH / NUM_COLUMNS)* (colNum)) + 5,// - (maxAtr * ((SCREEN_WIDTH / NUM_COLUMNS) / 2)),
						((SCREEN_HEIGHT / NUM_ROWS)* (rowNum)),// - (maxAtr * ((SCREEN_HEIGHT / NUM_ROWS) / 2)),
						(SCREEN_WIDTH / NUM_COLUMNS),// * (1.0 + (maxAtr / 2)),
						(SCREEN_HEIGHT / NUM_ROWS)// * (1.0 + (maxAtr / 2))
					);

					// Print class totals in each cell
					glPushMatrix();
					glMatrixMode(GL_MODELVIEW);
					glLoadIdentity();

					glColor3f(0.0, 0.0, 0.0);
					glRasterPos2f(-0.9, -0.9);

					
					//if (repCount < 0)
					//{
						// Extract neighborhood class ratio
						std::string s = *labelsIt;

						// Initialize font
						void* font = GLUT_BITMAP_9_BY_15;

						// Print neighborhhood class ratio
						for (std::string::iterator labelIt = s.begin(); labelIt != s.end(); ++labelIt)
						{	// Loop through string, displaying each character
							char c = *labelIt;
							glutBitmapCharacter(font, c);
						}
						glPopMatrix();
						++labelsIt;	// Increment label iterator
					//}
					
					

					glPushMatrix();		// Push new matrix
					glMatrixMode(GL_PROJECTION);
					glLineWidth(4.0);	// Line width = 4.0

					// Vertical translation factor
					GLfloat glyphTransVert = -(1.0 - axesSPC[1]) + 1.0;;
					// GHorizaontal translateion factor
					GLfloat glyphTransHort = -(1.0 - axesSPC[0]) + 0.1;
					GLfloat TRANS_DEF = 0.0;	// Default translation factor
					GLfloat SCALE_DEF = 1.0;	// Default scaling factor

					glPushMatrix();	// Push new matrix for glyph transformations
					// Translate glyph right a based on first SPC shift, and up a set value
					glTranslatef(glyphTransHort, glyphTransVert, TRANS_DEF);

					GLfloat glyphScaleFactor = 2.0;

					glEnable(GL_BLEND);		// Enable blending.
					glDepthMask(GL_FALSE);	// Disable depth masking
					// *********************** DRAW STICK FIGURE ***********************
					if (repCount < 2)
					{	// Draw representative glyphs in color
						glyph.drawGlyphSF(pos2, pos3, stickFig.begin(), *classIt, *turt, DYNAMIC_ANGLES, POS_ANGLE,
							GLYPH_SCALE_FACTOR, SF_SEGMENT_CONSTANT, SF_ANGLE_SCALE, ANGLE_FOCUS, true);
					}
					else
					{	// Draw other glyphs in grey, if selected
						glyph.drawGlyphSF(pos2, pos3, stickFig.begin(), *classIt, *turt, DYNAMIC_ANGLES, POS_ANGLE,
							GLYPH_SCALE_FACTOR, SF_SEGMENT_CONSTANT, SF_ANGLE_SCALE, ANGLE_FOCUS, BIRD_FOCUS);
					}

					// Disable blending and resume depth mask
					glDepthMask(GL_TRUE);
					glDisable(GL_BLEND);
					++repCount;
					// RESET THE CP AND CD
					turt->setCP(0.0, 0.0);

					// *********************** DRAW SPC AXES ***********************
					if (DRAW_AXES)
					{
						glLineWidth(2.0);	// Line width = 2.0
						glyph.drawAxesSPC(*pos2, *pos3, axesSPC.begin(),
							AXIS_LENGTH, GLYPH_SCALE_FACTOR, DOTTED_AXES);
					}

					glPopMatrix();	// Pop from gl matrix stack
					++colNum;		// Increment column index
				}

				if (repsIt == reps.end())
				{	// If at end of data, set loop exit conditions
					colNum = NUM_COLUMNS;
					rowNum = NUM_ROWS;
				}
				else
				{	// Otherwise increment iterators
					++repsIt;
					++classIt;
				}
			}

			// If exit conditions have not been met
			if (rowNum != 4)
			{
				colNum = 0;		// Reset column index
				++rowNum;		// Increment row index
			}
		}

		// Draw filled polygons
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		// Return to original modelview matrix
		glPopMatrix();
	}
	// ******************************** DISPLAY PARALLEL COORDINATES ************************************
	else
	{	// Draw PC
		glClearColor(1.0, 1.0, 1.0, 0.0);	// Clear window to white
		glClearDepth(1.0f);					// Clear depth buffer
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glMatrixMode(GL_PROJECTION);	// Load GL_PROJECTION matrix
		glLoadIdentity();				// Load identity matrix
		glPushMatrix();					// Push new matrix

		// Set orthographic projection to SCREEN_WIDTH x SCREEN_HEIGHT
		gluPerspective(0, float(SCREEN_WIDTH) / float(SCREEN_HEIGHT), 0.1, 100.0);
		glViewport(0.0, 0.0, SCREEN_WIDTH, SCREEN_HEIGHT);
		gluOrtho2D(0.0, SCREEN_WIDTH, 0.0, SCREEN_HEIGHT);

		GLfloat axisWidth = 4.0;
		int numDimensions = 8;		// Num dimensions - 1
		glColor3f(0.0, 0.0, 0.0);	// Black
		glLineWidth(axisWidth);		// Axis width = 4

		// Draw PC axes
		for (int i = 0; i <= numDimensions; i++)
		{
			glBegin(GL_LINES);
			glVertex2f(i * (SCREEN_WIDTH / numDimensions), 0.0);
			glVertex2f(i * (SCREEN_WIDTH / numDimensions), SCREEN_HEIGHT);
			glEnd();
		}

		// Display current neighborhood
		// Initialize data iterator
		std::vector<std::vector<GLfloat>>::iterator it = allData.begin();
		glColor3f(0.0, 1.0, 0.0);
		GLfloat lineWidth = 2.0;
		glLineWidth(lineWidth);		// Line width = 2

		// Loop through data points
		while (it != allData.end())
		{	// If within current neighborhood
			if (*threshold)
			{
				//Determine color by class
				if (*classVec) glColor4f(0.0, 0.0, 0.8, 0.7);
				else glColor4f(0.8, 0.0, 0.0, 0.7);

				// Draw PC graph of data point
				glBegin(GL_LINE_STRIP);
				for (int i = 0; i <= numDimensions; i++)
				{
					glVertex2f(i * (SCREEN_WIDTH / numDimensions), SCREEN_HEIGHT * (*it)[i]);
				}
				glEnd();
			}

			++it;	// Increment iterators
			++classVec;
			++threshold;
		}
		glPopMatrix();
	}

	glutSwapBuffers();	// Swap buffers
	glFlush();			// Flush buffer
}


int main(int argc, char** argv)
{
	// Load parameters from config file
	loadConfig();

	// Implement config struct
	glutInit(&argc, argv);

	// Enable color and depth buffer
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);

	glutInitWindowSize(SCREEN_WIDTH, SCREEN_HEIGHT);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("SPC-SF");

	glutDisplayFunc(myDisplay);
	glutIdleFunc(myIdle);
	glutKeyboardFunc(myKeyboard);
	glutSpecialFunc(keyboard_special);
	glutMouseFunc(mouse_button_callback);
	glutMainLoop();
}