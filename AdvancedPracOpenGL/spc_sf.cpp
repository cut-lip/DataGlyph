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
#include <random>
#include <queue>

int sumx1 = 0;
int sumy1 = 0;
int sumx2 = 0;
int sumy2 = 0;
int sumx3 = 0;
int sumy3 = 0;

int flockPosition = 1;
int sizeHB = 0;
std::vector<GLfloat> averagePoint(9);
float minX = 0.0;
float maxY = 0.0;
float maxX = 10.0;
float minY = 10.0;

std::vector<std::vector<GLfloat>> analyzeGlyphs{};

bool SIZE_VIEW = false;
int FLOCK_LENGTH = 2;
int FLOCK_HEIGHT = 2;

// NEW PARAMS
const int HYPERBLOCK_SIZE = 345;
const int HYPERBLOCK_DATA_SIZE = 8;

const int SEED_DATASET_SIZE = 210;
const int SEED_DATA_SIZE = 5;

const int STUDENT_DATASET_SIZE = 395;
std::vector<std::vector<GLfloat>> passStudents{};

bool STUDENT_HYPER_COLLECTED = false;
std::vector<std::vector<GLfloat>> studentHyperblocks{};
std::vector<std::string> studentLabels{};

int HB_CLASS = 1;

/************************* DATA AND DIMENSION CONSTANTS *******************************/
int SCREEN_WIDTH;						/* Screen Width */
int SCREEN_HEIGHT;						/* Screen Height */
int VIEWPORT_SCALE = 8;					/* viewport scaling constant */
float THRESHOLD_VALUE = 3.0;			/* threshold calculation value */
float MIN_THRESHOLD =2.0;
float ALLOWED_DIFFERENCES = 2;
float AXIS_LENGTH = 1.0;				/* SPC axis length constant 8 */
unsigned int DATA_SIZE = 683;			/* cardinality of data set */
unsigned int DATA_INDEX = 0;			/* current index of data set */
unsigned int MAX_SIG_INDEX = 9;			/* maximum significant data index */
const unsigned int HEIGHT_SCALE = 10;	/* height scaling constant */
const unsigned int WIDTH_SCALE = 4;		/* width scaling constant */
const float GRID_MARGIN =1.0;			/* glyph grid margin*/
const float MAR = 10.0;					/* general use margin */

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
std::vector<int> repsSize{};
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


// IMPORT LINCOLN'S HYPERBLOCK DATA FROM CSV FILES
void importHyperblockData(std::vector<std::vector<GLfloat>>* allData)
{
	// Read data 
	std::string line = "";
	std::ifstream myFile("hyperblocks/HB1.csv");

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

/* FUNCTION SIGNATURES */
void drawGridSPC(GLfloat originX, GLfloat originY, GLfloat endX, GLfloat endY, int dimension);

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
void drawLocatedGlyphs(std::vector<GLfloat>* normalData, bool classify, int size, int iteration, std::string hbLabel, std::string hbLabel2)
{
	// encode colors to bird glyph winds in located glyphs
	float colors[6];
	colors[0] = *(normalData->begin() + 2);
	colors[1] = *(normalData->begin() + 1);
	colors[2] = *(normalData->begin() + 0);
	colors[3] = *(normalData->begin() + 2);
	colors[4] = *(normalData->begin() + 7);
	colors[5] = *(normalData->begin() + 6);
	
	// If color values too high, round down so still visible
	for (unsigned int i = 0; i < 6; i++)
	{
		if (colors[i] > 0.8)
		{
			colors[i] = colors[i] - 0.1;
		}
	}
	
	glPushMatrix();

	// Encode custom attributes to SPC positions
	GLfloat ucsize = *(normalData->begin() + 1);

	// Uniformity of cell size
	GLfloat ucshape = *(normalData->begin() + 2);

	// Clump thickness
	GLfloat cl = *(normalData->begin());

	// Bland chromatin
	GLfloat bn = *(normalData->begin() + 5);

	// Marginal adhesion
	GLfloat bc = *(normalData->begin() + 6);

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
	std::vector<GLfloat> position{ cl, bn, ucsize, bc, bn, ucshape };
	//std::vector<GLfloat> position{ uc, bn, bc, cl, bn, mg };
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

	if (flockPosition == 1)	// draw flock based on position ID
	{
		// Place the viewport
		glViewport(
			x1,
			y1,
			SCREEN_WIDTH / VIEWPORT_SCALE,
			SCREEN_WIDTH / VIEWPORT_SCALE
		);

		
		//Glyph1
		glLineWidth(4.0);

		// Switch between glyphs and bars representing size of cluster
		if (!SIZE_VIEW)
		{
			glyph.drawGlyphSF(pos2, pos3, stickFig.begin(), classify, *turt, DYNAMIC_ANGLES, POS_ANGLE, GLYPH_SCALE_FACTOR, 0.1, 2.0, ANGLE_FOCUS, BIRD_FOCUS, colors);
			// RESET THE CP AND CD
			turt->setCP(0.0, 0.0);
			// *********************** DRAW SPC AXES ***********************
			if (DRAW_AXES)
			{
				glLineWidth(2.0);
				glyph.drawAxesSPC(*pos2, *pos3, axesSPC.begin(), AXIS_LENGTH, GLYPH_SCALE_FACTOR, DOTTED_AXES);
			}
		}

		/*
		// draw bar representing size of glyph
		if (classify)	// blue for benign
		{
			glColor3f(0.0, 0.0, 1.0);
			// If too small to be visible, round up
			if (size < 10) size += 10;
			// display bar
			glRectf(0, 0, 0.2, (GLfloat)size / 500.0);
		}
		if (!classify)	// red for malignant
		{
			glColor3f(1.0, 0.0, 0.0);
			// If too small, round size up
			if (size < 10) size += 20;
			// display bar
			glRectf(0, 0, 0.2, (GLfloat)size / 500.0);
		}
		*/
	}
	
	// Locate the viewport for glyph drawing
	GLfloat x2 = ((SCREEN_WIDTH / WIDTH_SCALE) * *++positionIt) + (SCREEN_WIDTH / 3) - ((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2));
	GLfloat y2 = ((SCREEN_HEIGHT - (SCREEN_HEIGHT / HEIGHT_SCALE)) * *++positionIt) - ((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2));

	// maintain record of minimum and maximum Y values in hyperblock,
	// for placement of average glyph
	if (y2 > maxY)
	{	// if current y value is greater than max
		maxY = y2;
		maxX = x2;
	}

	if (y2 < minY)
	{	// if current y value is less than max
		minY = y2;
		minX = x2;
	}

	if (flockPosition == 2)
	{

		// Locate the viewport for glyph drawing
		glViewport(
			x2,
			y2,
			SCREEN_WIDTH / VIEWPORT_SCALE,
			SCREEN_WIDTH / VIEWPORT_SCALE
		);
		
		//Glyph2
		glLineWidth(4.0);
		glyph.drawGlyphSF(pos2, pos3, stickFig.begin(), classify, *turt, DYNAMIC_ANGLES, POS_ANGLE, 0.25, 0.1, 2.0, ANGLE_FOCUS, BIRD_FOCUS, colors);
		// RESET THE CP AND CD
		turt->setCP(0.0, 0.0);

		// *********************** DRAW SPC AXES ***********************
		if (DRAW_AXES)
		{
			glLineWidth(2.0);
			glyph.drawAxesSPC(*pos2, *pos3, axesSPC.begin(), AXIS_LENGTH, 0.25, DOTTED_AXES);
		}


		/*
		if (classify)
		{
			glColor3f(0.0, 0.0, 1.0);
			if (size < 10) size += 10;
			glRectf(0, 0, 0.2, (GLfloat)size / 500.0);
		}
		if (!classify)
		{
			glColor3f(1.0, 0.0, 0.0);
			if (size < 10) size += 20;
			glRectf(0, 0, 0.2, (GLfloat)size / 500.0);
		}
		*/
	}


	// Locate the viewport for glyph drawing
	GLfloat x3 = ((SCREEN_WIDTH / WIDTH_SCALE) * *++positionIt) + (((2 * SCREEN_WIDTH) / 3)) - ((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2));
	GLfloat y3 = ((SCREEN_HEIGHT - (SCREEN_HEIGHT / HEIGHT_SCALE)) * *++positionIt) - ((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2));

	if (flockPosition == 3)
	{
		// Locate the viewport for glyph drawing
		glViewport(
			x3,
			y3,
			SCREEN_WIDTH / VIEWPORT_SCALE,
			SCREEN_WIDTH / VIEWPORT_SCALE
		);
		

		//Glyph3
		glLineWidth(4.0);
		glyph.drawGlyphSF(pos2, pos3, stickFig.begin(), classify, *turt, DYNAMIC_ANGLES, POS_ANGLE, 0.25, 0.1, 2.0, ANGLE_FOCUS, BIRD_FOCUS, colors);
		// RESET THE CP AND CD
		turt->setCP(0.0, 0.0);
		// *********************** DRAW SPC AXES ***********************
		if (DRAW_AXES)
		{
			glLineWidth(2.0);
			glyph.drawAxesSPC(*pos2, *pos3, axesSPC.begin(), AXIS_LENGTH, 0.25, DOTTED_AXES);
		}
		/*
		if (classify)	// blue for benign
		{
			glColor3f(0.0, 0.0, 1.0);
			if (size < 10) size += 10;
			glRectf(0, 0, 0.2, (GLfloat)size / 500.0);
		}
		if (!classify)	// red for malignant
		{
			glColor3f(1.0, 0.0, 0.0);
			// If too small, round size up
			if (size < 10) size += 20;
			// display bar
			glRectf(0, 0, 0.2, (GLfloat)size / 500.0);
		}
		*/

	}

	// Collect locations of edges as the function iterates,
	// to form an average pattern for all glyphs in the hyperblock
	sumx1 += x1 + (SCREEN_WIDTH / VIEWPORT_SCALE);
	sumy1 += y1 + (SCREEN_WIDTH / VIEWPORT_SCALE);
	sumx2 += x2 + (SCREEN_WIDTH / VIEWPORT_SCALE);
	sumy2 += y2 + (SCREEN_WIDTH / VIEWPORT_SCALE);
	sumx3 += x3 + (SCREEN_WIDTH / VIEWPORT_SCALE);
	sumy3 += y3 + (SCREEN_WIDTH / VIEWPORT_SCALE);


	/*
	if (iteration == sizeHB)	// if on last iteration of drawing hyperblock
	{
		glPushMatrix();
		// check for too high here
		// first iter just use high values
			// encode colors to bird glyph winds in located glyphs
		float avgColors[6];
		avgColors[0] = *(averagePoint.begin() + 2);
		avgColors[1] = *(averagePoint.begin() + 1);
		avgColors[2] = *(averagePoint.begin() + 0);
		avgColors[3] = *(averagePoint.begin() + 2);
		avgColors[4] = *(averagePoint.begin() + 7);
		avgColors[5] = *(averagePoint.begin() + 6);

		// Encode angles with most meaningful attributes
		std::vector<GLfloat> avgStickFig{};
		// Populate: CL (angle), UC (length), BN (angle), BC (length)
		avgStickFig.push_back(*(averagePoint.begin()));			// Angle 1
		avgStickFig.push_back(*(averagePoint.begin() + 1));	// Length 1
		avgStickFig.push_back(*(averagePoint.begin() + 5));	// Angle 2
		avgStickFig.push_back(*(averagePoint.begin() + 6));	// Length 2

		glPushMatrix();		// Push new Gl matrix
		glViewport(			// set glyph viewport to 1/3 of window
			(SCREEN_WIDTH) / 3,
			0,
			SCREEN_WIDTH / 3,
			SCREEN_HEIGHT
		);
		glTranslatef(-0.3, 0.3, 0.0);	// position glyph above graph
		glLineWidth(9.0);				// Line width 7

		// draw the glyph
		glyph.drawGlyphSF(pos2, pos3, avgStickFig.begin(), classify, *turt,
			DYNAMIC_ANGLES, POS_ANGLE, 0.25, 0.1, 2.0, ANGLE_FOCUS, BIRD_FOCUS, avgColors);
		glPopMatrix();

		// Draw average of all connections between glyphs in hyperblock
		sumx1 /= sizeHB;
		sumy1 /= sizeHB;
		sumx2 /= sizeHB;
		sumy2 /= sizeHB;
		sumx3 /= sizeHB;
		sumy3 /= sizeHB;

		// Set viewport to middle SPC axes
		glViewport(0.0, 0.0, SCREEN_WIDTH, SCREEN_HEIGHT);
		gluOrtho2D(0.0, SCREEN_WIDTH, 0.0, SCREEN_HEIGHT);
		if (classify)	glColor4f(0.0, 0.0, 1.0, 7.0);
		else			glColor4f(1.0, 0.0, 0.0, 7.0);

		glPushMatrix();		// Push new modelview matrix for translation
		glLineWidth(5.0);	// Line width 5
		// Enable line stipple
		glPushAttrib(GL_ENABLE_BIT);
		// glPushAttrib is done to return everything to normal after drawing

		glLineStipple(1, 0x00FF);  // Draw connection between SPC axes as dotted line
		glEnable(GL_LINE_STIPPLE);
		glPushMatrix();
		glTranslatef(-((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2)), -((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2)), 0);
		glTranslatef(-160, -160, 0);

		// Draw first edge
		glBegin(GL_LINES);
		glVertex2f(sumx1 + (SCREEN_WIDTH / VIEWPORT_SCALE), sumy1 + (SCREEN_WIDTH / VIEWPORT_SCALE));
		glVertex2f(sumx2 + (SCREEN_WIDTH / VIEWPORT_SCALE), sumy2 + (SCREEN_WIDTH / VIEWPORT_SCALE));
		glEnd();

		// Draw second edge
		drawArrow(Point2(sumx2 + (SCREEN_WIDTH / VIEWPORT_SCALE), sumy2 + (SCREEN_WIDTH / VIEWPORT_SCALE)),
			Point2(sumx3 + (SCREEN_WIDTH / VIEWPORT_SCALE), sumy3 + (SCREEN_WIDTH / VIEWPORT_SCALE)), (2 * (SCREEN_HEIGHT / 2)));
		glPopAttrib();
		GLfloat maxSizeHB = 400;
		GLfloat barShiftFactor = 40;

		glPopMatrix();

		// draw bar representing size of glyph
		if (classify)	// blue for benign
		{
			glColor3f(0.0, 0.0, 1.0);
			// If too small to be visible, round up
			if (sizeHB < 10) sizeHB += 10;
			// display bar
			glRectf(SCREEN_WIDTH - 55, barShiftFactor,
				SCREEN_WIDTH - barShiftFactor,
				barShiftFactor + (GLfloat)sizeHB * (SCREEN_HEIGHT/ maxSizeHB));
		}
		if (!classify)	// red for malignant
		{
			glColor3f(1.0, 0.0, 0.0);
			// If too small, round size up
			if (sizeHB < 10) sizeHB += 20;
			// display bar
			glRectf(SCREEN_WIDTH - 55,
				barShiftFactor, SCREEN_WIDTH - barShiftFactor,
				barShiftFactor + (GLfloat)sizeHB * (SCREEN_HEIGHT/ maxSizeHB));
		}
		
		// Prepare for printing labels
		glPushMatrix();
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glColor3f(0.0, 0.0, 0.0);	// Set color to black

		// Position font for label row1
		glRasterPos2f(SCREEN_WIDTH - 55, 25);
		// Initialize font
		void* font = GLUT_BITMAP_8_BY_13;

		// Print hyperblock class ratio
		for (std::string::iterator labelIt = hbLabel.begin(); labelIt != hbLabel.end(); ++labelIt)
		{	// Loop through string, displaying each character
			char c = *labelIt;
			glutBitmapCharacter(font, c);
		}

		// Position font for label row2
		glRasterPos2f(SCREEN_WIDTH - 55, 10);
		// Print hyperblock size
		for (std::string::iterator labelIt = hbLabel2.begin(); labelIt != hbLabel2.end(); ++labelIt)
		{	// Loop through string, displaying each character
			char c = *labelIt;
			glutBitmapCharacter(font, c);
		}

		glPopMatrix();
		
		glPopMatrix();
	}
	else
	{
		glPopMatrix();
	}

	*/
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glViewport(0.0, 0.0, SCREEN_WIDTH, SCREEN_HEIGHT);
	gluOrtho2D(0.0, SCREEN_WIDTH, 0.0, SCREEN_HEIGHT);
	glColor4f(0.0, 0.0, 0.0, 7.0);
	int gridMargin = 80;
	int heightMargin = 10;
	int numSPC = 3;
	int dimension = 5;
	glPushMatrix();
	glLineWidth(4.0);
	// Draw SPC axis frame * 3
	// Draw 1st SPC axes
	drawArrow(Point2(0, 0), Point2(0, SCREEN_HEIGHT - 10), SCREEN_HEIGHT);
	drawArrow(Point2(0, 0), Point2((SCREEN_WIDTH / 4) + (SCREEN_WIDTH / 24), 0), SCREEN_HEIGHT);
	// Draw grid inside axes
	drawGridSPC(0, 0,
		(SCREEN_WIDTH / 4) + (SCREEN_WIDTH / 24), SCREEN_HEIGHT - 10, dimension);
	glLineWidth(3.0);
	// Draw 2nd SPC axes
	drawArrow(Point2(((SCREEN_WIDTH) / numSPC), 0),
		Point2((SCREEN_WIDTH / numSPC), SCREEN_HEIGHT - heightMargin), SCREEN_HEIGHT);
	drawArrow(Point2((SCREEN_WIDTH / numSPC), 0),
		Point2(((7 * SCREEN_WIDTH) / 12) + (SCREEN_WIDTH / 24), 0), SCREEN_HEIGHT);
	// Draw grid inside axes
	drawGridSPC((SCREEN_WIDTH / numSPC) + gridMargin, 0,
		(SCREEN_WIDTH / 4) + (SCREEN_WIDTH / 24), SCREEN_HEIGHT - heightMargin, dimension);
	glLineWidth(3.0);
	// Draw 3rd SPC axes
	drawArrow(Point2(((2 * SCREEN_WIDTH) / numSPC), 0),
		Point2(((2 * SCREEN_WIDTH) / numSPC), SCREEN_HEIGHT - heightMargin), SCREEN_HEIGHT);
	drawArrow(Point2(((2 * SCREEN_WIDTH) / numSPC), 0),
		Point2(((11 * SCREEN_WIDTH) / 12) + (SCREEN_WIDTH / 24), 0), SCREEN_HEIGHT);
	// Draw grid inside axes
	drawGridSPC(((2 * SCREEN_WIDTH) / numSPC) + gridMargin, 0,
		(SCREEN_WIDTH / 4) + (SCREEN_WIDTH / 24), SCREEN_HEIGHT - heightMargin, dimension);

	// DON't POP MATRIX HERE
	
	// Draw edges between glyphs
	if (DRAW_EDGES)
	{	// Determine which class/color the edge belongs to
		if (classify)	glColor4f(0.0, 0.0, 1.0, 7.0);
		else			glColor4f(1.0, 0.0, 0.0, 7.0);
		glLineWidth(0.5);
		glPushMatrix();	// Push new modelview matrix for translation
		glTranslatef(-((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2)), -((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2)), 0);

		// Draw first edge
		glBegin(GL_LINES);
		glVertex2f(x1 + (SCREEN_WIDTH / VIEWPORT_SCALE), y1 + (SCREEN_WIDTH / VIEWPORT_SCALE));
		glVertex2f(x2 + (SCREEN_WIDTH / VIEWPORT_SCALE), y2 + (SCREEN_WIDTH / VIEWPORT_SCALE));
		glEnd();

		//drawArrow(Point2(x1 + (SCREEN_WIDTH / VIEWPORT_SCALE), y1 + (SCREEN_WIDTH / VIEWPORT_SCALE)),
			//Point2(x2 + (SCREEN_WIDTH / VIEWPORT_SCALE), y2 + (SCREEN_WIDTH / VIEWPORT_SCALE)), SCREEN_HEIGHT);

		// Draw second edge
		drawArrow(Point2(x2 + (SCREEN_WIDTH / VIEWPORT_SCALE), y2 + (SCREEN_WIDTH / VIEWPORT_SCALE)),
			Point2(x3 + (SCREEN_WIDTH / VIEWPORT_SCALE), y3 + (SCREEN_WIDTH / VIEWPORT_SCALE)), (2*(SCREEN_HEIGHT/3)));

		glPopMatrix(); // Pop modelview matrix used for translation
	}
	glPopMatrix();
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
	// modified HyClu hypercube clustering algorithm. algorithm now forms Hyperblocks.
	// Threshold value for attributes used in SPC shifts is MIN_THRESHOLD
	// Threshold value for attributes allowed to expand is THRESHOLD
	
	// Initialize iterators for vectors to be compared
	std::vector<float>::iterator it1 = vec1->begin();
	std::vector<float>::iterator it2 = vec2->begin();

	int count = 0;
	for (it1; it1 != vec1->end(); ++it1)
	{
		if (count == 0 || count == 1 || count == 2 || count == 5 || count == 6)
		{
			if (abs(*it1 - *it2) > MIN_THRESHOLD)
			{
				return false;
			}
		}
		else if (count == 3 || count == 4 || count == 7 || count == 8)
		{
			if (abs(*it1 - *it2) > THRESHOLD_VALUE)
			{
				return false;
			}
		}

		++count;
		++it2;
	}
	// Return true if all attributes are within threshold
	return true;

	
	/*
	std::vector<float>::iterator it1 = vec1->begin();
	std::vector<float>::iterator it2 = vec2->begin();
	int count = 0;
	for (it1; it1 != vec1->end(); ++it1)
	{
		if (abs(*it1 - *it2) > THRESHOLD_VALUE)
		{
			return false;
		}
		++it2;
	}
	// Return true if all attribtes are within threshold
	return true;
	*/
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

	else if (key == GLUT_KEY_F12)
	{
		SIZE_VIEW = !SIZE_VIEW;
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
	if (key == 'a') {
		// move flock to the left
		--flockPosition;
	}

	if (key == 's') {
		// move flock to the right
		++flockPosition;
	}

	// Redisplay with updated parameters
	glutPostRedisplay();
}

/* drawGridSPC - draw a grid for each SPC axis
*  @params -	originX - x coord of plane origin
*				originY - y coord of plane origin
*				endX	- max x value in plane
*				endY	- max y values in plane
*				dimension - dimensions of grid (dimension x dimension)
* */
void drawGridSPC(GLfloat originX, GLfloat originY, GLfloat endX, GLfloat endY, int dimension)
{
	// ************** DRAW SPC GRID (dimension x dimension) ****************
	glPushMatrix();
	glColor4f(0.0, 0.0, 0.0, 1.0);				// Draw in black
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);	// Draw outlines

	glLineWidth(1.0);		// Set line width to 1
	glPushMatrix();			// Push new matrix
	for (unsigned int rowNum = 0; rowNum < dimension; ++rowNum)		// Rows
	{
		for (unsigned int colNum = 0; colNum < dimension; ++colNum)	// Columns
		{
			glRectf(		// Draw outlined rectangles for SPC grid
				originX + ((endX / dimension) * colNum),
				originY + ((endY / dimension) * rowNum),
				((endX) / dimension) * (colNum + 1),
				((endY) / dimension) * (rowNum + 1)
			);
		}
	}
	// Reset OpenGL state
	glLineWidth(2.0);		// Set line width to 2
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);	// Fill shapes
	glPopMatrix();		// Pop a matrix from gl matrix stack
	glPopMatrix();		// Pop a matrix from gl matrix stack
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
				 (SCREEN_HEIGHT - (SCREEN_HEIGHT / NUM_ROWS) * rowNum) - 1.0,
				((SCREEN_WIDTH / NUM_COLUMNS) * (colNum + 1)) + GRID_MARGIN,
				  (SCREEN_HEIGHT - (SCREEN_HEIGHT / NUM_ROWS) * (rowNum + 1)) - 1.0
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

/*
analyzeGlyphShape
This function analyzes all data points in the set to form a subset of
points from both classes which are most similar to each other
(and thus at risk for misclassification)
@param		all_data	pointer to data set
			classify	pointer to vector of class labels
*/
void analyzeGlyphShape(std::vector<std::vector<GLfloat>>* all_Data,
	std::vector<bool>* classify)
{
	// Make deep copy of all data
	std::vector<std::vector<GLfloat>> allData(*all_Data);
	std::vector<bool> classVec(*classify);
	// Initialize iterator from the class vector
	std::vector<bool>::iterator classIt = classVec.begin();

	// Initialize temp variables to track wing size
	float leftWingMax = 0.0;
	float rightWingMax = 0.0;
	// Initalize temp vector to hold current max point
	std::vector<GLfloat> currMaxWingPoint;

	// Search for the point from class one with largest attributes for wing length
	for (std::vector<std::vector<GLfloat>>::iterator dataIt = allData.begin();
		dataIt != allData.end(); ++dataIt)
	{
		if (*classIt)	// if point belongs to the benign class
		{
			// analyze for max wing lengths
			std::vector<float>::iterator currIt = dataIt->begin();
			if ((*(currIt + 1) > leftWingMax) && (*(currIt + 6) > rightWingMax))
			{
				// save current max point to mariable
				currMaxWingPoint = *dataIt;

				// update max length variables
				leftWingMax = *(currIt + 1);
				rightWingMax = *(currIt + 6);
			}
		}
		++classIt;
	}

	// save selected point from class one for analysis
	analyzeGlyphs.push_back(currMaxWingPoint);

	classIt = classVec.begin();	// reset class iterator
		
	// Search for the most similar points from class 2, to the selected pointd
	// Search for the point from class one with largest attributes for wing length
	for (std::vector<std::vector<GLfloat>>::iterator dataIt = allData.begin();
		dataIt != allData.end(); ++dataIt)
	{
		if (!*classIt)	// if point belongs to the benign class
		{
			// initialize iterator to current point
			std::vector<float>::iterator currIt = dataIt->begin();
			// if current point is close to the selected point from class one
			if ((*(currIt + 1) <= leftWingMax) && (*(currIt + 6) <= rightWingMax))
			{
				// save point for analysis
				analyzeGlyphs.push_back(*dataIt);
			}

		}
		++classIt;	// Increment iterators
	}
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
		addThis = true;
		// Use threshold values to remove points
		while (dataIt != allData.end())
		{
			addThis = true;
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
								addThis = false;	// flag this data point as divergent from the center
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
		labels.push_back(pureLabel + "%, n=" + std::to_string(benCount + malCount) + "");
		repsSize.push_back(benCount + malCount);
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

	/*
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
	*/
	// Insert cluster representative glyph
	//mixedHood.insert(mixedHood.begin(), reps[CLUSTER + 2]);
	// Insert clas of representative glyph
	//mixedClass.insert(mixedClass.begin(), repsClass[CLUSTER + 2]);

	
	// Normalize data in REPS vector
	int i = 0;
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
void mergerHyperblock(std::vector<std::vector<GLfloat>>* all_Data, std::vector<bool>* classify)
{
	std::vector<std::vector<GLfloat>> allData(*all_Data);					// deep copy data
	std::vector<bool> classifyCopy(*classify);
	std::vector<bool> thresholds{};	// track hyperblock membership

	bool addThis = false;	// add curr point to hyperblock flag

	int passCount = 0;
	int failCount = 0;

	int count = 0;
	while (allData.size() != 0)	// continue til all points are clustered
	{
		// Compute points within threshold
		thresholds = computeAllDistances(&allData[0], &allData);

		std::vector< std::vector<GLfloat>> tempData{};							// temp vec for average glyph
		std::vector< std::vector<GLfloat>> tempDataPass{};
		std::vector< std::vector<GLfloat>> tempDataFail{};
		std::vector<std::vector<GLfloat>>::iterator dataIt = allData.begin();	// Initialize vector iterators
		std::vector<bool>::iterator classIt = classifyCopy.begin();
		std::vector<bool>::iterator threshIt = thresholds.begin();

		// Save first point in cluster for comparison
		std::vector<float> currPoint = (*dataIt);

		addThis = false;	// reset add curr point flag
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
					++passCount;
					tempDataPass.push_back(*dataIt);
					if (count == 2)
					{
						passStudents.push_back(*dataIt);
					}
				}
				else if (!*classIt)
				{
					++failCount;
					tempDataFail.push_back(*dataIt);
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

		++count;
		studentLabels.push_back(std::to_string(passCount) + " pass, " + std::to_string(failCount) + " fail");
		 
		// Initialize new rep vector
		std::vector<GLfloat> tempVec{};

		/* Compute average glyph of neighborhood by taking tempData,
		*  the set of vectors saved during clustering, and computing
		*  the average value for each attribute from the data points
		*  in the set.
		*/
		if (passCount > failCount)
		{
			for (unsigned int index = 0; index < 10; ++index)
			{	// For the current index of each vector
				GLfloat tempAttr = 0.0;
				for (auto& vec : tempDataPass)
				{	// For each vector
					tempAttr += vec[index];
				}
				// Compute average value from sum
				tempAttr /= tempDataPass.size();
				// Add to representative point
				tempVec.push_back(tempAttr);
			}
		}
		else
		{
			for (unsigned int index = 0; index < 10; ++index)
			{	// For the current index of each vector
				GLfloat tempAttr = 0.0;
				for (auto& vec : tempDataFail)
				{	// For each vector
					tempAttr += vec[index];
				}
				// Compute average value from sum
				tempAttr /= tempDataFail.size();
				// Add to representative point
				tempVec.push_back(tempAttr);
			}
		}

		// Add representative vector to set
		studentHyperblocks.push_back(tempVec);

		// Reset class counters
		passCount = 0;
		failCount = 0;
	}

	STUDENT_HYPER_COLLECTED = true;	// set collected flag
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

		// Duplicate data attributes x 4 to move from 9-D to 40-D
		//data.insert(data.end(), data.begin(), data.end());
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
	std::ifstream myFile("seeds_3.txt");

	// For each line in the hyperblock data file
	for (unsigned int i = 0; i < 27; i++)
	{
		getline(myFile, line);
		//int vecCount = SEED_DATA_SIZE;

		// Split string into a vector
		std::vector<int> dataFloat;
		std::stringstream ss(line);

		while (ss.good()) {
			std::string substr = "";
			getline(ss, substr, '\t');	// extract line from file
			dataFloat.push_back(stof(substr));	// convert to float and push to vector
			continue;
		}

		// Remove class attribute
		dataFloat.pop_back();

		std::vector<GLfloat> data(dataFloat.begin(), dataFloat.end());
		(*allData)[i] = data;	// add data point to hyperblock vector
	}
		myFile.close();
}

// IMPORT STUDENT DATA
void importStudentData(std::vector<std::vector<GLfloat>>* allData, std::vector<bool>* classify)
{
	std::string line = "";
	std::ifstream myFile("student_new_2.txt");


	// For each line in the hyperblock data file
	for (unsigned int i = 0; i < STUDENT_DATASET_SIZE; i++)
	{
		getline(myFile, line);
		//int vecCount = SEED_DATA_SIZE;

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

		// Set color determined by class
		if (*--(data.end()) == 0) classify->push_back(true);
		else
		{
			classify->push_back(false);
		}

		data.erase(--data.end());	// remove class label from datapoint

		(*allData)[i] = data;	// add data point to hyperblock vector
	}
	myFile.close();
}


// vector comparison function to determine the sum of differences between two vectors
GLfloat compareHyperblocks(const std::vector<GLfloat>& vec1, const std::vector<GLfloat>& vec2) {
	GLfloat sumDifference = 0.0;
	// sum all differences between attributes of the vectors together.
	for (std::size_t i = 0; i < vec1.size(); ++i) {
		sumDifference += std::abs(vec1[i] - vec2[i]);
	}
	// return the sum of the differences
	return sumDifference;
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
	std::vector<std::vector<GLfloat>> seeds(27);
	importSeedData(&seeds);

	// student container
	std::vector<bool> studentClass{};
	std::vector<std::vector<GLfloat>> students(STUDENT_DATASET_SIZE);
	importStudentData(&students, &studentClass);
	// Normalize student data
	unsigned int index = 0;
	for (auto& vec : students)
	{
		float studentNormalizeFactors[] = { 5, 5, 5, 5, 6, 6, 6, 6, 20, 20 };
		float studentNormalizeMins[] = { 0, 0, 1, 1, 1, 1, 1, 1, 0, 0 };

		std::vector<float> normalData;		// Normalize data to [0, 1]
		for (unsigned int i = 0; i < vec.size(); ++i)
		{
			normalData.push_back(((vec[i]) - studentNormalizeMins[i]) / (studentNormalizeFactors[i] - studentNormalizeMins[i]));
		}
		students[index] = normalData;
		++index;
	}
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
	std::vector<bool> hbClass = std::vector<bool>(HYPERBLOCK_SIZE, HB_CLASS);

	// SPLIT DATA VECTOR 90/10
	std::random_device rd;	// initialize random number generator
	std::mt19937 g(rd());	// ensure different seeds for different runs

	// Randomly shuffle the elements of the data vector
	std::shuffle(allData.begin(), allData.end(), g);

	// Determine what index to split the shuffled data vector
	std::size_t dataSize = allData.size();
	std::size_t splitIndex = (dataSize * 0.9);

	// Initialize training and testing vectors
	std::vector<std::vector<GLfloat>> trainingData(splitIndex);
	std::vector<std::vector<GLfloat>> testingData(dataSize - splitIndex);

	// Copy data from allData into training (90%) and testing (10%) vectors
	std::copy(allData.begin(), allData.begin() + splitIndex, trainingData.begin());
	std::copy(allData.begin() + splitIndex, allData.end(), testingData.begin());

	// ONE TIME OPERATIONS
	//if (!IDEAL_COLLECTED)	// Collect ideal class glyphs, if not done
		//getIdealGlyphs(&allData, &classify);

	if (!REPS_COLLECTED)	// Collect representative glyphs, if not already done
	{
		//getRepresentativeGlyphs(&allData, &classify);
		getRepresentativeGlyphs(&trainingData, &classify);
	}
	if (!STUDENT_HYPER_COLLECTED)
	{
		mergerHyperblock(&students, &studentClass);
	}

	analyzeGlyphShape(&allData, &classify);
	index = 0;
	for (auto& vec : analyzeGlyphs)
	{
		const unsigned int scaleFactor = 10;
		std::vector<float> normalData;		// Normalize data to [0, 1]
		for (std::vector<float>::iterator iter = vec.begin(); iter < vec.end(); iter++)
		{
			normalData.push_back(*iter / scaleFactor);
		}
		analyzeGlyphs[index] = normalData;
		++index;
	}

	// Compute points within threshold
	std::vector<bool> close = computeAllDistances(&allData[DATA_INDEX], &allData);
	// Normalize seed data
	index = 0;
	for (auto& vec : seeds)
	{
		float seedNormalizeFactors[] = { 21.2, 17.3, 6.7, 4.1, 6.6 };
		float seedMins[] = { 10.5, 12.4, 4.8, 2.6, 4.5 };

		std::vector<float> normalData;		// Normalize data to [0, 1]
		for (unsigned int i = 0; i < vec.size(); ++i)
		{
			normalData.push_back( (vec[i] - seedMins[i]) / (seedNormalizeFactors[i] - seedMins[i]));
		}
		seeds[index] = normalData;
		++index;
	}

	// Replicate seed data attributes
	for (auto& vec : seeds)
	{
		vec.insert(vec.end(), vec.begin(), vec.end());
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

	// Normalize data
	index = 0;
	for (auto& vec : testingData)
	{
		const unsigned int scaleFactor = 10;
		std::vector<float> normalData;		// Normalize data to [0, 1]
		for (std::vector<float>::iterator iter = vec.begin(); iter < vec.end(); iter++)
		{
			normalData.push_back(*iter / scaleFactor);
		}
		testingData[index] = normalData;
		++index;
	}

	// Normalize data
	index = 0;
	for (auto& vec : trainingData)
	{
		const unsigned int scaleFactor = 10;
		std::vector<float> normalData;		// Normalize data to [0, 1]
		for (std::vector<float>::iterator iter = vec.begin(); iter < vec.end(); iter++)
		{
			normalData.push_back(*iter / scaleFactor);
		}
		trainingData[index] = normalData;
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

	// Randomly pick a data point from the testing data
	std::uniform_int_distribution<std::size_t> distribution(0, testingData.size() - 1);
	std::size_t randomIndex = distribution(g);  // g is the random number generator

	// Save chosen point 
	std::vector<GLfloat> testingDataPoint = testingData[randomIndex];
	// Remove from testing data
	testingData.erase(testingData.begin() + randomIndex);

	// Create a priority queue to save the top five most similar points to the chosen testing point
	std::priority_queue<std::vector<GLfloat>, std::vector<std::vector<GLfloat>>, decltype(compareHyperblocks)*> similarVectors(compareHyperblocks);

	// struct VectorData allows sorting vectors and labels concurrently
	// while also storing sumDifference from chosen testing datapoint
	struct VectorData
	{
		std::vector<GLfloat> vector;
		std::string label;
		GLfloat sumDifference;

		// VectorData constructor
		VectorData(const std::vector<GLfloat>& v, const std::string& l, GLfloat sd)
			: vector(v), label(l), sumDifference(sd)
		{}
	};

	// Initialize vector of VectorData
	std::vector<VectorData> differences;

	// Loop through representative glyphs of hyperblocks
	for (std::size_t i = 0; i < reps.size(); ++i)
	{
		// Save index of current label and glyph
		const auto& vector = reps[i];
		const auto& label = labels[i];
		// Compute sum difference of current vector to testing data
		GLfloat sumDifference = compareHyperblocks(testingDataPoint, vector);
		// Save vector, label, and difference in differences
		differences.emplace_back(vector, label, sumDifference);
	}

	/*
	// Calculate the sumDifference for each vector in reps
	std::vector<std::pair<std::vector<GLfloat>, GLfloat>> differences;
	for (const auto& vector : reps)
	{
		// use the compareHyperblocks comparator function to determine the sum difference
		GLfloat sumDifference = compareHyperblocks(testingDataPoint, vector);
		// place into differences vector
		differences.emplace_back(vector, sumDifference);
	}
	*/

	/*
	// Sort the vectors based on the sumDifference in ascending order
	std::sort(differences.begin(), differences.end(),
		[](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
		*/

	std::sort(differences.begin(), differences.end(),
		[](const VectorData& lhs, const VectorData& rhs) { return lhs.sumDifference < rhs.sumDifference; });

	// Extract the top five most similar vectors
	std::vector<std::vector<GLfloat>> mostSimilarVectors;
	std::vector<std::string> mostSimilarLabels;

	/*
	// extract most similar vectors, checking if total is less than 5 for errors
	for (std::size_t i = 0; i < std::min<std::size_t>(5, differences.size()); ++i)
	{
		// push similar vector into mostSimilarVectors
		mostSimilarVectors.push_back(differences[i].first);
	}

	*/

	// Loop through differences, extractive 5 hyperblocks and labels which are the
	// least distance to the chosen testing data
	for (std::size_t i = 0; i < std::min<std::size_t>(5, differences.size()); ++i)
	{
		const auto& vector = differences[i].vector;
		const auto& label = differences[i].label;
		mostSimilarVectors.push_back(vector);
		mostSimilarLabels.push_back(label);
	}
	// Insert chosen training vector to front of most similar vectors for analysis
	mostSimilarVectors.insert(mostSimilarVectors.begin(), testingDataPoint);
	mostSimilarLabels.insert(mostSimilarLabels.begin(), "unlabeled");

	std::vector<bool> analyzeClass(analyzeGlyphs.size(), false);
	analyzeClass[0] = true;

	std::vector<std::string> analyzeLabels(analyzeGlyphs.size(), "mal");
	analyzeLabels[0] = "ben";
	
	// Find average of all points in current hyperblock, to show as a representative glyph
	std::vector<bool>::iterator currHB = close.begin();

	// Copy data vector
	std::vector<std::vector<GLfloat>>::iterator dataIt = (allData.begin());
	// Copy class vector
	std::vector<bool>::iterator hbClassIt = classify.begin();

	sizeHB = 0;	// Reset global values
	for (auto& attr : averagePoint)
	{
		attr = 0.0;
	}

	// Determine size and composition of hyperblock
	int benCount = 0;
	int malCount = 0;
	// iterate through all data points
	for (currHB; currHB != close.end(); ++currHB)
	{
		// if current point is in hyperblock
		if (*currHB)
		{
			// Add class to class count
			if (*hbClassIt)
			{
				++benCount;
			}
			else if (!*hbClassIt)
			{
				++malCount;
			}

			++sizeHB;	// increment size counter
			// point is in current hyperblock
			// add point to sum of all points
			for (unsigned int i = 0; i <9; ++i)
			{
				averagePoint[i] += (*dataIt)[i];
			}
		}
		++hbClassIt;
		++dataIt;	// increment iterator
	}

	// identify count of dominant class of cluster
	float domCount = std::max(benCount, malCount);
	// itentify total count of datapoints in cluster
	float totalCount = benCount + malCount;
	// calculate percent purity of cluster
	float percentPure = (domCount / totalCount) * 100.0;

	// round label to three decimal places
	std::string pureLabel1 = std::to_string((int)percentPure);
	std::string pureLabel = pureLabel1.substr(0, 5);

	// Colored HB label vector
	std::string hbLabel = pureLabel + "%";
	std::string hbLabel2 = "n=" + std::to_string(sizeHB);

	// divide averagePoint by size counter to retrieve average of hyperblock
	for (auto& attribute : averagePoint)
	{
		attribute /= sizeHB;
	}

	// reset hyperblock average edge sums
	sumx1 = 0;
	sumy1 = 0;
	sumx2 = 0;
	sumy2 = 0;
	sumx3 = 0;
	sumy3 = 0;

	sizeHB = 6;
	std::vector<bool>::iterator threshold = close.begin();		// initialize threshold iterator
	std::vector<bool>::iterator classVec = classify.begin();	// initialize class iterator
	std::vector<int> tempSize = std::vector<int>(allData.size(), 1);
	std::vector<int>::iterator sizeVec = tempSize.begin();
	//std::vector<int>::iterator sizeVec = repsSize.begin();

	// ************************************* DISPLAY SPC-SF GRAPH ***************************************
	if (DISPLAY_SELECTOR)
	{
		int iteration = 1;
		// Display points within threshold
		// Pass through twice: Draw classes sequentially depending on CLASS_SEPERATION_MODE flag
		for (std::vector<std::vector<GLfloat>>::iterator it = (mostSimilarVectors.begin());
			it < mostSimilarVectors.end(); ++it)
		{
			// Display the point if DISPLAY_ALL flag is set,
			// or if it is in the threshold of the current point
			if (true )//|| *threshold == true)
			// Use exclusive or here?
			//if ((!*classVec && CLASS_SEPERATION_MODE) && (DISPLAY_ALL || *threshold == true))
			//{
			//if (true)
			{

				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();
				gluPerspective(0, float(SCREEN_WIDTH) / float(SCREEN_HEIGHT), 0.1, 100.0);

				drawLocatedGlyphs(&(*it), *classVec, *sizeVec, iteration, hbLabel, hbLabel2);
				++iteration;
			}
			// don't iterate past end of vectors
			if (it != mostSimilarVectors.end())
			{
				++sizeVec;
				++classVec;
				++threshold;
			}
		}
		/*
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
	
		*/
		maxY = 0.0;
		maxX = 0.0;
		minX = 10.0;
		minY = 10.0;

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

					float colors[6];
					// *********************** DRAW STICK FIGURE ***********************
					glyph.drawGlyphSF(pos2, pos3, stickFig.begin(), *classVec, *turt, DYNAMIC_ANGLES, POS_ANGLE,
						GLYPH_SCALE_FACTOR, SF_SEGMENT_CONSTANT, SF_ANGLE_SCALE, ANGLE_FOCUS, BIRD_FOCUS, colors);

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
		std::vector<std::vector<GLfloat>>::iterator repsIt = mostSimilarVectors.begin();
		std::vector<bool>::iterator classIt = repsClass.begin();
		//std::vector<std::vector<GLfloat>>::iterator repsIt = mixedHood.begin();
		//std::vector<bool>::iterator classIt = mixedClass.begin();
		std::vector<std::string>::iterator labelsIt = mostSimilarLabels.begin();
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
				if (repsIt != mostSimilarVectors.end())
				{
					// Draw polygon outlines
					glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

					// Pointer to current data point
					std::vector<GLfloat> processedData = *repsIt;
					
					std::vector<GLfloat> axesSPC{};
					axesSPC.push_back(*(processedData.begin() + 3));	// X1
					axesSPC.push_back(*(processedData.begin() + 2));	// Y1
					axesSPC.push_back(*(processedData.begin() + 3));	// X2
					axesSPC.push_back(*(processedData.begin() + 4));	// Y2
					axesSPC.push_back(*(processedData.begin() + 8));	// X3
					axesSPC.push_back(*(processedData.begin() + 7));	// Y3

					// Encode angles with most meaningful attributes
					std::vector<GLfloat> stickFig{};
					// Populate: CL (angle), UC (length), BN (angle), BC (length)
					stickFig.push_back(*(processedData.begin()));			// Angle 1
					stickFig.push_back(*(processedData.begin() + 1));	// Length 1
					stickFig.push_back(*(processedData.begin() + 5));	// Angle 2
					stickFig.push_back(*(processedData.begin() + 6));	// Length 2

					/*
					// Vector holding Stick Figure attributes
					std::vector<GLfloat> stickFig{};
					// Populate: BN (angle), UCShape (length), UCSize (angle), NN (length)
					stickFig.push_back(*(processedData.begin()));	// Angle 1
					stickFig.push_back(*(processedData.begin() + 2));	// Length 1
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
					*/

					float colors[6];
					colors[0] = *(processedData.begin() + 2);
					colors[1] = *(processedData.begin() + 1);
					colors[2] = *(processedData.begin() + 0);
					colors[3] = *(processedData.begin() + 2);
					colors[4] = *(processedData.begin() + 7);
					colors[5] = *(processedData.begin() + 6);

					glViewport(		// (rowNum x colNum)
						// Encode shift based off of first SPC axis horizontal/vertical shift.
						// Shift calculation:
						((SCREEN_WIDTH / NUM_COLUMNS)* (colNum)) + 5,
						((SCREEN_HEIGHT / NUM_ROWS)* (rowNum)),
						(SCREEN_WIDTH / NUM_COLUMNS),
						(SCREEN_HEIGHT / NUM_ROWS)
					);

					// Print class totals in each cell
					glPushMatrix();
					glMatrixMode(GL_MODELVIEW);
					glLoadIdentity();
					glColor3f(1.0, 0.0, 0.0);
					if (*classIt)
					{
						glColor3f(0.0, 0.0, 1.0);
					}
					else
					{
						glColor3f(1.0, 0.0, 0.0);
					}
					glRasterPos2f(-0.9, -0.9);
	
					//if (repCount < 0)
					//{
					//glDisable(GL_LIGHTING);
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
						//glEnable(GL_LIGHTING);
						glPopMatrix();
						++labelsIt;	// Increment label iterator
					//}
					
					glPushMatrix();		// Push new matrix
					glMatrixMode(GL_PROJECTION);
					glLineWidth(4.0);	// Line width = 4.0
					// Vertical translation factor
					//GLfloat glyphTransVert = -(1.0 - axesSPC[1]) + 1.0;;
					// GHorizaontal translateion factor
					//GLfloat glyphTransHort = -(1.0 - axesSPC[0]) + 0.1;
					GLfloat TRANS_DEF = 0.0;	// Default translation factor
					GLfloat SCALE_DEF = 1.0;	// Default scaling factor
					glPushMatrix();	// Push new matrix for glyph transformations
					// Translate glyph right a based on first SPC shift, and up a set value
					glScalef(1.5, 1.5, 1.5);
					glTranslatef(-0.4, 0.3, 0.0);
					//glTranslatef(glyphTransHort, glyphTransVert - 0.1, TRANS_DEF);
					GLfloat glyphScaleFactor = 2.0;

					glEnable(GL_BLEND);		// Enable blending.
					glDepthMask(GL_FALSE);	// Disable depth masking
					// *********************** DRAW STICK FIGURE ***********************
					//if (repCount < 2)
					//{	// Draw representative glyphs in color
						//glyph.drawGlyphSF(pos2, pos3, stickFig.begin(), *classIt, *turt, DYNAMIC_ANGLES, POS_ANGLE,
							//GLYPH_SCALE_FACTOR, SF_SEGMENT_CONSTANT, SF_ANGLE_SCALE, ANGLE_FOCUS, true);
					//}
					//else
					//{	// Draw other glyphs in grey, if selected
						glyph.drawGlyphSF(pos2, pos3, stickFig.begin(), *classIt, *turt, DYNAMIC_ANGLES, POS_ANGLE,
							GLYPH_SCALE_FACTOR, SF_SEGMENT_CONSTANT, SF_ANGLE_SCALE, ANGLE_FOCUS, BIRD_FOCUS, colors);
					//}

					// Disable blending and resume depth mask
					glDepthMask(GL_TRUE);
					glDisable(GL_BLEND);
					++repCount;
					// RESET THE CP AND CD
					turt->setCP(0.0, 0.0);

					// *********************** DRAW SPC AXES ***********************
					//if (DRAW_AXES)
					//{
						//glLineWidth(2.0);	// Line width = 2.0
						//glyph.drawAxesSPC(*pos2, *pos3, axesSPC.begin(),
							//AXIS_LENGTH, GLYPH_SCALE_FACTOR, DOTTED_AXES);
					//}

					glPopMatrix();	// Pop from gl matrix stack
					++colNum;		// Increment column index
				}

				if (repsIt == mostSimilarVectors.end())
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