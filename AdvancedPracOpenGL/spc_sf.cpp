// spc_sf.cpp : This file implements the hybrid SPC-SF data visualization method
//				described in 'Visual Knowledge Discovery and Machine Learning'
//				by Dr. Boris Kovalerchuk, Figure 13.2
// AUTHOR:	Nicholas Cutlip (Central Washington University Computer Science)
// DATA:	29 December, 2022

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include "GL/glut.h"
#include "turtleg.h"
#include "spc_sf.h"


int SCREEN_WIDTH;
int SCREEN_HEIGHT;

int VIEWPORT_SCALE = 8;
const unsigned int HEIGHT_SCALE = 10;		// 5 instead of 4: extra space on bottom
const unsigned int WIDTH_SCALE = 4;		// 6 instead of 5: extra space on bottom
const float MARGIN = 2.0;		// Margin between window border and visualization
const float MAR = 20.0;
const float INNER_MARGIN = 5.0;		// Margin between glyphs
float AXIS_LENGTH = 1.0;
bool LINES_ON = false;
bool DRAW_EDGES = true;
bool DRAW_AXES = true;
bool DISPLAY_ALL = false;

// Switch to int later for multiple classes?
bool CLASS_SEPERATION_MODE = true;

unsigned int DATA_SIZE = 699;
unsigned int DATA_INDEX = 0;

//float GLYPH_SCALE_FACTOR = 0.6;

// Draw breast-cancer-wisconsin hypercubes
/*
void fillHypercubesDT()
{
	// Enable transparency
	// (non-transparent objects should already be drawn)
	glEnable(GL_BLEND); //Enable blending.
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //Set blending function.
	// Draw hypercubes based on DT paper
	glColor4f(1.0, 0.0, 0.0, 0.5);
	glRectf((SCREEN_HEIGHT / 20), (SCREEN_WIDTH / 10), (SCREEN_HEIGHT / 5), (SCREEN_WIDTH / 5));
	
}
*/

void drawLocatedGlyphs(std::vector<GLfloat>* normalData, bool classify)
{
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

	// Arrange greater SPC position attributes from optimal positioning described in Worland, Wagle, and Kovalerchuk
	std::vector<GLfloat> position{ bn, uc, cl, bc, bn, mg };
	std::vector<GLfloat>::iterator positionIt = position.begin();

	// Glyph 1
	std::vector<GLfloat>::iterator axes1 = position.begin();			// First six attributes used for SPC (6/10)
	std::vector<GLfloat>::iterator stick1 = normalData->begin() + 6;	// Last four attributes used for SF  (4/10)

	// Glyph 2
	std::vector<GLfloat>::iterator axes2 = position.begin();			// First six attributes used for SPC (6/10)
	std::vector<GLfloat>::iterator stick2 = normalData->begin() + 6;	// Last four attributes used for SF  (4/10)

	// Glyph 3
	std::vector<GLfloat>::iterator axes3 = position.begin();			// First six attributes used for SPC (6/10)
	std::vector<GLfloat>::iterator stick3 = normalData->begin() + 6;	// Last four attributes used for SF  (4/10)

	// *********************** DRAW SF GLYPHS ***********************
	// Note that we are currently converting 'radians to degrees'
	// Draw Stick Figure (save vertex positions)

	// Construct glyph tool
	SpcSfGlyph glyph = SpcSfGlyph(0.0, 0.0);
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
	glyph.drawGlyphSF(pos2, pos3, stick1, classify, *turt);
	// RESET THE CP AND CD
	turt->setCP(0.0, 0.0);
	// *********************** DRAW SPC AXES ***********************
	if (DRAW_AXES)
	{
		glLineWidth(2.0);
		glyph.drawAxesSPC(*pos2, *pos3, axes1, AXIS_LENGTH);
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
	glyph.drawGlyphSF(pos2, pos3, stick2, classify, *turt);
	// RESET THE CP AND CD
	turt->setCP(0.0, 0.0);

	// *********************** DRAW SPC AXES ***********************
	if (DRAW_AXES)
	{
		glLineWidth(2.0);
		glyph.drawAxesSPC(*pos2, *pos3, axes2, AXIS_LENGTH);
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
	glyph.drawGlyphSF(pos2, pos3, stick3, classify, *turt);
	// RESET THE CP AND CD
	turt->setCP(0.0, 0.0);
	// *********************** DRAW SPC AXES ***********************
	if (DRAW_AXES)
	{
		glLineWidth(2.0);
		glyph.drawAxesSPC(*pos2, *pos3, axes3, AXIS_LENGTH);
	}

	glPopMatrix();
	glViewport(0.0, 0.0, SCREEN_WIDTH, SCREEN_HEIGHT);
	gluOrtho2D(0.0, SCREEN_WIDTH, 0.0, SCREEN_HEIGHT);
	glColor3f(0.0, 0.0, 0.0);

	// Draw SPC axis frame * 3
	  
	drawArrow(Point2(MAR, MAR), Point2(MAR, SCREEN_HEIGHT - 10), SCREEN_HEIGHT);
	drawArrow(Point2(MAR, MAR), Point2((SCREEN_WIDTH / 4) + (SCREEN_WIDTH / 24), MAR), SCREEN_HEIGHT);

	drawArrow(Point2(((SCREEN_WIDTH) / 3) + MAR, MAR), Point2((SCREEN_WIDTH / 3) + MAR, SCREEN_HEIGHT - 10), SCREEN_HEIGHT);
	drawArrow(Point2((SCREEN_WIDTH / 3) + MAR, MAR), Point2(((7 * SCREEN_WIDTH) / 12) + (SCREEN_WIDTH / 24), MAR), SCREEN_HEIGHT);

	drawArrow(Point2(((2 * SCREEN_WIDTH) / 3) + MAR, MAR), Point2(((2 * SCREEN_WIDTH) / 3) + MAR, SCREEN_HEIGHT - 10), SCREEN_HEIGHT);
	drawArrow(Point2(((2 * SCREEN_WIDTH) / 3) + MAR, MAR), Point2(((11 * SCREEN_WIDTH) / 12) + (SCREEN_WIDTH / 24), MAR), SCREEN_HEIGHT);

	// Draw edges between glyphs
	if (DRAW_EDGES)
	{
		if (classify)	glColor3f(0.0, 0.0, 1.0);
		else			glColor3f(1.0, 0.0, 0.0);
		glLineWidth(1.0);
		glTranslatef(-((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2)), -((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2)), 0);

		drawArrow(Point2(x1 + (SCREEN_WIDTH / VIEWPORT_SCALE), y1 + (SCREEN_WIDTH / VIEWPORT_SCALE)),
			Point2(x2 + (SCREEN_WIDTH / VIEWPORT_SCALE), y2 + (SCREEN_WIDTH / VIEWPORT_SCALE)), SCREEN_HEIGHT);

		drawArrow(Point2(x2 + (SCREEN_WIDTH / VIEWPORT_SCALE), y2 + (SCREEN_WIDTH / VIEWPORT_SCALE)),
			Point2(x3 + (SCREEN_WIDTH / VIEWPORT_SCALE), y3 + (SCREEN_WIDTH / VIEWPORT_SCALE)), SCREEN_HEIGHT);
	}

	glPopMatrix();
	glFlush();
}

// Return the Euclidean distance between two vectors
float euclideanDistance(std::vector<float>* vec1, std::vector<float>* vec2)
{
	float sum = 0;
	std::vector<float>::iterator it2 = vec2->begin();

	// Sum the squares of the differences of the matching attributes between vectors
	for (auto& attr : *vec1)
	{
		sum += pow(attr - *it2, 2);
	}

	return sqrt(sum);
}

// Identify if two vectors are withn a certain threshold of each other
bool isClose(std::vector<float>* vec1, std::vector<float>* vec2)
{
	std::vector<float>::iterator it2 = vec2->begin();

	// Check each pair of elements for difference within threshold value
	for (auto& attr : *vec1)
	{
		if (abs(attr - *it2++) > 2)
			return false;
	}

	// Return true if all attribtes are within threshold
	return true;
}

std::vector<bool> computeAllDistances(std::vector<GLfloat>* curr, std::vector<std::vector<GLfloat>>* data)
{
	// This variable is for temp extracted data to compare
	//std::vector<float> data0();
	std::vector<bool> close{};

	// Populate boolean vector to determine
	// which data points are within threshold of current point
	for (std::vector<std::vector<GLfloat>>::iterator iter = data->begin(); iter < data->end(); ++iter)
	{
		close.push_back(isClose(curr, &(*iter)));
	}

	// Close file streams
	//output.close();

	return close;
}

// Load parameters from config file
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

	//getline(myFile, line);
	// GLYPH_SIZE
	//std::istringstream sin4(line.substr(line.find("=") + 1));
	//sin4 >> GLYPH_SCALE_FACTOR;
}

void myIdle() {}

// Mouse input callback function
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

	// Up arrow response:
	// toggle drawing edges between glyphs
	else if (key == GLUT_KEY_UP)
	{
		DRAW_EDGES = !DRAW_EDGES;
	}

	// Down arrow response:
	// toggle drawing dashed SPC axes for glyphs
	else if (key == GLUT_KEY_DOWN)
	{
		DRAW_AXES = !DRAW_AXES;
	}

	// Redisplay with updated parameters
	glutPostRedisplay();
}

void myKeyboard(unsigned char key, int x, int y)
{


	// Redisplay with updated DATA_INDEX
	//glutPostRedisplay();
	//glutIdleFunc(myIdle);
}

// *********************** Display SPC-SF Hybrid Visualization ***********************
// 
// Data from UCI Machine Learning Repository
// breast-cancer-wisconsin.DATA
void myDisplay()
{
	//glEnable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND); //Enable blending.
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //Set blending function.
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(0, float(SCREEN_WIDTH) / float(SCREEN_HEIGHT), 0.1, 100.0);
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClearDepth(1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Read data file
	std::string line = "";
	std::ifstream myFile("breast-cancer-wisconsin.DATA");

	// Process data into a 2-D vector for ease of use
	std::vector<std::vector<GLfloat>> allData(DATA_SIZE);
	std::vector<bool> classify{};
	// Generalize this by using push_back and getline in test
	// Also populate vector of classification here

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

		std::vector<float> data(dataInt.begin(), dataInt.end());

		// Set color determined by class
		if (*--(data.end()) == 4) classify.push_back(false);
		else					  classify.push_back(true);

		// Remove labels from data
		data.erase(data.begin());
		data.erase(--data.end());

		// Duplicate data attributes x 4 to move from 9-D to 36-D
		data.insert(data.end(), data.begin(), data.end());
		data.insert(data.end(), data.begin(), data.end());

		allData[i] = data;
	}
	myFile.close();

	// Compute points within threshold
	std::vector<bool> close = computeAllDistances(&allData[DATA_INDEX], &allData);
	//std::vector<bool>::iterator threshold = close.begin();
	//std::vector<bool>::iterator classVec = classify.begin();

	int i = 0;
	for (auto& vec : allData)
	{
		const unsigned int scaleFactor = 10;
		std::vector<float> normalData;		// Normalize data to [0, 1]
		for (std::vector<float>::iterator iter = vec.begin(); iter < vec.end(); iter++)
		{
			normalData.push_back(*iter / scaleFactor);
		}
		allData[i] = normalData;
		++i;
	}

	std::vector<bool>::iterator threshold = close.begin();
	std::vector<bool>::iterator classVec = classify.begin();

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

	//glPopMatrix();
	glViewport(0.0, 0.0, SCREEN_WIDTH, SCREEN_HEIGHT);
	//gluOrtho2D(0.0, SCREEN_WIDTH, 0.0, SCREEN_HEIGHT);
	
	// Draw hypercubes based on DT paper
	glColor4f(1.0, 0.0, 0.0, 0.5);
	glRectf(0.0, 0.0, 0.5, 0.5);
	//glRectf((SCREEN_HEIGHT / 20), (SCREEN_WIDTH / 10), (SCREEN_HEIGHT / 5), (SCREEN_WIDTH / 5));
	//fillHypercubesDT();
	glPopMatrix();
	glFlush();
}

int main(int argc, char** argv)
{
	// Load parameters from config file
	loadConfig();

	// Implement config struct
	glutInit(&argc, argv);
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