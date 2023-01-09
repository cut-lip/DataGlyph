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
int VIEWPORT_SCALE = 12;
const unsigned int HEIGHT_SCALE = 5;		// 5 instead of 4: extra space on bottom
const unsigned int WIDTH_SCALE = 6;		// 6 instead of 5: extra space on bottom
const float MARGIN = 2.0;		// Margin between window border and visualization
const float INNER_MARGIN = 5.0;		// Margin between glyphs
float AXIS_LENGTH = 1.0;
const float MAR = 20.0;

//float GLYPH_SCALE_FACTOR = 0.6;


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

	// Locate the viewport for glyph drawing
	GLfloat x1 = ((SCREEN_WIDTH / 3) * *positionIt) - ((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2));
	GLfloat y1 = (SCREEN_HEIGHT * *++positionIt) - ((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2));

	glPushMatrix();

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
	glLineWidth(2.0);
	glyph.drawAxesSPC(*pos2, *pos3, axes1, AXIS_LENGTH);

	// Locate the viewport for glyph drawing
	GLfloat x2 = (SCREEN_WIDTH / 3) + ((SCREEN_WIDTH / 3) * *++positionIt) - ((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2));
	GLfloat y2 = (SCREEN_HEIGHT * *++positionIt) - ((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2));

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
	glLineWidth(2.0);
	glyph.drawAxesSPC(*pos2, *pos3, axes2, AXIS_LENGTH);

	// Locate the viewport for glyph drawing
	GLfloat x3 = ((SCREEN_WIDTH / 3) * 2) + ((SCREEN_WIDTH / 3) * *++positionIt) - ((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2));
	GLfloat y3 = (SCREEN_HEIGHT * *++positionIt) - ((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE*2));

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
	glLineWidth(2.0);
	glyph.drawAxesSPC(*pos2, *pos3, axes3, AXIS_LENGTH);

	glPopMatrix();
	glViewport(0.0, 0.0, SCREEN_WIDTH, SCREEN_HEIGHT);
	gluOrtho2D(0.0, SCREEN_WIDTH, 0.0, SCREEN_HEIGHT);
	glColor3f(0.0, 0.0, 0.0);

	// Draw SPC axis frame * 3
	drawArrow(Point2(MAR, MAR), Point2(MAR, SCREEN_HEIGHT - 10), SCREEN_HEIGHT);
	drawArrow(Point2(MAR, MAR), Point2((SCREEN_WIDTH / 3) - MAR, MAR), SCREEN_HEIGHT);

	drawArrow(Point2((SCREEN_WIDTH / 3) + MAR, MAR), Point2((SCREEN_WIDTH / 3) + MAR, SCREEN_HEIGHT - 10), SCREEN_HEIGHT);
	drawArrow(Point2((SCREEN_WIDTH / 3) + MAR, MAR), Point2(((SCREEN_WIDTH / 3) * 2) - MAR, MAR), SCREEN_HEIGHT);

	drawArrow(Point2(((SCREEN_WIDTH / 3) * 2) + MAR, MAR), Point2(((SCREEN_WIDTH / 3) * 2) + MAR, SCREEN_HEIGHT - 10), SCREEN_HEIGHT);
	drawArrow(Point2(((SCREEN_WIDTH / 3) * 2) + MAR, MAR), Point2(SCREEN_WIDTH - MAR, MAR), SCREEN_HEIGHT);

	// Draw edges between glyphs
	if (classify)	glColor3f(0.0, 0.0, 1.0);
	else			glColor3f(1.0, 0.0, 0.0);
	glLineWidth(1.0);
	glTranslatef(-((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2)), -((GLint)SCREEN_WIDTH / (VIEWPORT_SCALE * 2)), 0);

	drawArrow(Point2(x1 + (SCREEN_WIDTH / VIEWPORT_SCALE), y1 + (SCREEN_WIDTH / VIEWPORT_SCALE)),
		Point2(x2 + (SCREEN_WIDTH / VIEWPORT_SCALE), y2 + (SCREEN_WIDTH / VIEWPORT_SCALE)), SCREEN_HEIGHT);

	drawArrow(Point2(x2 + (SCREEN_WIDTH / VIEWPORT_SCALE), y2 + (SCREEN_WIDTH / VIEWPORT_SCALE)),
		Point2(x3 + (SCREEN_WIDTH / VIEWPORT_SCALE), y3 + (SCREEN_WIDTH / VIEWPORT_SCALE)), SCREEN_HEIGHT);

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
		if (abs(attr - *it2++) > 1)
			return false;
	}

	// Return true if all attribtes are within threshold
	return true;
}

std::vector<float> computeAllDistances()
{

	// Open data file
	std::string line = "";
	std::ifstream myFile("breast-cancer-wisconsin0.DATA");

	// Create output file
	std::ofstream output;
	output.open("close.txt");

	// Extract first data point
	getline(myFile, line);

	// Split string into a vector
	std::vector<int> dataInt0;
	std::stringstream ss0(line);

	while (ss0.good()) {
		std::string substr;
		getline(ss0, substr, ',');
		dataInt0.push_back(stoi(substr));
		continue;
	}

	std::vector<float> data0(dataInt0.begin(), dataInt0.end());

	// Remove labels from data
	data0.erase(data0.begin());
	data0.erase(--data0.end());
	// Duplicate ninth attribute
	data0.push_back(*--data0.end());

	if (myFile.is_open())
	{
		while (getline(myFile, line))
		{
			int vecCount = 10;

			// Split string into a vector
			std::vector<int> dataInt;
			std::stringstream ss(line);

			while (ss.good()) {
				std::string substr;
				getline(ss, substr, ',');

				if (substr == "?") substr = std::to_string(*--dataInt.end());
				dataInt.push_back(stoi(substr));
				continue;
			}

			std::vector<float> data(dataInt.begin(), dataInt.end());

			// Remove labels from data
			data.erase(data.begin());
			data.erase(--data.end());
			// Duplicate ninth attribute
			data.push_back(*--data.end());

			// Determine if points are within threshold distance of each other

			if (isClose(&data0, &data))
			{
				output << "T" << "\n";
			}
			else
			{
				output << "F" << "\n";
			}

			// Push distance to output vector and write to output file
			// distances.push_back(dist);
		}
	}

	// Close file streams
	myFile.close();
	output.close();

	return std::vector<float>{};
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

// *********************** Display SPC-SF Hybrid Visualization ***********************
// 
// Data from UCI Machine Learning Repository
// breast-cancer-wisconsin.DATA
void myDisplay()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(0, float(SCREEN_WIDTH) / float(SCREEN_HEIGHT), 0.1, 100.0);

	// Clear background to white
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);

	// Determine points within threshold of current point
	computeAllDistances();

	// Read data file
	std::string line = "";
	std::ifstream myFile("breast-cancer-wisconsin.DATA");
	std::ifstream thresholds("close.txt");

	// Process data into a 2-D vector for ease of use
	std::vector<std::vector<GLfloat>>{};





	// Graph first point
	getline(myFile, line);

	int vecCount = 10;
	const unsigned int scaleFactor = 10;

	// Split string into a vector
	std::vector<int> dataInt;
	std::stringstream ss(line);

	// THIS SHOULD BE DIFFERENT FOR FIRST CASE
	while (ss.good()) {		// Replace "?" in data
		std::string substr;
		getline(ss, substr, ',');

		if (substr == "?") substr = "" + (*--dataInt.end());
		dataInt.push_back(stoi(substr));
		continue;
	}

	std::vector<float> data(dataInt.begin(), dataInt.end());

	// Remove labels from data
	data.erase(data.begin());
	data.erase(--data.end());

	// Duplicate data attributes x 4 to move from 9-D to 36-D
	data.insert(data.end(), data.begin(), data.end());
	data.insert(data.end(), data.begin(), data.end());

	std::vector<float> normalData;		// Normalize data to [0, 1]
	for (std::vector<float>::iterator iter = data.begin(); iter < data.end(); iter++)
	{
		normalData.push_back(*iter / scaleFactor);
	}

	// Draw LG (Located Glyphs)
	drawLocatedGlyphs(&normalData, true);

	if (myFile.is_open())
	{
		while (getline(myFile, line))
		{		
			std::string s = "";
			getline(thresholds, s);
			//if (s.find("T") != -1)
			//{
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluPerspective(0, float(SCREEN_WIDTH) / float(SCREEN_HEIGHT), 0.1, 100.0);

			int vecCount = 10;

			// Split string into a vector
			std::vector<int> dataInt;
			std::stringstream ss(line);

			while (ss.good()) {		// Replace "?" in data
				std::string substr;
				getline(ss, substr, ',');

				if (substr == "?") substr = std::to_string(*dataInt.begin());

				dataInt.push_back(stoi(substr));
				continue;
			}

			std::vector<float> data(dataInt.begin(), dataInt.end());

			// Set color based on tumor classification
			bool benign = true;
			// Set color determined by class
			if (*--(data.end()) == 4)
			{
				benign = false;
			}

			// Remove labels from data
			data.erase(data.begin());
			data.erase(--data.end());
			// Duplicate ninth attribute
			//data.push_back(*--data.end());

			// Duplicate data attributes to move from 10-D to 36-D
			data.insert(data.end(), data.begin(), data.end());
			data.insert(data.end(), data.begin(), data.end());

			std::vector<float> normalData;		// Normalize data to [0, 1]
			for (std::vector<float>::iterator iter = data.begin(); iter < data.end(); iter++)
			{
				normalData.push_back(*iter / scaleFactor);
			}

			drawLocatedGlyphs(&normalData, benign);
			//}
		}
	}
}

void myKeyboard(unsigned char key, int x, int y)
{
	// Move forwards in the list of data points
	if (key == GLUT_KEY_RIGHT)
	{
		glClearColor(1.0, 1.0, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT);

		computeAllDistances();

		// Redisplay the visualization with the next glyph and its neightbors
	}

	// Move backwards in the list of data points
	if (key == GLUT_KEY_LEFT)
	{
		glClearColor(1.0, 1.0, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT);

		computeAllDistances();

		// Redisplay the visualization with the last glyph and its neightbors
	}
}

int main(int argc, char** argv)
{
	// Load parameters from config file
	loadConfig();

	// Implement config struct
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(SCREEN_WIDTH, SCREEN_HEIGHT);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("SPC-SF");

	glutDisplayFunc(myDisplay);
	glutKeyboardFunc(myKeyboard);
	glutMainLoop();
}