#pragma once
#include "gl/glut.h"
#include "point2.h"
#include "vec2.h"


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> drawTileFrame >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void drawTileFrame(GLfloat margin, GLfloat height, GLfloat width, GLfloat heightScale, GLfloat widthScale)
{
    // (DISCOVERED: Outer brackets aren't necesary??)
    for (unsigned int rowNum = 0; rowNum < 4; rowNum++)		// Rows
    {
        for (unsigned int colNum = 0; colNum < 6; colNum++)	// Columns
        {
            glRectf(		// Draw black rectangle outline for tile frame
                margin,
                height - margin - (height / heightScale) * rowNum,
                (width / widthScale) * colNum,
                height - margin - (height / heightScale) * (rowNum + 1)
            );
        }
    }
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> drawArrow >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void drawArrow(Point2 start, Point2 end, const unsigned int winHeight)
{
    // Draw Line Segment
    glBegin(GL_LINE_STRIP);
    glVertex2i(start.getx(), start.gety());
    glVertex2i(end.getx(), end.gety());
    glEnd();

    // Draw trangle for arrow tip
    glBegin(GL_TRIANGLES);

    Vec2 p1 = { float(start.getx()), float(start.gety()) };
    Vec2 p2 = { float(end.getx()), float(end.gety()) };

    Vec2 direction = p2 - p1;
    direction.Normalize();

    Vec2 perpendicular = Vec2(direction);
    perpendicular.Perpendicu();

    Vec2 pt1 = p2 - direction * (winHeight / 80) * 3;
    Vec2 pt2 = pt1 - perpendicular * (winHeight / 80);
    Vec2 pt3 = pt1 + perpendicular * (winHeight / 80);

    glVertex2i(p2.getx(), p2.gety());
    glVertex2i(pt2.getx(), pt2.gety());
    glVertex2i(pt3.getx(), pt3.gety());

    glEnd();
};

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> drawDot >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void drawDot(GLint x, GLint y)
{
    glBegin(GL_POINTS);
    glPointSize(4.0);
    glVertex2i(x, y);
    glEnd();
};