/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// DEPRECATED
//
// Written: fmk 
// Created: 10/98
// Revision: A
//
// Description: This file contains the class interface for Renderer.
// Renderer is an abstract base class. An Renderer object is used
// to create an image of the domain.
//
#include "Renderer.h"
#include "ColorMap.h"
#include <string.h>

#include <Matrix.h>
#include <Vector.h>


Renderer::Renderer(ColorMap &_theMap)
  :theMap(&_theMap)
{

}


Renderer::Renderer(const char *title, ColorMap &_theMap)
  :theMap(&_theMap)
{
}

Renderer::~Renderer()
{

}

int
Renderer::saveImage(const char *fileName)
{
  return 0;
}


int
Renderer::saveImage(const char *rendererTitle, const char *fileName)
{
  return 0;
}

int
Renderer::drawVector(const Vector &position, const Vector &value, double factor, int tag)
{
  return 0;
}


void
Renderer::setColorMap(ColorMap &map)
{
  
}

int 
Renderer::drawCube(const Matrix &points, const Vector &values, int tag, int mode)
{
  return -1;
}

