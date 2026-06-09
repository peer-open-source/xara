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

#include <iostream>

int        Renderer::numRenderers(0);
char     **Renderer::theTitles =0;
Renderer **Renderer::theRenderers =0;

Renderer::Renderer(ColorMap &_theMap)
  :theMap(&_theMap)
{

}


Renderer::Renderer(const char *title, ColorMap &_theMap)
  :theMap(&_theMap)
{
  int loc = -1;

  // look for an empty slot
  for (int i=0; i<numRenderers; i++)
    if (theRenderers[i] == 0) {
      loc = i;
      i = numRenderers;
    }

  // if no space or not already there add
  if (loc == -1) {
    Renderer **theNewRenderers = new Renderer *[numRenderers+1];
    char **theNewTitles = new char *[numRenderers+1];

    for (int i=0; i<numRenderers; i++) {
      theNewRenderers[i] = theRenderers[i];
      theNewTitles[i] = theTitles[i];
    }

    loc = numRenderers;
    numRenderers++;
    
    if (theRenderers != 0) 
      delete [] theRenderers;
    if (theTitles != 0)
      delete [] theTitles;

    theRenderers = theNewRenderers;
    theTitles = theNewTitles;
  }

  // set this in current slot
  theRenderers[loc] = this;
  char *titleCopy = new char [strlen(title)+1];
  strcpy(titleCopy, title);
  theTitles[loc] = titleCopy;
}

Renderer::~Renderer()
{
  for (int i=0; i<numRenderers; i++)
    if (theRenderers[i] == this) {
      theRenderers[i] = 0;
      delete [] (theTitles[i]);
      theTitles[i] = 0;
    }
}

int
Renderer::saveImage(const char *fileName)
{
  return 0;
}


int
Renderer::saveImage(const char *rendererTitle, const char *fileName)
{
  for (int i=0; i<numRenderers; i++)
    if (theRenderers[i] != 0) 
      if (strcmp(rendererTitle, theTitles[i]) == 0)
	return theRenderers[i]->saveImage(fileName);

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
  theMap = &map;
}

int 
Renderer::drawCube(const Matrix &points, const Vector &values, int tag, int mode)
{
  return -1;
}

