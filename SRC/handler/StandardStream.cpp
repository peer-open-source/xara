/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.8 $
// $Date: 2009-01-09 00:55:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/StandardStream.cpp,v $

#include <StandardStream.h>
#include <Vector.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <string.h>

#define OPS_CONSOLE std::cout
using std::ios;

StandardStream::StandardStream(int indent, bool echo)
  :OPS_Stream(OPS_STREAM_TAGS_FileStream), 
   fileOpen(0), echoApplication(echo),  indentSize(indent), numIndent(-1)
{
  if (indentSize < 1)
    indentSize = 1;

  indentString = new char[indentSize+1];
  for (int i=0; i<indentSize; i++)
    strcpy(indentString, " ");

}
StandardStream::~StandardStream()
{
  if (fileOpen == 1)
    theFile.close();
  if (indentString != nullptr)
    delete [] indentString;
}

int 
StandardStream::setFile(const char *fileName, openMode mode, bool echo)
{
  if (fileOpen == 1) {
    theFile.close();
    fileOpen = 0;
  }

  if (mode == openMode::OVERWRITE) 
    theFile.open(fileName, ios::out);
  else
    theFile.open(fileName, ios::out| ios::app);

  if (theFile.bad()) {
    OPS_CONSOLE << "WARNING - StandardStream::setFile()";
    OPS_CONSOLE << " - could not open file " << fileName << std::endl;

    return -1;
  } else
    fileOpen = 1;

  echoApplication = echo;

  return 0;
}


int 
StandardStream::setPrecision(int prec)
{
  OPS_CONSOLE << std::setprecision(prec);

  if (fileOpen != false)
    theFile << std::setprecision(prec);

  return 0;
}

int 
StandardStream::setFloatField(OPS_Stream::Float field)
{
#ifndef _WIN32
  if (field == OPS_Stream::Float::Fixed) {
      OPS_CONSOLE << std::setiosflags(ios::fixed);
    if (fileOpen != false)
      theFile << std::setiosflags(ios::fixed);
  }
  else if (field == OPS_Stream::Float::Scientific) {
    OPS_CONSOLE << std::setiosflags(ios::scientific);
    if (fileOpen != 0)
      theFile << std::setiosflags(ios::scientific);
  }
#endif
  return 0;
}


int 
StandardStream::tag(const char *tagName)
{
  // output the xml for it to the file
  this->indent();
  (*this) << tagName << "\n";

  numIndent++;

  return 0;
}

int
StandardStream::tag(const char *tagName, const char *value)
{
  // output the xml for it to the file
  this->indent();
  (*this) << tagName << " " << value << "\n";


  numIndent++;

  return 0;
}

int 
StandardStream::endTag()
{
  numIndent--;

  return 0;
}

int 
StandardStream::attr(const char *name, int value)
{
  this->indent();
  (*this) << name << " = " << value << "\n";
  
  return 0;
}

int 
StandardStream::attr(const char *name, double value)
{
  this->indent();
  (*this) << name << " = " << value << "\n";

  return 0;
}

int 
StandardStream::attr(const char *name, const char *value)
{
  this->indent();
  (*this) << name << " = " << value << "\n";

  return 0;
}

int 
StandardStream::write(Vector &data)
{
  this->indent();
  (*this) << data;  

  return 0;
}


OPS_Stream& 
StandardStream::write(const char *s,int n)
{
  if (echoApplication == true)
    OPS_CONSOLE.write(s, n);

  if (fileOpen != 0)
    theFile.write(s, n);
  
  return *this;
}

OPS_Stream& 
StandardStream::write(const unsigned char*s, int n)
{
  if (echoApplication == true)
    OPS_CONSOLE.write((const char *) s, n);

  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
StandardStream::write(const signed char*s, int n)
{
  if (echoApplication == true)
    OPS_CONSOLE.write((const char *)s, n);

  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
StandardStream::write(const void *s, int n)
{
  if (echoApplication == true)
    OPS_CONSOLE.write((const char *)s, n);

  if (fileOpen != 0)
   theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
StandardStream::operator<<(char c)
{
  if (echoApplication == true)
    OPS_CONSOLE << c;

  if (fileOpen != false)
    theFile << c;

 return *this;
}
OPS_Stream& 
StandardStream::operator<<(unsigned char c)
{
  if (echoApplication == true)
    OPS_CONSOLE << c;

  if (fileOpen != 0)
    theFile << c;

 return *this;
}
OPS_Stream& 
StandardStream::operator<<(signed char c)
{
  if (echoApplication == true)
    OPS_CONSOLE << c;

  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
StandardStream::operator<<(const char *s)
{
  // note that we do the flush so that a "/n" before
  // a crash will cause a flush() - similar to what 
  if (echoApplication == true) {
    OPS_CONSOLE << s;
    OPS_CONSOLE.flush();
  }

  if (fileOpen != 0) {
    theFile << s;
    theFile.flush();
  }

  return *this;
}

OPS_Stream& 
StandardStream::operator<<(const unsigned char *s)
{
  if (echoApplication == true)
    OPS_CONSOLE << s;

  if (fileOpen != 0)
    theFile << s;

  return *this;
}


OPS_Stream& 
StandardStream::operator<<(const signed char *s)
{
  if (echoApplication == true)
    OPS_CONSOLE << s;

  if (fileOpen != 0)
    theFile << s;

  return *this;
}

OPS_Stream& 
StandardStream::operator<<(const void *p)
{
/*
//  cerr << p;

  if (fileOpen != 0)
    theFile << p;
*/
  return *this;
}
OPS_Stream& 
StandardStream::operator<<(int n)
{
  if (echoApplication == true)
    OPS_CONSOLE <<  n;

  if (fileOpen != 0)
    theFile << n;

  return *this;
}

OPS_Stream& 
StandardStream::operator<<(unsigned int n)
{
  if (echoApplication == true)
    OPS_CONSOLE << 1.0*n;

  if (fileOpen != 0)
    theFile << 1.0*n;

  return *this;
}
OPS_Stream& 
StandardStream::operator<<(long n)
{
/*
OPS_CONSOLE << n;

if (fileOpen != 0)
  theFile << n;
*/
  return *this;
}
OPS_Stream& 
StandardStream::operator<<(unsigned long n)
{
/*
  OPS_CONSOLE << n;

  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
StandardStream::operator<<(short n)
{
/*
  OPS_CONSOLE << n;

  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
StandardStream::operator<<(unsigned short n)
{
/*
  OPS_CONSOLE << n;

  if (fileOpen != 0)
    theFile << n;
*/
return *this;
}

OPS_Stream& 
StandardStream::operator<<(bool b)
{
/*
  OPS_CONSOLE << b;

  if (fileOpen != 0)
    theFile << b;
*/
 return *this;
}

OPS_Stream& 
StandardStream::operator<<(double n)
{
  if (echoApplication == true)
    OPS_CONSOLE << n;

  if (fileOpen != 0) {
    if (std::isnan(n)) {
      theFile << "NaN";
    } else if (std::isinf(n)) {
      theFile << "Inf";
    } else {
      theFile << n;
    }
  }

 return *this;
}

OPS_Stream& 
StandardStream::operator<<(float n)
{
  if (echoApplication == true)
    OPS_CONSOLE << n;

  if (fileOpen != 0)
    theFile << n;

  return *this;
}

OPS_Stream&
StandardStream::operator<<(std::string const&s)
{
  if (echoApplication == true)
    OPS_CONSOLE << s;

  if (fileOpen != 0)
    theFile << s;

  return *this;
}

void
StandardStream::indent()
{
  for (int i=0; i<numIndent; i++) {
    OPS_CONSOLE << indentString;
    if (fileOpen != 0)
      theFile << indentString;
  }
}

int 
StandardStream::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int 
StandardStream::recvSelf(int commitTag, Channel &theChannel, 
             FEM_ObjectBroker &theBroker)
{
  return 0;
}

int StandardStream::flush() {
  if (theFile.is_open() && theFile.good()) {
    theFile.flush();
  }
  return 0;
}
