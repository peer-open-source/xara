//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
//===----------------------------------------------------------------------===//
//

#include <StringStream.h>

#include <Vector.h>

#include <cmath>
#include <iomanip>

StringStream::StringStream(int indent)
  : OPS_Stream(OPS_STREAM_TAGS_StandardStream),
    nextWidth(0),
    indentSize(indent),
    numIndent(-1)
{
  if (indentSize < 1)
    indentSize = 1;

  indentString.assign(indentSize, ' ');
}

StringStream::StringStream(const std::string &value, int indent)
  : OPS_Stream(OPS_STREAM_TAGS_StandardStream),
    output(value),
    nextWidth(0),
    indentSize(indent),
    numIndent(-1)
{
  if (indentSize < 1)
    indentSize = 1;

  indentString.assign(indentSize, ' ');
}

StringStream::~StringStream()
{
}

int
StringStream::setFile(const char *fileName, openMode mode, bool echo)
{
  return 0;
}

int
StringStream::setPrecision(int precision)
{
  format << std::setprecision(precision);
  return 0;
}

int
StringStream::setFloatField(OPS_Stream::Float field)
{
  if (field == OPS_Stream::Float::Fixed)
    format << std::fixed;
  else if (field == OPS_Stream::Float::Scientific)
    format << std::scientific;

  return 0;
}

int
StringStream::precision(int precision)
{
  return this->setPrecision(precision);
}

int
StringStream::width(int width)
{
  nextWidth = width;
  return 0;
}

int
StringStream::flush()
{
  return 0;
}

const std::string &
StringStream::string() const
{
  return output;
}

std::string &
StringStream::string()
{
  return output;
}

std::string
StringStream::str() const
{
  return output;
}

void
StringStream::str(const std::string &value)
{
  output = value;
}

void
StringStream::clear()
{
  output.clear();
}

int
StringStream::tag(const char *tagName)
{
  this->indent();
  (*this) << tagName << "\n";

  numIndent++;

  return 0;
}

int
StringStream::tag(const char *tagName, const char *value)
{
  this->indent();
  (*this) << tagName << " " << value << "\n";

  numIndent++;

  return 0;
}

int
StringStream::endTag()
{
  numIndent--;
  return 0;
}

int
StringStream::attr(const char *name, int value)
{
  this->indent();
  (*this) << name << " = " << value << "\n";

  return 0;
}

int
StringStream::attr(const char *name, double value)
{
  this->indent();
  (*this) << name << " = " << value << "\n";

  return 0;
}

int
StringStream::attr(const char *name, const char *value)
{
  this->indent();
  (*this) << name << " = " << value << "\n";

  return 0;
}

int
StringStream::write(Vector &data)
{
  this->indent();
  (*this) << data;

  return 0;
}

OPS_Stream&
StringStream::write(const char *s, int n)
{
  if (s != nullptr && n > 0)
    output.append(s, n);

  return *this;
}

OPS_Stream&
StringStream::write(const unsigned char *s, int n)
{
  if (s != nullptr && n > 0)
    output.append(reinterpret_cast<const char *>(s), n);

  return *this;
}

OPS_Stream&
StringStream::write(const signed char *s, int n)
{
  if (s != nullptr && n > 0)
    output.append(reinterpret_cast<const char *>(s), n);

  return *this;
}

OPS_Stream&
StringStream::write(const void *s, int n)
{
  if (s != nullptr && n > 0)
    output.append(reinterpret_cast<const char *>(s), n);

  return *this;
}

OPS_Stream&
StringStream::write(const double *s, int n)
{
  if (n == 0)
    return *this;

  for (int i = 0; i < n; i++)
    (*this) << s[i] << " ";

  (*this) << "\n";

  return *this;
}

OPS_Stream&
StringStream::operator<<(char c)
{
  appendFormatted(c);

  return *this;
}

OPS_Stream&
StringStream::operator<<(unsigned char c)
{
  appendFormatted(c);

  return *this;
}

OPS_Stream&
StringStream::operator<<(signed char c)
{
  appendFormatted(c);

  return *this;
}

OPS_Stream&
StringStream::operator<<(const char *s)
{
  if (s != nullptr)
    appendFormatted(s);

  return *this;
}

OPS_Stream&
StringStream::operator<<(const unsigned char *s)
{
  if (s != nullptr)
    appendFormatted(s);

  return *this;
}

OPS_Stream&
StringStream::operator<<(const signed char *s)
{
  if (s != nullptr)
    appendFormatted(s);

  return *this;
}

OPS_Stream&
StringStream::operator<<(const void *p)
{
  appendFormatted(p);

  return *this;
}

OPS_Stream&
StringStream::operator<<(int n)
{
  appendFormatted(n);

  return *this;
}

OPS_Stream&
StringStream::operator<<(unsigned int n)
{
  appendFormatted(n);

  return *this;
}

OPS_Stream&
StringStream::operator<<(long n)
{
  appendFormatted(n);

  return *this;
}

OPS_Stream&
StringStream::operator<<(unsigned long n)
{
  appendFormatted(n);

  return *this;
}

OPS_Stream&
StringStream::operator<<(short n)
{
  appendFormatted(n);

  return *this;
}

OPS_Stream&
StringStream::operator<<(unsigned short n)
{
  appendFormatted(n);

  return *this;
}

OPS_Stream&
StringStream::operator<<(bool b)
{
  appendFormatted(b);

  return *this;
}

OPS_Stream&
StringStream::operator<<(double n)
{
  if (std::isnan(n))
    appendFormatted("NaN");
  else if (std::isinf(n))
    appendFormatted("Inf");
  else
    appendFormatted(n);

  return *this;
}

OPS_Stream&
StringStream::operator<<(float n)
{
  appendFormatted(n);

  return *this;
}

OPS_Stream&
StringStream::operator<<(std::string const &s)
{
  appendFormatted(s);

  return *this;
}

int
StringStream::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int
StringStream::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}

void
StringStream::indent()
{
  for (int i = 0; i < numIndent; i++)
    output += indentString;
}
