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
#pragma once

#include <OPS_Stream.h>

#include <ios>
#include <sstream>
#include <string>

class StringStream : public OPS_Stream
{
 public:
  StringStream(int indentSize=2);
  explicit StringStream(const std::string &value, int indentSize=2);
  ~StringStream();

  int setFile(const char *fileName, openMode mode = openMode::OVERWRITE, bool echo=false);
  int setPrecision(int precision);
  int setFloatField(OPS_Stream::Float);
  int precision(int precision);
  int width(int width);
  int flush();

  const std::string &string() const;
  std::string &string();
  std::string str() const;
  void str(const std::string &value);
  void clear();

  // xml stuff
  int tag(const char *);
  int tag(const char *, const char *);
  int endTag();
  int attr(const char *name, int value);
  int attr(const char *name, double value);
  int attr(const char *name, const char *value);
  int write(Vector &data);

  OPS_Stream& write(const char *s, int n);
  OPS_Stream& write(const unsigned char *s, int n);
  OPS_Stream& write(const signed char *s, int n);
  OPS_Stream& write(const void *s, int n);
  OPS_Stream& write(const double *s, int n);
  OPS_Stream& operator<<(char c);
  OPS_Stream& operator<<(unsigned char c);
  OPS_Stream& operator<<(signed char c);
  OPS_Stream& operator<<(const char *s);
  OPS_Stream& operator<<(const unsigned char *s);
  OPS_Stream& operator<<(const signed char *s);
  OPS_Stream& operator<<(const void *p);
  OPS_Stream& operator<<(int n);
  OPS_Stream& operator<<(unsigned int n);
  OPS_Stream& operator<<(long n);
  OPS_Stream& operator<<(unsigned long n);
  OPS_Stream& operator<<(short n);
  OPS_Stream& operator<<(unsigned short n);
  OPS_Stream& operator<<(bool b);
  OPS_Stream& operator<<(double n);
  OPS_Stream& operator<<(float n);
  OPS_Stream& operator<<(std::string const&);

  int sendSelf(int commitTag, Channel &);
  int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);

 private:
  template <class T>
  void appendFormatted(const T &value)
  {
    std::ostringstream stream;
    stream.copyfmt(format);

    if (nextWidth > 0)
      stream.width(nextWidth);

    stream << value;
    output += stream.str();
    nextWidth = 0;
  }

  void indent();

  std::string output;
  std::ostringstream format;
  int nextWidth;
  int indentSize;
  int numIndent;
  std::string indentString;
};

