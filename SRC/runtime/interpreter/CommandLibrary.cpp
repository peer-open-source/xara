#include <string>
#include <algorithm>

namespace  Xara::Internal {

std::string 
toLower(const std::string & s)
{
  std::string copy = s;
  transform( copy.begin( ), copy.end( ), copy.begin( ), 
      [](unsigned char c) { return std::tolower(c); });
  return copy;
}

bool 
equalsIgnoreCase( const std::string & lhs, const std::string & rhs )
{
  return toLower( lhs ) == toLower( rhs );
}

}