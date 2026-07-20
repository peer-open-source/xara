//
// TODO: add OPS_GetDoubleListInput here
//
#include <vector> 
#include <string>

bool
string_to_double(const std::string& text, double& num) 
{
  num = 0.0;
  try {
    num = std::stod(text);
    return true;
  }
  catch (...) {
    return false;
  }
}


// Adapted from ASDEA's string_to_list_of_doubles
bool 
SplitFloatList(const std::string& text, char sep, std::vector<double>& out) 
{
  if (out.size() > 0)
    out.clear();

  std::size_t start = 0, end = 0;

  double value;
  while (true) {
    end = text.find(sep, start);
    if (end == std::string::npos) {
      if (start < text.size()) {
        if (!string_to_double(text.substr(start), value))
          return false;
        out.push_back(value);
      }
      break;
    }
    std::string subs = text.substr(start, end - start);
    if (subs.size() > 0) {
      if (!string_to_double(subs, value))
        return false;
      out.push_back(value);
    }
    start = end + 1;
  }
  return true;
}
