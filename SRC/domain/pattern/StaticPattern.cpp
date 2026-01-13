#include "StaticPattern.h"

StaticPattern::StaticPattern(int tag, double fact)
  :LoadPattern(tag, PATTERN_TAG_StaticPattern, fact)
{

}

StaticPattern::~StaticPattern()
{
}

