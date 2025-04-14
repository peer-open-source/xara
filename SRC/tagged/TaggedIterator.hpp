#pragma once 
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>

template <typename T>
class TaggedIterator {
  public:
    TaggedIterator(TaggedObjectStorage* storage) : mIter(storage->getComponents()) {}
    
    void reset() {
        mIter.reset();
    }

    T* operator()() {
        TaggedObject* item = mIter();
        if (item == nullptr)
          return nullptr;
        return (T*)item;
    }

  private:
    TaggedObjectIter &mIter;
};