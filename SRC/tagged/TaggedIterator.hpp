#pragma once 
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>

template <typename T, typename StorageType=TaggedObjectStorage>
class TaggedIterator {
  public:
    // TaggedIterator(StorageType* storage) : mIter(storage->getIterRef()) {}
    TaggedIterator(StorageType* storage) : mIter(storage->getIter()) {}
    
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
    // StorageType::Iterator &mIter;
    typename StorageType::Iterator mIter;
};