# Contributing

## C++ Style

- DO NOT USE `exit()`. There is no reason to call this function. Validate all arguments to an object's constructor before calling it. Validation should not happen inside constructors.
- Dont use `void` as an argument type. This is good practice in C, but bad in C++ ([source](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#nl25-dont-use-void-as-an-argument-type))
- Dont write redundant comments. For example, dont do:
  ```c++
  class MyClass {
    // revert to last commit
    int revertToLastCommit();
  };
  ```
- Minimize use of static local variables. If the size of a vector or container is known at compile time, use a compile-time sized container that will be allocated on the stack instead of a statically-allocating a variable on the heap.
- All functions begining with `OPS_` are deprecated. Dont use them.
- Dont statically allocate global containers.
- If in doubt, refer to the official [C++ core guidelines](https://isocpp.github.io/CppCoreGuidelines/)
