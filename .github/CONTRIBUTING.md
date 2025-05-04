# Contributing

## C++ Style

- DO NOT USE `exit()`. There is no reason to call this function. Validate all arguments to an object's constructor before calling it. Validation should not happen inside constructors.
- Pull requests that increase the number of compiler warnings will not be accepted.
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
- Declare variables at first use, dont lump uninitialized declarations at the top of a routine. For example, do not do:
  ```c++
  int i,j,k; // BAD
  for (i=0; i<3; i++);
  ```
  Instead do
  ```c++
  for (int i=0; i<3; i++);
  ```
  The former misleadingly implies that the value of the loop index is used after the loop (this is only one of many reasons to prefer the latter style).
  
- If in doubt, refer to the official [C++ core guidelines](https://isocpp.github.io/CppCoreGuidelines/)
