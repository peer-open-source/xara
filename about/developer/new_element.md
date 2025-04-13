A new element is implemented in two steps:

1. The element class must be implemented. This is a class that inherits from the `Element` class and implements the basic state determination methods. These methods include:
   - `getResistingForce()`
   - `getTangentStiff()`
2. A parser must be implemented. This is a function that reads an input file and constructs an instance of the element. Once the parsing function is implemented, it  can be added to the element library in [`SRC/runtime/commands/modeling/element.hpp`](SRC/runtime/commands/modeling/element.hpp)


Once the `opensees` package is installed, you can have it use your locally compiled library by setting the environment variable `OPENSEESRT_LIB` to point to your local version of `libOpenSeesRT.so`