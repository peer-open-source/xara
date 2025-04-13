
3. Install the package:
   ``` shell
   python -m pip install -e .
   ```
   This will (1) carry out a full build, (2) copy the build artifacts into
   the project source directory, (3) install the `opensees` package
   into the current Python environment, and (4) remove the build directory. 

   If you do not expect to actively develop in C/C++, nothing else needs to
   be done. However, for users planning to actively develop in C/C++, it is 
   convenient to set up a *persistent* build directory that does not
   get cleaned out by `pip`. This is explained in the optional step 5.
