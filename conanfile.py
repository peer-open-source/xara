from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps

class OpenSeesConan(ConanFile):
    name = "opensees"
    version = "0.1.0"
    
    # Declare dependencies
    requires = [
        "tcl/8.6.11",
    ]
    
    generators = "CMakeDeps", "CMakeToolchain"
    default_options = {
        "tcl/*:shared": True, #
    }
    
#   def layout(self):
#       self.folders.source = "SRC"
#       self.folders.build = "build"
#       self.folders.generators = "build/generators"

    def generate(self):
        tc = CMakeToolchain(self)
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

