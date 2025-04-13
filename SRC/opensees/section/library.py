class FrameSection: pass

class FiberSection(FrameSection):
    def to_shear(self, material):
        pass

class ShearSection(FrameSection):
    pass


class PlateSection(ShellSection):
    pass


# Homogeneous constructor
<Type>Section(shape, material)

# Composite constructor
<Type>Section(shape)
