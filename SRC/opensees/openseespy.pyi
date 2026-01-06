from typing import overload, Literal


class _Materials:
    def material(self, type: str, tag: int, *args, **kwargs):
        """
        Create a material of the specified type.
        """

    def section(self, type: str, tag: int, *args, **kwargs):
        """
        Create a section of the specified type.
        """

class _Elements:

    @overload
    def element(self, type: str, tag: int, nodes: tuple, **kwargs) -> int:
        """
        Create an element of the specified type.
        """
    
    @overload
    def element(self,
                __elmType: Literal["Truss"],
                tag: int,
                nodes: tuple,
                section: int,
                *args
            ) -> int:
        """
        Create a truss element.
        
        :param tag: unique :ref:`element` tag
        :type tag: |integer|
        :param nodes: tuple of *two* integer :ref:`node` tags
        :type nodes: tuple
        :param section: tag of a previously-defined :ref:`section`
        :type section: |integer|
        """

    @overload
    def element(self,
                __elmType: Literal["PrismFrame"],
                tag: int,
                nodes: tuple,
                section: int,
                transform: int,
                *args
            ) -> int:
        """
        Create a linear elastic prismatic frame element.        
        
        :param tag: unique :ref:`element` tag
        :type tag: |integer|
        :param nodes: tuple of *two* integer :ref:`node` tags
        :type nodes: tuple
        :param section: tag of a previously-defined :ref:`section`
        :type section: |integer|
        :param transform: identifier for previously-defined coordinate-transformation
        :type transform: |integer|
        """

    @overload
    def element(self,
                __elmType: Literal["ForceFrame"],
                tag: int,
                nodes: tuple,
                nip: int,
                transform: int,
                *args
               ) -> int: ...

class _Algorithm:
    def numberer(self, type: str, *args, **kwargs):
        """
        Create a numberer of the specified type.
        """

    def constraints(self, type: str, *args, **kwargs):
        """
        Create a constraint handler of the specified type.
        """

    def system(self, type: str, *args, **kwargs):
        """
        Create a system of the specified type.
        """

    def integrator(self, type: str, *args, **kwargs):
        """
        Create an integrator of the specified type.
        """
    
    def algorithm(self, type: str, *args, **kwargs) -> int:
        """
        Create an algorithm of the specified type.
        """

    def algorithm(self, 
                  type: Literal["NewtonLineSearch"],
                  tol=0.8) -> int:
        """
        Create a specific algorithm for the analysis.
        """

    def analysis(self, type: str):
        """
        Create an analysis of the specified type.
        """

    def analyze(self, nsteps: int = 1, *args) -> int:
        """
        Perform the analysis for a specified number of steps.
        
        :param nsteps: number of analysis steps to perform
        :type nsteps: |integer|
        """

    def test(self, type: str, *args, **kwargs):
        """
        Create a convergence test of the specified type.
        """

    def numIter(self) -> int:
        """
        Get the number of iterations for the last analysis step.
        """

    def wipe(self)->None:
        """
        Clear the current model.
        """


class _Output:
    def nodeDisp(self, node: int, dof: int) -> float:
        """
        Get the displacement of a node in a specified degree of freedom.
        """

    def getNodeTags(self) -> list[int]:
        """
        Get a list of all node tags in the model.
        """

    def nodeCoord(self, node: int, index: int=None) -> list[float] | float:
        """
        Get the coordinates of a node.
        """

    @overload
    def nodeDisp(self, node: int) -> list[float]:
        """
        Return the displacements of a node.
        """

    @overload
    def nodeDisp(self, node: int, dof: int) -> float:
        """
        Return the displacement of a node.
        """

    @overload
    def nodeReaction(self, node: int) -> list[float]:
        """
        Get the reaction forces at a node.
        """
    @overload
    def nodeReaction(self, node: int, dof: int) -> float:
        """
        Get the reaction force(s) at a node.
        """

    def getEleTags(self) -> list[int]:
        """
        Get a list of all element tags in the model.
        """

    def eleResponse(self, ele: int, responseType: str) -> list[float]:
        """
        Get the response of an element.
        """

class Model(_Materials,_Elements,_Algorithm, _Output):
    def node(self, tag: int, coords: tuple, *args, **kwargs) -> int:
        """
        Create a node with the specified tag and coordinates.
        """

    def fix(self, node: int, *args) -> int:
        """
        Apply boundary conditions to a node.
        """

    def mass(self, node: int, *args):
        """
        Assign mass to a node.
        """

    def load(self, node: int, *args):
        """
        Apply loads to a node.
        """
