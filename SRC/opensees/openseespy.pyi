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
                __elmType: Literal["forceBeamColumn"],
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

class Model(_Elements,_Algorithm):
    pass
