from typing import overload, Literal

class _Elements:
    @overload
    def element(self,
                __elmType: Literal["PrismFrame"],
                tag: int,
                nodes: tuple,
                secTag: int,
                *args
               ) -> int: ...
    @overload
    def element(self,
                __elmType: Literal["forceBeamColumn"],
                tag: int,
                nodes: tuple,
                transform: int,
                nInts: int,
                pts: list,
                *args
               ) -> int: ...

    def element(self, *args, **kwargs) -> int: ...

class Model(_Elements):
    pass
