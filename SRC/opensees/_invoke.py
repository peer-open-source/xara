from functools import partial

class _Handle:
    def __init__(self, interpreter, type, tag, *args):
        self._interpreter = interpreter
        self._type = type 
        self._tag = tag

    # def __getattribute__(self, name):
    #     if name in ["getStress", "get_stress", "get_tangent"]:
    #         return partial(self._interpreter._invoke_proc, 
    #                    "invoke", self._type, self._tag, name)

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def set_strain(self, *args):
        pass

    def getStress(self, e, commit=False):
        self._interpreter._invoke_proc("invoke", self._type, self._tag, [f"strain {e};"])
        if commit:
            self._interpreter._invoke_proc("invoke", self._type, self._tag, "commit")
        return self._interpreter._invoke_proc("invoke", self._type, self._tag, "stress")

    def get_tangent(self, *args):
        pass
    
