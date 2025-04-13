from . import tcl

class Analysis:
    def __init__(self, model, strategy, patterns, recorders=None):
        #
        # Runtime
        #
        if not hasattr(model, "_rt") or model._rt is None:
            from . import tcl
            self.rt = tcl.TclRuntime()
            self.rt.send(model)
        else:
            self.rt = model

        #
        # Strategy
        #
        # Currently, all analysis constructors take the `strategy`
        # argument as a C++ std::vector<std::string>, so
        # parameters must be cast to strings.
        self._strategy = {
            k: [str(i) for i in v]
            for k,v in strategy.items()
        }

        #
        # Patterns
        #
        self.patterns = patterns
        if patterns is not None:
            for pattern in patterns.values():
                self.rt.eval(tcl.dumps(pattern))

        #
        # Recorders
        #
        if recorders is not None:
            self.recorders = {
                    self.rt.eval(tcl.dumps(recorder)): recorder
                    # recorder.link(runtime=self.rt): recorder
                for recorder in recorders
            }

    def system(self, *args): pass

    def integrator(self, *args): pass

    def handler(self, *args): pass

    def analyze(self, *args, **kwds):
        self._analysis.analyze(*args)



class StaticAnalysis(Analysis):
    def __init__(self, model, strategy=None, patterns=None):
        if strategy is None:
            strategy = {
                "algorithm": ["Newton"],
                "integrator": ["LoadControl"]
            }
        super().__init__(model, strategy=strategy, patterns=patterns)

        # Import C++ bindings and create an instance of the analysis
        from . import OpenSeesPyRT as libOpenSeesRT
        self._analysis = libOpenSeesRT._StaticAnalysis(self.rt._rt, self._strategy)





class DirectIntegrationAnalysis(Analysis):
    """
    DirectIntegrationAnalysis(patterns=[
        pattern.MultipleSupport(
            components = [
            #   node, dof,  history
                ( 1,   1,   pattern.ResponseComponent(accel, displ, veloc)),
                ( 2,   1,   pattern.ResponseComponent(accel, displ, veloc))
            ]
        ),

        pattern.UniformAcceleration([
        #    dof   history
            ( 1,   pattern.ResponseComponent(accel, displ, veloc))
          ]
        )

    ])
    """

    def __init__(self,
        model,
        patterns:  dict = None,
        strategy:  dict = None,
        recorders: list = None,
        inherit:   str  = None,
        gravity:   dict = None,
        time_step: float= None,
        steps:     int  = 1
    ):
        if strategy is None:
            strategy = {
                # "algorithm": ["Newton"],
                # "integrator": ["Newmark"]
            }
        super().__init__(model, strategy=strategy, patterns=patterns, recorders=recorders)
        # Import C++ bindings and create an instance of the analysis
        from . import OpenSeesPyRT as libOpenSeesRT
        self._analysis = libOpenSeesRT._DirectIntegrationAnalysis(self.rt._rt, self._strategy)

        # self.steps, self.time_step = min((
        #     (len(pat.series.values), pat.series.time_step)
        #         for pat in self.patterns.values() 
        #             if pat.series.time_step is not None
        #     ),
        #     key = lambda i: i[1]
        # ) if time_step is None else (steps, time_step)


    def analyze(self, steps=None, time_step=None, **kwds):
        time_step = time_step or kwds.get("dt", None) or self.time_step
        if time_step is None:
            raise ValueError("Unable to deduce time step size")

        steps = steps or self.steps
        if steps is None:
            raise ValueError("Unable to deduce time step size")
            while True:
                self._analysis.analyze(1, time_step)
        else:
            return self._analysis.analyze(steps, time_step)


