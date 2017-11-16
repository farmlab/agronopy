import logging
import sys
import pandas as pd
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)


class SimulatorMixin:
    default_params = {}
    default_variables = {}
    res = {}

    def __init__(self, inputs, variables={}, params={}):
        # dynamic variables
        self.inputs = inputs
        # Static variables
        self.variables = self.default_variables.copy()
        self.variables.update(variables)
        # parameters
        self.params = self.default_params.copy()
        self.params.update(params)

    def run(self):
        raise NotImplementedError("Subclasses should implement this!")

    def to_date(self, ts):
        return [pd.to_datetime(d + pd.datetime(1970, 1, 1)) for d in ts]
