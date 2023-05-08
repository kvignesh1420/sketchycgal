"""Run max-cut experiments through python bindings"""

import os
import time
from pprint import pprint
import pandas as pd
from sketchycgal import SketchyCGAL



class Runner:
    def __init__(self, datapath, max_iters, sketch_rank, tolerance):
        self.datapath = datapath
        self.max_iters = max_iters
        self.sketch_rank = sketch_rank
        self.tolerance = tolerance
        self.stats = {}
        self.initialize_stats()
        self.run()

    def initialize_stats(self):

        for dir_name in os.listdir(DATAPATH):
            if os.path.isdir(os.path.join(DATAPATH, dir_name)):
                filepath = os.path.join(DATAPATH, dir_name, dir_name, "{}.mtx".format(dir_name))
                self.stats[dir_name] = {"dataset": dir_name}
                self.stats[dir_name]["filepath"] = filepath
                self.stats[dir_name]["max_iters"] = self.max_iters
                self.stats[dir_name]["sketch_rank"] = self.sketch_rank
                self.stats[dir_name]["tolerance"] = self.tolerance

    def _helper(self, filepath):
        _runner = SketchyCGAL()
        _runner.setup(filepath, self.max_iters, self.sketch_rank, self.tolerance)
        start = time.time()
        _runner.run()
        end = time.time()
        return end - start

    def prepare_df(self):
        records = []
        for k, v in self.stats.items():
            records.append(v)
        df = pd.DataFrame(records).sort_values("dataset")
        print(df)

    def run(self):
        for dname in list(self.stats.keys()):
            print("Running SketchyCGAL on: {}".format(dname))
            filepath = self.stats[dname]["filepath"]
            duration = self._helper(filepath=filepath)
            self.stats[dname]["duration"] = duration
            self.stats[dname]["throughput"] = self.stats[dname]["max_iters"]/duration
        self.prepare_df()


if __name__ == "__main__":
    DATAPATH = "data/G"
    MAX_ITERS = 200
    TOLERANCE = 0.1
    SKETCH_RANK = 10
    r = Runner(datapath=DATAPATH, max_iters=MAX_ITERS, sketch_rank=SKETCH_RANK, tolerance=TOLERANCE)
