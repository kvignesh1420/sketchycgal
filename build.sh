
#!/bin/bash

bazel build //sketchycgal:all --config=optimization
bazel build //solver:all --config=optimization