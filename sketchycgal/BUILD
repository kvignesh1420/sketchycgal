load("@pybind11_bazel//:build_defs.bzl", "pybind_extension")

pybind_extension(
    name = "_libsketchycgal",
    srcs = ["py_api_binder.cc"],
    deps = [
        "//sketchycgal/cc:cgal",
        "//third_party/eigen3"
    ]
)