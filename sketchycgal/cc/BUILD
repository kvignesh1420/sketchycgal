cc_library(
    name = "tracer",
    hdrs = ["tracer.h"]
)

cc_library(
    name = "reader",
    srcs = ["reader.cc"],
    hdrs = ["reader.h"],
    visibility = ["//visibility:public"],
    deps = [
        "//third_party/eigen3"
    ]
)

cc_library(
    name = "nystrom",
    srcs = ["nystrom.cc"],
    hdrs = ["nystrom.h"],
    visibility = ["//visibility:public"],
    deps = [
        "//third_party/eigen3"
    ]
)

cc_library(
    name = "cgal",
    srcs = ["cgal.cc"],
    hdrs = ["cgal.h"],
    visibility = ["//visibility:public"],
    deps = [
        ":reader",
        ":tracer",
        ":nystrom",
        "//third_party/eigen3"
    ]
)