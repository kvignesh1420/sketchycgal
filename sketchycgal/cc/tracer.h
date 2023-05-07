#ifndef SKETCHYCGAL_CC_TRACER_H
#define SKETCHYCGAL_CC_TRACER_H

// TRACER macro for verbose tracing of intermediate actions
#define TRACER(fmt...) do { \
    printf(fmt); fflush(stdout);  \
} while(0)

#endif
