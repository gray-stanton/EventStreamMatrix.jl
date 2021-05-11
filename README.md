# EventStreamMatrix

This package implements non-materializing matrix operations for design matrices constructed via B-spline expansion of multi-source event-time data.
Additionally, benchmarking code against materializing matrix operartions is provided.

# Source Files
* IdentityEventStreamMatrix.jl - A discretized version of an event-time design matrix with no spline expansion
* BSplineEventStream.jl - Non-materializing multiplications an other operations for design matrices arising from B-spline basis expansion

# Script File
* core_op_benchmarks.jl - Benchmarking of relevant matrix operations.

# Example
Example call to `core_op_benchmarks.jl`:

```
core_op_benchmarks --mdensity 0.05 --n 5 --t 1e6 --fineness 0.5 --seed 5454 --seconds 10.0 --outfile ./out.csv --breakpoints 0.0 10.0 20.0 30.0 50.0 75.0 100.0
```