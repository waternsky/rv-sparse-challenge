# Problem Statement

RV-Sparse : Coding Challenge

Implement the `sparse_multiply` function:
* Scans a row-major matrix `A` and identifies its non-zero elements.
* Extracts them into Compressed Sparse Row (CSR) format using caller-provided buffers.
* Computes the matrix-vector product **y = A * x** using the extracted CSR data.
* Writes the result directly into a caller-provided output buffer.

Constraint: function must perform zero dynamic memory allocation. All memory is
pre-allocated by the caller
