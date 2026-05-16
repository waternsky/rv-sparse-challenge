# Problem Statement

RV-Sparse : Coding Challenge

Implement the `sparse_multiply` function:

- Scans a row-major matrix `A` and identifies its non-zero elements.
- Extracts them into Compressed Sparse Row (CSR) format using caller-provided buffers.
- Computes the matrix-vector product **y = A \* x** using the extracted CSR data.
- Writes the result directly into a caller-provided output buffer.

Constraint: function must perform zero dynamic memory allocation. All memory is
pre-allocated by the caller

# Implementation

- to check for zero in case of double precision numbers we need tolerance to
  check whether it is zero or not, in our case we set `tol = 1e-7`
- row_ptrs always start with 0, at index `i` (i > 0) gives number of non-zero elements
  till row (i - 1) (row 0 being first row)
- `row_ptrs[i + 1] - row_ptr[i]` gives number of element in ith row

To test the implementation run

```sh
gcc -lm -o run challenge.c
./run
```
