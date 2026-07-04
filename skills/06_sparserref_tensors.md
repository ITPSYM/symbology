# Skill 06 — SparseRREF tensor library essentials

This is a reference for the SparseRREF library types and functions used
throughout `symbology`. SparseRREF is not vendored — clone it separately
into `SparseRREF/` at the repository root.

## Core types

```cpp
using field_t = ...;  // field descriptor (FIELD_QQ for rationals)
using scalar_t = rat_t;  // rational number type from SparseRREF
using index_t = int32_t;

template <typename T, typename index_t, SPARSE_TYPE S>
struct sparse_tensor;  // S = SPARSE_CSR or SPARSE_COO

template <typename T, typename index_t>
struct sparse_vec;  // 1-d sparse vector

template <typename T, typename index_t>
struct sparse_mat;  // 2-d sparse matrix (row-major)
```

## CSR vs COO conventions

- **`bootstrap.cpp` and `bootstrap.hpp` use CSR** as the canonical format.
- **All tensor functions return CSR** (this is a hard project constraint).
- **COO is used only internally** for:
  - Reshaping (`reshape` is a COO-only method — CSR has no `reshape`).
  - The shuffle product (`tensor_shuffle_product_parallel` accumulates
    into an `unordered_map`, then builds a COO and converts to CSR).
- **WXF format is CSR-native** — write CSR tensors directly with
  `sparse_tensor_write_wxf`.

To reshape a CSR tensor:
```cpp
sparse_tensor<T, index_t, SPARSE_COO> coo(std::move(csr_tensor));
coo.reshape(new_dims);
auto back_to_csr = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(coo), pool);
```

## Key functions

### Tensor contraction

```cpp
auto C = tensor_contract(A_coo, B_coo, axis_a, axis_b, F, pool);
```
- Contracts `axis_a` of `A` with `axis_b` of `B`.
- The contracted axis of `B` moves to the **end** of the result.
- Inputs are COO; output is COO (convert to CSR if needed).

### Shuffle product

```cpp
auto C = tensor_shuffle_product_parallel<index_t, T>(A_coo, B_coo, F, pool);
```
- Generates all interleavings of positions from `A` and `B`.
- Result rank = `A.rank() + B.rank()`, all dims = `max(dim across A and B)`.
- **Pass `pool = nullptr` for the sequential variant** — the parallel
  variant produces incorrect results for boundary computation
  (`E1^2/2`).
- Implementation note: `B.dims()` returns by value (a temporary
  `std::vector<size_t>`); capture into a local before calling
  `begin()`/`end()` to avoid dangling iterators.

### Sparse RREF

```cpp
auto [rref, pivots, rank] = sparse_mat_rref_reconstruct(A, opt);
auto kernel = sparse_mat_rref_kernel(A, opt);
```
- `opt` is a `rref_option_t` (a 1-element array; pass as `opt`).
- `opt->method = 0` (right + left search), `opt->verbose = true`.
- `sparse_mat_rref_reconstruct` works over `Z / 2^61` and reconstructs
  to `rat_t`.

### Sparse matrix utilities

```cpp
auto joined = sparse_mat_join(A, B, pool);  // horizontal join
auto At = A.transpose();                    // transpose
auto rref = sparse_mat_rref_reconstruct(A, opt);
```

### WXF I/O

```cpp
auto u8arr = sparse_tensor_write_wxf(tensor_csr);  // CSR -> bytes
ofs.write(reinterpret_cast<const char*>(u8arr.data()), u8arr.size());

auto tensor_csr = sparse_tensor_read_wxf<T, index_t>(bytes, F, pool);
```

The project's `bootstrap.hpp` provides wrappers:
- `read_tensor(path, F, pool)` → CSR tensor (also prints CRC32).
- `write_tensor_file(path, std::move(tensor_csr))` → writes WXF file.
- `projection_read_tensor` / `projection_write_tensor` (in
  `projection.hpp`) — same but also `create_directories` on the parent
  path before writing.

## Tensor introspection

```cpp
tensor.rank();        // number of axes
tensor.dim(i);        // size of axis i
tensor.dims();        // std::vector<size_t> (by value — capture locally!)
tensor.nnz();         // non-zero count
tensor.alloc();       // allocated capacity
tensor.gen_perm();    // returns a permutation iterator for COO traversal
tensor.index_vector(i);  // the multi-index of entry i (std::vector<index_t>)
tensor.val(i);        // the value of entry i
```

## Push_back and accumulation

For COO tensors, the safe accumulation pattern (avoids the
`insert_add` corruption for large results):

```cpp
std::unordered_map<std::vector<index_t>, T, VectorHash> accumulated;
// ... fill accumulated ...
C.reserve(accumulated.size() + 1);
for (const auto& [idx, val] : accumulated) {
    if (val != T(0)) C.push_back(idx, val);
}
C.canonicalize();  // remove zeros from cancellation
C.sort_indices();  // lexicographic sort (CSR conversion requires this)
```

`push_back(idx, val)` on a COO tensor copies the index and appends.
`canonicalize()` removes zero entries created by cancellation.
`sort_indices()` sorts entries lexicographically (required before
converting to CSR or comparing tensors).

## Thread pool

```cpp
thread_pool pool(n_threads);  // BS::thread_pool
// or via rref_option:
rref_option_t opt;
opt->pool = thread_pool(n_threads);
thread_pool* pool = &(opt->pool);
```

- `pool->detach_blocks(start, end, lambda, n_blocks)` — parallel for.
- `pool->wait()` — block until all detached work finishes.
- Pass `nullptr` to functions that accept `thread_pool*` to force the
  sequential path.

## Common pitfalls

1. **`B.dims()` returns by value** — do not chain `begin()`/`end()` on
   two separate `B.dims()` calls. Capture locally.
2. **CSR has no `reshape`** — convert to COO first.
3. **`tensor_shuffle_product_parallel` parallel variant is incorrect**
   for boundary computation — pass `pool = nullptr`.
4. **`std::map<key, T>` overwrites duplicate keys** — when iterating a
   tensor with multiple entries per key, do not collapse into a map.
5. **WXF cannot serialize 0-nnz tensors** — skip writing empty
   projections.
6. **`sparse_tensor` move constructor** transfers ownership of the
   underlying buffers; the moved-from tensor is left empty (`_alloc = 0`).

## Reference

- Library: <https://github.com/munuxi/SparseRREF>
- Headers used: `SparseRREF/sparse_rref.h`, `SparseRREF/sparse_type.h`,
  `SparseRREF/sparse_mat.h`, `SparseRREF/scalar.h`,
  `SparseRREF/thread_pool.hpp`.
