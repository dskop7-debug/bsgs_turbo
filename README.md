# BSGS Turbo

BSGS Turbo is a high-performance C++ implementation of the Baby-step Giant-step (BSGS) algorithm for finding secp256k1 private keys from public keys.
It is heavily inspired by the BSGS mode of [keyhunt by AlbertoBSD](https://github.com/albertobsd/keyhunt), but completely rewritten from scratch to implement advanced mathematical and algorithmic optimizations.

## Features & Optimizations

- **Compact Hash Table:** Uses a memory-efficient Robin Hood hash table for the baby steps, completely eliminating Bloom filters and their false positives. Lookup time is O(1).
- **Batch Inversion (Montgomery's Trick):** Aggregates elliptic curve point conversions (Jacobian to Affine) using a single modular inversion, drastically speeding up baby-step precomputations.
- **Secp256k1 Endomorphism:** Implements the GLV endomorphism ($\lambda P = (\beta X, Y)$), effectively tripling (3x) the search speed without increasing memory usage.
- **Fast Point Addition:** Uses optimized Jacobian-Affine mixed point additions (`madd-2008-g`) during the Giant-step phase.
- **Multithreading:** Supports sequential and random search modes distributed across multiple CPU cores.
- **No external dependencies:** Pure C++17. Only uses standard libraries and custom 256-bit arithmetic tailored for the secp256k1 prime field.

## Compilation

The project uses a standard Makefile and requires a C++17 compatible compiler (like `g++` or `clang++`).

### Linux
```bash
make
```

### Windows (MSYS2 / MinGW-w64)
```bash
make windows
```

## Usage

The CLI arguments are designed to be familiar if you have used `keyhunt`.

```
Usage: ./bsgs_turbo [options]

Required:
  -f <file>      File with target public keys (one per line, hex format)
  -b <bits>      Bit range to search (e.g., 63 for puzzle #63)
  OR
  -r <from:to>   Explicit hex range (e.g., 4000000000000000:8000000000000000)

Optional:
  -t <threads>   Number of threads (default: CPU cores)
  -n <value>     Baby step count modifier (default: 0x100000 = 1,048,576)
  -k <factor>    K factor to multiply n (default: 1)
  -R             Random mode (threads pick random subranges instead of sequential)
  -e             Enable endomorphism (3x speed boost)
  -s <seconds>   Stats interval in seconds (default: 10)
  -q             Quiet mode (suppress per-thread target completion output)
  -o <file>      Output file for found keys (default: KEYFOUND.txt)
```

### Examples

Search for puzzle 63 sequentially using 8 threads and endomorphism:
```bash
./bsgs_turbo -f puzzle63.txt -b 63 -t 8 -e
```

Search a specific custom range in random mode (`-R`) with a larger baby-step table (to use more RAM but compute fewer giant steps):
```bash
./bsgs_turbo -f target.txt -r 2000000000000000:3FFFFFFFFFFFFFFF -t 4 -n 0x4000000 -k 4 -e -R
```

## Acknowledgements
- [AlbertoBSD's keyhunt](https://github.com/albertobsd/keyhunt) for the original idea and CLI parameter structures.
