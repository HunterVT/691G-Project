# ICS 691G Project: Benchmarking Underflow of Viterbi's algorithm for CpG island detection

This is an implementation of Viterbi's algorithm for CpG island detection and benchmarking. 

## Requirements:
- CMake 4.2.3 or newer
- C++17-compatible compiler (this project uses g++)
- GNU Time (for peak memory benchmarking, was done on Ubuntu)
- Python with pandas and Matplotlib

## Build from source

Building uses `cmake` and is as follows:

```bash
mkdir build
cd build
cmake ..
make
```
## Implementation and Usage

The executable is called `691G_Project`, and invoking it from the command-line
without any argument will display usage information. 
Use `./691G_Project <path_to_sequence> <data_type>`
The available `data_types` are `float`, `double`, and `long_double`.

Note that the `<path_to_sequence>` argument will only take a `.txt` file, and the input sequence
must be provided as a single line of uppercase nucleotides. If you have a FASTA file you would 
like to convert, use `preprocess_fasta.py` with the path to your FASTA file as an argument. 
