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
More specifically: `python3 preprocess_fasta.py <path_to_input_fasta> <path_to_output_txt>`

## Dataset & Preliminary Results
The dataset used for this project is courtesy of UCSC's genome library. As of right now, we are using a subsequence of nucleotides from
chromosome 12 of the human reference genome (hg38) for benchmarking. We are using nucleotides 6500001-6600000, and our sequence 
length is 100000 base pairs. This data can be found in the `data/` directory. 

Preliminary results can be found in `data_processing.ipynb`.

## Benchmarking
Execution time is measured internally using the C++ `chrono` library with `steady_clock` and is reported in milliseconds.
Peak memory usage is measured using GNU time. To reproduce benchmarks, use the following: 
```
/usr/bin/time -v ./691G_Project ../data/hg38_chr12.txt float
/usr/bin/time -v ./691G_Project ../data/hg38_chr12.txt double
/usr/bin/time -v ./691G_Project ../data/hg38_chr12.txt long_double
```
The peak memory is reported as the `Maximum Resident Set Set Size (kbytes)` in the terminal output. 
