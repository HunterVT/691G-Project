import sys
import argparse


def parse_fasta(filepath):
    """Parses a FASTA file and returns the sequence as a single string."""
    sequence = []
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence.append(line.strip().upper())
    return "".join(sequence)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("output", help="Output text file")
    args = parser.parse_args()

    sequence = parse_fasta(args.input)

    with open(args.output, 'w') as f:
        f.write(sequence)

    print(f"Sequence length: {len(sequence)}")


if __name__ == "__main__":
    main()