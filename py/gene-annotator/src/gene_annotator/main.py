import argparse
from gene_annotator.annotator import create_svg


def main():
    parser = argparse.ArgumentParser(description="Gene Annotator")
    parser.add_argument("--input", "-i", type=str, required=True, help="Input file")
    parser.add_argument("--output", "-o", type=str, required=True, help="Output file")
    args = parser.parse_args()

    create_svg(args.input, args.output)


if __name__ == "__main__":
    main()
