import argparse

from gene_annotator.annotator import GeneAnnotatorSVG
from gene_annotator.schema import InputSchema
from gene_annotator.sequence_data import SequenceData


def main():
    parser = argparse.ArgumentParser(description="Gene Annotator")
    parser.add_argument("--input", "-i", type=str, required=True, help="Input file")
    parser.add_argument("--output", "-o", type=str, required=True, help="Output file")
    args = parser.parse_args()

    input = InputSchema.parse_file(args.input)
    print(input.config)
    sequence_data_list = [SequenceData.from_schema(s) for s in input.sequences]

    annotator = GeneAnnotatorSVG(
        config=input.config, sequence_data_list=sequence_data_list
    )
    annotator.write_svg(args.output)


if __name__ == "__main__":
    main()
