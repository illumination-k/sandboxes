import json
from typing import Union, Optional, Literal
import xml.etree.ElementTree as ET
from Bio.Seq import Seq
import dataclasses


def compute_match(reference, sequences):
    match = []
    for i in range(len(reference)):
        if all(sequence[i] == reference[i] for sequence in sequences):
            match.append(True)
        else:
            match.append(False)
    return match


def translate_sequence(dna_sequence):
    seq = Seq(dna_sequence)
    protein = seq.translate()
    return str(protein)


@dataclasses.dataclass
class AminoAcidElement:
    text: str
    position: int  # corresponding to nucleiotide position
    number: int  # position in codon


@dataclasses.dataclass
class SequenceData:
    id: str
    nucleotide: str
    reference: bool = False
    frame: Literal[0, 1, 2] = 0

    def __post_init__(self):
        self.protein: list[AminoAcidElement] = []
        codon: list[tuple[str, i]] = []
        for i, nuc in enumerate(self.nucleotide[self.frame :]):
            if len(codon) == 3:
                print(codon)
                amino_acid = str(Seq("".join([c[0] for c in codon])).translate())
                self.protein.extend(
                    [
                        AminoAcidElement(amino_acid, position=codon[i][1], number=i)
                        for i in range(3)
                    ]
                )

                codon = []
            if nuc != "-":
                codon.append((nuc, i + self.frame))

            if len(codon) == 3:
                amino_acid = str(Seq("".join([c[0] for c in codon])).translate())
                self.protein.extend(
                    [
                        AminoAcidElement(amino_acid, position=codon[i][1], number=i)
                        for i in range(3)
                    ]
                )

    def get_amino_acid(self, i) -> str:
        return self.protein[i * 3].text


number = Union[int, float]


class GeneAnnotatorSVG:
    def __init__(
        self, sequence_data_list: list[SequenceData], viewBox: str = "0 0 1000 500"
    ) -> None:
        self.sequence_data_list = sequence_data_list
        # self.reference = next(
        #     sequence_data for sequence_data in sequence_data_list if sequence_data.reference
        # )

        self.current_x = 10
        self.current_y = 30
        self.fontsize = 20
        self.root = ET.Element(
            "svg", xmins="http://www.w3.org/2000/svg", viewBox=viewBox
        )

    def write_text(
        self, text: str, x: number, y: number, fill: str = "black", **kwargs
    ):
        text_element = ET.SubElement(
            self.root,
            "text",
            x=str(x),
            y=str(y),
            fill=fill,
            style=f"font-size: {self.fontsize}px",
            **{"text-anchor": "middle", "dominant-baseline": "central"},
            **kwargs,
        )
        text_element.text = text

    def write_ticks(self, start: int = 0, offset=10):
        length = len(self.sequence_data_list[0].nucleotide)

        for i in range(0, length, offset):
            x = self.current_x + i * self.fontsize
            self.write_text(str(start + i), x=x, y=self.current_y - 20)
            ET.SubElement(
                self.root,
                "line",
                x1=str(x),
                x2=str(x),
                y1=str(self.current_y - 5),
                y2=str(self.current_y),
                stroke="black",
            )

        self.current_y += self.fontsize

    def write_sequence_data(
        self, sequence_data: SequenceData, with_translation: bool = False
    ):
        for i, nuc in enumerate(sequence_data.nucleotide):
            # TODO! referenceとの差
            self.write_text(nuc, x=self.current_x + self.fontsize * i, y=self.current_y)

        if with_translation:
            self.current_y += self.fontsize
            for p in sequence_data.protein:
                top = self.current_y - self.fontsize / 2
                bottom = self.current_y + self.fontsize / 2
                if p.number == 0:
                    ET.SubElement(
                        self.root,
                        "polyline",
                        fill="none",
                        stroke="black",
                        points=" ".join(
                            [
                                f"{self.current_x + self.fontsize * p.position + self.fontsize / 2},{top}",
                                f"{self.current_x + self.fontsize * p.position - self.fontsize / 2},{top}",
                                f"{self.current_x + self.fontsize * p.position - self.fontsize / 2},{bottom}",
                                f"{self.current_x + self.fontsize * p.position + self.fontsize / 2},{bottom}",
                            ]
                        ),
                    )
                elif p.number == 1:
                    line_start = (
                        self.current_x + self.fontsize * p.position - self.fontsize / 2
                    )
                    line_end = (
                        self.current_x + self.fontsize * p.position + self.fontsize / 2
                    )
                    ET.SubElement(
                        self.root,
                        "line",
                        stroke="black",
                        x1=str(line_start),
                        x2=str(line_end),
                        y1=str(top),
                        y2=str(top),
                    )

                    ET.SubElement(
                        self.root,
                        "line",
                        stroke="black",
                        x1=str(line_start),
                        x2=str(line_end),
                        y1=str(bottom),
                        y2=str(bottom),
                    )

                    if p.text == "*":
                        self.write_text(
                            p.text,
                            x=self.current_x + self.fontsize * p.position,
                            y=self.current_y,
                            fill="purple",
                        )
                    else:
                        self.write_text(
                            p.text,
                            x=self.current_x + self.fontsize * p.position,
                            y=self.current_y,
                        )
                else:
                    ET.SubElement(
                        self.root,
                        "polyline",
                        fill="none",
                        stroke="black",
                        points=" ".join(
                            [
                                f"{self.current_x + self.fontsize * p.position - self.fontsize / 2},{top}",
                                f"{self.current_x + self.fontsize * p.position + self.fontsize / 2},{top}",
                                f"{self.current_x + self.fontsize * p.position + self.fontsize / 2},{bottom}",
                                f"{self.current_x + self.fontsize * p.position - self.fontsize / 2},{bottom}",
                            ]
                        ),
                    )

                    if p.text == "*":
                        break

    def write_svg(self, output_file: str):
        tree = ET.ElementTree(self.root)
        tree.write(output_file, encoding="utf-8", xml_declaration=True)


def create_svg(input_file, output_file):
    with open(input_file, "r") as f:
        data = json.load(f)
    print(data)

    sequence_data_list = [
        SequenceData(nucleotide=d["sequence"], id=d["id"], frame=2)
        for d in data["sequences"]
    ]
    svg = GeneAnnotatorSVG(sequence_data_list)
    svg.write_ticks()
    for sequence_data in svg.sequence_data_list:
        svg.write_sequence_data(sequence_data, with_translation=True)
        svg.current_y += svg.fontsize
        svg.current_y += svg.fontsize / 4

    svg.write_svg(output_file)


create_svg("input.json", "output.svg")
