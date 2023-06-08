import dataclasses
from typing import Literal, Optional, Sequence

from Bio.Seq import Seq  # type: ignore

from gene_annotator.schema import SequenceSchema, AnnotationSchema


def compute_match(reference: Sequence[str], sequence: Sequence[str]):
    match = []
    for i in range(len(reference)):
        match.append(reference[i] == sequence[i])
    return match


@dataclasses.dataclass
class AminoAcidElement:
    text: str
    position: int  # corresponding to nucleiotide position
    number: int  # position in codon

    def __str__(self) -> str:
        return self.text

    def __repr__(self) -> str:
        return f"{self.text} {self.position} {self.number}"


@dataclasses.dataclass(frozen=True)
class Annotation:
    start: int = 0
    end: int = 0
    color: str = "blue"
    text: Optional[str] = None
    position: Literal["top", "bottom"] = "top"

    @staticmethod
    def from_schema(schema: AnnotationSchema) -> "Annotation":
        return Annotation(
            start=schema.start,
            end=schema.end,
            color=schema.color,
            text=schema.text,
            position=schema.position,
        )


class SequenceData:
    def __init__(
        self,
        id: str,
        nucleotide: str,
        annotations: Optional[list[Annotation]] = None,
        reference: bool = False,
        frame: Literal[0, 1, 2] = 0,
        protein: Optional[list[AminoAcidElement]] = None,
    ) -> None:
        self.id = id
        self.reference = reference
        self.nucleotide = nucleotide
        self.frame = frame
        self.annotations = annotations if annotations is not None else []

        self.protein = protein if protein is not None else []
        if len(self.protein) == 0 and len(self.nucleotide) > 0:
            self.translate()

    def translate(self):
        codon: list[tuple[str, int]] = []
        for i, nuc in enumerate(self.nucleotide[self.frame :]):
            if len(codon) == 3:
                amino_acid = str(Seq("".join([c[0] for c in codon])).translate())
                self.protein.extend(
                    [
                        AminoAcidElement(
                            amino_acid, position=codon[number][1], number=number
                        )
                        for number in range(3)
                    ]
                )
                codon = []

                if amino_acid == "*":
                    break

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

    def chunk(self, n: int) -> list["SequenceData"]:
        size = n
        chunks = []
        for i in range(0, len(self.nucleotide), size):
            chunk_protein = [
                AminoAcidElement(
                    text=elem.text, position=elem.position - i, number=elem.number
                )
                for elem in self.protein
                if i <= elem.position < i + size
            ]

            chunk_annotation = [a for a in self.annotations if i <= a.start < i + size]
            chunks.append(
                SequenceData(
                    id=self.id,
                    nucleotide=self.nucleotide[i : i + size],
                    reference=self.reference,
                    frame=self.frame,
                    protein=chunk_protein,
                    annotations=chunk_annotation,
                )
            )
        return chunks

    def __repr__(self) -> str:
        return f"nucleotide: {self.nucleotide}\nprotein: {self.protein}"

    @staticmethod
    def from_schema(sequence: SequenceSchema) -> "SequenceData":
        annotations = None
        if sequence.annotations is not None:
            annotations = [
                Annotation.from_schema(annotation)
                for annotation in sequence.annotations
            ]

        return SequenceData(
            id=sequence.id,
            nucleotide=sequence.sequence,
            reference=sequence.reference,
            annotations=annotations,
        )
