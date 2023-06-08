import xml.etree.ElementTree as ET
from typing import Literal, Optional, Union

from gene_annotator.schema import ConfigSchema
from gene_annotator.sequence_data import SequenceData, compute_match

number = Union[int, float]


class GeneAnnotatorSVG:
    def __init__(
        self,
        sequence_data_list: list[SequenceData],
        config: ConfigSchema,
        viewBox: str = "0 0 1500 750",
    ) -> None:
        self.block_width = config.block_width
        self.blocks: list[list[SequenceData]] = list(
            map(
                list,
                zip(*[seq.chunk(config.block_width) for seq in sequence_data_list]),
            )
        )

        max_id_length = max(
            len(sequence_data.id) for sequence_data in sequence_data_list
        )

        # Font settings
        self.fontsize = config.fontsize
        self.fontfamily = config.fontfamily

        ## Calculate the position of the first text
        self.current_x = (max_id_length + 1) * self.fontsize
        self.current_y = self.fontsize * 1.5
        
        ## ticks settings
        self.offset = config.offset
        self.start = config.start
        self.tick_index = 0

        ## translation settings
        self.with_translation = config.with_translation
        
        self.root = ET.Element(
            "svg", xmins="http://www.w3.org/2000/svg", viewBox=viewBox
        )

    def write_text(
        self,
        text: Optional[str],
        x: number,
        y: number,
        fill: str = "black",
        text_acnhor: Literal["start", "middle", "end"] = "middle",
        dominant_baseline: Literal[
            "central",
            "hanging",
            "middle",
            "text-before-edge",
            "text-after-edge",
            "inherit",
        ] = "central",
        **kwargs,
    ):
        text_attribute: dict[str, str] = {
            "text-anchor": text_acnhor,
            "dominant-baseline": dominant_baseline,
            "font-family": self.fontfamily,
        }
        
        text_attribute.update(kwargs)
        text_element = ET.SubElement(
            self.root,
            "text",
            x=str(x),
            y=str(y),
            fill=fill,
            style=f"font-size: {self.fontsize}px",
            **text_attribute,
        )

        text_element.text = text

    def write_ticks(self, start: int):
        for i in range(0, self.block_width):
            self.tick_index += 1
            if self.tick_index % self.offset == 0:
                x = self.current_x + i * self.fontsize
                self.write_text(str(start + i + 1), x=x, y=self.current_y - 20)
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
        self,
        sequence_data: SequenceData,
        reference: SequenceData,
    ):
        top_annotations = [
            annotation
            for annotation in sequence_data.annotations
            if annotation.position == "top"
        ]
        # write top annotation
        if len(top_annotations) > 0:
            for annotation in top_annotations:
                self.write_text(
                    annotation.text,
                    x=self.current_x + self.fontsize * annotation.start * 1.5,
                    y=self.current_y,
                    fill="black",
                    text_acnhor="middle",
                )
                ET.SubElement(
                    self.root,
                    "line",
                    x1=str(self.current_x + self.fontsize * (annotation.start - 0.5)),
                    x2=str(self.current_x + self.fontsize * (annotation.end + 0.5)),
                    y1=str(self.current_y + self.fontsize / 2),
                    y2=str(self.current_y + self.fontsize / 2),
                    stroke=annotation.color,
                    **{"stroke-width": "5"},
                )

            self.current_y += self.fontsize * 1.3

        # write IDs
        self.write_text(
            sequence_data.id,
            x=self.current_x - self.fontsize * 1.5,
            y=self.current_y,
            text_acnhor="end",
        )

        match_nucleotide = compute_match(sequence_data.nucleotide, reference.nucleotide)

        # write nucleotide
        for i, nuc in enumerate(sequence_data.nucleotide):
            self.write_text(
                nuc,
                x=self.current_x + self.fontsize * i,
                y=self.current_y,
                fill="black" if match_nucleotide[i] else "red",
            )

        # write bottom annotation
        bottom_annotations = [
            annotation
            for annotation in sequence_data.annotations
            if annotation.position == "bottom"
        ]
        if len(bottom_annotations) > 0:
            self.current_y += self.fontsize * 0.3
            for annotation in bottom_annotations:
                self.write_text(
                    annotation.text,
                    x=self.current_x + self.fontsize * annotation.start * 1.5,
                    y=self.current_y + self.fontsize,
                    fill="black",
                    text_acnhor="middle",
                )
                ET.SubElement(
                    self.root,
                    "line",
                    x1=str(self.current_x + self.fontsize * (annotation.start - 0.5)),
                    x2=str(self.current_x + self.fontsize * (annotation.end + 0.5)),
                    y1=str(self.current_y + self.fontsize / 2),
                    y2=str(self.current_y + self.fontsize / 2),
                    stroke=annotation.color,
                    **{"stroke-width": "5"},
                )

            self.current_y += self.fontsize

        # write protein
        if self.with_translation:
            self.current_y += self.fontsize * 1.3
            match_protein = compute_match(
                [p.text for p in sequence_data.protein],
                [p.text for p in reference.protein],
            )

            for i, p in enumerate(sequence_data.protein):
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
                            fill="black" if match_protein[i] else "red",
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

    def write_block(self, i: int, block: list[SequenceData]):
        self.write_ticks(start=self.start + i * self.block_width)
        reference = next(
            sequence_data for sequence_data in block if sequence_data.reference
        )

        for sequence_data in block:
            self.write_sequence_data(
                sequence_data=sequence_data, reference=reference
            )
            self.current_y += self.fontsize * 1.5

    def write_blocks(self):
        for i, block in enumerate(self.blocks):
            self.write_block(i, block)
            self.current_y += self.fontsize * 3

    def write_svg(self, output_file: str):
        self.write_blocks()
        tree = ET.ElementTree(self.root)
        tree.write(output_file, encoding="utf-8", xml_declaration=True)
