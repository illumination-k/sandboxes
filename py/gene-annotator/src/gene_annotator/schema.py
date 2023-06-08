from typing import Literal, Optional

from pydantic import BaseModel


class ConfigSchema(BaseModel):
    fontsize: int = 20
    block_width: int = 100
    offset: int = 10
    start: int = 0


class AnnotationSchema(BaseModel):
    start: int
    end: int
    color: str = "blue"
    text: Optional[str] = None
    position: Literal["top", "bottom"] = "top"


class SequenceSchema(BaseModel):
    id: str
    sequence: str
    reference: bool = False
    annotations: list[AnnotationSchema] = []


class InputSchema(BaseModel):
    config: ConfigSchema
    sequences: list[SequenceSchema]
