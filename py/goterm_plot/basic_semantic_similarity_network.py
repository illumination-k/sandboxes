import abc
import numpy as np
import numpy.typing as npt

import networkx as nx

from typing import Optional
from goatools.obo_parser import GODag
from goatools.semantic import semantic_similarity, resnik_sim, lin_sim, TermCounts
from goatools.semsim.termwise.wang import SsWang

godag = GODag("go-basic.obo")

class IGoTermSemanticSimilarity:
    @abc.abstractmethod
    def get_semantic_similarity(self, term1: str, term2: str) -> float:
        raise NotImplementedError()
    
class BasicSemanticSimilarity(IGoTermSemanticSimilarity):
    def __init__(self, godag: GODag) -> None:
        super().__init__()

        self.godag = godag

    def get_semantic_similarity(self, term1: str, term2: str) -> float:
        return semantic_similarity(term1, term2, self.godag)
    
class ResnikSemanticSimilarity(IGoTermSemanticSimilarity):
    def __init__(self, godag: GODag, assosiactions: any) -> None:
        super().__init__()
        self.godag = godag
        self.term_count = TermCounts(godag, assosiactions)

    def get_semantic_similarity(self, term1: str, term2: str) -> float:
        return resnik_sim(term1, term2, self.godag, self.term_count)
    

class LinSemanticSimilarity(IGoTermSemanticSimilarity):
    def __init__(self, godag: GODag, assosiactions: any) -> None:
        super().__init__()
        self.godag = godag
        self.term_count = TermCounts(godag, assosiactions)

    def get_semantic_similarity(self, term1: str, term2: str) -> float:
        return lin_sim(term1, term2, self.godag, self.term_count)
    
class WangSeminaticSimilarity(IGoTermSemanticSimilarity):
    def __init__(self, godag: GODag, terms: list[str], relationships: Optional[dict] = None) -> None:
        super().__init__()
        self.godag = godag
        self.wang = SsWang(goids=terms, godag=godag, relationships=relationships)

    def get_semantic_similarity(self, term1: str, term2: str) -> float:
        return self.wang.get_sim(term1, term2)

def semantic_similarity_matrix(terms: list[str], similarity: IGoTermSemanticSimilarity) -> npt.NDArray[np.float64]:
    n = len(terms)
    matrix = np.zeros((n, n), dtype=np.float64)

    for i in range(n):
        for j in range(n):
            matrix[i][j] = similarity.get_semantic_similarity(terms[i], terms[j])

    return matrix


class SimilarityGraph:
    def __init__(self, terms: list[str], seed: int = 0) -> None:
            self.g = nx.Graph()
            self.sim_matrix = semantic_similarity_matrix(terms)
            n = len(terms)

            for i in range(n):
                for j in range(j):
                    self.g.add_edge(terms[i], terms[j], self.sim_matrix[i][j])

            self.layout = nx.layout.spring_layout(self.g)
            self.cluster = nx.community.louvain_communities(self.g, seed)