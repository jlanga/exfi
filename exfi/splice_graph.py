#!/usr/bin/env python3

"""
exfi.SpliceGraph: Class to store and manipulate a splice splice graph.


"""

import networkx as nx

class SpliceGraphComponent(object):
    """
    A SpliceGraphComponent is a nx.DiGraph where:
    - keys are the connected component name (a transcript or a gene) and
    - values are nx.DiGraph objects where
        - nodes are putative exons
        - nodes have an attribute called "coordinates" in bed3 notatio (seqid, start, end),
            0-indexed
        - edges are the connection between putative exons
        - edges have an attribute called "overlap" that shows how many nucleotides seem to be shared
            between the two connected exons. This number is an int:
            - if 0, then when one exon ends, the next starts (ideal case)
            - < 0 means that there is a gap that is unsolved by the build_splice_graph functions:
                - one or more consecutive exons, too short for the method to be detected
                - uncovered regions by the WGS experiment
                - too many snps at the end of a transcript
            - > 0 means that there is an overlap between the two exons. This happens when the
                intronic sequence and the next exons start with the same set of nucleotides
    """

    def __init__(self, node2coord=None, link2overlap=None):
        """Initialization with nodes, node attributes, edges and edge attributes"""
        self.graph = nx.DiGraph()
        if node2coord:
            self.graph.add_nodes_from(node2coord.keys())
            nx.set_node_attributes(
                G=self.graph, name="coordinates", values=node2coord
            )
        if link2overlap:
            self.graph.add_edges_from(link2overlap.keys())
            nx.set_edge_attributes(
                G=self.graph, name="overlaps", values=link2overlap
            )

    def add_nodes_from(self, node2coord=None):
        """(self, dict) -> None

        Add/overwirte nodes in node2coord = {node_id: (str, int, int)}
        """
        if node2coord:
            node2coord_generator = (
                (node_id, {"coordinates": coordinates})
                for node_id, coordinates in node2coord.items()
            )
            nx.DiGraph.graph.add_nodes_from(
                self.graph,
                node2coord_generator
            )

    def add_edges_from(self, edge2overlap=None):
        """(self, dict) -> None

        Add/overwrite edges in edge2overlap = {(str, str): int}. Add nodes if necessary.
        """
        if edge2overlap:
            nx.DiGraph.add_edges_from(
                self.graph,
                ebunch=edge2overlap.keys(),
                attr_dict={"overlap": edge2overlap}
            )

    #def remove_nodes_from()


class SpliceGraph(object):
    """
    A SpliceGraph is a dict of SplceGraphComponents
    """

    components = {}

    def __init__(self):
        """Init alone, no nodes, no nothing"""
        self.components = {}

    def add_component(self, component_name, splice_graph_component):
        """(str, nx.Digraph) ->
        Add/overwrite component in SpliceGraph
        """
        self.components[component_name] = splice_graph_component

    def delete_component(self, component_name):
        """(str) -> None
        Add/overwrite component in SpliceGraph
        """
        self.components.pop(component_name)

    @classmethod
    def get_component(cls, self, component_name):
        """(str) -> nx.DiGraph

        Get subgraph
        """
        return self.components[component_name]
