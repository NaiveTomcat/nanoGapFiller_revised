import networkx as nx
from Bio.Seq import Seq

from fastg_reader import parse_fastg


def find_site_locations(seq: Seq, site_seqs: list[Seq]) -> list[int]:
    locations: list[int] = []
    for site_seq in site_seqs:
        local_locations: list[int] = []
        if seq.find(site_seq) != -1:
            local_locations.append(seq.find(site_seq))
        while seq.find(site_seq, start=local_locations[-1]) != -1:
            local_locations.append(
                seq.find(site_seq, start=local_locations[-1]))
        locations += local_locations
    locations.sort()
    return locations


def assembly2site(assembly: nx.DiGraph, site_seqs: list[Seq]) -> nx.DiGraph:
    """Calculate site graph based on the assembly graph and site seqences given.

        Parameters
        ----------
        assembly : networkx.DiGraph
            A DiGraph which contains the assembly graph for calculation.
        site_seqs : list[Bio.Seq.Seq]
            A list of site sequences

        Returns
        -------
        networkx.DiGraph
            A DiGraph instance of calculated site graph for further uses.     
    """
    site_graph = nx.DiGraph()
    # Create site paths for each contig in the assembly
    for contig_name in assembly:
        seq: Seq = assembly.nodes(data='seq')[contig_name]
        locations = find_site_locations(seq, site_seqs)
        site_graph.add_node(contig_name+'_start',
                            placeholder=True, contig=contig_name)
        site_graph.add_node(contig_name+'_end',
                            placeholder=True, contig=contig_name)
        if locations == []:
            site_graph.add_edge(contig_name+'_start',
                                contig_name+'_end', distance=0)
            continue
        sites: list[str] = [contig_name + '_' + str(x) for x in locations]
        site_graph.add_nodes_from(sites, contig=contig_name)
        distances = [locations[0]] + [locations[i+1] - locations[i]
                                      for i in range(len(locations) - 1)] + [len(seq) - locations[-1]]
        site_graph.add_edge(contig_name+'_start',
                            sites[0], distance=distances[0])
        for i in range(1, len(sites)):
            site_graph.add_edge(sites[i-1], sites[i], distance=distances[i])
        site_graph.add_edge(sites[-1], contig_name +
                            '_end', distance=distances[-1])

    # Connect the end marker of each contig in site graph with the start marker
    # of the succeeding contig
    for contig_name in assembly:
        for successor_name in assembly.succ[contig_name]:
            site_graph.add_edge(contig_name+'_end',
                                successor_name+'_start', distance=0)

    # Get rig of start and end markers and connect the graph with re-calculated
    # distance
    # Connect the edges first
    for contig_name in assembly:
        last_site = next(iter(site_graph.pred[contig_name+'_end']))
        next_start_marker = next(iter(site_graph.succ[contig_name+'_end']))
        next_site = next(iter(site_graph.succ[next_start_marker]))
        distance = site_graph.edges[(last_site, contig_name+'_end')]['distance'] + \
            site_graph.edges[(next_start_marker, next_site)]['distance']
        site_graph.add_edge(last_site, next_site, distance=distance)
    # Then delete the marker edge
    site_graph.remove_nodes_from(
        [x for x, y in site_graph.nodes(date=True) if y['placeholder']])

    # Should the sitegraph been completed now
    return site_graph

def fastg2site(f: str, site_seqs: list[Seq]) -> nx.DiGraph:
    """Calculate site graph from given fastg file and site sequences.

    Parameters
    ----------
    f : str
        Path to SPADES generated fastg file
    site_seqs : list[Bio.Seq.Seq]
        List of site sequences

    Returns
    -------
    networkx.DiGraph
        A DiGraph instance of calculated site graph for further uses.
    """
    return assembly2site(parse_fastg(f), site_seqs)
