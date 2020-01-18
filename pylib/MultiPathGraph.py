import sys, os, re
import Splice_graph
import MultiPath
import MultiPathCounter
import networkx as nx
import Simple_path_utils as Simple_path_utils
from PASA_SALRAA_Globals import SPACER
from GenomeFeature import Intron, Exon
from MultiPathGraphNode import MultiPathGraphNode
from collections import defaultdict

import logging

logger = logging.getLogger(__name__)

    
class MultiPathGraph:


    def __init__(self, multiPathCounter, splice_graph, min_mpgn_read_count=1, allow_spacers=False): # //FIXME: enable interface to allow_spacers

        assert(type(multiPathCounter) == MultiPathCounter.MultiPathCounter)
        assert(type(splice_graph) == Splice_graph.Splice_graph)
        
        self._splice_graph = splice_graph
        
        mp_graph = nx.DiGraph()
        self._mp_graph = mp_graph

        self._mp_graph_nodes_list = list()

        self._incompatible_pairs = set() # store sets of pairs of incompatible nodes.
        
        multiPathCountPairs = multiPathCounter.get_all_MultiPathCountPairs()
        for mpCountPair in multiPathCountPairs:
            mp, count = mpCountPair.get_multipath_and_count()

            if count < min_mpgn_read_count:
                continue

            if (not allow_spacers) and SPACER in mp:
                continue
            
            first_node_id = mp[0]
            last_node_id = mp[-1]

            first_node_obj = splice_graph.get_node_obj_via_id(first_node_id)
            last_node_obj = splice_graph.get_node_obj_via_id(last_node_id)

            lend_coord = first_node_obj.get_coords()[0]
            rend_coord = last_node_obj.get_coords()[1]

            mp_graph_node = MultiPathGraphNode(mp, count, lend_coord, rend_coord, mpg=self)
            mp_graph.add_node(mp_graph_node)

            self._mp_graph_nodes_list.append(mp_graph_node)

            
        ## sort by rend
        self._mp_graph_nodes_list = sorted(self._mp_graph_nodes_list, key=lambda x: x._rend)

        ordered_nodes = self._mp_graph_nodes_list
        
        ## define edges, containments, and incompatibilities
        for i in range(0, len(ordered_nodes)):
            node_i = ordered_nodes[i]

            for j in range(i-1, -1, -1):
                node_j = ordered_nodes[j]

                #print("comparing {},{}".format(node_j, node_i))
                
                if node_j._rend < node_i._lend:
                    break # all earlier node j's will also be non-overlapping

                if node_i.contains_other_node(node_j):
                    # i contains j
                    node_i.add_containment(node_j)
                elif node_j.contains_other_node(node_i):
                    # j contains i
                    node_j.add_containment(node_i)
                elif node_i.compatible(node_j):
                    # draw edge between overlapping and compatible nodes.
                    self._mp_graph.add_edge(node_j, node_i)
                    logger.debug("adding edge: {},{}".format(node_j, node_i))
                else:
                    # incompatible pairs
                    incompatible_pair_token = MultiPathGraphNode.get_mpgn_pair_token(node_i, node_j)
                    self._incompatible_pairs.add(incompatible_pair_token)
                    
    
    def get_ordered_nodes(self):
        # these are sorted by rend
        return list(self._mp_graph_nodes_list)

    def init_mpgn_reweighting_flags(self):
        mpgn_nodes = self.get_ordered_nodes()
        for mpgn in mpgn_nodes:
            mpgn.set_reweighted_flag(False)
        return
    
    
    def has_edge(self, multiPath_before, multiPath_after):
        return self._mp_graph.has_edge(multiPath_before, multiPath_after)


    def successors(self, mpgn):
        return self._mp_graph.successors(mpgn)

    def predecessors(self, mpgn):
        return self._mp_graph.predecessors(mpgn)

    
    def incompatible_mpgn_pair(self, mpgn_A, mpgn_B):
        mpgn_pair_token = MultiPathGraphNode.get_mpgn_pair_token(mpgn_A, mpgn_B)
        if mpgn_pair_token in self._incompatible_pairs:
            return True
        else:
            return False

        
    def get_splice_graph(self):
        return self._splice_graph


    def define_disjoint_graph_components_via_graph_traversal(self):

        mpgn_list = self.get_ordered_nodes()
        mpgn_seen = set()

        component_list = list()
        
        while len(mpgn_list) > 0:

            queue = list()
            seed_node = mpgn_list.pop(0)

            if seed_node in mpgn_seen:
                continue

            # start a new component.
            component = list()
            queue.append(seed_node)

            while len(queue) > 0:
                node = queue.pop(0)
                if node not in mpgn_seen:
                    mpgn_seen.add(node)
                    component.append(node)
                    # add predecessors and successors to queue
                    if node.has_predecessors():
                        queue.extend(node.get_predecessors())
                    if node.has_successors():
                        queue.extend(node.get_successors())
                    if node.has_containments():
                        queue.extend(node.get_containments())
            
            if len(component) > 0:
                component_list.append(component)

        logger.info("identified {} disjoint graph components".format(len(component_list)))

        ## assign the component ids
        component_counter = 0
        for component in component_list:
            component_counter += 1
            #print("component: {}, counter: {}".format(component, component_counter))
            for mpgn in component:
                mpgn.set_component_id(component_counter)
        
        return component_list


    def define_disjoint_graph_components_via_shared_splice_graph_vertex(self):

        all_splice_graph_nodes = set()
        splice_graph_node_to_mpgns = defaultdict(set)

        # populate splice graph node to mpgn data structures
        mpgn_list = self.get_ordered_nodes()
        for mpgn in mpgn_list:
            splice_graph_nodes = mpgn.get_splice_graph_node_objs_for_path()
            for splice_graph_node in splice_graph_nodes:
                all_splice_graph_nodes.add(splice_graph_node)
                splice_graph_node_to_mpgns[splice_graph_node].add(mpgn)


        # build disjoint components
        splice_graph_nodes_visited = set()

        component_list = list()
        all_splice_graph_nodes = list(all_splice_graph_nodes) # convert earlier set to list
        while len(all_splice_graph_nodes) > 0:

            queue = list()
            seed_splice_graph_node = all_splice_graph_nodes.pop()

            if seed_splice_graph_node in splice_graph_nodes_visited:
                continue

            # start a new component
            component = list()
            queue.append(seed_splice_graph_node)

            while len(queue) > 0:
                splice_graph_node = queue.pop(0)
                mpgns = list(splice_graph_node_to_mpgns[splice_graph_node])
                splice_graph_nodes_visited.add(splice_graph_node)
                for mpgn in mpgns:
                    if mpgn not in component:
                        component.append(mpgn)
                        
                    mpgn_splice_graph_nodes = mpgn.get_splice_graph_node_objs_for_path()
                    for splice_graph_node in mpgn_splice_graph_nodes:
                        if ( (splice_graph_node not in splice_graph_nodes_visited)
                            and
                            (splice_graph_node not in queue) ):
                            queue.append(splice_graph_node)


            component_list.append(component)

        logger.info("identified {} disjoint graph components".format(len(component_list)))

        ## assign the component ids
        component_counter = 0
        for component in component_list:
            component_counter += 1
            #print("component: {}, counter: {}".format(component, component_counter))
            for mpgn in component:
                mpgn.set_component_id(component_counter)

            
        return component_list
                            
    
    
    def describe_graph(self, output_filename):

        ofh = open(output_filename, 'wt')

        mpgn_list = self.get_ordered_nodes()
        for mpgn in mpgn_list:
            ofh.write(str(mpgn) + "\n")
            for containment in mpgn.get_containments():
                ofh.write("\t" + str(containment) + "\n")
            ofh.write("\n")

        ofh.close()

        return
    

    def write_mp_graph_nodes_to_gtf(self, gtf_output_filename):

        contig_acc = self._splice_graph.get_contig_acc()
        
        ofh = open(gtf_output_filename, 'wt')

        mpgn_list = self.get_ordered_nodes()
        for mpgn in mpgn_list:
            splice_graph_nodes = mpgn.get_splice_graph_node_objs_for_path()

            mpgn_read_count = mpgn.get_count()
            mpgn_component_id = mpgn.get_component_id()
            
            mpgn_id = mpgn.get_id()
            trans_id = "t__count={}_Comp={}_.".format(mpgn_read_count, mpgn_component_id) + mpgn_id
            gene_id = "g__count={}_Comp={}_".format(mpgn_read_count, mpgn_component_id) + mpgn_id
            
            for splice_graph_node in splice_graph_nodes:
                if splice_graph_node is not None and type(splice_graph_node) == Exon:

                    coords = splice_graph_node.get_coords()
                    
                    ofh.write("\t".join([contig_acc,
                                        "MPGN",
                                        "exon",
                                         str(coords[0]),
                                         str(coords[1]),
                                         ".",
                                         "?",
                                         ".",
                                         "gene_id \"{}\"; transcript_id \"{}\";".format(gene_id, trans_id)]) + "\n")
                

        

        ofh.close()
        

        return


    def remove_small_components(self, mpg_components, min_transcript_length):


        surviving_components = list()
        
        for mpgn_list in mpg_components:
            max_seq_len = 0
            for mpgn in mpgn_list:
                max_seq_len += mpgn.get_seq_length()
            if max_seq_len < min_transcript_length:
                # component is too small to generate a sufficiently large transcript
                self._mp_graph.remove_nodes_from(mpgn_list)
            else:
                surviving_components.append(mpgn_list)

        return surviving_components
    
