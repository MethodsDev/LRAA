import sys, os, re
from GenomeFeature import GenomeFeature
import PASA_scored_path

class Transcript (GenomeFeature):

    trans_id_counter = 0
    

    def __init__(self, contig_acc, segment_coordinates_list, orient):

        segment_coordinates_list = sorted(segment_coordinates_list, key=lambda x: x[0])

        trans_lend = segment_coordinates_list[0][0]
        trans_rend = segment_coordinates_list[-1][1]

        super().__init__(contig_acc, trans_lend, trans_rend)

        self._orient = orient
        self._exon_segments = segment_coordinates_list

        Transcript.trans_id_counter += 1
        
        self._id = ":".join([contig_acc, str(trans_lend), str(trans_rend), "N:{}".format(Transcript.trans_id_counter), orient])

        self._gene_id = "g.{}".format(self._id)

        self._meta = None

        self._scored_path_obj = None #optional - useful if transcript obj was built based on a scored path

        self.read_names = list() # list of read names supporting the transcript structure.
        
        return



    def __repr__(self):

        return("Transcript: {} {}-{} [{}] segs: {}".format(self._contig_acc,
                                                           self._lend,
                                                           self._rend,
                                                           self._orient,
                                                           self._exon_segments) )
        
        
    def set_scored_path_obj(self, scored_path_obj):
        assert(type(scored_path_obj) == PASA_scored_path.PASA_scored_path)
        self._scored_path_obj = scored_path_obj
        
        
    def set_gene_id(self, gene_id):
        self._gene_id = gene_id


    def add_meta(self, meta_key, meta_val):
        if self._meta == None:
            self._meta = dict()

        self._meta[meta_key] = meta_val
        return


    def add_read_names(self, read_names):
        self.read_names = read_names
        
    

    def to_GTF_format(self):

        ## transcript line:

        gtf_text = ""
        
        if self.read_names:
            gtf_text = f"#{self._id}\t" + ",".join(self.read_names) + "\n"
        
        gtf_text += "\t".join([self._contig_acc,
                              "PASA-SALRAA",
                              "transcript",
                              str(self._lend),
                              str(self._rend),
                              ".",
                              self._orient,
                              ".",
                              "gene_id \"{}\"; transcript_id \"{}\";".format(self._gene_id, self._id)])

        if self._meta:
            for meta_key in sorted(self._meta.keys()):
                gtf_text += " {} \"{}\";".format(meta_key, self._meta[meta_key])

        gtf_text += "\n"

        for segment in self._exon_segments:
            gtf_text += "\t".join([self._contig_acc,
                                   "PASA-SALRAA",
                                   "exon",
                                   str(segment[0]),
                                   str(segment[1]),
                                   ".",
                                   self._orient,
                                   ".",
                                   "gene_id \"{}\"; transcript_id \"{}\";".format(self._gene_id, self._id)]) + "\n"


        if self._scored_path_obj:
            # include construction info as comments.
            gtf_text += "# derived from scored path obj:\n"

            mpgn_list = self._scored_path_obj.get_path_mpgn_list()
            for mpgn in mpgn_list:
                gtf_text += "# " + str(mpgn) + "\n"

            gtf_text += "\n"

        
        return gtf_text
        
        
        
