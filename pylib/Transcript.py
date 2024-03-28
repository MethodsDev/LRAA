import sys, os, re
from collections import defaultdict
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

        self._multipath = None # multipath obj

        self._simplepath = None
        
        return

    ## getters

    def get_exon_segments(self):
        return self._exon_segments.copy()

    def get_strand(self):
        return(self._orient)

    def get_orient(self):
        return(self.get_strand())
   

    def get_transcript_id(self):
        if 'transcript_id' in self._meta:
            return(self._meta['transcript_id'])
        else:
            return(self._id)

    def get_gene_id(self):
        if 'gene_id' in self._meta:
            return(self._meta['gene_id'])
        else:
            return(self._gene_id)
    
    
    def __repr__(self):

        text = "Transcript: {} {}-{} [{}] segs: {}".format(self._contig_acc,
                                                           self._lend,
                                                           self._rend,
                                                           self._orient,
                                                           self._exon_segments)

        if self._meta is not None:
            text += "\t" + str(self._meta)

        return text
        
        
    def set_scored_path_obj(self, scored_path_obj):
        assert(type(scored_path_obj) == PASA_scored_path.PASA_scored_path)
        self._scored_path_obj = scored_path_obj
        
        
    def set_gene_id(self, gene_id):
        self._gene_id = gene_id


    def add_meta(self, meta_key, meta_val=None):
        
        if self._meta == None:
            self._meta = dict()


        if meta_val is None and type(meta_key)==dict:
            self._meta = meta_key.copy()

        elif meta_val is not None:
            self._meta[meta_key] = meta_val
        else:
            raise RuntimeError("Error: not sure how to handle input params")
            
        return


    def add_read_names(self, read_names):
        if self.read_names == None:
            self.read_names = list()

        if type(read_names) == list:
            self.read_names.extend(read_names)
        else:
            self.read_names.append(read_names)
        
    

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
        
        
        
class GTF_contig_to_transcripts:

    @classmethod
    def parse_GTF_to_Transcripts(cls, gtf_filename):

        gene_id_to_meta = defaultdict(dict)
        transcript_id_to_meta = defaultdict(dict)
        transcript_id_to_genome_info = defaultdict(dict)
        
        with open(gtf_filename, "rt") as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                line = line.rstrip()
                vals = line.split("\t")
                contig = vals[0]
                feature_type = vals[2]
                lend = int(vals[3])
                rend = int(vals[4])
                strand = vals[6]
                info = vals[8]
                
                info_dict = cls._parse_info(info)

                if feature_type == 'gene':
                    gene_id = info_dict['gene_id']
                    gene_id_to_meta[gene_id] = info_dict
                    
                if 'transcript_id' not in info_dict:
                    continue

                if feature_type != 'exon':
                    continue

                
                transcript_id = info_dict['transcript_id']
                transcript_id_to_meta[transcript_id] = info_dict

                transcript_id_to_genome_info[transcript_id]['contig'] = contig
                transcript_id_to_genome_info[transcript_id]['strand'] = strand
                if ('coords' in transcript_id_to_genome_info[transcript_id].keys()):
                    transcript_id_to_genome_info[transcript_id]['coords'].append([lend, rend])
                else:
                    transcript_id_to_genome_info[transcript_id]['coords'] = [ [lend,rend] ]

        # convert to transcript objects

        contig_to_transcripts = defaultdict(list)

        for transcript_id in transcript_id_to_genome_info:
            transcript_info_dict = transcript_id_to_genome_info[transcript_id]
            contig = transcript_info_dict['contig']
            strand = transcript_info_dict['strand']
            coords_list = transcript_info_dict['coords']

            transcript_meta = transcript_id_to_meta[transcript_id]
            gene_id = transcript_meta['gene_id']
            gene_meta = gene_id_to_meta[gene_id]
            transcript_meta.update(gene_meta)
            
            transcript_obj = Transcript(contig, coords_list, strand)
            transcript_obj.add_meta(transcript_meta)

            contig_to_transcripts[contig].append(transcript_obj)

            
        return contig_to_transcripts


    
    # private
    @classmethod
    def _parse_info(cls, info):

        info_dict = dict()
        
        parts = info.split(";")
        for part in parts:
            part = part.strip()
            m = re.match('^(\S+) \\"([^\\"]+)\\"', part)
            if m:
                token = m.group(1)
                val = m.group(2)

                info_dict[token] = val

        return info_dict

        
if __name__=='__main__':

    # testing gtf parser
    usage = "usage: {} gtf_filename\n\n".format(sys.argv[0])

    if len(sys.argv) < 2:
        exit(usage)

    gtf_filename = sys.argv[1]
    
    contig_to_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(gtf_filename)

    for contig, transcript_list in contig_to_transcripts.items():
        for transcript_obj in transcript_list:
            print("\t".join([contig, str(transcript_obj)]))


    sys.exit(0)


