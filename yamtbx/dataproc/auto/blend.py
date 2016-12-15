"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

"""
CLUSTERS.txt

to ccp4-6.5 looks like:
 Cluster     Number of         Cluster         LCV      aLCV      Datasets
  Number      Datasets          Height                            ID
                               
     001             2           0.061        0.15      0.13      81 203
....

from ccp4-7.0 looks like:
 Cluster     Number of         Cluster         LCV      aLCV   Furthest    Datasets
  Number      Datasets          Height                         Datasets          ID
                               
     001             2           0.061        0.15      0.13    81  203    81 203
....
"""

import os
import numpy
import shutil
from yamtbx.dataproc.xds import xds_ascii
from yamtbx.dataproc.xds import integrate_hkl_as_flex
from yamtbx import util
from cctbx.array_family import flex
from cctbx import miller
from libtbx.utils import null_out

blend_comm = "blend"
sg_to_cent = {1:"P",2:"P",3:"P",4:"P",5:"C",6:"P",7:"P",8:"C",9:"C",10:"P",\
11:"P",12:"C",13:"P",14:"P",15:"C",16:"P",17:"P",18:"P",19:"P",20:"C",\
21:"C",22:"F",23:"I",24:"I",25:"P",26:"P",27:"P",28:"P",29:"P",30:"P",\
31:"P",32:"P",33:"P",34:"P",35:"C",36:"C",37:"C",38:"A",39:"A",40:"A",\
41:"A",42:"F",43:"F",44:"I",45:"I",46:"I",47:"P",48:"P",49:"P",50:"P",\
51:"P",52:"P",53:"P",54:"P",55:"P",56:"P",57:"P",58:"P",59:"P",60:"P",\
61:"P",62:"P",63:"C",64:"C",65:"C",66:"C",67:"C",68:"C",69:"F",70:"F",\
71:"I",72:"I",73:"I",74:"I",75:"P",76:"P",77:"P",78:"P",79:"I",80:"I",\
81:"P",82:"I",83:"P",84:"P",85:"P",86:"P",87:"I",88:"I",89:"P",90:"P",\
91:"P",92:"P",93:"P",94:"P",95:"P",96:"P",97:"I",98:"I",99:"P",100:"P",\
101:"P",102:"P",103:"P",104:"P",105:"P",106:"P",107:"I",108:"I",109:"I",110:"I",\
111:"P",112:"P",113:"P",114:"P",115:"P",116:"P",117:"P",118:"P",119:"I",120:"I",\
121:"I",122:"I",123:"P",124:"P",125:"P",126:"P",127:"P",128:"P",129:"P",130:"P",\
131:"P",132:"P",133:"P",134:"P",135:"P",136:"P",137:"P",138:"P",139:"I",140:"I",\
141:"I",142:"I",143:"P",144:"P",145:"P",146:"R",147:"P",148:"R",149:"P",150:"P",\
151:"P",152:"P",153:"P",154:"P",155:"R",156:"P",157:"P",158:"P",159:"P",160:"R",\
161:"R",162:"P",163:"P",164:"P",165:"P",166:"R",167:"R",168:"P",169:"P",170:"P",\
171:"P",172:"P",173:"P",174:"P",175:"P",176:"P",177:"P",178:"P",179:"P",180:"P",\
181:"P",182:"P",183:"P",184:"P",185:"P",186:"P",187:"P",188:"P",189:"P",190:"P",\
191:"P",192:"P",193:"P",194:"P",195:"P",196:"F",197:"I",198:"P",199:"I",200:"P",\
201:"P",202:"F",203:"F",204:"I",205:"P",206:"I",207:"P",208:"P",209:"F",210:"F",\
211:"I",212:"P",213:"P",214:"I",215:"P",216:"F",217:"I",218:"P",219:"F",220:"I",\
221:"P",222:"P",223:"P",224:"P",225:"F",226:"F",227:"F",228:"F",229:"I",230:"I"}

def run_blend(wdir, xds_ascii_files, logout="blend_a.log"):
    ofs_lst = open(os.path.join(wdir, "files.lst"), "w")
    map(lambda x: ofs_lst.write(os.path.relpath(x, wdir)+"\n"), xds_ascii_files)
    ofs_lst.close()
    util.call(blend_comm, "-aDO files.lst", stdin="tolerance 10\n",
              wdir=wdir,
              stdout=open(os.path.join(wdir, logout), "w"))
# run_blend()

def run_blend0R(wdir, xds_ascii_files, logout="blend0.log"):
    ofs_cell = open(os.path.join(wdir, "forR_macropar.dat"), "w")
    ofs_lst = open(os.path.join(wdir, "xds_lookup_table.txt"), "w")
    ofs_files = open(os.path.join(wdir, "NEW_list_of_files.dat"), "w") # Just to avoid error
    for i, xac in enumerate(xds_ascii_files):
        symm = xds_ascii.XDS_ASCII(xac, read_data=False).symm
        cell = " ".join(map(lambda x: "%7.3f"%x, symm.unit_cell().parameters()))
        sgnum =   int(symm.space_group().type().number())
        ofs_lst.write("dataset_%.3d %s\n" % (i, xac))
        ofs_cell.write("%4d %s 0 0 0 \"%s\"\n" % (i+1, cell, sg_to_cent[sgnum]))
        ofs_files.write("%s\n"%xac)

    ofs_lst.close()
    ofs_cell.close()
    ofs_files.close()

    shutil.copyfile(os.path.join(wdir, "forR_macropar.dat"),
                    os.path.join(wdir, "forR_macropar.dat.bak"))

    pathblend = "modules/yamtbx/blend/R/blend0.R"
    useNC = os.getenv("BLEND_USE_NCDIST","no")
    if useNC != "no":
        pathblend = "modules/yamtbx/blend/R/blend0NC.R"

    util.call("Rscript", os.path.join(os.environ["PHENIX"], pathblend),
              wdir=wdir,
              stdout=open(os.path.join(wdir, logout), "w"))

    open(os.path.join(wdir, "hctojson.R"), "w").write("""\
# reference: http://www.coppelia.io/2014/07/converting-an-r-hclust-object-into-a-d3-js-dendrogram/
library(rjson)
HCtoJSON<-function(hc){
  labels<-hc$labels
  merge<-data.frame(hc$merge)
  for (i in (1:nrow(merge))) {
    if (merge[i,1]<0 & merge[i,2]<0) {eval(parse(text=paste0("node", i, "<-list(name=\\"", i, "\\", children=list(list(name=labels[-merge[i,1]]),list(name=labels[-merge[i,2]])))")))}
    else if (merge[i,1]>0 & merge[i,2]<0) {eval(parse(text=paste0("node", i, "<-list(name=\\"", i, "\\", children=list(node", merge[i,1], ", list(name=labels[-merge[i,2]])))")))}
    else if (merge[i,1]<0 & merge[i,2]>0) {eval(parse(text=paste0("node", i, "<-list(name=\\"", i, "\\", children=list(list(name=labels[-merge[i,1]]), node", merge[i,2],"))")))}
    else if (merge[i,1]>0 & merge[i,2]>0) {eval(parse(text=paste0("node", i, "<-list(name=\\"", i, "\\", children=list(node",merge[i,1] , ", node" , merge[i,2]," ))")))}
  }
  eval(parse(text=paste0("JSON<-toJSON(node",nrow(merge), ")")))
  return(JSON)
}

load("BLEND0.RData",.GlobalEnv)
JSON<-HCtoJSON(npar.hc_ward)
cat(JSON, file="dendro.json")
""")
    util.call("Rscript", "hctojson.R",
              wdir=wdir,
              stdout=open(os.path.join(wdir, logout), "a"))

# run_blend0R()

def load_xds_data_only_indices(xac_files, d_min=None):
    miller_sets = {}
    for f in xac_files:
        if xds_ascii.is_xds_ascii(f):
            print "Loading", f
            ms = xds_ascii.XDS_ASCII(f, i_only=True).as_miller_set()
            miller_sets[f] = ms.resolution_filter(d_min=d_min)
        elif integrate_hkl_as_flex.is_integrate_hkl(f):
            print "Sorry, skipping", f
        else:
            print "Skipping unrecognized:", f
    return miller_sets
# load_xds_data_only_indices()

class BlendClusters:
    def __init__(self, workdir=None, d_min=None, load_results=True):
        self.workdir = workdir
        self.d_min = d_min

        self.clusters = {} # clno -> (cluster_height, LCV, aLCV, IDs)
        self.files = None

        if load_results:
            self.read_results()
            self.miller_sets = load_xds_data_only_indices(xac_files=self.files, d_min=self.d_min)
    # __init__()

    def read_results(self):
        self.clusters = {}

        lookup_table_txt = os.path.join(self.workdir, "xds_lookup_table.txt")
        self.files = map(lambda l: l.split()[1], open(lookup_table_txt))
        
        clusters_txt = os.path.join(self.workdir, "CLUSTERS.txt")
        ifs = open(clusters_txt)
        ifs.readline() # skip first (blank)
        firstline = ifs.readline()
        if firstline.startswith(" Cluster     Number of         Cluster         LCV      aLCV      Datasets"):
            with_furhest = False
        elif firstline.startswith(" Cluster     Number of         Cluster         LCV      aLCV   Furthest    Datasets"):
            with_furhest = True
        else:
            raise Exception("Unexpected header of CLUSTERS.txt")
        for i in xrange(2): ifs.readline()
        for l in ifs:
            sp = l.split()
            if with_furhest:
                clno, num, clh, lcv, alcv, f1, f2 = sp[:7]
                ids = map(int, sp[7:])
            else:
                clno, num, clh, lcv, alcv = sp[:5]
                ids = map(int, sp[5:])
            assert int(num) == len(ids)
            assert max(ids) <= len(self.files)
            self.clusters[int(clno)] = (float(clh), float(lcv), float(alcv), ids)
    # read_results()

    def cluster_completeness(self, clno, anomalous_flag, calc_redundancy=True):
        if clno not in self.clusters:
            print "Cluster No. %d not found" % clno
            return

        cls = self.clusters[clno][3]
        msets = map(lambda x: self.miller_sets[self.files[x-1]], cls)
        num_idx = sum(map(lambda x: x.size(), msets))
        all_idx = flex.miller_index()
        all_idx.reserve(num_idx)
        for mset in msets: all_idx.extend(mset.indices())

        # Calc median cell
        cells = numpy.array(map(lambda x: x.unit_cell().parameters(), msets))
        median_cell = map(lambda i: numpy.median(cells[:,i]), xrange(6))
        symm = msets[0].customized_copy(unit_cell=median_cell)

        assert anomalous_flag is not None
        # XXX all must belong to the same Laue group and appropriately reindexed..
        all_set = miller.set(indices=all_idx,
                             crystal_symmetry=symm, anomalous_flag=anomalous_flag)
        all_set = all_set.resolution_filter(d_min=self.d_min)

        # dummy for redundancy calculation. dirty way..
        if calc_redundancy:
            dummy_array = miller.array(miller_set=all_set, data=flex.int(all_set.size()))
            merge = dummy_array.merge_equivalents()
            cmpl = merge.array().completeness()
            redun = merge.redundancies().as_double().mean()
            return cmpl, redun
        else:
            cmpl = all_set.unique_under_symmetry().completeness()
            return cmpl
    # cluster_completeness()

    def show_cluster_summary(self, out=null_out()):
        tmp = []
        for clno in self.clusters:
            cluster_height, LCV, aLCV, IDs = self.clusters[clno]
            cmpl, redun = self.cluster_completeness(clno, anomalous_flag=False)
            acmpl, aredun = self.cluster_completeness(clno, anomalous_flag=True)
            tmp.append((clno, IDs, cluster_height, cmpl*100., redun, acmpl*100., aredun, LCV, aLCV))

        tmp.sort(key=lambda x: (-x[4], -x[3])) # redundancy & completeness
        out.write("# d_min= %.3f\n" % (self.d_min))
        out.write("# Sorted by redundancy and completeness\n")
        out.write("Cluster Number   CLh   Cmpl Redun  ACmpl ARedun  LCV  aLCV\n")
        for clno, IDs, clh, cmpl, redun, acmpl, aredun, LCV, aLCV in tmp:
            out.write("%7d %6d %5.1f %6.2f %5.1f %6.2f %5.1f %5.1f %5.1f\n" % (clno, len(IDs), clh, cmpl, redun,
                                                                               acmpl, aredun, LCV, aLCV))
        return tmp
    # show_cluster_summary()

# class BlendClusters


if __name__ == "__main__":
    import sys

    bc = BlendClusters(sys.argv[1], d_min=3.)
    bc.show_cluster_summary(sys.stdout)
