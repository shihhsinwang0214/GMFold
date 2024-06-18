import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi
from forgi.graph.bulge_graph import BulgeGraph
from typing import Any, Dict, Tuple
import sys
import os
import copy
from tqdm import tqdm
from operator import itemgetter
from itertools import combinations
import networkx as nx
from networkx.algorithms.clique import find_cliques
from concurrent.futures import ProcessPoolExecutor
import math
from typing import List, Tuple


# Example path you want to add
new_path = r'../src/'
# Check if the path already exists to avoid duplication
if new_path not in sys.path:
    sys.path.append(new_path)
    
from Types import Comp, MultiBranch, BpEnergy, LoopEnergy, Energies
from dna import DNA_ENERGIES
from substructures import substructure_ID
from energy_functions import _hairpin, _stack, _bulge, _internal_loop, _multi_branch, compute_energy
emap = DNA_ENERGIES




class Aptamer_Fold():
    def __init__(self, sequence='ATA', n_tmpl=4, l_fix=4, structures_DB = None, parallel = False):
        self.l_fix= None # number of fixed base pairs in the 5' 3' bonds.
        self.n_tmpl = None
        self.n_wrld = None
        self.sequence = sequence # orginal sequnce in capital letters
        self.sqnc_num = None # original sequence converted in numbers
        self.B_g = None #matrix  storing all posible bonds between base-pairs
        self.dlt_g  = None #matrix sotring backbone bonds
        self.C_dct = {}#candidate dictionary, storing for each node in the tmpl its associated condadidates in wrld
        self.motifs = []  #list of all tmpl graphs found
        self.structures = None
        self.structures_DB =  structures_DB #secondary structures in Dot brachets  notation
        self.energies = []
        self.min_energy_structure_DB = None
        self.min_energy_structures = []
        self.maximal_sets = None
        self.parallel = parallel
        self.min_energy = None
        self.graph = None
        self.graph_m = None
        self.graph_i_m = None
        self.min_motifs_structures = None
        self.bps = []
        self.M = None
        
    def all_combinations(self, lst):
        all_combos = []
        for r in range(2, len(lst) + 1):
            combos = combinations(lst, r)
            all_combos.extend(combos)
        return all_combos
                 
    

    
     
    def init_candidate_dict(self):
        """
        Initializes a dictionary of candidate vertices for motifs.
        Each key represents a motif vertex, and the value is a list of candidate vertices from the world graph.
        """
        for i in range(self.n_tmpl):
            self.C_dct[f'{i}'] = list(range(self.n_wrld))
        return 

    def convert_sequence_to_numeric(self, sequence):
        """
        Converts a nucleotide string into a list of numeric representations.
        """
        mapping = {'A': 0, 'a': 0, 'T': 3, 't': 3, 'C': 1, 'c':1, 'G': 2, 'g':2}
        self.sqnc_num = []
        for char in self.sequence:
            if char in mapping:
                self.sqnc_num.append(mapping[char])
            else:
                print('Invalid character. Only allowed nucleotides are ACTG')
                break
        return 

    def build_BG_matrix(self,):
        """
        Builds a binary adjacency matrix representing the BG motif structure.
        """
        self.B_g = np.zeros((self.n_wrld, self.n_wrld), dtype=int)
        for i in range(self.n_wrld):
            for j in range(i+self.n_tmpl//2 +1,self.n_wrld):
                if (i <= self.l_fix -1 or j <= self.l_fix-1) and (i + j == self.n_wrld - 1):
                    self.B_g[i, j] = 1
                    self.B_g[j, i] = 1
                if self.n_wrld - self.l_fix -1 >= i > self.l_fix -1 and self.n_wrld - self.l_fix -1>= j > self.l_fix -1 and self.sqnc_num[i] + self.sqnc_num[j] in [3] and np.abs(i - j) >= 4:
                    self.B_g[i, j] = 1
                    self.B_g[j, i] = 1
        
        return


    def motif_ngh_dct(self,):
        """
        Defines the neighborhood of each motif node.
        """
        graph = {}
        for i in range(self.n_tmpl):
            if i == 0:
                graph[str(i)] = [1, self.n_tmpl - 1]
            elif i == self.n_tmpl - 1:
                graph[str(i)] = [0, self.n_tmpl - 2]
            elif i == self.n_tmpl // 2 or i == self.n_tmpl // 2 - 1:
                graph[str(i)] = [i - 1, i + 1]
            else:
                graph[str(i)] = [i - 1, i + 1, self.n_tmpl - i - 1]
        return graph

    def backbone_bonds_matrix(self,):
        """
        Constructs the backbone bonds matrix.
        """
        array = np.zeros((self.n_wrld, self.n_wrld), dtype=int)
        array[np.diag_indices(self.n_wrld - 1)[1], np.diag_indices(self.n_wrld - 1)[0] + 1] = 1
        self.dlt_g = array + array.T
        
        return 

    def stats_filter(self,):
        """
        Filters candidate vertices based on statistical criteria.
        """
        ngh_mtf = np.ones(self.n_tmpl, dtype=int) if self.n_tmpl == 2 else 3 * np.ones(self.n_tmpl, dtype=int)
        ngh_mtf[[0, self.n_tmpl - 1, self.n_tmpl // 2, self.n_tmpl // 2 - 1]] = 2
        self.backbone_bonds_matrix()
        ngh_wrld = np.sum(self.dlt_g + self.B_g, axis=1)
        for i in range(self.n_tmpl):
            if i not in [0, self.n_tmpl - 1, self.n_tmpl // 2, self.n_tmpl // 2 - 1]:
                condition = np.where(ngh_wrld < ngh_mtf[i])
                self.C_dct[str(i)] = list(np.delete(np.array(self.C_dct[str(i)]), condition))
        return 

    def topology_filter(self, ):
        """
        Filters candidate vertices based on topological criteria.
        """
        nhg_mtf = self.motif_ngh_dct() 
        D_g = self.dlt_g + self.B_g
        for i in range(self.n_tmpl):
            for j in nhg_mtf[str(i)]:
                indx_remove = np.where(np.sum(D_g[self.C_dct[str(i)]][:,self.C_dct[str(j)]], axis=1) == 0)
                self.C_dct[str(i)] = list(np.delete(np.array(self.C_dct[str(i)]), indx_remove))
        return 

    def apt_filter(self, ):
        """
        Filters candidate vertices based on aptamer criteria.
        """
        self.M= np.zeros((len(self.sequence), len(self.sequence)), dtype = int)
        nhg_mtf = self.motif_ngh_dct()
        range_loop = copy.deepcopy(self.C_dct[str(0)])
        range_loop_2 = copy.deepcopy(self.C_dct[str(self.n_tmpl // 2)])
        for i in range_loop[:]:
            if self.l_fix - 1 <= i <= self.n_wrld - self.l_fix - 1:
                flag = False
                if np.any([i + count not in self.C_dct[str(self.n_tmpl // 2 - k)] for count, k in enumerate(range(self.n_tmpl // 2 - 1, 0, -1), 1)]):
                    for k in range(self.n_tmpl // 2):
                        try:
                            self.C_dct[str(k)].remove(i + k)
                        except ValueError:
                            pass
                else:
                    for j in range_loop_2:
                        if np.any([j + i not in  self.C_dct[str(self.n_tmpl // 2 + i)] for i in range(1, self.n_tmpl // 2)]):
                            for k in range(self.n_tmpl // 2, self.n_tmpl):
                                try:
                                    self.C_dct[str(k)].remove(j + k - self.n_tmpl // 2)
                                except ValueError:
                                    pass
                        elif (j >= i + self.n_tmpl//2 +1 and j <= self.n_wrld - self.l_fix -1 and
                              np.all(np.diag(self.B_g[list(range(i, i + self.n_tmpl // 2))][:, list(range(j + self.n_tmpl // 2 - 1, j - 1, -1))]) == np.ones(self.n_tmpl // 2, dtype=int))):
                            flag = True
                            self.motifs.append([(i + k, j + self.n_tmpl // 2 - 1 - k) for k in range(self.n_tmpl // 2)])
                            for m in self.motifs[-1]:
                                self.M[m[0], m[1]]= 1
                                self.M[m[1], m[0]]= 1
                    if not flag:
                        for k in range(self.n_tmpl // 2):
                            try:
                                self.C_dct[str(k)].remove(i + k)
                            except ValueError:
                                pass
        return 

    def check_segments(self,t1, t2):


        # Extract the endpoints
        t1_start, t1_end = t1
        t2_start, t2_end = t2
      
        
        # Check if one segment is contained in the other
        if (t1_start <= t2_start and t1_end >= t2_end) or (t2_start <= t1_start and t2_end >= t1_end):
            return True
        
        # Check if the segments do not intersect
        if t1_end < t2_start or t2_end < t1_start:
            return True
        
        # Segments intersect
        return False

    def compatible(self, motif_1, motif_2):
            '''
            defines compatibility between two motifs with the same verteces and edges but different
            vertex labels
            
            '''
            '''
            if upper left corner of motif1 is smaller than that of motif2 than we must have that lower right corner of 
            motif1 is bigger than motif2 for compatibility otherwise we have edge crossing bewteen motifs ATTENTION THIS IS NOT
            TRUe FOR EX IN CASE OF MULTI-BRANCHES
            '''
            #if motif_1[0][0] < motif_2[0][0] and motif_1[-1][1] < motif_2[-1][1]:
             #   return False
            #if motif_2[0][0] < motif_1[0][0] and motif_2[-1][1] < motif_1[-1][1]:
             #   return False
            '''
            store all touples in motifs in list check if pairs of motifs contradict each other, if node is
            bonded with two different nodes
            '''
            if motif_1 ==[] or motif_2 ==[]:
                return True
            bp = [('A','T'), ('T','A'), ('G','C'), ('C','G')]
            tuple_list = motif_1 + motif_2
            
          
            compatible = True
            for i in range(len(tuple_list)):
                for j in range(i + 1, len(tuple_list)):
                    if 0 < len(set(tuple_list[i]).intersection(set(tuple_list[j]))) < 2:
                        compatible = False
                        break
                    
                    
            for pair1 in motif_1:
                for pair2 in motif_2:
                    if self.check_segments(pair1, pair2) == False:
                        compatible = False
                        break
            
      
            return compatible

    def compatible_bp(self, bp1, bp2):
        
          
            bp = [('A','T'), ('T','A'), ('G','C'), ('C','G')]
         
            compatible = True
            if 0 < len(set(bp1).intersection(set(bp2))) < 2:
                compatible = False
                
            if self.check_segments(bp1, bp2) == False:
                compatible =  False
            
            return compatible
         
  

# Example usage in your class context
# Assuming `self` is your class instance and the methods/attributes exist in your class
# self.parallel_processing()


    def find_maximal_sets(self, graph):
            '''
            given graph with i-th entry associated with list of indeces of motifs
            compatible to the i-th motif, find the list with maximal compatible motif,
            excludding overlapping components. Thus find the list with maximal compatible bonds.
            '''
            self.maximal_sets = [[0]]#[[0]] # start with list of motifs compatible to the 0-th motif
            graph_m = {str(0): graph[str(0)]} 
         
            # graph storing lists of motifs compatible with all the motifs in the i-th maximal set of motifs ( this case 0-th)
            
            for i in range(1, len(graph)): # loop over all motifs skipping the first one
                flag = 0 
                for j, m in enumerate(self.maximal_sets): # loop over all lists of maximal sets
                    if i in graph_m[str(j)]: # if i-th is compatible with all the elements in the j-th maximal set
                        m.append(i) #add the i-th motif to the j-th maximal set
                        graph_m[str(j)] = list(set(graph_m[str(j)]) & set(graph[str(i)])) # intersect list of of all motifs compatible with the j-th maximal set with the list of motifs compatile to the i-th motif
                        flag = 1
                if flag == 0: # in case the i-th motifs does not fit in any maximal set, that is, it is not bcompatible with at list one motiff in each maximal set, create a new maximal set with initialized with the i-th motif 
                    new_m = [i] 
                    graph_m[f'{len(self.maximal_sets)}'] = graph[str(i)] # lists of motifs compatible with all the elements in the new maximal set ( at the beginning only with the elements of the i-th motif)
                    for j in range(i): # loop over all the motifs previously considered and check wethere they fit in the new maximal set
                        if j in graph_m[f'{len(self.maximal_sets)}']: # if j-th motif is in list of compatible elements of new maximal set
                            new_m.append(j) # include j-th motif in list of new maximal sets
                            graph_m[f'{len(self.maximal_sets)}'] = list(set(graph_m[f'{len(self.maximal_sets)}']) & set(graph[str(j)])) # interesect list of motifs compatible to the new maximal set with that of the motifs compatible to the j-th motifs
                    self.maximal_sets.append(new_m)# add new list, associated with i-th motif,  to the list of lists of maximal set
            
            return  
                

    
    def motifs_ordered_clusters(self, ):
        self.graph_m = {}
        self.graph_i_m = {}
        for count, mot in enumerate(self.motifs):
            if mot ==[]:
                self.graph_m['empty']= [[]]
                self.graph_i_m['empty']= [len(self.motifs)-1]
            else:
                try:
                    self.graph_m[str(mot[0][0])].append(mot)
                    self.graph_i_m[str(mot[0][0])].append(count)
                except KeyError:
                    self.graph_m[str(mot[0][0])] = []
                    self.graph_i_m[str(mot[0][0])] = []
                    self.graph_m[str(mot[0][0])].append(mot)
                    self.graph_i_m[str(mot[0][0])].append(count)
        for key in list(self.graph_m.keys()):
            self.graph_m[key].append([])
            self.graph_i_m[key].append(len(self.motifs))
        return  
    
    def create_dict_d(self,):
        '''
        for each possible length d = j- i, with (i,j) bp, store al bp with length d in the respective entry of the dictionary
        
        '''
        
        self.dict_d = {}
        

        for bp in self.bps:
                    try:
                        self.dict_d[str(bp[1]-bp[0])].add(bp)
                        
                    except KeyError:
                        self.dict_d[str(bp[1]-bp[0])] = set([])
                        
                        self.dict_d[str(bp[1]-bp[0])].add(bp)
        return  
    
    def create_dict_Sij(self, ):
        
        self.dict_Sij = {}
     
        for count, bp1 in enumerate(self.bps):
                if '({}, {})'.format(*bp1) not in list(self.dict_Sij.keys()):
                     self.dict_Sij[str(bp1)]  =[]
                for bp2 in self.bps[count+1:]:
                    
                    if bp1[0]<bp2[0]<bp2[1]<bp1[1]:
                        try:
                            self.dict_Sij[str(bp1)].append(bp2)
                            
                        except KeyError:
                            
                            self.dict_Sij[str(bp1)]  =[]
                            self.dict_Sij[str(bp1)].append(bp2)
                            
                    elif bp2[0]<bp1[0]<bp1[1]<bp2[1]:
                        try:
                            self.dict_Sij[str(bp2)].append(bp1)
                            
                        except KeyError:
                            
                            self.dict_Sij[str(bp2)]  =[]
                            self.dict_Sij[str(bp2)].append(bp1)
   
        return  
        


    
    def bp_ordered_clusters(self, ):
        self.graph_m = {}
        self.graph_i_m = {}
      
        
        for count, bp in enumerate(self.bps):
      
            try:
                self.graph_m[str(bp[0])].add(bp)
                self.graph_i_m[str(bp[0])].add(count)
                #self.graph_m[str(mot[1][0])].append(mot[1])
                #self.graph_i_m[str(mot[0][0])].append(count)
            except KeyError:
                self.graph_m[str(bp[0])] = set([])
                self.graph_i_m[str(bp[0])] = set([])
                self.graph_m[str(bp[0])].add(bp)
                self.graph_i_m[str(bp[0])].add(count)
        
        for key in list(self.graph_m.keys()):
            
            self.graph_m[key] =  list(self.graph_m[key] )
            self.graph_i_m[key] = list(self.graph_i_m[key])
 
        return  

    def greedy_bp4_selection(self,):
        cache = []
        bps = []
        for key in sorted(list(self.dict_d.keys()), key = int):
                bps.extend(self.dict_d[key])
        while len(bps)>0:
        
            bps_old = copy.deepcopy(bps[1:])
            l = [bps[0]]
            to_delete = []
            for count, bp in enumerate(bps_old):
                if  (bp[0]==l[-1][0]-1 and bp[1] ==l[-1][1] +1): 
                    l.append(bp)
                    to_delete.append(count)
            for ele in sorted(to_delete, reverse = True):        
                            del bps_old[ele]
            bps = copy.deepcopy(bps_old)         
            cache.append(l)
            
        #print(cache)
        return cache, 0

    
    def greedy_bp3_selection(self,):
        cache = [[]]
        cache_comp = [self.bps]
        i_cache =  [[]]
        e_cache = [np.inf]
    
        for key in sorted(list(self.dict_d.keys()), key = int):
            for i_bp, bp in enumerate(self.dict_d[key]):
            
                l_optimal = []
                i_l_optimal = []
                l_comp_optimal = [] 
                e_opt = []
                flag = 0
                remove_optimal=[]
                for count, l  in enumerate(cache):
                        
                        l_m = copy.deepcopy(l)
                        l_m.append(bp)
                        e_l = copy.deepcopy(e_cache[count])
                        
                        if  len(l) <1 or (bp[0]==l[-1][0]-1 and bp[1] ==l[-1][1] +1):
                            if bp in cache_comp[count]:

                                try:
                                    e_l_m = compute_energy(self.sequence, self.from_bp_to_DB( l_m , self.sequence,))
                                except Exception:
                                    e_l_m = np.inf
                                    
                                if e_l > e_l_m:
                                    flag = 1
            
                                    remove_optimal.append(count)
                                    l_optimal.append(l_m)  # Use the modified copy
                                    i_m  = copy.deepcopy(i_cache[count])
                                    i_m.append(i_bp)
                                    i_l_optimal.append(i_m)
                                    m_comp =  copy.deepcopy(cache_comp[count])
                                    l_comp_optimal.append(list(set(self.graph_bp['({}, {})'.format(*bp)]).intersection(set(m_comp))))
                                    e_opt.append(e_l_m)
                       
                if flag ==1 :

                    for ele in sorted(remove_optimal, reverse = True): 
                        if ele!=0 :
                            del cache[ele]
                            del i_cache[ele]
                            del cache_comp[ele]
                            del e_cache[ele]
                    cache.extend(l_optimal)
                    i_cache.extend(i_l_optimal) 
                    cache_comp.extend(l_comp_optimal)
                    e_cache.extend(e_opt)
                    
                
        #print(cache)
        return cache, e_cache
    
    
    def greedy_bp2_selection(self,):
        cache = [[]]
        cache_comp = [self.bps]
        i_cache =  [[]]
        e_cache = [np.inf]
    
        for key in sorted(list(self.dict_d.keys()), key = int):
            for i_bp, bp in enumerate(self.dict_d[key]):
            
                l_optimal = []
                i_l_optimal = []
                l_comp_optimal = [] 
                e_opt = []
                flag = 0
                
                for count, l  in enumerate(cache):
                        l_m = copy.deepcopy(l)
                        l_m.append(bp)
                        e_l = copy.deepcopy(e_cache[count])
                        
                        if l != [] and bp[0]!=l[-1][0]-1 and bp[1] !=l[-1][1] +1:
                                continue
                        elif l ==[]:
                            try:
                                e_l_m = compute_energy(self.sequence, self.from_bp_to_DB( l_m , self.sequence,))
                            except Exception:
                                e_l_m = np.inf
                            
                            flag = 1
                            l_optimal.append(l_m)  # Use the modified copy
                            i_m  = copy.deepcopy(i_cache[count])
                            i_m.append(i_bp)
                            i_l_optimal.append(i_m)
                            m_comp =  copy.deepcopy(cache_comp[count])
                            l_comp_optimal.append(list(set(self.graph_bp['({}, {})'.format(*bp)])))
                            e_opt.append(e_l_m)
                        elif bp in cache_comp[count]:

                            try:
                                e_l_m = compute_energy(self.sequence, self.from_bp_to_DB( l_m , self.sequence,))
                            except Exception:
                                e_l_m = np.inf
                                
                            if e_l > e_l_m:
                                flag = 1
                                l_optimal.append(l_m)  # Use the modified copy
                                i_m  = copy.deepcopy(i_cache[count])
                                i_m.append(i_bp)
                                i_l_optimal.append(i_m)
                                m_comp =  copy.deepcopy(cache_comp[count])
                                l_comp_optimal.append(list(set(self.graph_bp['({}, {})'.format(*bp)]).intersection(set(m_comp))))
                                e_opt.append(e_l_m)
                            
                if flag == 1:
                    cache.extend(l_optimal)
                    i_cache.extend(i_l_optimal) 
                    cache_comp.extend(l_comp_optimal)
                    e_cache.extend(e_opt)
        #print(cache)
        return cache, e_cache
    
    
    
    def from_motifs_to_DB(self,structures= None, lst= None, start = None, end = None):
        
            n_lst = len(lst)    
            db_string = '('*self.l_fix + '.'*(n_lst-(self.l_fix*2)) + ')'*self.l_fix
            string_list = list(db_string)
            for structure in structures:
                    for bond in structure:
                        string_list[bond[0]] = '('
                        string_list[bond[1]] = ')' 
            if start is not None and end is not None: 
                return ''.join(string_list[start:end])
            else: 
                return ''.join(string_list)
            
    def from_bp_to_DB(self,structures= None, lst= None, start = None, end = None):
        
            n_lst = len(lst)    
            db_string = '('*self.l_fix + '.'*(n_lst-(self.l_fix*2)) + ')'*self.l_fix
            string_list = list(db_string)
            for bond in structures:
                        string_list[bond[0]] = '('
                        string_list[bond[1]] = ')' 
            if start is not None and end is not None: 
                return ''.join(string_list[start:end])
            else: 
                return ''.join(string_list)
            
            
            
    def lowest_n_percent(self, numbers, n=5):
        # Calculate the 5th percentile value
        percentile_value = np.percentile(numbers, n)
        
        lowest_indices = [i for i, x in enumerate(numbers) if x < percentile_value]
    
        # Sort the indices based on the corresponding values
        lowest_indices_sorted = sorted(lowest_indices, key=lambda i: numbers[i])
            
        return lowest_indices_sorted 
    def find_maximal_cliques(self,G):
        
        cliques = list(find_cliques(G))
        max_size = max(len(clique) for clique in cliques)
        largest_cliques = [clique for clique in cliques if len(clique) == max_size]
        return largest_cliques
    
    def fit_fold(self,sequence=[] , n_tmpl=4, l_fix=4, return_min = True, filters = False):
        
            self.l_fix= l_fix # number of fixed base pairs in the 5' 3' bonds.
            self.n_tmpl = n_tmpl
            self.n_wrld = len(sequence)
            self.sequence = sequence # orginal sequnce in capital letters
            self.convert_sequence_to_numeric(self.sequence)
            self.build_BG_matrix()
            self.init_candidate_dict()
            if filters: 
                 self.stats_filter()
                 self.topology_filter()
            self.apt_filter()
            '''
            '''
            self.bps = set([])
            for mot in self.motifs:
                for bp in mot:
                    self.bps.add(bp)
            self.bps = list(self.bps)
            '''
            self.graph_bp = {str(bp): [] for bp in self.bps} # for each motif stores indices of associated compatible motifs
            for i, bp1 in enumerate(self.bps):
                for j, bp2 in enumerate(self.bps[i + 1:], i + 1):
                        if self.compatible_bp(bp1,bp2):
                            self.graph_bp['({}, {})'.format(*bp2)].append(bp1)
                            self.graph_bp['({}, {})'.format(*bp1)].append(bp2)  
            '''
            self.create_dict_d()
            self.create_dict_Sij()
            
            
            #cache, e_cache = self.greedy_bp4_selection()
            
            #self.create_dict_Sij()
            #self.create_dict_d()
            return #cache, e_cache# self.M,  self.dict_Sij
            '''
            self.graph_m = {}
            
            for  bp in self.bps:
                    try:
                        self.graph_m[str(bp[0])].append(bp)
            
                    except KeyError:
                        self.graph_m[str(bp[0])] = []
                
                        self.graph_m[str(bp[0])].append(bp)
            self.sorted_bps = []            
            for key in sorted(list(self.graph_m.keys()), key=int):
                self.graph_m[key] = sorted(self.graph_m[key] , key=lambda x: x[1])[::-1]
                self.sorted_bps.extend(self.graph_m[key])
            '''
            '''
            self.graph_bp = {str(i): [] for i in range(len(self.sorted_bps))} # for each motif stores indices of associated compatible motifs
            for i, bp1 in enumerate(self.sorted_bps):
                for j, bp2 in enumerate(self.sorted_bps[i + 1:], i + 1):
                        if self.compatible_bp(bp1,bp2):
                            self.graph_bp[str(i)].append(j)
                            self.graph_bp[str(j)].append(i)    
            create_dict_Sij      
                     
            self.graph = {str(i): [i] for i in range(len(self.motifs))} # for each motif stores indices of associated compatible motifs
            for i, motifs1 in enumerate(self.motifs):
                for j, motifs2 in enumerate(self.motifs[i + 1:], i + 1):
                    if self.compatible(motifs1, motifs2):
                        self.graph[str(i)].append(j)
                        self.graph[str(j)].append(i)
            a, b = self.greedy_bp2_selection()
            '''   
           
            '''
            self.graph_m = {}
            
            for  bp in list(self.bps):
                    try:
                        self.graph_m[str(bp[0])].append(bp)
            
                    except KeyError:
                        self.graph_m[str(bp[0])] = []
                
                        self.graph_m[str(bp[0])].append(bp)
            sorted_bps = []            
            for key in sorted(list(self.graph_m.keys())):
                self.graph_m[key] = sorted(self.graph_m[key] , key=lambda x: x[1])[::-1]
                sorted_bps += self.graph_m[key]
                
            #self.bps =  list(dict.fromkeys(self.bps))
        
     
            #self.maximal_sets =  self.find_maximal_cliques(G)
            return  self.graph_m
            '''
            
            G = nx.Graph()
            n = len(self.motifs)
            for i, motifs1 in enumerate(self.motifs):
                for j, motifs2 in enumerate(self.motifs[i + 1:], i + 1):
                    if self.compatible(motifs1, motifs2):
                        G.add_edge(i, j)
                        G.add_edge(j, i)
            self.maximal_sets =  self.find_maximal_cliques(G)
            '''

            if len(self.motifs) == 0:
                if return_min: # if not motif is found retorun only stem
                    return '('*self.l_fix + '.'*(self.n_wrld-(self.l_fix*2)) + ')'*self.l_fix
                else: # if no retrun is requested inform no motif was found
                    print('secondary structure consists only of given initial stem')
                    stem = []
                    for i in range(self.l_fix):
                        stem.append((i, self.n_wrld-1-i))
                    self.motifs.append(stem)
            
            self.bps =set([])
            for count, mot in enumerate(self.motifs):
                for m in mot:
                    self.bps.add(m)
            self.bps = list(self.bps)
            '''
            '''
            self.motifs.append([])  
            self.graph = {str(i): [i] for i in range(len(self.motifs))} # for each motif stores indices of associated compatible motifs
            for i, motifs1 in enumerate(self.motifs):
                for j, motifs2 in enumerate(self.motifs[i + 1:], i + 1):
                    if self.compatible(motifs1, motifs2):
                        self.graph[str(i)].append(j)
                        self.graph[str(j)].append(i)
            
            self.find_maximal_sets(self.graph)
            self.graph_m = {}
            for  mot in self.motifs:
                for  bp in mot:
                    try:
                        self.graph_m[str(bp[0])].add(bp)
            
                    except KeyError:
                        self.graph_m[str(bp[0])] = set([])
                
                        self.graph_m[str(bp[0])].add(bp)
            return self.graph_m
            '''
            return self.maximal_sets
            #self.bps = sorted( self.bps, key=lambda x: x[0])       
            '''
            self.motifs.append([])  
            self.graph = {str(i): [i] for i in range(len(self.motifs))} # for each motif stores indices of associated compatible motifs
            for i, motifs1 in enumerate(self.motifs):
                for j, motifs2 in enumerate(self.motifs[i + 1:], i + 1):
                    if self.compatible(motifs1, motifs2):
                        self.graph[str(i)].append(j)
                        self.graph[str(j)].append(i)
            
            
            
            self.graph_bp = {str(i): [i] for i in range(len(self.bps))} # for each motif stores indices of associated compatible motifs
            for i, bp1 in enumerate(self.bps):
                for j, bp2 in enumerate(self.bps[i + 1:], i + 1):
                    if bp1[0] != bp2[0]:
                        if self.compatible_bp(bp1,bp2):
                            self.graph_bp[str(i)].append(j)
                            self.graph_bp[str(j)].append(i)          
                        
            self.bp_ordered_clusters()
            self.min_motifs_structures, self.energies = self.greedy_bp2_selection()       
            #self.find_largest_compatible_sets(self.motifs)
            #self.structures_DB = self.from_motifs_to_DB(self.structures, self.sequence) 
            #print(np.argmin(self.energies))
            #self.min_energy = self.lowest_n_percent(self.energies) #self.energies[np.argmin(self.energies)]
   
            #self.min_energy_structure_DB = [self.from_motifs_to_DB(  self.min_motifs_structures[e], self.sequence)  for e in self.min_energy]
            
            
            #return self.bps
          
           
        
    def plot_structures(self, threshold = 'min_energy'):
        print('Number of computed structures is', len(self.structures) )
        if threshold == 'min_energy':
                   print('Free energy', np.min(self.energies))
                   plt.figure(figsize=(5,5))
                   bg = BulgeGraph.from_dotbracket(self.min_energy_structure_DB,self.sequence)
                   fvm.plot_rna(bg, text_kwargs={"fontweight":"black"}, lighten=0.7,backbone_kwargs={"linewidth":3})
                   plt.show()
        
        elif threshold.isdigit(): 
            self.structures_DB = self.from_motifs_to_DB(self.structures, self.sequence) 
            count_bonds  = []
            for st in self.structures_DB:
                count_bonds.append(Counter(st)['('])   
            max_bonds = np.max(count_bonds)
            for i, st in enumerate(self.structures_DB):
                if count_bonds[i] >=  max_bonds - threshold:
                   plt.figure(figsize=(5,5))
                   bg = BulgeGraph.from_dotbracket(st,self.sequence)
                   fvm.plot_rna(bg, text_kwargs={"fontweight":"black"}, lighten=0.7,backbone_kwargs={"linewidth":3})
                   plt.show()
        else: 
                self.structures_DB = self.from_motifs_to_DB(self.structures, self.sequence) 
                for  st in self.structures_DB:
                        plt.figure(figsize=(5,5))
                        bg = BulgeGraph.from_dotbracket(st,self.sequence)
                        fvm.plot_rna(bg, text_kwargs={"fontweight":"black"}, lighten=0.7,backbone_kwargs={"linewidth":3})
                        plt.show()
            '''
