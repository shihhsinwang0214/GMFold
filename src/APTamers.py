import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi
from forgi.graph.bulge_graph import BulgeGraph

class Aptamer_Fold():
    def __init__(self, sequence='ATA', n_tmpl=4, l_fix=4, structures_DB = None):
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
                if self.n_wrld - self.l_fix -1 >= i > self.l_fix -1 and self.n_wrld - self.l_fix -1>= j > self.l_fix -1 and self.sqnc_num[i] + self.sqnc_num[j] == 3 and np.abs(i - j) >= 7:
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
        nhg_mtf = self.motif_ngh_dct()
        range_loop = self.C_dct[str(0)].copy()
        range_loop_2 = self.C_dct[str(self.n_tmpl // 2)].copy()
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
                            
                    if not flag:
                        for k in range(self.n_tmpl // 2):
                            try:
                                self.C_dct[str(k)].remove(i + k)
                            except ValueError:
                                pass
        return 

    
    def compatible(self, motif_1, motif_2):
            '''
            defines compatibility between two motifs with the same verteces and edges but different
            vertex labels
            '''
            if motif_1[0][0] < motif_2[0][0] and motif_1[-1][1] < motif_2[-1][1]:
                return False
            if motif_2[0][0] < motif_1[0][0] and motif_2[-1][1] < motif_1[-1][1]:
                return False
            tuple_list = motif_1 + motif_2
            compatible = True
            for i in range(len(tuple_list)):
                for j in range(i + 1, len(tuple_list)):
                    if 0 < len(set(tuple_list[i]).intersection(set(tuple_list[j]))) < 2:
                        compatible = False
                        break
            return compatible
            
            
    def find_maximal_sets(self, graph):
            '''
            given graph with i-th entry associated with list of indeces of motifs
            compatible to the i-th motif, find the list with maximal compatible motif,
            excludding overlapping components. Thus find the list with maximal compatible bonds.
            '''
            self.maximal_sets = [[0]] # start with list of motifs compatible to the 0-th motif
            graph_m = {str(0): graph[str(0)]}  # graph storing lists of motifs compatible with all the motifs in the i-th maximal set of motifs ( this case 0-th)
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
            return self.maximal_sets    
            
            
              
    def find_largest_compatible_sets(self,motifs ):
        """
        Finds the largest compatible sets of motifs.
        """

        graph = {str(i): [i] for i in range(len(motifs))} # for each motif stores indices of associated compatible motifs
        for i, motifs1 in enumerate(motifs):
            for j, motifs2 in enumerate(motifs[i + 1:], i + 1):
                if self.compatible(motifs1, motifs2):
                    graph[str(i)].append(j)
                    graph[str(j)].append(i)
        maximal_sets = self.find_maximal_sets(graph)
        self.structures = []
        for m in maximal_sets:
            strc = []  # list of all bonds provided by motifs in maximal set
            energy = 0  # compute energy associated with given structure by summing the bond energy of bonds in that structure
            for i in m:
                strc += motifs[i]
            unique_strc = list(set(strc))  # remove duplicates
            self.structures.append(unique_strc)  
            #print(len(self.structures))   
            for j in unique_strc:  # TO DO: adapt to additional type of bonds GA, TG, etc...
                if j[0] in ['A', 'T']:
                    energy -= 1.5
                else:
                    energy -= 3.6
            self.energies.append(energy)
            #print(len(self.energies))   
        return
            
    def from_motifs_to_DB(self,structures= None, lst= None):
            n_lst = len(lst)
            list_strings = []
            for structure in structures:
                    db_string = '('*self.l_fix + '.'*(n_lst-(self.l_fix*2)) + ')'*self.l_fix
                    string_list = list(db_string)
                    for bond in structure:
                        string_list[bond[0]] = '('
                        string_list[bond[1]] = ')'
                        #new_string = db_string[:bond[0]] + '(' + db_string[bond[0]+1:bond[1]]+ ')'  +db_string[bond[1]+1:]
                        #db_string = new_string
                    db_string = ''.join(string_list)
                    list_strings.append( db_string)    
            return list_strings
            
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
            
            if len(self.motifs) == 0:
                if return_min: # if not motif is found retorun only stem
                    return '('*self.l_fix + '.'*(self.n_wrld-(self.l_fix*2)) + ')'*self.l_fix
                else: # if no retrun is requested inform no motif was found
                    print('secondary structure consists only of given initial stem')
                    stem = []
                    for i in range(self.l_fix):
                        stem.append((i, self.n_wrld-1-i))
                    self.motifs.append(stem)
                        
            self.find_largest_compatible_sets(self.motifs)
            #self.structures_DB = self.from_motifs_to_DB(self.structures, self.sequence) 
            self.min_energy_structure_DB =  self.from_motifs_to_DB([self.structures[np.argmin(self.energies)]], self.sequence)[0]   
            if return_min:
                return self.min_energy_structure_DB
            else:
                return
        
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

