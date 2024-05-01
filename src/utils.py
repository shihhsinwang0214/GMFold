import pandas as pd
import random
import re

def remove_spaces_and_convert_to_upper(input_string, remove_first = 5):
    #input_string = input_string.upper()
    input_string = input_string.replace(" ", "")
    input_string = input_string.replace("\u200b", "")
    return input_string[remove_first:]

def remove_spaces(input_string):
    #input_string = input_string.upper()
    input_string = input_string.replace(" ", "")
    input_string = input_string.replace("\u200b", "")
    return input_string



def substitute_random_chars(input_string, d = 4 , l= 8):
    
    '''
    d is the number of nucleotides we want to swap in  the string
    
    l i s the nuber of intial/final nucleatides that are fixed
    
    '''
    
    if len(input_string) <= 2:  # If string length is less than or equal to 2, cannot perform substitution
        return input_string

    # Number of elements to exclude from the start and end
    length = len(input_string)
    if length <= 2 * l:  # If string length is too small to perform substitution
        return input_string

    # Selecting characters between the first l and last l elements
    substr_start = l
    substr_end = length - l-1
    substr = input_string[substr_start:substr_end]

    # Randomly selecting d characters to substitute
    substitute_indices = random.sample(range(len(substr)), min(d, len(substr)))

    # Substituting selected characters with random characters from set (A, C, T, G)
    substitution_set = ['A', 'C', 'T', 'G']
    substr_list = list(substr)
    for index in substitute_indices:
        substr_list[index] = random.choice(substitution_set)

    # Reconstructing the string with substitutions
    output_string = input_string[:substr_start] + ''.join(substr_list) + input_string[substr_end:]
    return output_string


def change_name( input_sequence, seed ):
    
    return input_sequence + f'_seed_{seed}'

def create_batches(csv_file_path, times):
    
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_file_path)
    df['sequence']= df['sequence'].apply(remove_spaces_and_convert_to_upper)
    sd = df

    for i in range( times ):
        random.seed(i)
        cd = df.copy()
        cd['sequence']= cd['sequence'].apply(substitute_random_chars)
        cd['name'] = cd['name'].apply(change_name, args=(i,))
        cd['times'] = 0
        sd = pd.concat([sd, cd], ignore_index=True)
     #dropp duplicates   
    sd.drop_duplicates(subset=['sequence'])

    return  sd

def check_char(sequence):
    mapping = {'A': 0, 'a': 0, 'T': 3, 't': 3, 'C': 1, 'c':1, 'G': 2, 'g':2}
    flag = True
    for char in sequence:
        if char not in mapping:
             flag = False

    return flag

def check_if_string(sequence):
    if isinstance(sequence, str):
        return True
    else:
        return False

# Define the check_structure function
def check_structure(input_string):
    pattern = r'^.*?____(.*?)____.*?$'
    if re.match(pattern, input_string):
        return True
    else:
        return False
    
def check_empty_or_short(input_string, threshold = 20): 
    if len(input_string) >=threshold:
        return True
    else:
        return False
    
def find_substrings(input_string):
    return re.findall(r'____(.*?)____', input_string)[0]

def attach_stems(main_string):
    return 'GGGACGAC' + main_string + 'GTCGTCCC'