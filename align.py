#!/usr/bin/python3

"""
DESCRIPTION:
    Template code for the Dynamic Programming assignment in the Algorithms in Sequence Analysis course at the VU.
    
INSTRUCTIONS:
    Complete the code (compatible with Python 3!) upload to CodeGrade via corresponding Canvas assignment.

AUTHOR:
    Name: Ville Lehtonen
    VUnet id: vln490
    Student number: 2658063

"""



import argparse
import pickle



def parse_args():
    "Parses inputs from commandline and returns them as a Namespace object."

    parser = argparse.ArgumentParser(prog = 'python3 align.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Aligns the first two sequences in a specified FASTA\n'
        '  file with a chosen strategy and parameters.\n'
        '\n'
        'defaults:\n'
        '  strategy = global\n'
        '  substitution matrix = pam250\n'
        '  gap penalty = 2')
        
    parser.add_argument('fasta', help='path to a FASTA formatted input file')
    parser.add_argument('output', nargs='*', 
        help='path to an output file where the alignment is saved\n'
             '  (if a second output file is given,\n'
             '   save the score matrix in there)')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
        help='print the score matrix and alignment on screen', default=False)
    parser.add_argument('-s', '--strategy', dest='strategy',
        choices=['global','semiglobal','local'], default="global")
    parser.add_argument('-m', '--matrix', dest='substitution_matrix',
        choices=['pam250','blosum62','identity'], default='pam250')
    parser.add_argument('-g', '--gap_penalty', dest='gap_penalty', type=int,
        help='must be a positive integer', default=2)

    args = parser.parse_args()

    args.align_out = args.output[0] if args.output else False
    args.matrix_out = args.output[1] if len(args.output) >= 2 else False
                      # Fancy inline if-else statements. Use cautiously!
                      
    if args.gap_penalty <= 0:
        parser.error('gap penalty must be a positive integer')

    return args



def load_substitution_matrix(name):
    "Loads and returns the specified substitution matrix from a pickle (.pkl) file."
    # Substitution matrices have been prepared as nested dictionaries:
    # the score of substituting A for Z can be found with subst['A']['Z']
    # NOTE: Only works if working directory contains the correct folder and file!
    
    with open('substitution_matrices/%s.pkl' % name, 'rb') as f:
        subst = pickle.load(f)
    return subst
    

def load_sequences(filepath):
    "Reads a FASTA file and returns the first two sequences it contains."
    
    seq1 = []
    seq2 = []
    with open(filepath,'r') as f:
        for line in f:
            if line.startswith('>'):
                if not seq1:
                    current_seq = seq1
                elif not seq2:
                    current_seq = seq2
                else:
                    break # Stop if a 3rd sequence is encountered
            else:
                current_seq.append(line.strip())
    
    if not seq2:
        raise Exception('Error: Not enough sequences in specified FASTA file.')
    
    seq1 = ''.join(seq1)
    seq2 = ''.join(seq2)
    return seq1, seq2


def find_starting_index_for_global_alignment(score_matrix):
    i = len(score_matrix) - 1
    j = len(score_matrix[0]) - 1
    return [i, j]


def find_starting_index_for_local_alignment(score_matrix):
    """Returns the index [i,j] of the cell with the highest value in a score matrix"""
    result = [0,0]
    for i in range(len(score_matrix)):
        for j in range(len(score_matrix[i])):
            if score_matrix[i][j] > score_matrix[result[0]][result[1]]:
                result = [i, j]
    return result


def find_starting_index_for_semiglobal_alignment(score_matrix):
    """Returns the index [i, j] of the starting cell for tracebacking semiglobal alignment"""
    ## First, find the maximum value (and its index) from the bottom row of the matrix
    bottom_row_max_value = max(score_matrix[-1])
    index_of_bottom_row_max_value = [len(score_matrix) - 1, score_matrix[-1].index(bottom_row_max_value)]

    ## Then, find the maximum value (and its index) from the rightmost column of the matrix
    max_rightmost_column_values = list(map(max, score_matrix))
    max_rightmost_column_value = max(max_rightmost_column_values)
    index_of_max_rightmost_column_value = [max_rightmost_column_values.index(max_rightmost_column_value), len(score_matrix[0]) - 1]

    ## Return the index of the max value in a cell either from the bottom row or from the most rightmost column
    if (max_rightmost_column_value >= bottom_row_max_value):
        return index_of_max_rightmost_column_value
    else:
        return index_of_bottom_row_max_value


def find_starting_point(strategy, score_matrix):
    if strategy == "global":
        return find_starting_index_for_global_alignment(score_matrix)

    elif strategy == "semiglobal":
        return find_starting_index_for_semiglobal_alignment(score_matrix)

    elif strategy == "local":
        return find_starting_index_for_local_alignment(score_matrix)


def dp_function(seq1, seq2, score_matrix, i, j, strategy, gap_penalty, substitution_matrix):
    """Calculates the correct value for a cell in score matrix"""
    #Score for match / mismatch
    score = substitution_matrix[seq1[i-1]][seq2[j-1]]

    #Scores for moving vertically, horizontally, or diagonally
    vertical = score_matrix[i-1][j] - gap_penalty
    horizontal = score_matrix[i][j-1] - gap_penalty
    diagonal = score_matrix[i-1][j-1] + score

    #for global & semiglobal, return the max value gained by moving into any of the three directions
    #for local alignment, add a zero to the formula
    if strategy == 'global' or strategy == 'semiglobal':
        return max(vertical, horizontal, diagonal)
    elif strategy == 'local':
        return max(vertical, horizontal, diagonal, 0)


def tb_function(i, j, score_matrix, gap_penalty, strategy):
    """Returns the direction (up, left, up-left) from which the cell (i,j)'s score was derived"""

    # the score in the cell (i,j) of the score matrix
    cell_score = score_matrix[i][j]

    #check if the vertical cell (the cell above) lead to the cell (i,j)
    cell_above = score_matrix[i][j-1]
    if cell_above - gap_penalty == cell_score:
        return "up"
    
    # check the negation of the horizontal cell (-> diagonal / horizontal)
    cell_left = score_matrix[i-1][j]
    if cell_left - gap_penalty != cell_score:
        return "up-left"
    else:
        return "left"


def fix_sequence_beyond_the_starting_point(seq1, seq2, seq1_final, seq2_final, i, j):
    """ Fixes the sequence when the starting point doesn't cover all the letters in the sequences """
    counter_i = len(seq1)
    counter_j = len(seq2)

    # (If needed) Add missing letters to seq1 and gaps to seq2    
    while counter_i > i:
        seq1_final += seq1[counter_i-1]
        seq2_final += "-"
        counter_i -= 1

    # (If needed) Add missing letters to seq2 and gaps to seq1
    counter_j = len(seq2)
    while counter_j > j:
        seq2_final += seq2[counter_j-1]
        seq1_final += "-"
        counter_j -= 1

    return seq1_final, seq2_final


def align(seq1, seq2, strategy, substitution_matrix, gap_penalty):
    "Do pairwise alignment using the specified strategy and parameters."
    # This function consists of 3 parts:
    #
    #   1) Initialize a score matrix as a "list of lists" of the appropriate length.
    #      Fill in the correct values for the first row and column given the strategy.
    #        (local / semiglobal = 0  --  global = stacking gap penalties)
    #   2) Fill in the rest of the score matrix using Dynamic Programming, accounting
    #      for the selected alignment strategy, substitution matrix and gap penalty.
    #   3) Perform the correct traceback routine on your filled in score matrix.
    #
    # Both the resulting alignment (sequences with gaps and the corresponding score)
    # and the filled in score matrix are returned as outputs.
    #
    # NOTE: You are strongly encouraged to think about how you can reuse (parts of)
    #       your code between steps 2 and 3 for the different strategies!
    
    
    ### 1: Initialize
    M = len(seq1)+1
    N = len(seq2)+1
    score_matrix = []
    for i in range(M):
        row = []
        score_matrix.append(row)
        for j in range(N):
            row.append(0)

    if strategy == 'global':
        #####################
        # START CODING HERE #
        #####################

        # Change the zeros in the first row to the correct value
        for j in range(len(score_matrix[0])):
            if j > 0:
                score_matrix[0][j] = score_matrix[0][j-1] - gap_penalty

        # Change the zeros in the first column to the correct value
        for i in range(len(score_matrix)):
            if i > 0:
                score_matrix[i][0] = score_matrix[i-1][0] - gap_penalty        

        #####################
        #  END CODING HERE  #
        #####################
    
    ### 2: Fill in Score Matrix
 
    #####################
    # START CODING HERE #
    #####################

    for i in range(1, M):
        for j in range(1,N):
            score_matrix[i][j] = dp_function(seq1, seq2, score_matrix, i, j, strategy, gap_penalty, substitution_matrix)
    
    #####################
    #  END CODING HERE  #
    #####################   
    
    ### 3: Traceback
    #####################
    # START CODING HERE #
    #####################
    
    ## Initialize the final sequence strings
    seq1_final = ""
    seq2_final = ""

    ## Find the coordinates (i, j) of the cell where the traceback starts
    start_cell = find_starting_point(strategy, score_matrix)
    i = start_cell[0]
    j = start_cell[1]

    ## For semiglobal alignments, fix the sequences with cells beyond the starting point
    if strategy == "semiglobal":
        seq1_final, seq2_final = fix_sequence_beyond_the_starting_point(seq1, seq2, seq1_final, seq2_final, i, j)

    #Set alignment score
    align_score = score_matrix[i][j]        # with gaps inserted at the appropriate positions.

    while True:
        if(i > 0 and j > 0):
            direction = tb_function(i, j, score_matrix, gap_penalty, strategy)
            if direction == "up": #vertical
                seq1_final += "-"
                seq2_final += seq2[j-1]
                j -= 1

            elif direction == "up-left": #diagonal
                seq1_final += seq1[i-1]
                seq2_final += seq2[j-1]
                j -= 1
                i -= 1

            elif direction == "left": #horizontal
                seq1_final += seq1[i-1]
                seq2_final += "-"
                i -= 1
        else:
            break

    ## Reverse the final alignments to get them in the correct form
    aligned_seq1 = seq1_final[::-1]
    aligned_seq2 = seq2_final[::-1]

    #####################
    #  END CODING HERE  #
    #####################   
    alignment = (aligned_seq1, aligned_seq2, align_score)
    return (alignment, score_matrix)


def print_score_matrix(s1,s2,mat):
    "Pretty print function for a score matrix."
    
    # Prepend filler characters to seq1 and seq2
    s1 = '-' + s1
    s2 = ' -' + s2
    
    # Print them around the score matrix, in columns of 5 characters
    print(''.join(['%5s' % aa for aa in s2])) # Convert s2 to a list of length 5 strings, then join it back into a string
    for i,row in enumerate(mat):               # Iterate through the rows of your score matrix (and keep count with 'i').
        vals = ['%5i' % val for val in row]    # Convert this row's scores to a list of strings.
        vals.insert(0,'%5s' % s1[i])           # Add this row's character from s2 to the front of the list
        print(''.join(vals))                   # Join the list elements into a single string, and print the line.



def print_alignment(a):
    "Pretty print function for an alignment (and alignment score)."
    
    # Unpack the alignment tuple
    seq1 = a[0]
    seq2 = a[1]
    score = a[2]
    
    # Check which positions are identical
    match = ''
    for i in range(len(seq1)): # Remember: Aligned sequences have the same length!
        match += '|' if seq1[i] == seq2[i] else ' ' # Fancy inline if-else statement. Use cautiously!
            
    # Concatenate lines into a list, and join them together with newline characters.
    print('\n'.join([seq1,match,seq2,'','Score = %i' % score]))



def save_alignment(a,f):
    "Saves two aligned sequences and their alignment score to a file."
    with open(f,'w') as out:
        out.write(a[0] + '\n') # Aligned sequence 1
        out.write(a[1] + '\n') # Aligned sequence 2
        out.write('Score: %i' % a[2]) # Alignment score


    
def save_score_matrix(m,f):
    "Saves a score matrix to a file in tab-separated format."
    with open(f,'w') as out:
        for row in m:
            vals = [str(val) for val in row]
            out.write('\t'.join(vals)+'\n')
    


def main(args = False):
    # Process arguments and load required data
    if not args: args = parse_args()
    
    sub_mat = load_substitution_matrix(args.substitution_matrix)
    seq1, seq2 = load_sequences(args.fasta)

    # Perform specified alignment
    strat = args.strategy
    gp = args.gap_penalty
    alignment, score_matrix = align(seq1, seq2, strat, sub_mat, gp)

    # If running in "verbose" mode, print additional output
    if args.verbose:
        print_score_matrix(seq1,seq2,score_matrix)
        print('') # Insert a blank line in between
        print_alignment(alignment)
    
    # Save results
    if args.align_out: save_alignment(alignment, args.align_out)
    if args.matrix_out: save_score_matrix(score_matrix, args.matrix_out)



if __name__ == '__main__':
    main()