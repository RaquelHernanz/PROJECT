###################################################### TO DO ##################################################################################
#document data classes and function
        # explaining what they mean 
        # what they dO 
        # defining preconditions for their usage
        # use type hints in every function

# implement the different features in the same order as presented in Section 4

#RULES 
"""
     - the user will always be correct and consistent with the command definition
     - every FASTA file provided by the user exists and is correctly formatted
     - All commands can be invoked at any time, in particular, the user should be able to load another file at any time. 

"""

# every function written with the help of an AI assistant (e.g., ChatGPT, Copilot, Claude) must explicitly include the statement " AI-assisted function ""

# !! Positions in a sequence are 1-based, as is standard in biology: the first base of a sequence is at position 1 (careful: in Python, it is at index 0!).


# Check the entire code with AI, ask for corrections, and see how we can improve its appearance


############################################## MAIN GOAL ###########################################################################################

# How does it work?
""" 
  1. loads DNA sequences from files; 
  2. supports a set of operations over individual sequences or over the whole collection of loaded sequences; 
  3. Find where enzymes cut and simulating a digest.
""" 

#Essential concepts

"""
  - NULCEOTIDES- carries one of four bases: A (adenine), C (cytosine), G (guanine) and T (thymine).


Rules for holding the two strands together by hydrogen bonds between base pairs: 
- COMPLEMENTARY - Each base only pairs with its complementary base: A-T, C-G
      Each such pair is a base pair (bp)- the length of DNA is measured in base pairs.

- ANTIPARALLELISM - Each strand has a direction. The two strands run in opposite directions: 5' → 3'; 3' → 5'
       A strand is always read and written in the 5' → 3' (direction in which cells read and copy DNA)


"""


######################################################### INFORMATION #############################################################
"""
#What is a restriction digest?
A restriction digest, consists in cutting a DNA sample with one or more enzymes and then running the resulting fragments on an agarose gel, 
which separates them by size.
"""

# Program that simulates RESTRICTION DIGEST experiments in silico. 
"""This program helps the user PLANNING and INTERPRETING restriction digests. """


################################################################ PROJECTS ###########################################################
""" PAB Project 2026-2027

Programme description: 
A command-line program that simulates restriction enzyme digests in silico. 
It loads and validates DNA sequences from FASTA files, reports basic statistics (length, base counts, GC content), 
finds enzyme recognition sites (including overlapping ones) and lists the fragments produced by digesting a sequence 
with one or more enzymes, ordered as they would appear on an agarose gel (Western Blot).

@author: Beatriz Fernansdes and Raquel Hernanz

"""

################################################# CONSTANTS  #############################################################
DATA_FOLDER: str = "data"
ENZYMES_FILE: str = "enzymes.txt"
PROMPT: str = "Command to execute:"
FASTA_HEADER_PREFIX: str = ">"
GOODBYE_MESSAGE: str = "Goodbye!"
UNKNOWN_COMMAND: str = "Unknown command."
REPEATED_FASTA: str = ""


################################################# DATACLASSES AND THEIR METHODS   ######################################
from dataclasses import dataclass

"""
A DNA sequence loaded from a FASTA record.
    
    Attributes:
    - seq_id: Sequence identifier: the text after ">" up to the first space.
    - bases: Full sequence, with all FASTA lines joined.
"""
@dataclass
class Sequence:
    seq_id: str
    bases: str


"""
A restriction enzyme, as described in the enzymes file.

    Attributes:
    - name: Enzyme name (e.g. "EcoRI"); case-sensitive.
    - site: Recognition site, as a string over A, C, G, T.
    - cut_offset: Number of bases of the site that come before the cut
            on the top strand (e.g. 1 for G^AATTC).
"""
@dataclass
class Enzyme:
    name: str
    site: str
    cut_offset: int 


"""
A fragment produced by a digested Sequence.

    Attributes:
    - start: First position of the fragment (1-based, inclusive).
    - end: Last position of the fragment (1-based, inclusive).

    Preconditions:
        1 <= start <= end <= length of the original sequence.
"""
@dataclass
class Fragment:
    start: int
    end: int


##########################################  INPUT/OUTPUT and FILE READING FUNCTIONS  ################################
"""
1. PROVID INFORMATION:- command EXIT terminates the execution
2. COMMAND NOT FOUND? 
     - print ("Command to execute")
"""

import os

"""
readString return a String from the input() function
    Arguments:
    - message: A string message used for the input function
"""
def readString(message: str)-> str:
    return input(message)


##################################### STEP 3. check the sequence ID ##################################################

# In thi step the program going to read the enzymes file, enzymes.txt (always valid) in the order they were given, 
#   reporting only the first unknown one

#Attributes: 
#    - file_path: the enzyme.txt file path



def readEnzymeFile (file: str, data_folder: str)-> dict[str, Enzyme]:
        enzymes: dict[str,Enzyme] = {}
        file_path = os.path.join(data_folder, file)
        print(file_path)
        
        with open(file_path) as enzyme_txt:
                for line in enzyme_txt:
                        complete_line = line.strip() #Read each line 
                        name, site, cutoff_site = complete_line.split() #Split into the following variables: name, site and cutoff_site
                        enzymes[name] = Enzyme(name,site,int(cutoff_site))
        return enzymes


###################################  STEP 4. (4.2) LOADING AND VERIFYING A FASTA FILE ###############################################

# The pram going to read and add a FASTA file to a dict. 
# 1.  Create a dict of the available enzymes inside the file. 
#     Each line follows the next pattern: name + site + cutoff_site
#     BamHI GGATCC 1



def readFastaFile(file:str, data_folder: str, prefix: str)-> list[Sequence]:
    sequences: list[Sequence] = []
    file_path = os.path.join(data_folder, file)
    name_seq: str = ""
    nucl_seq: str = ""
    with open(file_path) as fasta:
        for line in fasta:
            complete_line = line.strip()
            if (line == ""):
                continue
            else:
                if (complete_line[0] == prefix): #If equal to >
                    complete_line = complete_line.replace(prefix,"")# > must be remove
                    complete_line = complete_line.split()
                    current_seq_name = complete_line[0] #Obtain the name sequence
                    
                    #If it wasnt then we store the old one in the list and add the new opne
                    if (current_seq_name not in sequences):
                        sequences.append(Sequence(previous_seq_name,nucl_seq))
                        previous_seq_name = current_seq_name
                        nucl_seq = ""
                    
                    else:
                        nucl_seq = nucl_seq + complete_line


###################################  STEP 5. (4.3) LIST ALL LOADED SEQUENCES #################################################

# command LOAD_SEQUENCES (used to load the sequences of a FASTA file)




###################################  STEP 6. (4.4) SEQUENCE INFORMATION  #####################################################



###################################  STEP 7. (4.5) FINDING RECOGNITIONM SITES  ###############################################



######################################  STEP 8. (4.6) SIMULATING A DIGEST  ##################################################
# The prorgam going to FIND_SITES AND DIGEST FUNCTIONS




##############################################   MAIN Function  ##########################################################
def main() -> None:

    enzymes: dict[str, Enzyme] = readEnzymeFile(ENZYMES_FILE, DATA_FOLDER) #Checked
    loaded_fasta_files: set [str] = set()
    sequences: dict [str, Sequence] = {}
    #This is useful to avoid uploading the same file again.

    """
    Considering the commmands follow this structur:
    - SEQUENCE_INFO demo
    - DIGEST pGFP_insert EcoRI HindIII

    We must separate the whole command line into 
    smaller strings by the split()

    Also take into account if the user accidently introduce an 
    enter without introducing a letter (in that case, unkwon message)
    """

    # We should do a io function for this part
    # Specially if we want to control the input of 
    # the user
    command = readString(PROMPT).split()
    command_name: str = ""
    if (command):
        command_name = command[0]

   #The command EXIT terminates the execution saying "Goodbye".
    while command_name != "EXIT":

        if command_name == "LOAD_SEQUENCES":
            print("Function complete")
        elif command_name == "LIST_SEQUENCES":
            print("Function complete")
        elif command_name == "SEQUENCE_INFO":
            print("Function complete")
        elif command_name == "FIND_SITES":
            print("Function complete")
        elif command_name == "DIGEST":
            print("Function complete")

        else:
            print("Goodbye!")                                          #4.1. Terminating the Program 

        command = readString(PROMPT).split()
        command_name = ""
        if (command):
            command_name = command[0]

    print(GOODBYE_MESSAGE)


if __name__ == "__main__":
    main()
