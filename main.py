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

# AT THE END: 
# Check the entire code with AI, ask for corrections, and see how we can improve its appearance
# The folder named data should be in the same folder as your main.py file.


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
PROMPT: str = "Command to execute: "
FASTA_HEADER_PREFIX: str = ">"
GOODBYE_MESSAGE: str = "Goodbye!"
UNKNOWN_COMMAND: str = "Unknown command" 
REPEATED_FASTA: str = ""

# HACIAS ASÍ? List of valid commands in the system (Case-Sensitive)
VALID_COMMANDS = ["LOAD_SEQUENCES", "LIST_SEQUENCES", "SEQUENCE_INFO", "FIND_SITES", "DIGEST", "EXIT"]


################################################# DATACLASSES AND THEIR METHODS   ######################################
from dataclasses import dataclass

"""
A DNA sequence loaded from a FASTA record.
    
    Attributes:
    - seq_id: Sequence identifier: the text after ">" up to the first space.
    - bases: Full sequence, with all FASTA lines joined.
"""

# It will be used in step 4.3. LIST ALL LOADED SEQUENCES
@dataclass
class Sequence: 
    # It represents a charged biological sequence.
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


########################################## FILE READING FUNCTIONS  #################################################################################

# To avoid problems with files across different Operating Systems, we use the os.path.join function to create each file path. 
# WHy? Because each operating system joins folders and files differently to create a path, some with \ and others with / in the file path.
# What does the function do? It detects which operating system the program is running on and joins the words using the correct slash for that system
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

#if ENZIME not in ENZIME_NAME:  -print (f"The enzyme {ENZIME_NAME} is not available.")


def readEnzymeFile (file: str, data_folder: str)-> dict[str, Enzyme]:
        enzymes: dict[str,Enzyme] = {}
        
        #The file_path variable holds the path that we provide to the open function
        file_path = os.path.join(DATA_FOLDER, file)                            
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
# command LOAD_SEQUENCES (used to load the sequences of a FASTA file)

# if sequence not in LOAD_SEQUENCES: -print(f"The sequence {SEQUENCE_ID} is not loaded.")


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

# loading the files demo.fasta and constructs.fasta


###################################  STEP 5. (4.3) LIST ALL LOADED SEQUENCES #################################################

# 1. In this part the command LIST_SEQUENCES needs to list all currently loaded sequences, sorted alphabetically by ID, 
#    in ascending order, so we going to use the function "sorted()"

def get_sorted_sequences(sequences: dict[str, Sequence] -> list[sequence]: 
     #Returns a list of sequences ordered alphabetically by ID.
        return sorted(sequences.values(), extrair_seq.id(Sequence): 

                      
# 2. prints the number of sequences loaded and lists all sequences, with their length
#This function will execute the LIST_SEQUENCES command, displaying the sequences in the terminal.
def list_sequences_io(sequences: dict[str, Sequence]) -> None:

        if not sequences


@dataclass
class Sequence: 
    # It represents a charged biological sequence.
    seq_id: str
    bases: str

###################################  STEP 6. (4.4) SEQUENCE INFORMATION  #####################################################



###################################  STEP 7. (4.5) FINDING RECOGNITIONM SITES  ###############################################



######################################  STEP 8. (4.6) SIMULATING A DIGEST  ##################################################
# The prorgam going to FIND_SITES AND DIGEST FUNCTIONS
# The command "DIGEST" lists the The command DIGEST fragments in this same order, from the top to the bottom of the 
# gel


# TO SIMULATED A DIGEST OF A SEQUENCE s in sequence: 
#       1. we find all occurrences of the recognition sites of the enzymes in s.

#       2. Each time it finds a pattern, it records the starting position of that pattern as p (starting the count at 1).

#       Enzymes do not necessarily cut at the beginning of the pattern; they have an offset (referred to as offset/ or o) 

#       The cut is made exactly between position (p + o - 1) and position p + o)

#NOTE: A linear sequence with n distinct cuts is thus split into n + 1 fragments
#     the lengths of the fragments always add up to the length of the sequence.
#In the end, if we have time, we can use the NEBcutter online tool to find out exactly what the actual digestion results would be.


#COMMENT: SHOULD WE ASK FOR THE OFFSET (o)?


##############################################   MAIN Function  ##########################################################
def main() -> None:        
        #STEP 1. Program reads the enzymes files. 
        enzymes: dict[str, Enzyme] = readEnzymeFile(ENZYMES_FILE, DATA_FOLDER) #Checked
        loaded_fasta_files: set [str] = set()
        sequences: dict [str, Sequence] = {}
        
        #Loop control variable used to determine whether the loop should continue executing or stop
        command_name: str = ""

        #STEP2. repeatedly asks the user for a command (using while), executes it, and prints the corresponding messages, until the user types EXIT
                 # - WE SHOULD PROVID INFORMATION:- command EXIT terminates the execution

                 #NOTE: Considering the commmands follow this structur: 
                 #          - SEQUENCE_INFO demo
                 #          - DIGEST pGFP_insert EcoRI HindIII
                 # We must separate the whole command line into smaller strings by the split()- Also take into account if the user accidently introduce an enter 
                 # without introducing a letter (in that case, unkwon message)
    
        while command_name != "EXIT":
                
                #1. Read user input and split it by spaces. 
                     # - readString reads all characters received on the serial port and groups them into a single string
                command_parts = readString(PROMPT).split()

                #2. WE NEED TO CHECK THIS PART! 
                #This part is useful to avoid uploading the same file again.
                if command_parts:
                        command_name = command_parts[0]
                
                #3. Check if the user provides an unknown command: 
                if command_name not in VALID_COMMANDS: 
                        print ("UNKNOWN_COMMAND")
                        continue 

                #4. Execution of valid commands
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
                        #The command EXIT terminates the execution saying "Goodbye".
                        print(GOODBYE_MESSAGE)                                          



#- We should do a io function for this part. Specially if we want to control the input of  the user
       
# 4.3. LIST_SEQUENCE 
        #1. Ask user for sequence 
        #2. Count number of sequences loaded
        #count_number_of_sequence_loaded 
        #if count_number_of_sequence_loaded = 0 : 
                #print("No loaded sequences.")

        #3. Print number of sequences loaded 
        #else: 
                #print(f"Number of sequences loaded: {count_number_of_sequence_loaded}")
        
        #4. Call the function to give us the list of uploaded sequences
    
    print(GOODBYE_MESSAGE)


if __name__ == "__main__":
    main()

#-----------------------------------------------------
# INFORMATIONS THAT CAN BE USEFULL
#Program reads the enzyme 
    # - CUT LOCATION: always palindromic (it is read backwards)
    # - WE give the CUT OFFSET (0< int > length of its site | number of bases of the site that come before the cut)
    #- WE give  the CLEAVAGE SEQUENCE  and the site (file enzimes.text)

# TOP STRAND (forward DNA sequence)
# Create BOTTOM SEQUENCE (reverse sequence): 
         # REVERSE COMPLEMENT SEQUENCE: e.g ATGAATTCGG is CCGAATTCAT
                  #1. Complementary
                  #2. read backwards

# Restriction Digest 
#    - Find the sites and cut (CHECK)
#    - Fragments
#    - Agarose Gel - the fragments are separated by size
#          - command DIGEST lists:  order, from the top (larger the fragment) to the bottom of the gel (smallest one).

#----------------------------------
