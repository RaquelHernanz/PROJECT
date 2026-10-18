#################################################### INFORMATION ########################################################
"""
#What is a restriction digest?
A restriction digest, consists in cutting a DNA sample with one or more enzymes and then running the resulting fragments 
on an agarose gel, which separates them by size.
"""

###################################################### PROJECT ##########################################################
""" PAB Project 2026-2027

Purpose of the program: 
This program helps the user PLANNING and INTERPRETING restriction digests. 

Programme description: 
A command-line program that simulates restriction enzyme digests in silico. 
It loads and validates DNA sequences from FASTA files, reports basic statistics (length, base counts, GC content),  
finds enzyme recognition sites (including overlapping ones) and lists the fragments produced by digesting a sequence 
with one or more enzymes, ordered as they would appear on an agarose gel (Western Blot).

# How does it work?
  1. loads DNA sequences from files; 
  2. supports a set of operations over individual sequences or over the whole collection of loaded sequences; 
  3. Find where enzymes cut and simulating a digest.


@author: Beatriz Fernansdes and Raquel Hernanz 

"""
############################################# IMPORTS ###################################################################

import os
#To avoid problems with files across different Operating Systems (Windows, MacOS, Linux), we use the os.path.join function to create
#each file path. 
#Each operating system joins folders and files differently to create a path, some with \ and others with / in the file path.

#What does the function do? It detects which operating system the program is running on and joins the words using the correct slash 
#for that system


from dataclasses import dataclass

#################################################  CONSTANTS  #############################################
#CONSTANTS 
DATA_FOLDER: str = "data"
ENZYMES_FILE: str = "enzymes.txt"
PROMPT: str = "Command to execute: "
FASTA_HEADER_PREFIX: str = ">"
NUCLEOTIDES: set[str] = {'A','C','G','T'}

#não entendi 
REPEATED_FASTA: str = ""

# List of valid commands in the system (Case-Sensitive)
VALID_COMMANDS: lits[str]= ["LOAD_SEQUENCES", "LIST_SEQUENCES", "SEQUENCE_INFO", "FIND_SITES", "DIGEST", "EXIT"]


######################################### DATACLASSES AND THEIR METHODS ################################################################

#A @dataclass is a simplified and automated way to create classes in Python focused exclusively on storing 
#data, allowing you to save information for a specific field.In this project, we have 3 fields:

"""1. A DNA sequence loaded from a FASTA record (It will be used in step 4.3. LIST ALL LOADED SEQUENCES)."""
    #Attributes:
      # - seq_id: Sequence identifier: the text after ">" up to the first space.
      # - bases: Full sequence, with all FASTA lines joined.

@dataclass
class Sequence: 
    seq_id: str
    bases: str



""" 2. This dataclass represents restriction enzymes as described in the enzymes file, which function 
    as molecular 'scissors' for the sequence undergoing digestion"""
    #Attributes:
       #- name: Enzyme name (e.g. "EcoRI"); case-sensitive.
       #- site: Recognition site, as a string over A, C, G, T (e.g. GAATTC);
       #- cut_offset: Number of bases of the site that come before the cut on the top strand (e.g. 1 for G^AATTC).
             # - Precondition: 0 < cut_offset < len(site).

@dataclass
class Enzyme:
    name: str
    site: str
    cut_offset: int 



"""3. A fragment produced by a digested Sequence."""
    #Attributes:
    #- start: First position of the fragment (1-based, inclusive).
    #- end: Last position of the fragment (1-based, inclusive).
    #Preconditions:
        # 1 <= start <= end (length of the original sequence).

@dataclass
class Fragment:
    start: int
    end: int
    


################################################# DOMAIN FUNCTIONS #################################################################

"""--------------------------------
     READ THE ENZYMES (start-up) 
----------------------------------- 
This function allows to read the enzyme restriction file, enzymes.txt (always valid) """
#load the file and returns a dictionary where the key is the enzyme name:  dict[str, Enzyme]
# - Arguments: 
       #- file: the name of the file
       #- data_folder: the folder where the file is contained

def read_Enzyme_File (file: str, data_folder: str)-> dict[str, Enzyme]:
    enzymes: dict[str,Enzyme] = {}
        
    #The file_path variable holds the path that we provide to the open function
    file_path = os.path.join(data_folder, file)                            
        
    with open(file_path) as enzyme_txt:
        for line in enzyme_txt:
            complete_line = line.strip()                          #Read each line 
            name, site, cutoff_site = complete_line.split()       #Split into the following variables: name, site and cutoff_site
            enzymes[name] = Enzyme(name,site,int(cutoff_site))

    return enzymes

"""--------------------------------------------------------
                      READ FASTA FILE 
     STEP(4.2) LOADING SEQUENCE (VERIFYING A FASTA FILE) 
-----------------------------------------------------------
The program going to read and add the sequences from a FASTA file to a list[] in the order they appear. """
# Return: list[Sequence]
# Arguments:
     #- file: the name of the file
     #- data_folder: the folder where the file is contained
     #- prefix: helps to identify the fasta file
 
def read_Fasta_File(file:str, data_folder: str, prefix: str)-> list[Sequence]:
  # AI-assited function: In this case, we ask for some minor corrections structure
    sequences: list[Sequence] = []
    file_path = os.path.join(data_folder, file)

    with open(file_path) as fasta:
        for line in fasta:
            complete_line = line.strip()

            #If there is no words in that line
            if (complete_line == ""):
                continue                                                       #In this case we use continue 
              
            #If we have a new sequence
            if (complete_line[0] == prefix):                                                  #If equal to >
                complete_line = complete_line[1:]                   #We can use the scling to remove the ">"
                words = complete_line.split()                    #We are separating the words from each line
                seq_id = words[0]                                               #The name/id of the sequence
                sequences.append(Sequence(seq_id, ""))                     #Add to the list the new sequence

            #If we have bases from a sequence
            else:
                #We are using the last sequence added to the list
                current_sequence = sequences[-1]
                # And we modify the bases from that sequence
                current_sequence.bases = current_sequence.bases + complete_line
    
    return sequences


"""
readString return a String from the input() function

    Arguments: - message: A string message used for the input function
"""
def readString(message: str)-> str:
    return input(message)
"""
Identifies if there is a sequence from the fasta file that contains invalid nucleotides

        Return: str (name of the sequence). If empty, there is no sequence that has invalid nuclotides

        Arguments:
        - sequences: the list of sequences obtained from a fasta file
        - nucleotides: the verified nucleotides that we can used for our verification
"""

def nucleotideVerification (sequences: list[Sequence], nucleotides: set[str])-> str:
    sequence_id: str = ""
    for sequence in sequences:
        for i in range(len(sequence.bases)):
            if (sequence.bases[i] not in nucleotides):
                sequence_id = sequence.seq_id
                break

    return sequence_id

"""
Identifies if there is a sequence from the fasta file that is already in our loaded sequences

        Return: str (name of the sequence). If empty, there is no sequence that is repeated

        Arguments:
        - sequences: the list of sequences obtained from a fasta file
        - loaded_sequences: we have to see if they already uploaded
"""

def sequenceRepeated (sequences: list[Sequence], loaded_sequences: dict[str,Sequence]):
    sequence_id: str = ""
    for sequence in sequences:
            if(sequence.seq_id in loaded_sequences): #We are comparing the keys, which is the name/id of the sequences
                sequence_id = sequence.seq_id
                break
    
    return sequence_id

"""
Add the valid sequences from a fasta file to the dictionary  

        Return: dict[str,Sequence]

        Arguments:
        - sequences: the list of sequences obtained from a fasta file
        - loaded_sequences: we have to upload the new sequences to the dictionary
"""

def addSequencesToDic (sequences: list[Sequence], loaded_sequences: dict[str,Sequence])-> dict[str,Sequence]:
    for sequence in sequences:
        loaded_sequences[sequence.seq_id] = sequence

    return loaded_sequences



###################################  STEP (4.3) LIST ALL LOADED SEQUENCES (inspect) #################################################

# LIST_SEQUENCES COMMAND 

"""1. Returns the ID of the first sequence with invalid nucleotides, or "" if all are valid."""
# Argument: sequences list with objects type Sequence (list[Sequence])

def findInvalidSequence(sequences: list[Sequence]) -> str:
    #AI-assisted function. 
    for sequence in sequences:
        for base in sequence.bases:
            if base not in NUCLEOTIDES:
                return sequence.seq_id
    return ""


# 1. list all currently loaded sequences, sorted alphabetically by ID, in ascending order using function "sorted()"
# 2. prints the number of total sequences loaded 
# 3. prints each sequence ordered by bases pairs lenght

# Arguments: 
#      - sequences: the sequences that will be ordered and printed later
#      - the dictionary key is the seq_id

def list_Sequences(sequences: dict[str, Sequence])-> None:
  # ensure it is created as a dict, not as a list
  sequences = dict(sorted(sequences.items()))
  
  print (f"Number of sequences loaded: {len(sequences)}")      #2. 
  
  for sequence in sequences.values():                          #3. 
    print(f" {sequence.seq_id} ({len(sequence.bases)} bp)")



# LOADED SEQUENCES (with 3 validation rules)
#This function will execute the LIST_SEQUENCES command, displaying the sequences in the terminal
# 1. Count number of sequences loaded
        #if count_number_of_sequence_loaded = 0 : 
                #print("No loaded sequences.")
        

###################################  STEP (4.4) SEQUENCE INFORMATION (inspect) #####################################################

"""----------------------------
          GC PERCENTAGE
-------------------------------
Returns the percentage of G and C in the sequence, rounded to 2 decimal places."""
#Precondition: the sequence has at least one base."""

def gc_content(sequence: Sequence) -> str: 
    count_gc = sequence.bases.count("G") + sequence.bases.count("C")
    porcentage= (count_gc / len(sequence.bases)) * 100
    
    # Retorna uma string formatada
    return f"{porcentagem:.2f}"


"""----------------------------
          FRAGMENT LENGTH
-------------------------------
Returns the length (in bp) of a DNA fragment of a fragment."""
#Preconditions: 1 <= start <= end (length of the original sequence).
# Argument: fragment; 
# Classe (type hint): Fragment
# Return: integer number

def fragment_length(fragment: Fragment) -> int:
  return fragment.end - fragment.start + 1 

# If a fragment begins at position 1 and ends at position 3, it contains 3 bases (bases 1, 2, and 3). By adding 1 (3 - 1 + 1 = 3), 
# we ensure that the initial base is also included in the total count.


###################################  STEP (4.5) FINDING RECOGNITIONM SITES (analyse) ###############################################

#FIND_SITES COMMAND 
"""--------------------------------------
          FINDING TE SITE POSITION 
-----------------------------------------
The function will return a list containing the position numbers where the enzyme recognition site was found."""
#IMPORTANT: In biology, positions are counted starting from 1 (the first base is position 1). Since lists in Python 
#           begin at index 0, so we added 1 (positions.append(i + 1)) to convert the Python index into the correct biological position.

#Arguments:
# - sequence: Sequence(dataclass) which contains the DNA strand where we will perform the operations
# - enzyme: Enzyme (dataclass) which contains information about the enzyme, specifically its recognition site (enzyme.site)

def find_site_positions(sequence: Sequence, enzyme: Enzyme) -> list[int]:
    positions: list[int] = []
    size = len(enzyme.site)
  
    for i in range(len(sequence.bases) - size + 1):
        if sequence.bases[i: i + size] == enzyme.site:
            positions.append(i + 1)                                           # i starts at 0 but the positions start at 1
    return positions


######################################  STEP (4.6) SIMULATING A DIGEST  ##################################################
# The prorgam going to FIND_SITES AND DIGEST FUNCTIONS
# The command "DIGEST" lists the The command DIGEST fragments in this same order, from the top to the bottom of the 
# gel
#DIGEST COMMAND (We must pay special attention to the custom sorting logic in case of a tie)

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
# Main function that acts as a command interpreter.
        #Features are organized in 3 groups:
          #1. Load - read FASTA files (LOAD_SEQUENCES)
          #2. Inspect - list the loaded file and show information about each one (LIST_SEQUENCES, SEQUENCE_INFO)
          #3. Analyse - use enzyme to recognition sites and simulate digests (FIND_SITES, DIGEST)


def main() -> None:        
        
    #STEP 1. Program reads the enzymes files. 
    # if ENZIME not in ENZIME_NAME:  -print (f"The enzyme {ENZIME_NAME} is not available.")
    enzymes: dict[str, Enzyme] = readEnzymeFile(ENZYMES_FILE, DATA_FOLDER)  #Checked
        
    loaded_fasta_files: set [str] = set()
    sequences: dict [str, Sequence] = {}
        
    #Loop control variable used to determine whether the loop should continue executing or stop
    command_name: str = ""
    command_parts: list[str]

    #STEP2. repeatedly asks the user for a command (using while), executes it, and prints the corresponding messages, until the user types EXIT
    # - WE SHOULD PROVIDE INFORMATION- command EXIT terminates the execution

                
    while command_name != "EXIT":

    #The commands entered by the user do not come alone, they bring additional arguments. 
    # - eg. If the user wants to upload a FASTA file:  LOAD_SEQUENCES demo.fasta ; SEQUENCE_INFO demo ; DIGEST pGFP_insert EcoRI HindIII
                            
    #1. Read user input and split it by spaces-  We must separate the whole command line into smaller strings by the split()
    # - readString reads all characters received on the serial port and groups them into a single string
    # - The command_parts variable becomes a list with two elements: ["LOAD_SEQUENCES", "demo.fasta"]
    
    #2. stores the first element of the list (the main command) in the command_name variable and check if the user has not entered anything

        command_parts = readString(PROMPT).split()

        if command_parts:
            command_name = command_parts[0]    
        else:
            command_name = "" 
            continue 
 
                                
        #3. Check if the user provides an unknown command: 
        if command_name not in VALID_COMMANDS: 
            print ("Unknown command" )
            continue 

        if command_name == "LOAD_SEQUENCES":

            #To avoid issuss with command_parts, put the value inside file variable
            file = command_parts[1]

            # If the file was already uploaded previously 
            if (file not in loaded_fasta_files):

                #Read all of the sequences in the file
                sequences_in_fasta: list[Sequence] = readFastaFile(file, DATA_FOLDER, FASTA_HEADER_PREFIX)

                #Verify if the nuclotides in the sequences are correct
                if(nucleotideVerification(sequences_in_fasta, NUCLEOTIDES) != ""):
                    invalid_sequence = nucleotideVerification(sequences_in_fasta, NUCLEOTIDES)
                    print(f"File {file} could not be loaded. Sequence {invalid_sequence} contains invalid nucleotides.")

                #Also verify if they are repeated sequences   
                elif (sequenceRepeated(sequences_in_fasta,sequences) != ""):
                    repeated_sequence = sequenceRepeated(sequences_in_fasta,sequences)
                    print (f"File {file} could not be loaded. Sequence ID {repeated_sequence} is already in use.")

                #If everything is correct, then we can upload the sequences to the dictionary
                else:
                    loaded_fasta_files.add(file)
                    sequences = addSequencesToDic(sequences_in_fasta,sequences)
                    print(f"File {file} validated and loaded!")

            else:
                print(f"File {file} already loaded!")
                         
        elif command_name == "LIST_SEQUENCES":
            # It the dictionarry of sequences is empty
            if (len(sequences) == 0):
                print("No sequences loaded.")
            else:
                list_Sequences(sequences)

        elif command_name == "SEQUENCE_INFO":
          # This command presents basic information about a loaded sequence: its length, the number of occurrences of each base, and
          # its GC content
            print("Function complete")
          
        elif command_name == "FIND_SITES":
            print("Function complete")

        elif command_name == "DIGEST":
            if len(command_parts[2:]) < 3: 
                print("Error: The DIGEST command requires a sequence ID and at least one enzyme.")        
                print("Function complete")
                #print(f"Digestion of the sequence {sequence.ID?} with {enzime_used} produced the following fragments: {DIGEST}") 
                #- WE NEED TO DEFINE sequence.ID
                #- We need to define the enzime_used
        else:
            #The command EXIT terminates the execution saying "Goodbye".
            print("Goodbye!")                                          

#The rest of the list (command_parts[1:]) contains the arguments that we will need to pass to the functions       

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
#!!!!!!!!!! order: imports, constants, data classes, domain functions, _io functions, the main function, and the if __name__ == "__main__": guard.

#Essential concepts

"""
  - NULCEOTIDES- carries one of four bases: A (adenine), C (cytosine), G (guanine) and T (thymine).


Rules for holding the two strands together by hydrogen bonds between base pairs: 
- COMPLEMENTARY - Each base only pairs with its complementary base: A-T, C-G
      Each such pair is a base pair (bp)- the length of DNA is measured in base pairs.

- ANTIPARALLELISM - Each strand has a direction. The two strands run in opposite directions: 5' → 3'; 3' → 5'
       A strand is always read and written in the 5' → 3' (direction in which cells read and copy DNA)


"""
#----------------------------------
