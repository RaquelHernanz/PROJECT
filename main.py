"""################################ INFORMATION ###########################
#What is a restriction digest?
#A restriction digest, consists in cutting a DNA sample with one or more enzymes and then running the resulting fragments 
#on an agarose gel, which separates them by size.


################################## PROJECT ###################################
PAB Project 2026-2027 | @author: Beatriz Fernansdes and Raquel Hernanz

Purpose of the program: 
#This program helps the user PLANNING and INTERPRETING restriction digests. 

Programme description:
#A command-line program that simulates restriction enzyme digests in silico. 
#It loads and validates DNA sequences from FASTA files, reports basic statistics (length, base counts, GC content),  
#finds enzyme recognition sites (including overlapping ones) and lists the fragments produced by digesting a sequence 
#with one or more enzymes, ordered as they would appear on an agarose gel (Western Blot).

How does it work?
  # 1. loads DNA sequences from files; 
  # 2. supports a set of operations over individual sequences or over the whole collection of loaded sequences; 
  # 3. Find where enzymes cut and simulating a digest.

"""


####################################################### IMPORTS #####################################################

import os
#To avoid problems with files across different Operating Systems (Windows, MacOS, Linux), we use the os.path.join function to create
#each file path. Each operating system joins folders and files differently to create a path, some with \ and others with / in the file path.
#What does the function do? It detects which operating system the program is running on and joins the words using the correct slash 
#for that system


from dataclasses import dataclass


###################################################  CONSTANTS  ##################################################

#CONSTANTS- fixed values ​​that should not be changed during program execution

DATA_FOLDER: str = "data"
ENZYMES_FILE: str = "enzymes.txt"
PROMPT: str = "Command to execute: "
FASTA_HEADER_PREFIX: str = ">"
NUCLEOTIDES: set[str] = {'A','C','G','T'}


# List of valid commands in the system (Case-Sensitive)
VALID_COMMANDS: list[str]= ["LOAD_SEQUENCES", "LIST_SEQUENCES", "SEQUENCE_INFO", "FIND_SITES", "DIGEST", "EXIT"]



############################################### DATACLASSES AND THEIR METHODS ############################################

#A @dataclass is a simplified and automated way to create classes in Python focused exclusively on storing 
#data, allowing you to save information for a specific field.In this project, we have 3 fields:

"""
-------------------------------------------------SEQUENCE DATA CLASSE--------------------------------------------------------
    A DNA sequence loaded from a FASTA record (It will be used in step 4.3. LIST ALL LOADED SEQUENCES).
        Attributes:
        - seq_id: Sequence identifier: the text after ">" up to the first space.
        - bases: Full sequence, with all FASTA lines joined.
"""

@dataclass
class Sequence:
    seq_id: str
    bases: str


"""
-------------------------------------------------ENZYME DATA CLASSE--------------------------------------------------------
    This dataclass represents restriction enzymes as described in the enzymes file, which function as molecular 'scissors' for the 
    sequence undergoing digestion
        Attributes:
        - name: Enzyme name (e.g. "EcoRI"); case-sensitive.
        - site: Recognition site, as a string over A, C, G, T (e.g. GAATTC);
        - cut_offset: Number of bases of the site that come before the cut on the top strand (e.g. 1 for G^AATTC).
        Precondition: 0 < cut_offset < len(site).
"""

@dataclass
class Enzyme:
    name: str
    site: str
    cut_offset: int 


"""-------------------------------------------------FRAGMENT DATA CLASSE--------------------------------------------------------
    A fragment produced by a digested Sequence.
        Attributes:
        - start: First position of the fragment (1-based, inclusive).
        - end: Last position of the fragment (1-based, inclusive).
        Preconditions: 1 <= start <= end (length of the original sequence).

"""

@dataclass
class Fragment:
    start: int
    end: int


################################################# DOMAIN FUNCTIONS #################################################################

# DOMAIN FUCTIONS- Functions that serve strictly to compute logical or mathematical results. 
# They are responsible for the program's core logic (such as simulating enzyme cuts or calculating DNA sequence statistics).


################ READ THE ENZYMES (start-up) ##############

def read_enzyme_file (file: str, data_folder: str)-> dict[str, Enzyme]:
    """
    This function allows to read the enzyme restriction file (enzymes.txt), load the file and returns a dictionary where the 
    key is the enzyme name:  dict[str, Enzyme]
        return: dict[str, Enzyme]
        Arguments: 
        - file: the name of the file
        - data_folder: the folder where the file is contained
    """
    enzymes: dict[str,Enzyme] = {}
        
    #The file_path variable holds the path that we provide to the open function
    file_path = os.path.join(data_folder, file)                            
        
    with open(file_path) as enzyme_txt:
        for line in enzyme_txt:
            complete_line = line.strip()                          
            name, site, cutoff_site = complete_line.split() #Split into the following enzyme attributes: name, site and cutoff_site
            enzymes[name] = Enzyme(name,site,int(cutoff_site))

    return enzymes
  

############## LOADING SEQUENCE STEP (4.2), INCLUDING READING FASTA FILE, VERIFICATION FUNCTIONS AND ADDITION OF SEQUENCES ##############

def read_fasta_file(file:str, data_folder: str, prefix: str)-> list[Sequence]:
    """
    The program going to read and add the sequences from a FASTA file to a list[] in the order they appear.
    AI-assited function: In this case, we ask for some minor corrections in the code structure
    
        return: list[Sequence] 
        Arguments:
        - file: the name of the file
        - data_folder: the folder where the file is contained
        - prefix: helps to identify the fasta file
    """
    sequences: list[Sequence] = []
    file_path = os.path.join(data_folder, file)
    
    with open(file_path) as fasta:
        for line in fasta:
            complete_line = line.strip()
            
            #If there is no words in that line
            if (complete_line == ""):
                continue                                                       

            #If we have a new sequence
            if (complete_line[0] == prefix): 
                complete_line = complete_line[1:]  #We can use the scling to remove the ">"
                
                words = complete_line.split()          
                
                #The name/id of the new sequence and we later add it to the list
                seq_id = words[0]      
                sequences.append(Sequence(seq_id, ""))                     

            #If we have bases from a sequence    
            else:
                current_sequence = sequences[-1] #We are using the last sequence added to the list
                current_sequence.bases = current_sequence.bases + complete_line # And we modify the bases from that sequence
    
    return sequences


def nucleotide_verification (sequences: list[Sequence], nucleotides: set[str])-> str:
    """
    Identifies if there is a sequence from the fasta file that contains invalid nucleotides

        return: str (name of the sequence). If empty, there is no sequence that has invalid nuclotides
        Arguments:
        - sequences: the list of sequences obtained from a fasta file
        - nucleotides: the verified nucleotides that we can used for our verification
    """

    for sequence in sequences:
        for base in sequence.bases:
            if (base not in nucleotides):
                return sequence.seq_id
    return ""


def sequence_repeated (sequences: list[Sequence], loaded_sequences: dict[str,Sequence])-> str:
    """
    Identifies if there is a sequence from the fasta file that is already in our loaded sequences from other files 
    or if in the same file we have duplicated sequences 
            return: str (name of the sequence). 
            Arguments:
            - sequences: the list of sequences obtained from a fasta file
            - loaded_sequences: we have to see if they already uploaded
    """

    #create an empty set where we will store the identifiers of the 
    #sequences that are in the fasta file. The idea is to check if the sequence
    # is equal to the previous sequences (if there are duplicates in the same file)
    check_id: set[str] = set() 

    for sequence in sequences:
        if(sequence.seq_id in loaded_sequences or sequence.seq_id in check_id): 
            return sequence.seq_id
        
        check_id.add(sequence.seq_id)
      
    #If empty, there is no sequence that is repeated
    return ""


def add_sequences_to_dic (sequences: list[Sequence], loaded_sequences: dict[str,Sequence])-> dict[str,Sequence]:
    """
    Add the valid sequences from a fasta file to the dictionary  
        Return: dict[str,Sequence]
        Arguments:
        - sequences: the list of sequences obtained from a fasta file
        - loaded_sequences: we have to upload the new sequences to the dictionary
    """
    for sequence in sequences:
        loaded_sequences[sequence.seq_id] = sequence

    return loaded_sequences



############## LIST ALL LOADED SEQUENCES STEP (4.3) ##############

def list_sequences_io(sequences: dict[str, Sequence])-> None:
    """
     Prints the number of loeaded sequences and list them. The list is sorted alphabetically by ID, in ascending order using function "sorted()
        Arguments: 
        - sequences: the sequences that will be ordered and printed later
        - the dictionary key is the seq_id
    """

    sequence: Sequence

    print(f"Number of sequences loaded: {len(sequences)}")    
    
    #We use 'sorted()' function to sorted alphabetically by ID, keep in mind we need to use dict()
    # because sorted() returns by default a list, not a dictionary
    
    for id_key in dict(sorted(sequences)):  
        sequence = sequences[id_key]
        print(f" - {sequence.seq_id} ({len(sequence.bases)} bp)")
      


############## SEQUENCE INFORMATION (inspect) STEP (4.4) ##############

def count_bases(sequence: Sequence) -> dict[str, int]:
    """
     In this part the count_bases going to count how many base occurances are in the sequence loaded.
        return: dict[str, int] (the nucleotide letter as key and the count as the integer)
        arguments: sequence

    """
    counts: dict[str, int] = {
       "A": sequence.bases.count("A"),
       "C": sequence.bases.count("C"),
       "G": sequence.bases.count("G"),
       "T": sequence.bases.count("T")
    }
    return counts


def gc_content(sequence: Sequence) -> str:
    """
        Shows the percentage of G and C in the sequence, rounded to 2 decimal places.
        return: str
        arguments: sequence
        Precondition: the sequence has at least one base.
    """
    count_gc: int = sequence.bases.count("G") + sequence.bases.count("C")
    percentage: float = (count_gc / len(sequence.bases)) * 100
    
    return f"{percentage:.2f}" 

  
############## FINDING TE SITE POSITION STEP (4.5)  ##############

"""
#IMPORTANT: In biology, positions are counted starting from 1 (the first base is position 1). Since lists in Python 
#           begin at index 0, so we added 1 (positions.append(i + 1)) to convert the Python index into the correct biological position.

#NOTE: A linear sequence with n distinct cuts is thus split into n + 1 fragments
#     the lengths of the fragments always add up to the length of the sequence.
#In the end, if we have time, we can use the NEBcutter online tool to find out exactly what the actual digestion results would be.

"""

def find_site_positions(sequence: Sequence, enzyme: Enzyme) -> list[int]:
    """
    The function will return a list containing the position numbers where the enzyme recognition site was found.

        Arguments:
        - sequence: Sequence(dataclass) which contains the DNA strand where we will perform the operations
        - enzyme: Enzyme (dataclass) which contains information about the enzyme, specifically its recognition site (enzyme.site)

    """
    positions: list[int] = []
    size = len(enzyme.site)

    # Take into account that we are searching from i to i + size positions
    for i in range(len(sequence.bases) - size + 1): 
        if sequence.bases[i: i + size] == enzyme.site:
            positions.append(i + 1) # i starts at 0, but the user will see from 1-n                                         
    return positions


############## SEQUENCE DIGEST STEP (4.6) ##############

"""
The command "DIGEST" lists the The command DIGEST fragments in this same order, from the top to the bottom of the gel
DIGEST COMMAND (We must pay special attention to the custom sorting logic in case of a tie)
"""

def digest_sequence(sequence: Sequence, enzymes_used: list[Enzyme]) -> list[Fragment]:
    """
    This function will return the fragments of the digest, in order that they appear in the gel or an empty list if there are no cuts.
    AI-assisted function, in the final part of the function related to creating the fragments

        return: list[Fragment] (the fragments of the digest)
        Arguments: 
        - sequence: Sequence (dataclass) contains the linear DNA sequence that will be digested.
        - enzymes_used: objetcs list within Enzyme (list[Enzyme]), contains the restriction enzymes selected by the user to perform cuts in the sequence
    """

    cuts: list[int] = []
    fragments: list[Fragment] = []
    start: int = 1

    for enzyme in enzymes_used:
        for position in find_site_positions(sequence, enzyme):
          #calculate the exact position (considering we start at i = 0 in Python) where the DNA strand will be cut, 
          # based on the site where the enzyme bound.
            cut = position + enzyme.cut_offset - 1

            # the same cut position is only counted once
            if cut not in cuts:                                                                         
                cuts.append(cut)

    # If the selected enzyme didnt find any spots for cutting the sequences            
    if len(cuts) == 0:
        return []

    #Sort the cleavage sites in ascending order so that we can cleave the DNA sequentially.
    cuts.sort()                                        

    #AI-assisted function in this part

    for cut in cuts:
        #Creates a new fragment that extends from the current start point to the cut position.
        fragments.append(Fragment(start, cut))  
         # Updated the start position of the next fragment, considering the previous fragment cut
        start = cut + 1

    # The final piece of DNA—the segment between the last cut and the absolute end of the sequence—still needs to be 
    # processed. This line ensures that this last fragment, which extends to the full length of the sequence (len(sequence.bases)), is also saved.
    fragments.append(Fragment(start, len(sequence.bases)))                                                                      
    
    return fragments


def fragment_length(fragment: Fragment) -> int:
  """
    Returns the length (in bp) of a DNA fragment of a fragment.
        return: int
        Argument: fragment
        - fragment: Fragment Class
        Preconditions: 1 <= start <= end (length of the original sequence).
  """
  return fragment.end - fragment.start + 1 

# If a fragment begins at position 1 and ends at position 3, it contains 3 bases (bases 1, 2, and 3). By adding 1 (3 - 1 + 1 = 3), 
# we ensure that the initial base is also included in the total count.



""" ORGANIZE FRAGMENTS IN GEL
    Creates a list of tuples that are sorted according to the bp
        return: list [tuple[int, Fragment]]
        Arguments:
        - fragments: the previous fragments we created with the sequence_digest
"""
# Creates a list of tuples that are sorted according to the bp
# - return: list [tuple[int, Fragment]]
# Arguments:
#      - fragments: the previous fragments we created with the sequence_digest

def organize_fragments_in_gel(fragments: list[Fragment])->list [tuple[int, Fragment]]:
   
    organized_fragments : list[tuple[int, Fragment]] = []
    for fragment in fragments:
        organized_fragments.append((fragment_length(fragment), fragment))

    organized_fragments = sorted((organized_fragments), reverse = True)

    return organized_fragments



def print_fragments_in_gel (fragments: list[tuple[int, Fragment]])-> None:
    """
    prints the ordered fragments (from higher to lower bp) from the requested digest command
            Arguments:
            - fragments: a list of tupples, includes the lenght of the fragment and the fragment itself
    """
    for length_fragment in fragments:
        length, fragment = length_fragment
        print(f"- {length} bp (positions {fragment.start}-{fragment.end})")




################################ INPUT/OUTPUT FUNCTIONS (DOES NOT INCLUDE FILE READING FUNCTIONS) ###########################################

############## READ STRING and PRINT MESSAGE ##############

def read_string(message: str)-> str:
    """
    This function going to return a String from the input() function
        return: str
        Arguments: 
        - message: A string message used for the input function   
    """ 

    return input(message)


def print_message(message: str)-> None:
    """
    This function prints string messages
        return: none
        Arguments:
        - message: string printed
    """
    print(message)



#########################################################   MAIN FUNCTION  #######################################################
# Main function that acts as a command interpreter.
    #Features are organized in 3 groups:
        #1. Load - read FASTA files (LOAD_SEQUENCES)
        #2. Inspect - list the loaded file and show information about each one (LIST_SEQUENCES, SEQUENCE_INFO)
        #3. Analyse - use enzyme to recognition sites and simulate digests (FIND_SITES, DIGEST)


def main() -> None:        
        
    #STEP 1. Program reads the enzymes files. 
    # if ENZIME not in ENZIME_NAME:  -print (f"The enzyme {ENZIME_NAME} is not available.")
    enzymes: dict[str, Enzyme] = read_enzyme_file(ENZYMES_FILE, DATA_FOLDER)  #Checked
        
    loaded_fasta_files: set [str] = set()
    sequences: dict [str, Sequence] = {}
        
    #Loop control variable used to determine whether the loop should continue executing or stop
    command_name: str = ""


    #STEP2. repeatedly asks the user for a command (using while), executes it, and prints the corresponding messages, until the user types EXIT
    # - WE SHOULD PROVIDE INFORMATION- command EXIT terminates the execution

                
    while command_name != "EXIT":
      #The commands entered by the user do not come alone, they bring additional arguments. 
      # - eg. If the user wants to upload a FASTA file:  LOAD_SEQUENCES demo.fasta ; SEQUENCE_INFO demo ; DIGEST pGFP_insert EcoRI HindIII
                       
      #1. Read user input and split it by spaces-  We must separate the whole command line into smaller strings by the split()
      # - read_string reads all characters received on the serial port and groups them into a single string
      # - The command_parts variable becomes a list with two elements: ["LOAD_SEQUENCES", "demo.fasta"]
    
      command_parts: list[str] = read_string(PROMPT).split()

      #2. check whether the user wants to end the program or has entered something invalid
      #2.1. check whether the user has not entered any command. In that case, return to the beginning of the cycle without executing anything else.
      if not command_parts:   
          continue

      #2.2. stores the first element of the list (the main command) in the command_name variable 
      command_name = command_parts[0]    
                            
      #2.3. Check if the user provides an unknown command: 
      if command_name not in VALID_COMMANDS: 
        print ("Unknown command")
        continue 

      #2.4.The command EXIT terminates the execution saying "Goodbye".
      elif command_name == "EXIT": 
        print("Goodbye!") 
      
      elif command_name == "LOAD_SEQUENCES":
        #To avoid issuss with command_parts, put the value inside file variable
        file: str = command_parts[1]
        
        # 3.1 Check if the file was already uploaded previously 
        if file in loaded_fasta_files:
            print(f"File {file} already loaded!")
            continue

        # 3.2. Check if the nuclotides in the sequences are correct
        # We call the read_fasta_file function to read the FASTA file, and then we call the nucleotide_verification function 
        # to validate whether the nucleotides are valid.
        sequences_in_fasta: list[Sequence] = read_fasta_file(file, DATA_FOLDER, FASTA_HEADER_PREFIX)

        #If empty, there is no sequence that has invalid nuclotides
        invalid_sequence: str = nucleotide_verification(sequences_in_fasta, NUCLEOTIDES)
        if (invalid_sequence != ""):
          print(f"File {file} could not be loaded. Sequence {invalid_sequence} contains invalid nucleotides.")
          continue 

        #3.2. Check if the sequence ID is already in use
        repeated_sequence : str = sequence_repeated(sequences_in_fasta,sequences)
        if (repeated_sequence != ""):
          print (f"File {file} could not be loaded. Sequence ID {repeated_sequence} is already in use.")
          continue 
          
        #3.3. Finally, once all the checks have been completed and everything has passed, load the sequence.   
        add_sequences_to_dic(sequences_in_fasta, sequences)
        loaded_fasta_files.add(file)        
        print(f"File {file} validated and loaded!")


      
      elif command_name == "LIST_SEQUENCES":
            # It the dictionarry of sequences is empty
            if (len(sequences) == 0):
                print("No sequences loaded.")
            else:
                list_sequences_io(sequences)

      
      
      elif command_name == "SEQUENCE_INFO":
        
        sequence_id: str = command_parts[1]

        # 1. We check whether the sequence has been loaded
        if sequence_id not in sequences: 
              print(f"The sequence {sequence_id} is not loaded.")
              # The 'continue' command terminates this command because this condition was true
              continue


        # 2. If the condition above was not met, we will retrieve a specific sequence from the sequences dictionary and store 
        # # it in the sequence variable.
        sequence: Sequence = sequences[sequence_id]

        #3. we call the function count_bases and store the result in the variable 'count' as a dictionary
        counts: dict[str, int] = count_bases(sequence) 
      
        print(f"Information about the sequence {sequence_id}:") 
      
        #4. Its length
        print(f"- Length: {len(sequence.bases)} bp.") 
            
        # The number of occurences of each base
        print(f"- Base counts: A={counts['A']}, C={counts['C']}, G={counts['G']}, T={counts['T']}.")

        #GC content
        print(f"- GC content: {gc_content(sequence)}%.")



      elif command_name == "FIND_SITES":
        sequence_id = command_parts[1] 
        enzyme_name: str = command_parts[2] 

        #The command fails in 3 situations: 
        #1. The sequence is not loaded
        if sequence_id not in sequences: 
          print(f"The sequence{sequence_id} is not loaded.")
          continue 
       
        #2. the enzyme is not available 
        if enzyme_name not in enzymes: 
          print(f"The enzyme {enzyme_name} is not available.") 
          continue

        #3. the enzyme has no recognition site in sequence
        enzyme: Enzyme = enzyme[enzyme_name] 
        position: list[int] = find_site_positions(sequences[sequence_id], enzyme)

        if not positions: 
          print(f"The enzyme {enzyme_name} does not cut the sequence {sequence_id}.")
          continue 


        # If all the checks above have been passed, this command will find all the recognition sites of
        # a given enzyme in a given sequence. 
        positions_output: str = ", ".join(map(str, positions))
        # AI assistent function in this part: map(str, positions) converts each position into text, and join() combines them with commas.
        print(f"Recognition sites of {enzyme_name} ({enzyme.site}) in the sequence {sequence_id}:") 
        print (f"- Number of sites: {len(positions)}.")
        print(f"- Positions: {positions_text}.")


      elif command_name == "DIGEST":
         sequence_id = command_parts[1]
         #reads from position 2 up to an undefined final position (can be 1 or more enzymes)
         enzyme_names: list[str] = command_parts[2:]

         #The command fails in 3 situations: 
         #1. The sequence is not loaded:
         if sequence_id not in sequences:
                print(f"The sequence {sequence_id} is not loaded.")
                continue

         #2. one of the enzymes is not available 
         for enzyme_name in enzyme_names:
           if enzyme_name not in enzymes:
             print(f"The enzyme {enzyme} is not available.")
             break
             
         #3. none of the enzymes cuts the sequence 

         #TENEMOS DE VER ESTA PARTE BIEN 
          if not fragments:
            print(f"The sequence {sequence_id} is not cut by the ""provided enzymes.")
            continue
        
        
        enzymes_output: str = "," ".join(enzyme_names) 
        print(f"Digestion of the sequence {sequence_id} with {enzyme_output} produced the following fragments:") 

        #TENEMOS DE VER ESTA PARTE 
        print_fragments_in_gel


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
