# f8ForeignStuff
This is the code to find the foreign 15mers. It is specifically written around the MLOF database. (File not included due to PII)

## Disclaimer 
I don't think this will work as is but it should at least explain the process.

## the files

- get_foreign_exons.py 
    - this takes the file main.csv and reads column 44 or 41 to see what the mutation is. |There are hard coded rules for I22 and I1 inversions. otherwise it uses functions located in the utility_scripts_for_get_foreign.py to mutate the nucleotide sequence and recode the protein. the helper functions return a list of foreign 15mers
  
- utility_scripts_for_get_foreign_exons.py 
  - get_affinitiy this uses a file of affinities and alleles I pregenerated. from memory it is the netMHCIIpan3.1 affinities
  - the add non-clinical mutation references a file I will email
  - the turn_f8_string_to_foreign_15mers function just matches what is in the wild-type f8 with the subjects new protein and returns what doesn't match
  - the rest of the functions are mutation specific. I feel like one or two mutation types might not be accounted for

- save_for_affinities_mlof.py is how I was saving the affinities. It will take some work to recreate this but on the plus side, my username and password are in the file. the database may still be on the workstation.

- add_foreign_to_all.py is an older version of get_foreign_exons.py I believe that the main difference is that it does not add HLA. The newer file does.

## General workflow

get_foreign exons reads the mlof data files. looks at the mutations, figures out what is foreign and appends that to the data, calculates the patient specific affinities for DRB1 and writes MLOF_final_drb_dpb_dqb_and_foreign.csv
