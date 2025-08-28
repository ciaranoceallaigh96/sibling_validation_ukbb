##SCRIPT FOR REMOVING RELATED INDIVIDUALS FROM UKBB DATA (WILL LEAVE ONE PERSON PER FAMILY CLUSTER)

import random

related_file = "ukb_rel_0884.txt"  #List of related pairs
exclude_file = "related_exclude_list2.txt"  #Individuals to exclude

#Set to track excluded individuals
excluded = set()
kept = set()

#Process related pairs
with open(related_file, 'r') as infile, open(exclude_file, 'w') as outfile:
    for line in infile:
        person1, person2 = line.strip().split()[:2]

        #Check if one of them is already excluded
        if person1 in excluded or person2 in excluded:
            #Ensure the already excluded individual remains excluded
            to_exclude = person1 if person1 in excluded else person2
        else:
            #Randomly exclude one if neither is excluded
            to_exclude = random.choice([person1, person2])

        #Update sets
        excluded.add(to_exclude)
        kept.add(person1 if to_exclude == person2 else person2)

        #Write to exclude file
        outfile.write(f"{to_exclude}\n")



