#%%
#!/usr/bin/env python3
from pybtex.database import parse_file
input_file = 'all_reduced.bib'
bib_database = parse_file(input_file)

#%%
for k in bib_database.entries.keys():
    item = bib_database.entries[k]
    persons = item.persons
    print(k)
    for i in range(len(persons['author'])):
        oldFirstNames = persons['author'][i].first_names
        newFirstNames = []
        for j in oldFirstNames:
            if j[-1] == '.':
                newFirstNames.append(j)
            else:
                newFirstNames.append(j[0].upper() + '.')
        if ' '.join(oldFirstNames) != ' '.join(newFirstNames):
            print('{:s} => {:s}'.format(' '.join(oldFirstNames), ' '.join(newFirstNames)))
            persons['author'][i].first_names = newFirstNames

#%%
authors = [ ' and '.join(['{:s}, {:s} {:s}'.format(' '.join(p.last_names),  ' '.join(p.first_names), ' '.join(p.middle_names))
                      for p in v.persons['author']]) for v in bib_database.entries.values() ]

#%%
# I don't like the default writer of pybtex
# It prefers double quote instead of curly bracket
import bibtexparser
from bibtexparser.bwriter import BibTexWriter
with open(input_file,'r', encoding='utf-8') as bibtex_file:
    bibtex_str = bibtex_file.read()
bib = bibtexparser.loads(bibtex_str)
for k in bib.entries_dict:
    bib.entries_dict[k]['author'] = ' and '.join(['{:s}, {:s} {:s}'.format(' '.join(p.last_names),  ' '.join(p.first_names), ' '.join(p.middle_names))
                      for p in bib_database.entries[k].persons['author']])
writer = BibTexWriter()
with open('all_reduced_abbrv.bib', 'w',encoding='utf-8') as bibfile:
    bibfile.write(writer.write(bib))


# %%
