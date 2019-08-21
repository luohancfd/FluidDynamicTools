#!/usr/bin/env python

import bibtexparser,re,os

from bibtexparser.bwriter import BibTexWriter

os.chdir('/home/Draft')
with open('DSMCMF2.bib','r', encoding='utf-8') as bibtex_file:
    bibtex_str = bibtex_file.read()

bib_database = bibtexparser.loads(bibtex_str)

AIAA_conference= [('2018','Aerospace Sciences Meeting','Kissimmee','Florida'),
                  ('2018','Joint Thermophysics', 'Atlanta', 'Georgia'),
                  ('6th', 'Joint Thermophysics', 'Colorado Springs','CO'),
                  ('47th','Aerospace Sciences Meeting','Orlando', 'Florida',),
                  ('47th','Thermophysics','Denver', 'Colorado'),
                  ('43rd','Thermophysics','New Orleans', 'Louisiana'),
                  ('52nd','Aerospace Sciences Meeting','National Harbor', 'Maryland'),
                  ('55th','Aerospace Sciences Meeting','Grapevine','Texas'),
                  ('54th','Aerospace Sciences Meeting','San Diego','California'),
                  ('None','None','None','None')]
bib_post = bib_database
for i in bib_post.entries:
    if i["ENTRYTYPE"] == "inproceedings" :
        print(i)
        confname = i["booktitle"]
        for conf in AIAA_conference:
            if conf[1] in confname and conf[0] in confname:
                break
        if conf[0] == 'None':
            print('No conference data found for the following:')
            print(i)
        else:
            doi = i["doi"]
            if doi:
                w =  re.findall("\d\.(\d{4})-(\d{4})",doi)
                if w:
                    year = w[0][0]
                    number = w[0][1]
                    if "address" in i:
                        i.pop('address',None)
                    if "publisher" in i:
                        i.pop('publisher',None)
                    
                    i["series"] = "AIAA Paper %s-%s"%(year,number)
                    i["year"] = year
                    print(i)
                else:
                    print("Doi number is not correct for the following:")
                    print(i)
            else:
                print("No doi found for the following:")
   # i["title"] = i["title"].replace("O2","O$_2$")
   # i["title"] = i["title"].replace("N2","N$_2$")
   # i["title"] = i["title"].replace("\mathrm{_2}","_2")


writer = BibTexWriter()
with open('DSMCMF2.bib', 'w',encoding='utf-8') as bibfile:
    bibfile.write(writer.write(bib_post))               
        