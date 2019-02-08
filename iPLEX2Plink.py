#!/usr/bin/env python3

import argparse
import requests
import json

##############################################
# args

parser = argparse.ArgumentParser(description='Converts an iPLEX MassARRAY data file to Plink .ped and .map files')
parser.add_argument('--iplex',required=True,help="iPLEX file")
parser.add_argument('--output',required=True,help="output .ped and .map base name")

args = parser.parse_args()
# print(args)


###############################################
# functions

def writemap(line,mapfile):
    elements = line.strip().split()[1:]
    outmap = open(mapfile,"w")
    for e in elements:
        info = ["0","0","0","0"]
        if e.startswith("rs"):
            # print(e[2:])
            info = getdbsnp(e[2:])
        outmap.write("\t".join([info[0],e,"0",info[1]])+"\n")
    outmap.close()

def getdbsnp(rsid):
    url = "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/"+rsid
    web = requests.get(url)
    rs_obj = json.loads(web.text)
    if 'primary_snapshot_data' in rs_obj and \
        rs_obj['primary_snapshot_data']['variant_type'] == "snv":
        info = placements(rs_obj['primary_snapshot_data']['placements_with_allele'])
        return(info)

def placements(info):
    '''
    rs genomic positions
    '''
    for alleleinfo in info:
        if len(alleleinfo['placement_annot']['seq_id_traits_by_assembly']) > 0:
            assembly_name = alleleinfo['placement_annot'] \
                                      ['seq_id_traits_by_assembly'] \
                                      [0]['assembly_name']
            if(assembly_name == "GRCh37.p13"):
                for a in alleleinfo['alleles']:
                    spdi = a['allele']['spdi']
                    if spdi['inserted_sequence'] != spdi['deleted_sequence']:
                        (ref, alt, pos, seq_id) = (spdi['deleted_sequence'],
                                                   spdi['inserted_sequence'],
                                                   spdi['position'],
                                                   spdi['seq_id'])
                        chr = int(seq_id[3:9])
                        return([str(chr), str(pos+1), ref, alt])
                        break


###############################################
# parse

ped = args.output+".ped"
map = args.output+".map"

outped = open(ped,"w")
with open(args.iplex) as f:
    rsids = f.readline()
    writemap(rsids,map)
    for line in f:
        outl = ""
        elements = line.strip().split()
        id = elements[0]
        outl = outl + id + "\t" + id + "\t" + "0" + "\t" + "0" + "\t" + "0" + "\t" + "0"
        elements = elements[1:]
        for e in elements:
            if e == "NA":
                e = "0/0"
            alleles = e.split("/")
            outl = outl + "\t" + alleles[0] + "\t" + alleles[1]
        outped.write(outl+"\n")
outped.close()
