#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 15:58:18 2020

@author: shanto

structure:
    >> distinguish between non muon and muon events (X)
    >> append event num and num mu generated into data lines
    >> get rid of E lines
    >> add the header line to the processed text file 
    >>> send the file to pandas code
    >>>> analysis of events
"""
#libraries:
import sys

#string parser for data line formatting
def string2num(st):
    st = st.split("\n")
    st = st[0]
    st = st.split("=")
    st = st[1]
    return st

def read(inp, out):
    event = "1,1"
    len_of_ev_line = 89
    
    #code to distinguish between non muon and muon events
    event_id = []
    no_mu_event = []
    mu_events = []
    
    o = open(out, "w")
    f = open(inp, "r")
    o.write("mu_id, num_mu, tk, step, part, ke, edep, x, y, z, depth, loc, mother, status,\n")
    lines = f.readlines()
    for line in enumerate(lines):
       # print(line)
        if line[1][0] == "E":
           # print(line)
            event_id.append(line[0])
            next_line = (line[0]+1)
        if (line[0] == (line[0]+1)) and (line[1] == "\n"):
            no_mu_event.append(next_line-1)
    for ev in event_id:
        if ev  not in no_mu_event:
            mu_events.append(ev) 
            
    #append event num and num mu generated into data lines
    for line in enumerate(lines):
        if line[0] not in no_mu_event:
            if line[0] in mu_events:
            #    event = line[1][2:5]  #OG
                event = string2num(line[1])
              #  print(event)
            else:
             #   print(event+","+line[1])
                if len(event+","+line[1]) >= len_of_ev_line:
                    o.write(event+","+line[1])
                
    f.close()
    o.close()

inp = sys.argv[1]
outp = sys.argv[2]

read(inp, outp)

    