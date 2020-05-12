#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 15:27:35 2020

@author: sshanto
structure:
    > import data files into dataframe
    > for each data file
        > only keep events that had mu-
        > only keep events that had a 4/4 coincidence
    > compare events from control case and wtp case
        > for each of these events
            > if event was in control and not in wtp keep them in Bin1
                > check the loc information of such event
                > if last loc is in wtp/tank/tube
                    > the mu- in this event was absorbed
                > else
                    > the mu- in this was scattered
            > if event was in wtp and not in control keep them in Bin2
                > these events were scattered in 
    > record the numbers of scattered_in, scattered_out and absorbed events
    > calculate the required numbers
    
@Note: 
    >> need to fix data format for name of trays
    >> need to rethink analysis if no common events found
        >> possible tracking of each muons not gone through the tower 
    
#results:
    > more events with wtp than air
    > only one common event!?
    > check for mu+ included case
    > 12% absorption, 88% scattering
    > differences in the same mu event
"""

#============================== Libraries ==============================
import pandas as pd

control = pd.read_csv('processed_control_50k.txt')
control.columns = control.columns.str.strip() 

wtp = pd.read_csv('processed_wtp_50k.txt')
wtp.columns = wtp.columns.str.strip() 

#============================== Constants, Variables and Inputs ==============================
arr_name_mv = ['SC8 0','Station1 0', 'Tank1 0', 'Tank2 0', 'Tank3 0', 'Tower 0', 'Tray1 0', 'Tray2 1', 'Tray3 2', 'Tray4 3','Tray1 1','Tray1 2', 'Tray1 3','World 0', 'tube 0']                    
arr_name_det_control = ['Tray1 0', 'Tray2 1', 'Tray3 2', 'Tray4 3']
arr_name_det_wtp = ['Tray1 0','Tray1 1','Tray1 2', 'Tray1 3']
arr_name_loc_tower = ['Tower 0','Tank1 0','Tank2 0','Tank3 0','WatTank1 0','WatTank3 0','tube 0','tubeWater 0',' Tank1 0',' Tank2 0',' Tank3 0',' Tower 0',' WatTank1 0',' WatTank3 0',' tube 0',' tubeWater 0']
#============================== Functions ==============================

def filter_muon_minus_only(df):  #function only allows events involving mu-
    df = df[(df.part=='mu-')].reset_index()
    return df

def debug_print(param1, param2):   #function to display and debug
    print(param1)
    print("\n\n")
    print(param2)
            
def group_ana(df):
    for name, group in df:
        print("muon id: " + str(name) + ": " + str(len(group)))
      #  print(str(group.mother) in arr_name_det)
        print(group.mother[0])
        print("\n")

def fbyf_wtp(R1):
    R1_four = R1.query('mother=="Tray1 0" or mother=="Tray1 1" or mother=="Tray1 2" or mother=="Tray1 3"').reset_index()
    mu_ids = list(R1_four['mu_id'])
    mvs = list(R1_four['mother'])
    arr = list(set(zip(mu_ids,mvs)))
    arr_new = [x[0] for x in arr ]
    freqs = pd.DataFrame(pd.value_counts(arr_new), columns = ['freq']).reset_index()
    true_muons = freqs.query('freq==4') #this is the data frame
    evs = list(true_muons['index'])
    return evs

def fbyf_control(R1):
    R1_four = R1.query('mother=="Tray1 0" or mother=="Tray2 1" or mother=="Tray3 2" or mother=="Tray4 3"').reset_index()
    mu_ids = list(R1_four['mu_id'])
    mvs = list(R1_four['mother'])
    arr = list(set(zip(mu_ids,mvs)))
    arr_new = [x[0] for x in arr ]
    freqs = pd.DataFrame(pd.value_counts(arr_new), columns = ['freq']).reset_index()
    true_muons = freqs.query('freq==4') #this is the data frame
    evs = list(true_muons['index'])
    return evs

def common(a,b): 
    c = [value for value in a if value in b] 
    return c

def delete__by_values(lst, values):
    values_as_set = set(values)
    return [ x for x in lst if x not in values_as_set ]

def det_miss_ana(R1_det, R2_det, R1_miss, R2_miss):
    print("Detected: \n")
    debug_print(R1_det, R2_det)
    print()
    print("Missed: \n")
    debug_print(R1_miss, R2_miss)

def return_last_mu_index(R1_four):
    indexes = []
    for row in range(len(R1_four['mu_id'])):
        j = row+1
        if j < len(R1_four['mu_id']):
            if R1_four['mu_id'][row] != R1_four['mu_id'][j]:    
                indexes.append(R1_four['index'][row])
        if row ==(len(R1_four['mu_id'])-1):
            indexes.append(R1_four['index'][row])
    return indexes

def scatter_abs_ana(R1_ana):
    absmu = 0
    scatmu = 0
    for i in R1_ana['loc']:
        if (i in arr_name_loc_tower):
            absmu+=1
        else:
            scatmu+=1
    print("absorbed muon events: " + str(absmu * 100. / (scatmu+absmu)))
    print("scattered muon events: " + str(scatmu * 100. / (scatmu+absmu)))
    return (absmu, scatmu)

#==============================  Start of main code ==============================
#only mu-
wtp = filter_muon_minus_only(wtp)   
control = filter_muon_minus_only(control)
#reducing data fields
R1 = wtp.filter(items=['mu_id', 'loc', 'mother'])   
R2 = control.filter(items=['mu_id', 'loc', 'mother'])
#only 4/4 events index list
R1_ev = fbyf_wtp(R1)
R2_ev = fbyf_control(R2)
#removing common events in the two cases
common_ele = common(R1_ev,R2_ev)
R1_ev = delete__by_values( R1_ev, common_ele)
R2_ev = delete__by_values( R2_ev, common_ele)
#events with 4/4 coincidence 
R1_det = R1.loc[R1['mu_id'].isin(R1_ev)].reset_index() 
R2_det = R2.loc[R2['mu_id'].isin(R2_ev)].reset_index()
#events that didn't have a 4/4
R1_miss = R1.loc[~R1['mu_id'].isin(R1_ev)].reset_index()
R2_miss = R2.loc[~R2['mu_id'].isin(R2_ev)].reset_index() 
#Scatt1ring-absorption analysis:
ind_r1 = return_last_mu_index(R1_miss)
ind_r2 = return_last_mu_index(R2_miss)
#To do:
R1_ana = R1_miss.loc[R1_miss['index'].isin(ind_r1)].reset_index() 
R2_ana = R2_miss.loc[R2_miss['index'].isin(ind_r2)].reset_index()    
#results in tuple format
results_r1 = scatter_abs_ana(R1_ana)
results_r2 = scatter_abs_ana(R2_ana)


#print(R1_ana.info())
#print(R1_ana)
