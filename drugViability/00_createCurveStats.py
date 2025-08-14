'''
pull data and run curve fitting code
'''

import pandas as pd
import synapseclient as sc
import os
import subprocess
from datetime import date
#import runpy
syn = sc.login()
##this is a pain, need to move to file

#comboFilelist = syn.tableQuery("select entityId,individualID,specimenID from syn65473033 where dataType='drug screen' AND dataSubtype='processed' AND assay='cell viability assay'").asDataFrame()
# raw files are mistakenly marked as 'processed' so check 'dataSubtype' within file
comboFiles = pd.read_csv(syn.get("syn66330284").path) # previous: syn66330284 but points to combo
comboFiles = comboFiles[comboFiles['dataSubtype']=="processed"]

singleFiles = pd.read_csv(syn.get("syn65473034").path) # previous: syn66330284 but points to combo
singleFiles = singleFiles[singleFiles['dataSubtype']=="processed"]
# remove data from 2025-01-22 due to different method where media was refreshed
singleFiles = singleFiles[~singleFiles.Filename.str.contains("250122.csv")]

###pull files

singledose = []
multidose = []
comboMulti = []
for index,row in comboFilelist.iterrows():
  #print(row['id'])
  dfile = pd.read_csv(syn.get(row['entityId']).path)
  if all(dfile['dataSubtype'] == 'processed'):
      dfile['improve_sample_id'] = row['specimenID']
      dfile = dfile.reset_index()
      
      # replace nan with drug 1 or 2 names as appropriate
      drug1 = dfile['drugOneName'].dropna().unique()[0]
      dfile['drugOneName'].fillna(drug1, inplace=True)
      drug2 = dfile['drugTwoName'].dropna().unique()[0]
      dfile['drugTwoName'].fillna(drug2, inplace=True)
      if all(dfile['drugOneName'] == drug1) & all(dfile['drugTwoName'] == drug2):
          comboMulti.append(dfile)
              
for index,row in singleFiles.iterrows():
  #print(row['id'])
  dfile = pd.read_csv(syn.get(row['entityId']).path)
  if all(dfile['dataSubtype'] == 'processed'):
      dfile['improve_sample_id'] = row['specimenID']
      dfile = dfile.reset_index()
      
      ##get counts of drugs to only keep drugs with >1 dose for curve fitting
      if "drugName" in dfile.columns:
          dcounts = dfile.groupby("drugName").count().reset_index()
          #print(dcounts)
          if any(dcounts['concentration']>1):
              #print("multi")
              more = dcounts[dcounts['concentration']>1]['drugName']
              #print(list(set(more)))
              tempMulti = dfile[dfile['drugName'].isin(list(set(more)))]
              #print(tempMulti)
              multidose.append(tempMulti)
              #print(len(multidose))
          
          ## also get single data points
          if any(dcounts['concentration']==1):
              #print("single")
              sings = dcounts[dcounts['concentration']==1]['drugName']
              #print(list(set(sings)))
              tempSingle = dfile[dfile['drugName'].isin(list(set(sings)))]
              #print(tempSingle)
              singledose.append(tempSingle)
              #print(len(singledose))

####first fit multidose curves...
#os.chdir("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/MPNST-PDX-MT/drugViability/")
if len(multidose) > 0:
    print("compiling multi: ", len(multidose))
    fulltab = pd.concat(multidose)
    #fulltab['DOSE']=fulltab.concentration#+0.0001 # check if should add 0.0001? perhaps cNF data had 0 values
    fulltab = fulltab.rename(columns={"concentration": "DOSE", "percentViability": "GROWTH",
                            "timePoint": "time", "timePointUnit": "time_unit",
                            "drugName": "Drug"})
    fulltab['study']='mpnstPDXMT'
    fulltab['source']='synapse'
    ##mutate the values create new columns
    ncols=['DOSE','GROWTH','study','source','improve_sample_id','Drug','time','time_unit']
    fulltab = fulltab[ncols]
    # there was a percentViability == "147.128*" - seems like it might be 147.1288 but double check
    # if any(fulltab['GROWTH'].str.contains("\*")):
    #     print("replacing asterisk")
    #     fulltab["GROWTH"] = fulltab["GROWTH"].replace(r'\*','', regex=True)
    fulltab['GROWTH'] = pd.to_numeric(fulltab['GROWTH'])
    ##change file headers to DOSE/RESPONSE values needed by other script
    fulltab.to_csv('mpnst_drug_response.tsv',sep='\t')
    
    ##fit curve
    #if not os.path.exists("fit_curve.py"):
    script='https://raw.githubusercontent.com/PNNL-CompBio/coderdata/refs/heads/main/build/utils/fit_curve.py'
    subprocess.run(['wget',script])
    #subprocess.run(['python3','-m','pip','install','matplotlib']) # due to matplotlib error
    subprocess.run(['python3','fit_curve.py','--input','mpnst_drug_response.tsv','--output','mpnstDrugOutput'])
    #runpy.run_path("fit_curve.py") # replaces subprocess call
    #os.getcwd()
    otab = pd.read_csv('mpnstDrugOutput.0',sep='\t')
else:
    otab = pd.DataFrame()

#####now we can take single drug points and format those
if len(singledose) > 0:
    print("compiling single: ", len(singledose))
    stab = pd.concat(singledose)
    stab = stab.rename(columns={"drugName": "improve_drug_id", 
                            "timePoint": "time", "timePointUnit": "time_unit"})
    stab['study']='mpnstPDXMT'
    stab['source']='synapse'
    stab['percentViability'] = pd.to_numeric(stab['percentViability'])
    stab['dose_response_value'] = stab.percentViability/100.00
    stab['dose_response_metric'] = 'uM_viability' # not sure, but concentrationUnit='uM' and assay='cell viability assay', platform='3D CellTiter-Glo'

    curve_cols = ['source','improve_sample_id','improve_drug_id','study','time','time_unit',
              'dose_response_metric','dose_response_value']
    stab = stab[curve_cols]
else:
    stab = pd.DataFrame()
    
    
newtab =pd.concat([otab,stab])
##rename file
newtab.to_csv('viability_curves.tsv',sep='\t',index=False)

##now add in single-point drug measurements

#store on synapse
syn.store(sc.File('mpnst_drug_response.tsv',parentId='syn52369034'))
syn.store(sc.File('viability_curves.tsv',parentId='syn52369034'))
    
####first fit multidose curves...
if len(comboMulti) > 0:
    print("compiling multi combos: ", len(comboMulti))
    fulltab = pd.concat(comboMulti)
    fulltab = fulltab[fulltab['drugThreeName'].isna()]
    
    # replace blanks with 0 values
    fulltab['drugOneConcentration'].fillna(0, inplace=True)
    fulltab['drugTwoConcentration'].fillna(0, inplace=True)
    
    #fulltab['DOSE']=fulltab.concentration#+0.0001 # check if should add 0.0001? perhaps cNF data had 0 values
    # combine cell line and time into "sample" column
    fulltab['sample'] = fulltab['sampleName'] + "_" + fulltab['timePoint'].astype(str) + fulltab['timePointUnit'] # if need replicate, then use specimenID
    fulltab['drug1.units'] = fulltab['concentrationUnit']
    fulltab['drug2.units'] = fulltab['concentrationUnit']
    fulltab = fulltab.rename(columns={"percentViability": "effect",
                            "drugOneName": "drug1", "drugTwoName": "drug2",
                            "drugOneConcentration": "drug1.conc", "drugTwoConcentration": "drug2.conc"})
    # there was a percentViability == "147.128*" - seems like it might be 147.1288 but double check
    # if any(fulltab['effect'].str.contains("\*")):
    #     print("replacing asterisk")
    #     fulltab["effect"] = fulltab["effect"].replace(r'\*','', regex=True)
    fulltab['effect'] = pd.to_numeric(fulltab['effect'])/100.00 # convert to fraction e.g., 0.9 instead of 90%
    ##change file headers to DOSE/RESPONSE values needed by other script
    fulltab['expt.date'] = date.today()
    ncols=['expt.date','drug1.conc','drug2.conc','effect','sample','drug1','drug2','drug1.units','drug2.units']
    fulltab2 = fulltab[ncols]
    fulltab2.to_csv('mpnst_combo_drug_response.csv', index=False)
    
    # also create redacted version for upload to MuSyC web app
    fulltab3 = fulltab
    fulltab3['sampleName'] = fulltab3['sampleName'].str[-3:]
    fulltab3['sample'] = fulltab3['sampleName'] + "_" + fulltab3['timePoint'].astype(str) + fulltab3['timePointUnit'] # if need replicate, then use specimenID
    fulltab3['drug1'] = fulltab3['drug1'].str[1:4]
    fulltab3['drug2'] = fulltab3['drug2'].str[1:4]
    fulltab3 = fulltab3[ncols]
    fulltab3.to_csv('combo_drug_response.csv', index=False) # MuSyC expects CSV file

##now add in single-point drug measurements

#store on synapse
syn.store(sc.File('mpnst_combo_drug_response.csv',parentId='syn52369040'))

