'''
pull data and run curve fitting code
'''

import pandas as pd
import synapseclient as sc
import os
import subprocess
syn = sc.login()
##this is a pain, need to move to file

filelist = syn.tableQuery("select entityId,individualID,specimenID from syn65473033 where dataType='drugScreen' AND dataSubtype='processed' AND assay='cell viability assay'").asDataFrame()
# raw files are mistakenly marked as 'processed' so check 'dataSubtype' within file

###pull files

singledose = []
multidose = []
for index,row in filelist.iterrows():
  #print(row['id'])
  dfile = pd.read_csv(syn.get(row['entityId']).path)
  if all(dfile['dataSubtype'] == 'processed'):
      dfile['improve_sample_id'] = row['specimenID']
      dfile = dfile.reset_index()
      
      ##get counts of drugs to only keep drugs with >1 dose for curve fitting
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
    if any(fulltab['GROWTH'].str.contains("\*")):
        print("replacing asterisk")
        fulltab["GROWTH"] = fulltab["GROWTH"].replace(r'\*','', regex=True)
    fulltab['GROWTH'] = pd.to_numeric(fulltab['GROWTH'])
    ##change file headers to DOSE/RESPONSE values needed by other script
    fulltab.to_csv('mpnst_drug_response.tsv',sep='\t')
    
    ##fit curve
    script='https://raw.githubusercontent.com/PNNL-CompBio/coderdata/refs/heads/main/build/utils/fit_curve.py'
    subprocess.run(['wget',script])
    subprocess.run(['python','fit_curve.py','--input','mpnst_drug_response.tsv','--output','mpnstDrugOutput'])
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
syn.store(sc.File('viability_curves.tsv',parentId='syn52369034'))

