'''
pull data and run curve fitting code
'''

from dataclasses import dataclass

import pandas as pd
import synapseclient as sc
import os
import subprocess
from datetime import date
import sys
#import runpy


@dataclass
class Study:
    """Data class to capture crutial study information"""
    name: str
    data_folder_synapse_id: str
    result_folder_synapse_id: str
    combination_agent_trial: bool
    

def main():
    """
    Main script calling the individual curve fitting & drug response
    metric calculations.
    """

    """
    The lines below setup the studies with the respective synapse ids
    pointing to the folders that contain the data files, as well as the
    synapse folder where results should be stored.

    To include another seperate study, add another `Study` object
    (including the defined class variables) into the dict(k,v) with the
    study name as key k and the `Study` object as value v 
    """
    studies = {
        'combo': Study(
            name='combo',
            data_folder_synapse_id='syn66330226',
            result_folder_synapse_id='syn52369040',
            combination_agent_trial=True,
        ),
        'mek+egfr': Study(
            name='mek+egfr',
            data_folder_synapse_id='syn71857386',
            result_folder_synapse_id='syn71857385',
            combination_agent_trial=True,
        ),
        'single': Study(
            name='single',
            data_folder_synapse_id='syn65473019',
            result_folder_synapse_id='syn52369034',
            combination_agent_trial=False,
        ),    
    }

    """
    Define for which of the studies (see the `studies` dictionary above)
    the analysis should be run. This is a simple list of keys / names.
    This can be as little as `studies_to_run = ['study_name'] to run
    only one study analysis. 
    """
    studies_to_run = ['combo', 'mek+egfr', 'single']

    """
    Calls to the functions that contain the analysis logic.
    """
    for study_name in studies_to_run:
        if study_name in studies:
            study = studies[study_name]
        else:
            print(f'Study "{study_name}" not defined!', file=sys.stderr)
            continue
        
        if study.combination_agent_trial:
            combination_agent_trial(
                study.data_folder_synapse_id,
                study.result_folder_synapse_id
            )
        elif not study.combination_agent_trial:
            single_agent_trial(
                study.data_folder_synapse_id,
                study.result_folder_synapse_id
            )


def single_agent_trial(in_data_folder_id: str, result_folder_id: str) -> None:
    """
    Helper function collecting the processing logic to single agent 
    studies

    Parameters
    ----------
    in_data_folder_id : str
        synapse id for folder that contains all datafiles for ingestion
    result_folder_id : str
        synapse id for folder that resul files will be stored to

    Returns
    -------
    None
    """

    syn = sc.login()

    ##pull files
    
    #singleFiles = pd.read_csv(syn.get("syn65473034").path) # previous: syn66330284 but points to combo
    #singleFiles = singleFiles[singleFiles['dataSubtype']=="processed"]

    singleFiles = syn.tableQuery(f"SELECT * FROM syn52369043 where parentId='{in_data_folder_id}'")
    singleFiles = singleFiles.asDataFrame()
    # remove data from 2025-01-22 due to different method where media was refreshed
    singleFiles = singleFiles[~singleFiles.name.str.contains("250122.csv")]

    singledose = []
    multidose = []

    for index,row in singleFiles.iterrows():
        #print(row['id'])
        dfile_raw = pd.read_csv(syn.get(row['id']).path)
        dfile = dfile_raw.drop(dfile_raw.index[dfile_raw.isna().all(axis=1)])
        if all(dfile['dataSubtype'] == 'processed'):
            dfile['improve_sample_id'] = row['specimenID']
            if pd.isna(row['specimenID']):
                dfile['improve_sample_id'] = dfile['sampleName']
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

    ##first fit multidose curves...
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
        #if any(fulltab['GROWTH'].str.contains("\*")):
        #     print("replacing asterisk")
        #     fulltab["GROWTH"] = fulltab["GROWTH"].replace(r'\*','', regex=True)
        fulltab['GROWTH'] = pd.to_numeric(fulltab['GROWTH'])
        ##change file headers to DOSE/RESPONSE values needed by other script
        fulltab.to_csv('mpnst_drug_response.tsv',sep='\t')
        
        ##fit curve
        script='https://raw.githubusercontent.com/PNNL-CompBio/coderdata/refs/heads/main/coderbuild/utils/fit_curve.py'
        if not os.path.exists("fit_curve.py"):
            subprocess.run(['wget',script, '-O', 'fit_curve.py'])
            #subprocess.run(['python3','-m','pip','install','matplotlib']) # due to matplotlib error
        subprocess.run([
            'python', 'fit_curve.py',
            '--input', 'mpnst_drug_response.tsv',
            '--output', 'mpnstDrugOutput'
        ])
        #runpy.run_path("fit_curve.py") # replaces subprocess call
        #os.getcwd()
        otab = pd.read_csv('mpnstDrugOutput.0',sep='\t')
    else:
        otab = pd.DataFrame()

    #####now we can take single drug points and format those
    ##there should be none for this dataset
    # if len(singledose) > 0:
    #     print("compiling single: ", len(singledose))
    #     stab = pd.concat(singledose)
    #     stab = stab.rename(columns={"drugName": "improve_drug_id", 
    #                             "timePoint": "time", "timePointUnit": "time_unit"})
    #     stab['study']='mpnstPDXMT'
    #     stab['source']='synapse'
    #     stab['percentViability'] = pd.to_numeric(stab['percentViability'])
    #     stab['dose_response_value'] = stab.percentViability#/100.00
    #     stab['dose_response_metric'] = 'uM_viability' # not sure, but concentrationUnit='uM' and assay='cell viability assay', platform='3D CellTiter-Glo'

    #     curve_cols = ['source','improve_sample_id','improve_drug_id','study','time','time_unit',
    #               'dose_response_metric','dose_response_value']
    #     stab = stab[curve_cols]
    # else:
    #     stab = pd.DataFrame()
        
        
    newtab = otab #pd.concat([otab,stab])
    ##rename file
    newtab.to_csv('viability_curves.tsv',sep='\t',index=False)

    ##now add in single-point drug measurements

    #store on synapse
    # syn.store(sc.File('mpnst_drug_response.tsv',parentId=result_folder_id))
    # syn.store(sc.File('viability_curves.tsv',parentId=result_folder_id))

    return None


def combination_agent_trial(in_data_folder_id: str, result_folder_id: str) -> None:
    """
    Helper function collecting the processing logic to multi agent / 
    drug combination studies

    Parameters
    ----------
    in_data_folder_id : str
        synapse id for folder that contains all datafiles for ingestion
    result_folder_id : str
        synapse id for folder that resul files will be stored to

    Returns
    -------
    None
    """

    syn = sc.login()

    #comboFilelist = syn.tableQuery("select entityId,individualID,specimenID from syn65473033 where dataType='drug screen' AND dataSubtype='processed' AND assay='cell viability assay'").asDataFrame()
    # raw files are mistakenly marked as 'processed' so check 'dataSubtype' within file
    #comboFiles = pd.read_csv(syn.get("syn66330284").path) # previous: syn66330284 but points to combo
    #comboFiles = comboFiles[comboFiles['dataSubtype']=="processed"]

    comboFiles = syn.tableQuery(f"SELECT * FROM syn52369043 where parentId='{in_data_folder_id}'")
    comboFiles = comboFiles.asDataFrame()

    comboMulti = []
    for index,row in comboFiles.iterrows():
        #print(row['id'])
        dfile_raw = pd.read_csv(syn.get(row['id']).path)
        dfile = dfile_raw.drop(dfile_raw.index[dfile_raw.isna().all(axis=1)])
        if all(dfile['dataSubtype'] == 'processed'):
            dfile['improve_sample_id'] = row['specimenID']
            if pd.isna(row['specimenID']):
                dfile['improve_sample_id'] = dfile['sampleName']
            dfile = dfile.reset_index()
            
            # replace nan with drug 1 or 2 names as appropriate
            drug1 = dfile['drugOneName'].dropna().unique()[0]
            dfile['drugOneName'] = dfile['drugOneName'].fillna(drug1)#.fillna(drug1, inplace=True)
            drug2 = dfile['drugTwoName'].dropna().unique()[0]
            dfile['drugTwoName'] = dfile['drugTwoName'].fillna(drug2)#.fillna(drug2, inplace=True)

            if all(dfile['drugOneName'] == drug1) & all(dfile['drugTwoName'] == drug2):
                comboMulti.append(dfile)
                



        
    ####now we move to the drug combinations
    if len(comboMulti) > 0:
        print("compiling multi combos: ", len(comboMulti))
        fulltab = pd.concat(comboMulti)
        # fulltab = fulltab[fulltab['drugThreeName'].isna()]
        
        # replace blanks with 0 values
        fulltab['drugOneConcentration'].fillna(0, inplace=True)
        fulltab['drugTwoConcentration'].fillna(0, inplace=True)
        
        #fulltab['DOSE']=fulltab.concentration#+0.0001 # check if should add 0.0001? perhaps cNF data had 0 values
        # combine cell line and time into "sample" column
        fulltab['sample'] = fulltab['sampleName'].astype(str) + "_" + fulltab['timePoint'].astype(str) + fulltab['timePointUnit'].astype(str) # if need replicate, then use specimenID
        fulltab['drug1.units'] = fulltab['concentrationUnit']
        fulltab['drug2.units'] = fulltab['concentrationUnit']
        fulltab = fulltab.rename(columns={"percentViability": "effect",
                                "drugOneName": "drug1", "drugTwoName": "drug2",
                                "drugOneConcentration": "drug1.conc", "drugTwoConcentration": "drug2.conc"})
        # there was a percentViability == "147.128*" - seems like it might be 147.1288 but double check
        #if any(fulltab['effect'].str.contains("\*")):
        #     print("replacing asterisk")
        #     fulltab["effect"] = fulltab["effect"].replace(r'\*','', regex=True)
        #we need percentages (SG):
        fulltab['effect'] = pd.to_numeric(fulltab['effect'])#/100.00 # convert to fraction e.g., 0.9 instead of 90%
        ##change file headers to DOSE/RESPONSE values needed by other script
        fulltab['expt.date'] = date.today()
        ncols=['expt.date','drug1.conc','drug2.conc','effect','sample','drug1','drug2','drug1.units','drug2.units']
        fulltab2 = fulltab[ncols]
        fulltab2.to_csv('mpnst_combo_drug_response.csv', index=False)

        
        # also calculate auc values for the combinations
        fulltab['improve_sample_id'] = fulltab['improve_sample_id']+'_'+fulltab['drug1']+'+'+fulltab['drug2']
            
        tabs1 = fulltab[(fulltab['drug1.conc']>0) & (fulltab['drug2.conc']==0)]
        tabs1 = tabs1.rename(columns={"drug1.conc": "DOSE", "effect": "GROWTH",
                                "timePoint": "time", "timePointUnit": "time_unit",
                                "drug1": "Drug"})
        tabs2 = fulltab[(fulltab['drug2.conc']>0) & (fulltab['drug1.conc']==0)]
        tabs2 = tabs2.rename(columns={"drug2.conc": "DOSE", "effect": "GROWTH",
                                "timePoint": "time", "timePointUnit": "time_unit",
                                "drug2": "Drug"})
        tabsCombo1 = fulltab[(fulltab['drug2.conc']>0) & (fulltab['drug1.conc']>0)]
        tabsCombo1['Drug'] = tabsCombo1['drug1']+'+'+tabsCombo1['drug2']#+'('+round(tabsCombo1['drug1.conc']/tabsCombo1['drug2.conc'],2).astype(str)+')1'
        tabsCombo1 = tabsCombo1.rename(columns={"drug1.conc": "DOSE", "effect": "GROWTH",
                                "timePoint": "time", "timePointUnit": "time_unit"})
        tabsCombo2 = fulltab[(fulltab['drug2.conc']>0) & (fulltab['drug1.conc']>0)]
        tabsCombo2['Drug'] = tabsCombo2['drug1']+'+'+tabsCombo2['drug2']#+'('+round(tabsCombo2['drug1.conc']/tabsCombo2['drug2.conc'],2).astype(str)+')2'
        tabsCombo2 = tabsCombo2.rename(columns={"drug2.conc": "DOSE", "effect": "GROWTH",
                                "timePoint": "time", "timePointUnit": "time_unit"})


        
        ncols0=['DOSE','GROWTH','improve_sample_id','Drug','time','time_unit']
        fulltab = pd.concat([tabs1[ncols0], tabs2[ncols0], tabsCombo1[ncols0], tabsCombo2[ncols0]], ignore_index=True)


        #fulltab = pd.concat([tabs1, tabs2, tabsCombo1, tabsCombo2])
        fulltab['study']='mpnstPDXMT'
        fulltab['source']='synapse'
        ##mutate the values create new columns
        ncols=['DOSE','GROWTH','study','source','improve_sample_id','Drug','time','time_unit']
        fulltab = fulltab[ncols]
        #there was a percentViability == "147.128*" - seems like it might be 147.1288 but double check
        #if any(fulltab['GROWTH'].str.contains("\*")):
        #     print("replacing asterisk")
        #     fulltab["GROWTH"] = fulltab["GROWTH"].replace(r'\*','', regex=True)
        #SG: script requires percentages, so converting growth values to percentages here
        fulltab['GROWTH'] = pd.to_numeric(fulltab['GROWTH'])#*100.00
        ##change file headers to DOSE/RESPONSE values needed by other script
        fulltab.to_csv('mpnst_combo_drug_response_forCurves.tsv',sep='\t')
        
        ##fit curve
        if not os.path.exists("fit_curve.py"):
            script='https://raw.githubusercontent.com/PNNL-CompBio/coderdata/refs/heads/main/coderbuild/utils/fit_curve.py'
            subprocess.run(['wget',script, '-O', 'fit_curve.py'])
            #subprocess.run(['python3','-m','pip','install','matplotlib']) # due to matplotlib error
        subprocess.run([
            'python', 'fit_curve.py',
            '--input', 'mpnst_combo_drug_response_forCurves.tsv',
            '--output', 'mpnstDrugComboOutput'
            ])
        #runpy.run_path("fit_curve.py") # replaces subprocess call
        #os.getcwd()
        otab = pd.read_csv('mpnstDrugComboOutput.0',sep='\t')
        
        # # also create redacted version for upload to MuSyC web app
        # fulltab3 = fulltab
        # fulltab['expt.date'] = date.today()
        # fulltab3 = fulltab3.rename(columns={"improve_sample_id": "sampleName"})
        # fulltab3['sampleName'] = fulltab3['sampleName'].str[-3:]
        # fulltab3['sample'] = fulltab3['sampleName'] + "_" + fulltab3['timePoint'].astype(str) + fulltab3['timePointUnit']
        # if need replicate, then use specimenID
        # fulltab3['drug1'] = fulltab3['drug1'].str[1:4]
        # fulltab3['drug2'] = fulltab3['drug2'].str[1:4]
        # ncols=['expt.date','drug1.conc','drug2.conc','effect','sample','drug1','drug2','drug1.units','drug2.units']
        # fulltab3 = fulltab3[ncols]
        # fulltab3.to_csv('combo_drug_response.csv', index=False) # MuSyC expects CSV file
        ##SG wrote results to files
        otab.to_csv('mpnstDrugComboOutput.txt',sep='\t')
        wtab = otab.pivot(index=['source','improve_sample_id','improve_drug_id','study','time','time_unit'],columns='dose_response_metric',values='dose_response_value')
    #    print(wtab)
        wtab.to_csv('mpnstDrugComboMatrix.tsv',sep='\t')
        wtab = pd.read_csv('mpnstDrugComboMatrix.tsv',sep='\t')

        single = wtab[['+' not in a for a in wtab.improve_drug_id]]

        combo = wtab[['+' in a for a in wtab.improve_drug_id]]
        combo[['drug1','drug2']]=combo['improve_drug_id'].str.split('+',expand=True)

        s1 = single[['improve_drug_id','improve_sample_id','time','auc']]
        s1 = s1.rename(columns={'improve_drug_id':'drug1','auc':'drug1_auc'})
        s2 = s1.rename(columns={'drug1':'drug2','drug1_auc':'drug2_auc'})


        combo = combo.merge(s1)
        combo = combo.merge(s2)
        combo[['sample','experiment']] = combo['improve_sample_id'].str.split('_',expand=True)
        combo = combo[['sample','experiment','drug1','drug2','time','time_unit','aac','auc','dss','fit_auc','fit_ec50','fit_ec50se','fit_einf','fit_hs','fit_r2','drug1_auc','drug2_auc']]
        combo.to_csv('mpnstDrugComboMatrix.tsv',sep='\t')
    ##now add in single-point drug measurements

    """
    store 2 agent combination results to synapse
    """
    # syn.store(sc.File('mpnst_combo_drug_response.csv', parentId=f'{result_folder_id}'))
    # syn.store(sc.File('mpnstDrugComboOutput.txt', parentId=f'{result_folder_id}'))
    # syn.store(sc.File('mpnstDrugComboMatrix.tsv', parentId=f'{result_folder_id}'))

    #syn.store(sc.File('mpnst_combo_drug_response_forCurves.tsv',parentId=result_folder_id))

    return None

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
