#!/usr/bin/env python
# coding: utf-8

from neuroHarmonize import harmonizationLearn
from neuroHarmonize import harmonizationApply, loadHarmonizationModel
import neuroHarmonize as nh
import pandas as pd
import numpy as np
import argparse


# Creazione di un parser per gli argomenti
parser = argparse.ArgumentParser(description='Subj_ID')

# Definisci un argomento
parser.add_argument('--ID', type=str, help='Subject ID')

# Analizza gli argomenti dalla linea di comando
args = parser.parse_args()

my_data  = pd.read_csv('IMAGING_RR_QC_LIGHT_VERSION_for_Suicide_Squad_per_armonizzazione_TOT_only_CLM.csv', index_col='Sub-id')
columns_name = my_data.columns
my_data.replace(['.', 'nan', 'Nan', ' '], np.nan, inplace = True)
my_data = my_data.astype(float)
covars = pd.read_csv('covariates_TOT_PROVA_only_CLM.csv', index_col='Sub-id')
new_subject = pd.read_csv(f'{args.ID}_Output_TOTALE_IMAGING_selected_features.csv', index_col=False)

# Definisci la lista delle colonne in comune
features_selected = ['Left-Inf-Lat-Vent','CC_Mid_Anterior','lh_caudalmiddlefrontal_thickness','rh_caudalmiddlefrontal_thickness','lh_hippocampal-fissure', 
                  'MD_Avg_Weight_mcp_avg16_syn_bbr']  

# Filtra i DataFrame per mantenere solo le colonne in comune
my_data_filtered = my_data[features_selected]
new_subject_filtered = new_subject[features_selected]

# Aggiungi le righe del secondo DataFrame al primo
df_for_harmonization = pd.concat([my_data_filtered, new_subject_filtered], ignore_index=True)

# Trasforma in NP array
df_for_harmonization = np.array(df_for_harmonization)

# aggiungo nuovo soggetto a db delle covariate
index_name = 'new_subj'
site_value = 'site_B'

# Aggiungo la nuova riga
covars.loc[index_name] = site_value

#Inizializza un DataFrame vuoto per i risultati armonizzati
harmonized_results = pd.DataFrame(index=range(df_for_harmonization.shape[0]), columns=range(df_for_harmonization.shape[1]))

#Loop attraverso le colonne 
for col_idx in range(df_for_harmonization.shape[1]):
    #Seleziona la colonna corrente 
    col_data = df_for_harmonization[:, col_idx]
    
    #Trova le righe senza NaN nella colonna corrente 
    valid_rows = ~np.isnan(col_data)
    
    #Filtra per le righe valide
    data_valid = col_data[valid_rows].reshape(-1, 1)
    covars_valid = covars[valid_rows]
    print(f"Column: {features_selected[col_idx]}: data_valid shape: {data_valid.shape}, covars_valid shape: {covars_valid.shape}")
    
    #Applica l'armonizzazione ai dati validi
    my_model_col_idx = loadHarmonizationModel(f'MY_MODEL_{features_selected[col_idx]}')
    data_harmonized_1_subj = harmonizationApply(data_valid, covars_valid, my_model_col_idx)
    
    #Aggiorna harmonized_results con i valori armonizzati
    harmonized_results.loc[valid_rows, col_idx] = data_harmonized_1_subj.flatten()

#Visualizza i risultati
#print(harmonized_results)
#harmonized_results.to_excel(f'harmonized_results_new_subj_PROVA.xlsx') 

# Estrai l'ultima riga 
harmonized_new_subject = harmonized_results.tail(1)
#Aggiungi colonne
harmonized_new_subject.columns = features_selected
harmonized_new_subject.to_excel('New_subject_harmonized.xlsx',index=False)
#print(harmonized_new_subject)







