import pandas as pd
import pingouin as pg

def pw_corr(data_path="data/TMA36_project/Radiomics/processed/rad_healthmyne.csv", 
            cde_path="data/TMA36_project/CDE/CDE_TMA36_2020FEB25_SA_MF.csv"):
    rad_hm = pd.read_csv(data_path, index_col=0)
    cde = pd.read_csv(cde_path, index_col=1)
    cde_sila = pd.DataFrame(cde['SILA'])
    rad_hm_sila = pd.merge(rad_hm, cde_sila, how='left', left_index=True, right_index=True)
    pairwise = rad_hm_sila.pairwise_corr(method='spearman',padjust='holm', columns=['SILA'])
    pairwise_sig = pairwise[pairwise['p-corr']<0.05]

    return pairwise_sig