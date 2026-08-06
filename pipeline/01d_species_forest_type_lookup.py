"""01d_species_forest_type_lookup.py — map the inventory species onto forest-type groups.

The imputation matches plots on their field forest type, which means the species recorded at a plot
have to be expressed in the same ~15 groups the forest-type map uses. This builds that crosswalk.

Assignment is by trait and taxonomy rather than by a stored species list, so a species new to the
inventory is placed by rule instead of falling through. The traits come from the curated species
table: needle-leaved against broad-leaved, and evergreen against deciduous.

  conifers   Pinus densiflora 14 · P. thunbergii 11 · P. rigida, taeda, banksiana 13 ·
             P. koraiensis and other nut pines 15 · Larix 12 · Chamaecyparis 16 · other conifer 10
             The two regional races of P. densiflora are merged into 14.
  broadleaf  Castanea 30 · Quercus mongolica 23, variabilis 21, acutissima 22, evergreen oaks 40,
             other deciduous oaks 20 · all other broadleaf: evergreen 40, deciduous 50

Growth form (tree, small tree, shrub, vine) is carried through as a reference column. It does not
affect the mapping, because plot dominance is computed from canopy trees only and a shrub cannot
be dominant.

In : data/raw/nfi_curated/species_NFI.xlsx
     optionally a second table to cross-check against, via SPECIES_VG_CSV
Out: data/species_forest_type_lookup.csv        (Name_Korean, Genus, Species, FT)
     data/interim/species_ft_audit.csv          every species with the rule that placed it
     data/interim/species_ft_crosscheck_disagree.csv   only when a cross-check table is supplied
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
XLSX = os.path.join(ROOT, 'data', 'raw', 'nfi_curated', 'species_NFI.xlsx')
# Optional second table to cross-check against. Point SPECIES_VG_CSV at one to enable the comparison;
# without it the lookup is built from species_NFI.xlsx alone, which is the authoritative source.
LEG = os.environ.get('SPECIES_VG_CSV', '')

FT_NAME = {10: '기타침엽수', 11: '곰솔', 12: '낙엽송', 13: '리기다소나무', 14: '소나무', 15: '잣나무',
           16: '편백나무', 20: '기타참나무', 21: '굴참나무', 22: '상수리나무', 23: '신갈나무', 30: '밤나무',
           40: '상록활엽수', 50: '기타활엽수', 60: '침활혼효림'}
EVER = {'상록수', '반상록수'}
PINE15 = {'koraiensis', 'parviflora', 'pumila', 'strobus', 'densa'}
PINE13 = {'rigida', 'taeda', 'banksiana', 'rigitaeda'}
PINE14_KO = {'소나무', '강원지방소나무', '중부지방소나무'}   # the two regional races and the species itself
# The evergreen oaks, fixed by Latin species because the trait table records Quercus salicina as
# deciduous.
QUERCUS_EVERGREEN = {'salicina', 'acuta', 'glauca', 'myrsinifolia', 'gilva', 'phillyraeoides'}

def resolve(ko, genus, sp, grp, phen):
    ko, genus, sp = str(ko).strip(), str(genus).strip(), str(sp).strip()
    grp, phen = str(grp).strip(), str(phen).strip()
    if genus == 'Ginkgo':                                      return 10, 'ginkgo10'   # in group 10 by its definition, though the trait table calls it broad-leaved
    if grp == '침엽수':
        if genus == 'Pinus':
            if sp.startswith('densiflora') or ko in PINE14_KO: return 14, 'pine_densiflora'
            if sp == 'thunbergii':                             return 11, 'pine'
            if sp in PINE13:                                   return 13, 'pine'
            if sp in PINE15:                                   return 15, 'pine'
            return 10, 'pine_other10'
        if genus == 'Larix':                                   return 12, 'conifer'
        if genus == 'Chamaecyparis':                           return 16, 'conifer'
        return 10, 'conifer_other10'
    # broadleaf
    # Borderline semi-evergreens go to 40. Ligustrum ovalifolium arrives there already, through the
    # 반상록수 trait; Lindera glauca and Sageretia thea are recorded deciduous but hold their leaves
    # through the southern winter, so they are placed by name.
    if ko in ('감태나무', '상동나무'):                          return 40, 'semi_evergreen_user'
    if genus == 'Castanea':                                    return 30, 'castanea'
    if genus == 'Quercus':
        if sp == 'mongolica':                                  return 23, 'oak'
        if sp == 'variabilis':                                 return 21, 'oak'
        if sp == 'acutissima':                                 return 22, 'oak'
        if sp in QUERCUS_EVERGREEN or phen in EVER:            return 40, 'oak_evergreen'
        return 20, 'oak_other20'
    if phen in EVER:                                           return 40, 'evergreen_broadleaf'
    return 50, 'deciduous_broadleaf'

x = pd.read_excel(XLSX)
if LEG and os.path.exists(LEG):
    leg = pd.read_csv(LEG, encoding='utf-8'); leg.columns = ['Species', 'VG']
    leg_map = dict(zip(leg['Species'].astype(str).str.strip(), leg['VG'].astype(int)))
else:
    leg_map = {}                       # cross-check skipped; the trait-based resolve() is authoritative
    print('SPECIES_VG_CSV not set or absent -> cross-check skipped')

rows = []
for r in x.itertuples():
    ft, src = resolve(r.Name_Korean, r.Genus, r.Species, r.Species_group, r.Leaf_phenology)
    lv = leg_map.get(str(r.Name_Korean).strip())
    agree = (lv == ft) if lv is not None else None
    rows.append(dict(Name_Korean=r.Name_Korean, Genus=r.Genus, Species=r.Species, FT=ft, FT_name=FT_NAME[ft],
                     Leaf_phenology=r.Leaf_phenology, Tree_type=r.Tree_type, source=src,
                     ref_VG=lv, agree=agree, review=(agree is False) or (str(r.Leaf_phenology).strip() == '반상록수')))
df = pd.DataFrame(rows).sort_values(['FT', 'Genus', 'Name_Korean'])

AUDIT = os.path.join(ROOT, 'data', 'interim')
os.makedirs(AUDIT, exist_ok=True)
df[['Name_Korean', 'Genus', 'Species', 'FT']].to_csv(os.path.join(ROOT, 'data', 'species_forest_type_lookup.csv'), index=False, encoding='utf-8')
df.to_csv(os.path.join(AUDIT, 'species_ft_audit.csv'), index=False, encoding='utf-8')
if leg_map:
    dis = df[df.agree == False]
    dis[['Name_Korean', 'Genus', 'Species', 'Leaf_phenology', 'Tree_type', 'FT', 'FT_name', 'ref_VG']].to_csv(os.path.join(AUDIT, 'species_ft_crosscheck_disagree.csv'), index=False, encoding='utf-8')

print(f'{len(df)} species -> data/species_forest_type_lookup.csv  (trait-based; the second table only cross-checks)')
print('FT distribution:', {f'{k[0]} {k[1]}': v for k, v in df.groupby(['FT', 'FT_name']).size().to_dict().items()})
print('source:', df['source'].value_counts().to_dict())
if leg_map:
    print(f'cross-checked: {int(df.ref_VG.notna().sum())}/{len(df)} | AGREE {int((df.agree==True).sum())} | DISAGREE {int((df.agree==False).sum())} -> species_ft_crosscheck_disagree.csv')
print(f'evergreen broadleaf (FT40): {int((df.FT==40).sum())} | 반상록수 flagged for review')
print('pine check:', df[df.Genus=='Pinus'][['Name_Korean','Species','FT']].to_dict('records'))
