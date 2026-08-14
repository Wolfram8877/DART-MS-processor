"""
DART-MS Core Processing Engine.
Handles the heavy lifting of parsing raw MS spectra, performing accurate mass matching
against a database, applying heuristic chemical scoring, and exporting the final report.
"""

import pandas as pd
import numpy as np
import math
import os
import re
import hashlib
from PIL import Image, ImageOps

from rdkit import Chem
from rdkit.Chem import Draw

# Import our custom chemistry module
from chemistry_utils import get_mol_props

def generate_safe_filename(cas, name):
    """
    Generates a secure and universally valid filename for a compound.
    Uses the CAS number if available. Otherwise, falls back to a deterministic MD5 hash of the name.
    
    Args:
        cas (str): The CAS registry number.
        name (str): The IUPAC or common name of the compound.
        
    Returns:
        str: A safe filename without extension.
    """
    cas_clean = str(cas).strip() if pd.notna(cas) else ""
    if cas_clean and cas_clean != "nan":
        return cas_clean.replace("/", "-")
        
    # Deterministic hash generation for compounds without CAS
    seed = str(name).strip()
    hash_object = hashlib.md5(seed.encode('utf-8'))
    short_hash = hash_object.hexdigest()[:10]
    return f"mol_{short_hash}"

def process_dart_ms(exp_file, db_file, sample_matrix, min_intensity, max_ppm, output_dir, max_candidates=5, contamination_targets=None, sort_by="m/z", generate_png=True, generate_mol=False):
    """
    Main pipeline for processing a single DART-MS experimental file against the database.
    """
    if contamination_targets is None:
        contamination_targets = []
        
    # Parse the input filename to generate matching output reports dynamically
    base_name = os.path.basename(exp_file)
    name_no_ext = os.path.splitext(base_name)[0]
    
    match = re.match(r"(?i)^ms\s+spectra\s+(.+)$", name_no_ext)
    if match:
        sample_id = match.group(1).strip()
        output_file = os.path.join(output_dir, f"DART_Results_{sample_id}.xlsx")
    else:
        output_file = os.path.join(output_dir, f"DART_Results_{name_no_ext}.xlsx")

    # ==========================================
    # STEP 1: DATABASE LOADING & PREPARATION
    # ==========================================
    print("--- Loading Database ---")
    
    db = pd.read_excel(db_file, header=None, engine='openpyxl')

    # Dynamically locate the actual header row to handle badly formatted Excel files
    header_row = 0
    for i in range(min(10, len(db))):
        row_vals = [str(x).strip().lower() for x in db.iloc[i].values]
        if 'monoisotopic_mass' in row_vals or 'exact_mass' in row_vals or 'canonical_smiles' in row_vals or 'smiles' in row_vals:
            header_row = i
            break

    db.columns = db.iloc[header_row].astype(str).str.strip()
    db = db.iloc[header_row + 1:].reset_index(drop=True)
    
    # Drop duplicated columns to prevent Pandas UserWarnings when converting to dictionary later
    db = db.loc[:, ~db.columns.duplicated()]

    mass_col_name_db = 'monoisotopic_mass' if 'monoisotopic_mass' in db.columns else 'exact_mass'
    
    possible_smiles = ['canonical_smiles', 'smiles', 'SMILES', 'Canonical_SMILES']
    smiles_col = next((c for c in db.columns if c in possible_smiles), None)
    if not smiles_col:
        raise ValueError(f"ERROR: SMILES column not found. Detected columns: {list(db.columns)}")

    name_col = next((c for c in db.columns if c in ['iupac_name', 'name', 'Name', 'Compound Name']), db.columns[0])
    cas_col = 'cas' if 'cas' in db.columns else None
    group_col = 'Groups' if 'Groups' in db.columns else None
    func_col = 'Harmonized_functions' if 'Harmonized_functions' in db.columns else None
    list_col = 'PlastChem_lists' if 'PlastChem_lists' in db.columns else None
    
    def clean_mass(val):
        try: return float(str(val).replace(',', '.'))
        except: return np.nan
        
    db['monoisotopic_mass'] = db[mass_col_name_db].apply(clean_mass)

    print("Filtering out multi-component mixtures and salts...")
    initial_len = len(db)
    mask_dot = db[smiles_col].astype(str).str.contains(r'\.', na=False, regex=True)
    mask_semi = db[name_col].astype(str).str.contains(r';', na=False, regex=True)
        
    db = db[~(mask_dot | mask_semi)].reset_index(drop=True)
    print(f"Removed {initial_len - len(db)} mixture/disconnected entries from the database.")

    print("Pre-calculating molecular properties (RDKit) to optimize speed...")
    mols, num_cs, num_cls, num_brs, num_sis, num_ss, num_ns = [], [], [], [], [], [], []
    has_g1s, fams, formulas = [], [], []
    
    for s in db[smiles_col]:
        m, c, cl, br, si, s_atom, n, g1, fam, form = get_mol_props(s)
        mols.append(m); num_cs.append(c); num_cls.append(cl); num_brs.append(br); num_sis.append(si); num_ss.append(s_atom); num_ns.append(n)
        has_g1s.append(g1); fams.append(fam); formulas.append(form)

    # Calculate exact mass column securely
    calc_mass_col = pd.Series(np.nan, index=db.index)
    if 'monoisotopic_mass' in db.columns:
        calc_mass_col = pd.to_numeric(db['monoisotopic_mass'].astype(str).str.replace(',', '.'), errors='coerce')
    if 'exact_mass' in db.columns:
        exact_clean = pd.to_numeric(db['exact_mass'].astype(str).str.replace(',', '.'), errors='coerce')
        calc_mass_col = calc_mass_col.fillna(exact_clean)

    # Bulk concatenation of new columns to prevent DataFrame fragmentation warnings
    new_cols = pd.DataFrame({
        'RDKit_Mol': mols, 'num_C': num_cs, 'num_Cl': num_cls, 'num_Br': num_brs,
        'num_Si': num_sis, 'num_S': num_ss, 'num_N': num_ns,
        'has_group1': has_g1s, 'dart_families': fams, 'formula_str': formulas,
        'calc_mass': calc_mass_col
    })
    db = pd.concat([db, new_cols], axis=1)

    # Chemical constants for exact mass physics
    MASS_PROTON = 1.007276
    MASS_NH4 = 18.033823
    MASS_C13_DIFF = 1.003355
    MASS_M1_B_DIFF = 0.9970   
    MASS_M2_DIFF = 1.99705     

    mols_dir = os.path.join(output_dir, "mols")
    images_dir = os.path.join(output_dir, "images")
    if generate_png: os.makedirs(images_dir, exist_ok=True)
    if generate_mol: os.makedirs(mols_dir, exist_ok=True)

    # ==========================================
    # STEP 2: EXPERIMENTAL SPECTRUM LOADING
    # ==========================================
    print("Loading Spectrum Data...")
    
    if exp_file.endswith('.xlsx') or exp_file.endswith('.xls'):
        exp = pd.read_excel(exp_file)
    else:
        exp = pd.read_csv(exp_file, on_bad_lines='skip', sep=None, engine='python')

    exp.columns = exp.columns.str.strip()
    
    # Safely identify m/z and Intensity columns regardless of slight naming variations
    mz_col_name = 'm/z' if 'm/z' in exp.columns else (exp.columns[1] if len(exp.columns) > 1 else exp.columns[0])
    int_col_name = 'I' if 'I' in exp.columns else ('Intensity' if 'Intensity' in exp.columns else (exp.columns[2] if len(exp.columns) > 2 else exp.columns[-1]))

    print("Cleaning numeric data...")
    exp[mz_col_name] = pd.to_numeric(exp[mz_col_name].astype(str).str.replace(',', '.'), errors='coerce')
    exp[int_col_name] = pd.to_numeric(exp[int_col_name].astype(str).str.replace(',', '.'), errors='coerce')

    # Apply global threshold filtering
    exp_filtered = exp[exp[int_col_name] >= min_intensity].reset_index(drop=True)
    
    # Apply user-defined sorting
    if sort_by == "Intensity":
        print("Sorting peaks by Intensity (Descending)...")
        exp_filtered = exp_filtered.sort_values(by=int_col_name, ascending=False).reset_index(drop=True)
    else:
        print("Sorting peaks by m/z (Ascending)...")
        exp_filtered = exp_filtered.sort_values(by=mz_col_name).reset_index(drop=True)

    # ==========================================
    # STEP 3: ISOTOPE DETECTION (M+1, M+2)
    # ==========================================
    identified_peaks = {} 
    parent_formulas = {}  

    print("Detecting and flagging auto-isotopes...")
    isotope_flags = {} 
    
    for idx_curr, row_curr in exp_filtered.iterrows():
        mz_curr = row_curr[mz_col_name]
        int_curr = row_curr[int_col_name]
        if pd.isna(mz_curr): continue
        
        # Check for Carbon-13 M+1 isotope
        target_parent_13c = mz_curr - MASS_C13_DIFF
        tol_13c = target_parent_13c * (max_ppm / 1000000)
        parents_13c = exp_filtered[
            (exp_filtered[mz_col_name] >= target_parent_13c - tol_13c) & 
            (exp_filtered[mz_col_name] <= target_parent_13c + tol_13c)
        ]
        
        parent_found = False
        
        if not parents_13c.empty:
            p_row = parents_13c.loc[parents_13c[int_col_name].idxmax()]
            parent_mz = p_row[mz_col_name]
            observed_ratio = int_curr / p_row[int_col_name]
            
            # Theoretical max limits based on rough maximum carbon counts for a given mass
            max_carbons = max(1, int(parent_mz / 13))
            max_theoretical_ratio = max_carbons * 0.011 * 1.5 
            
            if p_row[int_col_name] > int_curr and observed_ratio <= max_theoretical_ratio: 
                isotope_flags[idx_curr] = {
                    "parent_peak_id": parents_13c[int_col_name].idxmax() + 1,
                    "parent_mz": parent_mz,
                    "ratio": observed_ratio,
                    "type": "+1 Da (C13)"
                }
                parent_found = True
                continue 

        # Check for Silicon/Nitrogen M+1 isotope
        target_parent_m1b = mz_curr - MASS_M1_B_DIFF
        tol_m1b = target_parent_m1b * (max_ppm / 1000000)
        parents_m1b = exp_filtered[
            (exp_filtered[mz_col_name] >= target_parent_m1b - tol_m1b) & 
            (exp_filtered[mz_col_name] <= target_parent_m1b + tol_m1b)
        ]
        
        if not parent_found and not parents_m1b.empty:
            p_row = parents_m1b.loc[parents_m1b[int_col_name].idxmax()]
            observed_ratio = int_curr / p_row[int_col_name]
            if p_row[int_col_name] > int_curr and observed_ratio <= 0.80:
                isotope_flags[idx_curr] = {
                    "parent_peak_id": parents_m1b[int_col_name].idxmax() + 1,
                    "parent_mz": p_row[mz_col_name],
                    "ratio": observed_ratio,
                    "type": "+1 Da (Si/N)"
                }
                parent_found = True
                continue

        # Check for Cl/Br/Si/S M+2 isotope
        target_parent_m2 = mz_curr - MASS_M2_DIFF
        tol_m2 = target_parent_m2 * (max(10, max_ppm) / 1000000)
        parents_m2 = exp_filtered[
            (exp_filtered[mz_col_name] >= target_parent_m2 - tol_m2) & 
            (exp_filtered[mz_col_name] <= target_parent_m2 + tol_m2)
        ]
        
        if not parents_m2.empty:
            p_row = parents_m2.loc[parents_m2[int_col_name].idxmax()]
            ratio_m2 = int_curr / p_row[int_col_name]
            
            if 0.01 <= ratio_m2 <= 2.50: 
                isotope_flags[idx_curr] = {
                    "parent_peak_id": parents_m2[int_col_name].idxmax() + 1,
                    "parent_mz": p_row[mz_col_name],
                    "ratio": ratio_m2,
                    "type": "+2 Da (Cl/Br/Si/S)"
                }

    mass_col = 'calc_mass'
    results = []

    # ==========================================
    # STEP 4: DB SEARCH AND SCORING (HEURISTICS)
    # ==========================================
    print("Analyzing peaks and ranking Database vs Isotopic hypotheses...")
    
    # OPTIMIZATION: Convert Pandas DataFrame to a native Python dictionary array.
    # This prevents the massive overhead of Pandas `.iterrows()` and speeds up the search loop by 10x-50x.
    db_records = db.to_dict('records') 
    
    for idx, row in exp_filtered.iterrows():
        mz_exp = row[mz_col_name]
        intensity = row[int_col_name]

        if pd.isna(mz_exp): continue

        valid_matches = []

        # Evaluate if current peak is already identified as an isotope
        if idx in isotope_flags:
            parent_info = isotope_flags[idx]
            parent_id = parent_info['parent_peak_id']
            iso_type = parent_info['type']
            obs_ratio = parent_info['ratio']
            
            is_proven = False
            is_possible = False
            iso_score = 0
            iso_explanation = ""
            
            # If the parent peak was successfully matched to a database formula, cross-validate the isotope mathematically
            if identified_peaks.get(parent_id, False) == True:
                p_form = parent_formulas.get(parent_id, {'C': 0, 'Cl': 0, 'Br': 0, 'Si': 0, 'S': 0, 'N': 0})
                
                if "2 Da" in iso_type:
                    theo_m2 = p_form['Cl'] * 31.98 + p_form['Br'] * 97.28 + p_form['Si'] * 3.36 + p_form['S'] * 4.47
                    if theo_m2 > 0 and abs(obs_ratio*100 - theo_m2) / theo_m2 <= 0.50:
                        is_proven = True
                        terms = []
                        if p_form['Cl'] > 0: terms.append(f"Cl({p_form['Cl']})")
                        if p_form['Br'] > 0: terms.append(f"Br({p_form['Br']})")
                        if p_form['Si'] > 0: terms.append(f"Si({p_form['Si']})")
                        if p_form['S'] > 0: terms.append(f"S({p_form['S']})")
                        formula_str = "+".join(terms)
                        iso_explanation = f"Validated M+2 Contribution -> Theory: [{formula_str}] = {round(theo_m2,1)}% | Exp: {round(obs_ratio*100,1)}%"
                        
                elif "1 Da" in iso_type:
                    theo_m1 = p_form['C'] * 1.08 + p_form['Si'] * 5.07 + p_form['S'] * 0.80 + p_form['N'] * 0.37
                    if theo_m1 > 0 and abs(obs_ratio*100 - theo_m1) / theo_m1 <= 0.50:
                        is_proven = True
                        terms = []
                        if p_form['C'] > 0: terms.append(f"C({p_form['C']})")
                        if p_form['Si'] > 0: terms.append(f"Si({p_form['Si']})")
                        if p_form['S'] > 0: terms.append(f"S({p_form['S']})")
                        if p_form['N'] > 0: terms.append(f"N({p_form['N']})")
                        formula_str = "+".join(terms)
                        iso_explanation = f"Validated M+1 Contribution -> Theory: [{formula_str}] = {round(theo_m1,1)}% | Exp: {round(obs_ratio*100,1)}%"
                
                if is_proven:
                    iso_score = 150 
            else:
                is_possible = True
                iso_score = 45 
                iso_explanation = f"Potential isotope of Unknown Peak {parent_id} (Exp Ratio: {round(obs_ratio*100,1)}%)"
                
            if is_proven or is_possible:
                valid_matches.append({
                    "is_special_isotope": True,
                    "score": iso_score,
                    "compound_name": f"Isotope ({iso_type}) of Peak {parent_id}",
                    "adduct": "Isotopic adduct",
                    "ppm": "",
                    "groups": "[ISOTOPE DETECTED]",
                    "calc_details": iso_explanation,
                    "db_row": None,
                    "matrix_status": "-",
                    "contamination": "/"
                })

        matches = []
        for db_row in db_records:
            neutral_mass = db_row[mass_col]
            if pd.isna(neutral_mass): continue

            mz_H = neutral_mass + MASS_PROTON
            mz_NH4 = neutral_mass + MASS_NH4

            # Calculate precise PPM errors
            ppm_H = (abs(mz_exp - mz_H) / mz_H) * 1000000
            ppm_NH4 = (abs(mz_exp - mz_NH4) / mz_NH4) * 1000000

            if ppm_H <= max_ppm: 
                matches.append({"db_row": db_row, "adduct": "[M+H]+", "ppm": ppm_H, "neutral_mass": neutral_mass})
            if ppm_NH4 <= max_ppm: 
                matches.append({"db_row": db_row, "adduct": "[M+NH4]+", "ppm": ppm_NH4, "neutral_mass": neutral_mass})
        
        for m in matches:
            db_row = m["db_row"]
            ppm_val = m["ppm"]
            
            num_C = db_row['num_C']
            num_Cl = db_row['num_Cl']
            num_Br = db_row['num_Br']
            num_Si = db_row['num_Si']
            num_S = db_row['num_S']
            num_N = db_row['num_N']
            has_group1 = db_row['has_group1']
            dart_families_txt = db_row['dart_families']

            # Gaussian decay score based on PPM error
            sigma = max_ppm / 3.0
            base_score = 100 * math.exp(-0.5 * (ppm_val / sigma)**2)

            # M+1 Penalty Logic
            m1_penalty = 0
            m1_details = ""
            theo_m1_ratio = num_C * 0.0108 + num_Si * 0.0507 + num_S * 0.0080 + num_N * 0.0037
            expected_m1_intensity = intensity * theo_m1_ratio
            
            if expected_m1_intensity >= min_intensity: 
                target_m1_mz = mz_exp + MASS_C13_DIFF
                tol_m1 = target_m1_mz * (max(10, max_ppm) / 1000000) 
                m1_peaks = exp[(exp[mz_col_name] >= target_m1_mz - tol_m1) & (exp[mz_col_name] <= target_m1_mz + tol_m1)]
                
                if m1_peaks.empty:
                    m1_penalty = -30
                    m1_details = f" | Penalty: -30 (Missing M+1, expected Int: {int(expected_m1_intensity)})"
                else:
                    m1_max_int = float(m1_peaks[int_col_name].max())
                    obs_m1_ratio = m1_max_int / intensity
                    
                    if obs_m1_ratio < (theo_m1_ratio * 0.25):
                        m1_penalty = -30
                        m1_details = f" | Penalty: -30 (M+1 too small: {round(obs_m1_ratio*100,1)}% vs exp {round(theo_m1_ratio*100,1)}%)"
                    else:
                        m1_details = f" | M+1 Verified" 

            # Halogen M+2 Bonus and Penalty Logic
            halogen_bonus = 0
            halogen_penalty = 0
            halogen_details = ""
            
            if num_Cl > 0 or num_Br > 0 or num_Si > 0 or num_S > 0:
                theo_m2 = num_Cl * 31.98 + num_Br * 97.28 + num_Si * 3.36 + num_S * 4.47
                expected_m2_intensity = intensity * (theo_m2 / 100.0)
                
                if expected_m2_intensity >= min_intensity:
                    target_m2_mz = mz_exp + MASS_M2_DIFF
                    tol_m2 = target_m2_mz * (max(10, max_ppm) / 1000000)
                    m2_peaks = exp[(exp[mz_col_name] >= target_m2_mz - tol_m2) & (exp[mz_col_name] <= target_m2_mz + tol_m2)]
                    
                    if m2_peaks.empty: 
                        halogen_penalty = -30
                        halogen_details = f" | Penalty: -30 (Missing M+2, expected Int: {int(expected_m2_intensity)})"
                    else:
                        m2_int = float(m2_peaks[int_col_name].max())
                        observed_m2_ratio = (m2_int / intensity) * 100
                        
                        m2_error = abs(observed_m2_ratio - theo_m2) / theo_m2
                        if m2_error <= 0.50: 
                            terms = []
                            if num_Cl > 0: terms.append(f"Cl({num_Cl})")
                            if num_Br > 0: terms.append(f"Br({num_Br})")
                            if num_Si > 0: terms.append(f"Si({num_Si})")
                            if num_S > 0: terms.append(f"S({num_S})")
                            formula_str = "+".join(terms)
                            
                            halogen_bonus = 30
                            halogen_details = f" | M+2 Verified ({formula_str}): +30 ({round(observed_m2_ratio,1)}% vs exp {round(theo_m2,1)}%)"
                        elif observed_m2_ratio < (theo_m2 * 0.25):
                            halogen_penalty = -30
                            halogen_details = f" | Penalty: -30 (M+2 Ratio Mismatch: {round(observed_m2_ratio,1)}% vs exp {round(theo_m2,1)}%)"

            matrix_bonus = 0
            matrix_details = ""
            matrix_status = "N/A"
            if sample_matrix in db.columns and sample_matrix != "Unknown":
                try:
                    val_raw = str(db_row[sample_matrix]).replace(',', '.')
                    val_float = float(val_raw)
                    if val_float == 2.0: 
                        matrix_bonus = 30
                        matrix_details = f" | Matrix (2.0): +30"
                        matrix_status = f"{sample_matrix} (Released)"
                    elif val_float == 1.0: 
                        matrix_bonus = 30
                        matrix_details = f" | Matrix (1.0): +30"
                        matrix_status = f"{sample_matrix} (Present)"
                    elif val_float == 0.5:
                        matrix_bonus = 20
                        matrix_details = f" | Matrix (0.5): +20"
                        matrix_status = f"{sample_matrix} (Used)"
                    elif val_float == 0.25:
                        matrix_bonus = 10
                        matrix_details = f" | Matrix (0.25): +10"
                        matrix_status = f"{sample_matrix} (Inconclusive)"
                    elif val_float == 0:
                        matrix_status = f"{sample_matrix} (Not Detected)"
                except (ValueError, TypeError):
                    matrix_status = "Unknown Data"
            else:
                matrix_status = "N/A (Matrix Unknown)"
                
            m["matrix_status"] = matrix_status

            contaminants = []
            if contamination_targets:
                matrix_groups = {
                    "PE": ["PE_combined", "PE", "HDPE", "LDPE"],
                    "PS": ["PS_combined", "PS", "EPS", "HIPS"],
                    "PP": ["PP"], "PET": ["PET"], "PVC": ["PVC"], 
                    "PUR": ["PUR"], "PA": ["PA"], "EVA": ["EVA"],
                    "All Known Plastics": ["PE_combined", "PE", "HDPE", "LDPE", "PS_combined", "PS", "EPS", "HIPS", "PP", "PET", "PVC", "PUR", "PA", "EVA"]
                }
                
                target_cols = []
                for target in contamination_targets:
                    target_cols.extend(matrix_groups.get(target, []))
                    
                target_cols = list(set(target_cols)) # Remove duplicates safely
                target_cols = [col for col in target_cols if col != sample_matrix]
                
                for other_mat in target_cols:
                    if other_mat in db.columns:
                        try:
                            val_raw = str(db_row[other_mat]).replace(',', '.')
                            val_float = float(val_raw)
                            if val_float == 2.0:
                                contaminants.append(f"{other_mat} (Released)")
                            elif val_float == 1.0:
                                contaminants.append(f"{other_mat} (Present)")
                            elif val_float == 0.5:
                                contaminants.append(f"{other_mat} (Used)")
                            elif val_float == 0.25:
                                contaminants.append(f"{other_mat} (Inconclusive)")
                        except (ValueError, TypeError):
                            pass
                            
            if contaminants:
                m["contamination"] = f"Contamination possible from: {', '.join(sorted(contaminants))}"
            else:
                m["contamination"] = "/"

            if m["adduct"] == "[M+NH4]+":
                if not has_group1:
                    chem_bonus = -50
                    chem_text = "NH4: -50"
                else:
                    chem_bonus = 0
                    chem_text = ""
            else:
                chem_bonus = 0
                chem_text = ""

            # Co-adduct detection bonus
            co_adduct_bonus = 0
            co_adduct_details = ""
            
            if has_group1:
                target_other_mass = m["neutral_mass"] + MASS_NH4 if m["adduct"] == "[M+H]+" else m["neutral_mass"] + MASS_PROTON
                other_adduct_name = "[M+NH4]+" if m["adduct"] == "[M+H]+" else "[M+H]+"
                tol_mz = target_other_mass * (max_ppm / 1000000)

                co_peaks = exp[(exp[mz_col_name] >= target_other_mass - tol_mz) & (exp[mz_col_name] <= target_other_mass + tol_mz)]
                if not co_peaks.empty: 
                    co_peak_max_int = float(co_peaks[int_col_name].max())
                    if co_peak_max_int >= (0.05 * intensity):
                        co_adduct_bonus = 30
                    else:
                        co_adduct_bonus = 10
                        
                    co_adduct_details = f" | Co-adduct {other_adduct_name}: +{co_adduct_bonus}"

            raw_total = base_score + m1_penalty + halogen_bonus + halogen_penalty + matrix_bonus + co_adduct_bonus + chem_bonus
            m["score"] = round(max(0, raw_total), 1) 
            m["dart_families"] = dart_families_txt
            
            details = f"ppm: {round(base_score,1)}"
            details += matrix_details
            if chem_text != "": details += f" | {chem_text}"
            details += co_adduct_details
            details += m1_details
            details += halogen_details
            
            m["calc_details"] = details + f" => Total: {m['score']}"
            m["is_special_isotope"] = False
            
            valid_matches.append(m)

        # ==========================================
        # STEP 5: FINAL RANKING & EXPORT PREPARATION
        # ==========================================
        if not valid_matches:
            results.append({
                "Peak ID": idx + 1, "m/z": mz_exp, "Intensity": intensity, "Identified Compound": "Unknown",
                "Ionization & Chemistry": "", "Error (ppm)": "", "Score": "", "Score comment": f"No valid match found",
                "Matrix Status": "", "Function": "", "Possible contamination": "", "PlastChem List": "", "CAS Number": "",
                "Structure_Image_Path": "", "Chemical Formula": "", "Chemical Explanation": "", "Is_Isotope_Row": "No"
            })
            identified_peaks[idx + 1] = False
        else:
            valid_matches = sorted(valid_matches, key=lambda x: x["score"], reverse=True)
            
            top_db_match = next((match for match in valid_matches if not match.get("is_special_isotope")), None)
            if top_db_match:
                identified_peaks[idx + 1] = True
                dr = top_db_match["db_row"]
                parent_formulas[idx + 1] = {
                    'C': dr['num_C'], 'Cl': dr['num_Cl'], 'Br': dr['num_Br'],
                    'Si': dr['num_Si'], 'S': dr['num_S'], 'N': dr['num_N']
                }
            else:
                identified_peaks[idx + 1] = False

            for i, m in enumerate(valid_matches[:max_candidates]):
                if m.get("is_special_isotope"):
                    results.append({
                        "Peak ID": idx + 1 if i == 0 else "", "m/z": mz_exp if i == 0 else "",
                        "Intensity": intensity if i == 0 else "", 
                        "Identified Compound": m["compound_name"],
                        "Ionization & Chemistry": f"Isotopic adduct | {m['groups']}", "Error (ppm)": m["ppm"], 
                        "Score": m["score"], 
                        "Score comment": f"Candidate {i+1}/{min(len(valid_matches), max_candidates)}.",
                        "Matrix Status": m["matrix_status"], "Function": "", 
                        "Possible contamination": m.get("contamination", "/"), 
                        "PlastChem List": "", "CAS Number": "",
                        "Structure_Image_Path": "", "Chemical Formula": "", 
                        "Chemical Explanation": m["calc_details"], "Is_Isotope_Row": "Yes"
                    })
                else:
                    db_row = m["db_row"]
                    compound_name = db_row.get(name_col, "N/A")
                    if pd.isna(compound_name): compound_name = "N/A"
                    cas_number = cas_col and db_row.get(cas_col, "") or ""
                    
                    formula_str = db_row['formula_str']
                    mol = db_row['RDKit_Mol']
                    img_path = ""

                    if mol is not None:
                        try:
                            # Use deterministic hash function for missing CAS to ensure Cross-Screening compatibility
                            safe_name = generate_safe_filename(cas_number, compound_name)
                            img_path = os.path.join(images_dir, f"{safe_name}.png")
                            
                            if generate_png and not os.path.exists(img_path): 
                                options = Draw.MolDrawOptions()
                                options.bondLineWidth = 3.0
                                options.minFontSize = 16
                                options.clearBackground = True
                                Draw.MolToFile(mol, img_path, size=(290, 290), options=options)
                                img = Image.open(img_path)
                                img_with_border = ImageOps.expand(img, border=5, fill='#bdc3c7')
                                img_with_border.save(img_path)
                                
                            if generate_mol:
                                mol_path = os.path.join(mols_dir, f"{safe_name}.mol")
                                if not os.path.exists(mol_path):
                                    Chem.MolToMolFile(mol, mol_path)
                        except: pass

                    comment = f"Candidate {i+1}/{min(len(valid_matches), max_candidates)}."
                    if i == 0: comment += " Top score DB match."
                    
                    groupes_db = group_col and str(db_row.get(group_col, "")) or ""
                    dart_families = m.get("dart_families", "")
                    final_group = f"{groupes_db} [DART: {dart_families}]" if dart_families else groupes_db
                    
                    ionization_and_chem = f"{m['adduct']} | {final_group}"

                    func_val = db_row.get(func_col, "") if func_col else ""
                    list_val = db_row.get(list_col, "") if list_col else ""

                    results.append({
                        "Peak ID": idx + 1 if i == 0 else "", "m/z": mz_exp if i == 0 else "",
                        "Intensity": intensity if i == 0 else "", "Identified Compound": compound_name,
                        "Ionization & Chemistry": ionization_and_chem, "Error (ppm)": round(m["ppm"], 2) if m["ppm"] != "" else "",
                        "Score": m["score"], "Score comment": comment, "Matrix Status": m["matrix_status"],
                        "Function": func_val, "Possible contamination": m["contamination"], "PlastChem List": list_val,
                        "CAS Number": cas_number if not pd.isna(cas_number) else "",
                        "Structure_Image_Path": img_path, 
                        "Chemical Formula": formula_str,
                        "Chemical Explanation": m.get("calc_details", ""),
                        "Is_Isotope_Row": "No"
                    })

    print("Creating Excel report...")
    df_res = pd.DataFrame(results)
    df_res = df_res.replace([np.inf, -np.inf], np.nan).fillna("")

    try:
        writer = pd.ExcelWriter(output_file, engine='xlsxwriter')
    except ImportError:
        print("ERROR: xlsxwriter module is missing. Please run: pip install xlsxwriter")
        return

    df_excel = df_res.drop(['Structure_Image_Path', 'Is_Isotope_Row'], axis=1, errors='ignore')
    
    # Insert visual column dynamically before CAS Number
    cas_idx = df_excel.columns.get_loc('CAS Number') if 'CAS Number' in df_excel.columns else len(df_excel.columns) - 1
    df_excel.insert(cas_idx, "Structure (Visualization)", "")
    
    df_excel.to_excel(writer, sheet_name='DART Results', index=False, header=False, startrow=1)

    workbook  = writer.book
    worksheet = writer.sheets['DART Results']
    
    # Establish unified professional font size (12pt) applied to all formats
    header_format = workbook.add_format({'bold': True, 'text_wrap': True, 'valign': 'vcenter', 'align': 'center', 'bg_color': '#4F81BD', 'font_color': 'white', 'border': 1, 'font_size': 12})
    cell_format = workbook.add_format({'text_wrap': True, 'valign': 'vcenter', 'align': 'center', 'border': 1, 'border_color': '#D3D3D3', 'font_size': 12})
    peak_id_format = workbook.add_format({'bold': True, 'valign': 'vcenter', 'align': 'center', 'border': 1, 'border_color': '#D3D3D3', 'font_size': 14})
    num_format = workbook.add_format({'num_format': '0.00', 'valign': 'vcenter', 'align': 'center', 'border': 1, 'border_color': '#D3D3D3', 'font_size': 12})
    score_format = workbook.add_format({'bold': True, 'font_color': '#006100', 'bg_color': '#C6EFCE', 'valign': 'vcenter', 'align': 'center', 'border': 1, 'font_size': 12})
    score_bronze_format = workbook.add_format({'bold': True, 'font_color': '#5C4033', 'bg_color': '#E5D3B3', 'valign': 'vcenter', 'align': 'center', 'border': 1, 'font_size': 12})

    for col_num, value in enumerate(df_excel.columns.values): worksheet.write(0, col_num, value, header_format)
        
    # Dynamic column widths mapping
    for col_num, col_name in enumerate(df_excel.columns):
        if col_name == "Chemical Explanation":
            worksheet.set_column(col_num, col_num, 80, cell_format)
        elif col_name in ["Structure (Visualization)", "Ionization & Chemistry", "Identified Compound"]:
            worksheet.set_column(col_num, col_num, 35, cell_format)
        elif col_name in ["Score comment", "Possible contamination", "Function", "PlastChem List"]:
            worksheet.set_column(col_num, col_num, 30, cell_format)
        elif col_name in ["m/z", "Intensity", "Error (ppm)", "Peak ID", "Score"]:
            worksheet.set_column(col_num, col_num, 12, cell_format)
        else:
            worksheet.set_column(col_num, col_num, 20, cell_format)

    worksheet.set_default_row(150)
    worksheet.set_row(0, 30)

    for row_idx, row in df_res.iterrows():
        excel_row = row_idx + 1
        peak_val = row.get("Peak ID")
        is_isotope_row = row.get("Is_Isotope_Row") == "Yes" or "[ISOTOPE DETECTED]" in str(row.get("Ionization & Chemistry", ""))
        
        for col_idx, col_name in enumerate(df_excel.columns):
            val_to_write = row.get(col_name, "") if col_name != "Structure (Visualization)" else ""
            
            if col_idx == 0 and peak_val != "":
                worksheet.write(excel_row, 0, peak_val, peak_id_format)
            elif col_name == "Score": 
                if is_isotope_row:
                    worksheet.write(excel_row, col_idx, row.get("Score", ""), score_bronze_format)
                else:
                    worksheet.write(excel_row, col_idx, row.get("Score", ""), score_format)
            elif col_name == "Error (ppm)":
                worksheet.write(excel_row, col_idx, row.get("Error (ppm)", ""), num_format)
            else:
                worksheet.write(excel_row, col_idx, val_to_write, cell_format)

        img_p = row.get("Structure_Image_Path")
        if isinstance(img_p, str) and img_p != "" and os.path.exists(img_p):
            visu_idx = df_excel.columns.get_loc('Structure (Visualization)')
            worksheet.insert_image(excel_row, visu_idx, img_p, {'x_scale': 0.6, 'y_scale': 0.6, 'x_offset': 35, 'y_offset': 10, 'positioning': 1})
            
    writer.close()
    print(f"\nSUCCESS! Results generated: {output_file}")