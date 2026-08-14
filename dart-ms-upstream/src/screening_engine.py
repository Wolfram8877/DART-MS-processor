"""
Cross-Screening Engine.
Aggregates and compares DART-MS result files to highlight compound presence
across multiple sample environments based on highly precise m/z matching.
"""

import pandas as pd
import numpy as np
import os
import hashlib

def generate_safe_filename(cas, name):
    """
    Generates a secure and universally valid filename for a compound.
    Matches the exact hashing logic found in the processing_engine to correctly 
    re-link structural images during the screening report generation.
    
    Args:
        cas (str): The CAS registry number.
        name (str): The IUPAC or common name.
        
    Returns:
        str: A safe filename without extension.
    """
    cas_clean = str(cas).strip() if pd.notna(cas) else ""
    if cas_clean and cas_clean != "nan":
        return cas_clean.replace("/", "-")
        
    seed = str(name).strip()
    hash_object = hashlib.md5(seed.encode('utf-8'))
    short_hash = hash_object.hexdigest()[:10]
    return f"mol_{short_hash}"

def run_cross_screening(file_paths, output_dir, mz_tol, multi_polymer_mode):
    """
    Core engine aggregating data from multiple DART-MS runs.
    
    Args:
        file_paths (list): List of paths to the Excel result files.
        output_dir (str): Directory where the compiled Master Report will be saved.
        mz_tol (float): m/z deviation tolerance for grouping identical peaks across runs.
        multi_polymer_mode (bool): If True, consolidates primary matrix and contamination lists.
        
    Returns:
        tuple: (success_boolean, output_path_or_error_message)
    """
    if not file_paths:
        return False, "No files provided."
        
    print("\n========================================")
    print("STARTING CROSS-SCREENING ANALYSIS")
    print("========================================")
    print(f"Grouping m/z tolerance: +/- {mz_tol} Da")
    if multi_polymer_mode:
        print("Mode: Multi-Polymer (Consolidating Matrix info)")
    
    master_compounds = {}
    sample_names = []
    
    # ==========================================
    # STEP 1: READ AND CONSOLIDATE ALL FILES
    # ==========================================
    for filepath in file_paths:
        base_name = os.path.basename(filepath)
        sample_name = base_name.replace("DART_Results_", "").replace(".xlsx", "")
        sample_names.append(sample_name)
        
        print(f"Reading: {sample_name}...", flush=True)
        try:
            df = pd.read_excel(filepath)
            
            if 'Identified Compound' not in df.columns or 'm/z' not in df.columns:
                print(f"  -> Warning: Skipped {sample_name}, missing required columns.", flush=True)
                continue
                
            # Forward fill structural variables derived from merged cells
            df['m/z'] = pd.to_numeric(df['m/z'], errors='coerce').ffill()
            if 'Intensity' in df.columns:
                df['Intensity'] = pd.to_numeric(df['Intensity'], errors='coerce').ffill()
                
            # Filter non-applicable data
            df_valid = df[~df['Identified Compound'].isin(['Unknown', '', np.nan])]
            df_valid = df_valid[~df_valid['Identified Compound'].astype(str).str.contains("Isotope", na=False, case=False)]
            
            if 'Score' in df_valid.columns:
                df_valid['Score'] = pd.to_numeric(df_valid['Score'], errors='coerce').fillna(0)
            else:
                df_valid['Score'] = 0
                
            # Keep only the top scoring instance per unique compound per file
            df_valid = df_valid.sort_values(by=['Score'], ascending=[False])
            best_hits = df_valid.drop_duplicates(subset=['Identified Compound'], keep='first')

            for _, row in best_hits.iterrows():
                comp_name = row['Identified Compound']
                score = float(row.get('Score', 0))
                intensity = float(row.get('Intensity', 0))
                func = str(row.get('Function', '')).strip()
                
                adduct_raw = str(row.get('Ionization & Chemistry', '')).split('|')[0].strip()
                
                try:
                    ppm_raw = row.get('Error (ppm)', "")
                    if pd.isna(ppm_raw) or str(ppm_raw).strip() == "":
                        ppm_val = ""
                    else:
                        ppm_val = round(float(ppm_raw), 2)
                except:
                    ppm_val = ""
                
                mat_status = str(row.get('Matrix Status', '')).strip()
                contam = str(row.get('Possible contamination', '')).replace('Contamination possible from:', '').strip()
                
                # Register compound globally if encountered for the first time
                if comp_name not in master_compounds:
                    cas_clean = str(row.get('CAS Number', '')).strip()
                    
                    # Identical hash resolution strategy across both processor and screener
                    safe_name = generate_safe_filename(cas_clean, comp_name)
                        
                    # Graceful image location fallback (Source dir vs Output dir)
                    img_path_src = os.path.join(os.path.dirname(filepath), "images", f"{safe_name}.png")
                    img_path_dest = os.path.join(output_dir, "images", f"{safe_name}.png")
                    
                    if os.path.exists(img_path_src):
                        img_path = img_path_src
                    elif os.path.exists(img_path_dest):
                        img_path = img_path_dest
                    else:
                        print(f"  -> [Image not found] for: {comp_name} (Expected path: {safe_name}.png)")
                        img_path = ""
                        
                    try:
                        exp_mz = round(float(row.get('m/z', 0)), 4)
                    except (ValueError, TypeError):
                        exp_mz = 999999.0
                        
                    master_compounds[comp_name] = {
                        'Name': comp_name,
                        'Exp_mz': exp_mz,
                        'Best_Score': score,
                        'Max_Intensity': intensity,
                        'CAS Number': cas_clean,
                        'Chemical Formula': str(row.get('Chemical Formula', '')).strip(),
                        'Function': func if func != "nan" else "",
                        'Adducts': set(),
                        'Matrix_Statuses': set(),
                        'Contaminations': set(),
                        'Consolidated_Polymers': set(),
                        'Structure_Image_Path': img_path,
                        'Samples': {}
                    }
                else:
                    # Update optimal representative values if current run provides a stronger hit
                    if score > master_compounds[comp_name]['Best_Score']:
                        try:
                            master_compounds[comp_name]['Exp_mz'] = round(float(row.get('m/z', 0)), 4)
                            master_compounds[comp_name]['Best_Score'] = score
                        except (ValueError, TypeError):
                            pass
                            
                    if intensity > master_compounds[comp_name]['Max_Intensity']:
                        master_compounds[comp_name]['Max_Intensity'] = intensity
                            
                if adduct_raw and adduct_raw != "nan": 
                    master_compounds[comp_name]['Adducts'].add(adduct_raw)
                    
                if multi_polymer_mode:
                    if mat_status and mat_status not in ["N/A", "N/A (Matrix Unknown)", "/", "-", "nan", ""]: 
                        master_compounds[comp_name]['Consolidated_Polymers'].add(mat_status)
                    if contam and contam not in ["/", "-", "nan", ""]:
                        for p in contam.split(','):
                            p = p.strip()
                            if p: master_compounds[comp_name]['Consolidated_Polymers'].add(p)
                else:
                    if mat_status and mat_status not in ["N/A", "N/A (Matrix Unknown)", "/", "-", "nan", ""]: 
                        master_compounds[comp_name]['Matrix_Statuses'].add(mat_status)
                    if contam and contam not in ["/", "-", "nan", ""]: 
                        master_compounds[comp_name]['Contaminations'].add(contam)
                
                # Flag presence in this specific sample run
                master_compounds[comp_name]['Samples'][sample_name] = {
                    'score': score,
                    'ppm': ppm_val
                }
                
        except Exception as e:
            print(f"  -> Error processing {base_name}: {e}", flush=True)

    # ==========================================
    # STEP 2: GROUPING & SORTING 
    # ==========================================
    print("\nGrouping by m/z and sorting peaks by Max Intensity...", flush=True)
    
    comp_list = list(master_compounds.values())
    comp_list.sort(key=lambda x: x['Exp_mz'])
    
    peaks = []
    current_peak_comps = []
    last_mass = -1.0
    
    for comp in comp_list:
        mass = comp['Exp_mz']
        if current_peak_comps and abs(mass - last_mass) > mz_tol:
            peaks.append(current_peak_comps)
            current_peak_comps = []
            
        current_peak_comps.append(comp)
        last_mass = mass
        
    if current_peak_comps:
        peaks.append(current_peak_comps)
        
    peak_data = []
    for peak_comps in peaks:
        peak_comps.sort(key=lambda x: x['Max_Intensity'], reverse=True) 
        peak_max_int = peak_comps[0]['Max_Intensity']
        peak_data.append({
            'peak_max_int': peak_max_int,
            'comps': peak_comps
        })
        
    peak_data.sort(key=lambda x: x['peak_max_int'], reverse=True)

    # ==========================================
    # STEP 3: EXCEL REPORT GENERATION
    # ==========================================
    aggregated_results = []
    current_peak_id = 1
    
    for pd_dict in peak_data:
        comps = pd_dict['comps']
        
        for i, comp in enumerate(comps):
            is_first_in_peak = (i == 0)
            
            disp_peak_id = f"Peak {current_peak_id}" if is_first_in_peak else ""
            disp_mass = comp['Exp_mz'] if is_first_in_peak and comp['Exp_mz'] != 999999.0 else ""
            disp_int = f"{int(comp['Max_Intensity']):,}".replace(',', ' ')
            
            row_dict = {
                'Peak ID': disp_peak_id,
                'm/z': disp_mass,
                'Peak Max Int.': disp_int,
                'Identified Compound': comp['Name'],
                'CAS Number': comp['CAS Number'],
                'Chemical Formula': comp['Chemical Formula'],
                'Adduct': ", ".join(sorted(comp['Adducts'])),
                'Function': comp['Function']
            }
            
            if multi_polymer_mode:
                poly_list = list(comp['Consolidated_Polymers'])
                row_dict['Polymer Presence (DB)'] = " | ".join(sorted(poly_list)) if poly_list else "/"
            else:
                row_dict['Matrix Status'] = " | ".join(sorted(comp['Matrix_Statuses'])) if comp['Matrix_Statuses'] else "N/A"
                row_dict['Possible Contamination'] = " | ".join(sorted(comp['Contaminations'])) if comp['Contaminations'] else "/"
                
            row_dict['Structure_Path'] = comp['Structure_Image_Path']
            
            for s in sample_names:
                s_data = comp['Samples'].get(s, None)
                if s_data is not None:
                    ppm_str = f" | {s_data['ppm']} ppm" if s_data['ppm'] != "" else ""
                    row_dict[s] = f"✔ : {s_data['score']}{ppm_str}"
                else:
                    row_dict[s] = ""
                    
            aggregated_results.append(row_dict)
            
        current_peak_id += 1
            
    print("\nExporting Master Report...", flush=True)
        
    df_export = pd.DataFrame(aggregated_results)
    if df_export.empty:
        return False, "No valid compounds found across the selected files."
        
    report_path = os.path.join(output_dir, "Master_Screening_Report.xlsx")
    
    try:
        try:
            writer = pd.ExcelWriter(report_path, engine='xlsxwriter', engine_kwargs={'options': {'nan_inf_to_errors': True}})
        except TypeError:
            writer = pd.ExcelWriter(report_path, engine='xlsxwriter', options={'nan_inf_to_errors': True})
            
        df_clean = df_export.drop(['Structure_Path'], axis=1)
        
        if multi_polymer_mode:
            idx_samples_start = df_clean.columns.get_loc('Polymer Presence (DB)') + 1
        else:
            idx_samples_start = df_clean.columns.get_loc('Possible Contamination') + 1
            
        df_clean.insert(idx_samples_start, "Structure", "")
        
        df_clean.to_excel(writer, sheet_name='Cross-Screening', index=False, header=False, startrow=1)
        
        workbook = writer.book
        worksheet = writer.sheets['Cross-Screening']
        
        header_format = workbook.add_format({'bold': True, 'text_wrap': True, 'valign': 'vcenter', 'align': 'center', 'bg_color': '#4F81BD', 'font_color': 'white', 'border': 1, 'font_size': 11})
        cell_format = workbook.add_format({'text_wrap': True, 'valign': 'vcenter', 'align': 'center', 'border': 1, 'border_color': '#D3D3D3', 'font_size': 11})
        peak_id_format = workbook.add_format({'bold': True, 'valign': 'vcenter', 'align': 'center', 'border': 1, 'border_color': '#D3D3D3', 'font_size': 12})
        check_format = workbook.add_format({'bold': True, 'valign': 'vcenter', 'align': 'center', 'font_color': '#006100', 'bg_color': '#C6EFCE', 'border': 1, 'border_color': '#D3D3D3', 'font_size': 11})
        
        for col_num, value in enumerate(df_clean.columns.values): 
            worksheet.write(0, col_num, value, header_format)
            
        for col_num, col_name in enumerate(df_clean.columns):
            if col_name in ['Peak ID', 'm/z']:
                worksheet.set_column(col_num, col_num, 12, cell_format)
            elif col_name == 'Peak Max Int.':
                worksheet.set_column(col_num, col_num, 16, cell_format)
            elif col_name == 'Identified Compound':
                worksheet.set_column(col_num, col_num, 40, cell_format)
            elif col_name in ['CAS Number', 'Chemical Formula', 'Adduct']:
                worksheet.set_column(col_num, col_num, 15, cell_format)
            elif col_name in ['Function', 'Matrix Status', 'Possible Contamination']:
                worksheet.set_column(col_num, col_num, 25, cell_format)
            elif col_name == 'Polymer Presence (DB)':
                worksheet.set_column(col_num, col_num, 35, cell_format)
            elif col_name == 'Structure':
                worksheet.set_column(col_num, col_num, 35, cell_format)
            elif col_name in sample_names:
                worksheet.set_column(col_num, col_num, 22, cell_format)
        
        worksheet.set_default_row(150)
        worksheet.set_row(0, 40)
        
        for row_idx, row in df_export.reset_index(drop=True).iterrows():
            excel_row = row_idx + 1
            
            for col_idx, col_name in enumerate(df_clean.columns):
                val = df_clean.iloc[row_idx][col_name]
                
                if col_name in ['Peak ID', 'm/z'] and val != "":
                    worksheet.write(excel_row, col_idx, val, peak_id_format)
                elif col_name in sample_names and "✔" in str(val):
                    worksheet.write(excel_row, col_idx, val, check_format)
                elif col_name != 'Structure':
                    worksheet.write(excel_row, col_idx, val, cell_format)
                
            struct_p = row.get("Structure_Path")
            if isinstance(struct_p, str) and struct_p != "" and os.path.exists(struct_p):
                idx_struct = df_clean.columns.get_loc('Structure')
                worksheet.insert_image(excel_row, idx_struct, struct_p, {'x_scale': 0.6, 'y_scale': 0.6, 'x_offset': 20, 'y_offset': 10, 'positioning': 1})
                
        writer.close()
        print(f"\nSUCCESS! Master report generated: {report_path}", flush=True)
        return True, report_path
        
    except Exception as e:
        print(f"\n[CRITICAL ERROR] Failed to save Excel file: {e}", flush=True)
        # Throw specific exception to be caught and displayed by the isolated GUI tier
        raise Exception(f"Failed to save the Excel file.\n\nError: {e}\n\nPlease verify that the file is not currently open in Excel!")