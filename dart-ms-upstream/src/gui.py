"""
Graphical User Interface for the DART-MS Data Processor.
Built using Tkinter, providing a user-friendly way to supply inputs
to the underlying processing engine without exposing the codebase.
"""

import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os

# Separation of concerns: Import the processing logic from the engine file
from processing_engine import process_dart_ms

def launch_gui():
    """
    Constructs and executes the Tkinter main loop for the Data Processor.
    """
    root = tk.Tk()
    root.title("DART-MS Data Processor")
    root.geometry("600x680")
    root.configure(padx=20, pady=20)
    
    style = ttk.Style()
    style.theme_use('clam')
    
    # State variables
    selected_specs = []  
    spec_file_var = tk.StringVar()
    db_file_var = tk.StringVar(value="database.xlsx")
    out_dir_var = tk.StringVar(value=os.getcwd())
    
    all_matrices = ["Unknown", "PET", "PP", "PE", "HDPE", "LDPE", "PE_combined", "PVC", "PUR", "PA", "EVA", "PS", "EPS", "HIPS", "PS_combined"]
    matrix_var = tk.StringVar(value="PET")
    
    contam_options = ["PE", "PS", "PP", "PET", "PVC", "PUR", "PA", "EVA", "All Known Plastics"]
    contam_vars = {opt: tk.BooleanVar(value=False) for opt in contam_options}
    
    min_int_var = tk.StringVar(value="5000")
    max_ppm_var = tk.StringVar(value="5")
    max_cand_var = tk.StringVar(value="5")
    sort_by_var = tk.StringVar(value="Intensity") 
    
    # Export options specifically designed to save computation time
    export_png_var = tk.BooleanVar(value=True)
    export_mol_var = tk.BooleanVar(value=False)

    def browse_spec():
        nonlocal selected_specs
        filenames = filedialog.askopenfilenames(title="Select Spectrum Files", filetypes=(("Excel/CSV Files", "*.xlsx *.xls *.csv"), ("All Files", "*.*")))
        if filenames:
            selected_specs = list(filenames)
            if len(filenames) == 1:
                spec_file_var.set(os.path.basename(filenames[0]))
            else:
                spec_file_var.set(f"{len(filenames)} files selected")

    def browse_db():
        filename = filedialog.askopenfilename(title="Select Database File", filetypes=(("Excel Files", "*.xlsx *.xls"), ("All Files", "*.*")))
        if filename: db_file_var.set(filename)

    def browse_out_dir():
        directory = filedialog.askdirectory(title="Select Output Folder")
        if directory: out_dir_var.set(directory)

    def run_analysis():
        """
        Validates user inputs and triggers the processing engine
        for each selected spectral file. Catches and displays potential errors.
        """
        if not selected_specs:
            messagebox.showerror("Error", "Please select at least one Spectrum file.")
            return
            
        db = db_file_var.get()
        out_dir = out_dir_var.get()
        
        if not db or not os.path.exists(db):
            messagebox.showerror("Error", "Please select a valid Database file.")
            return
        if not os.path.exists(out_dir):
            messagebox.showerror("Error", "Please select a valid Output Folder.")
            return
            
        try:
            min_int = int(min_int_var.get())
            ppm = float(max_ppm_var.get())
            max_cand = int(max_cand_var.get())
        except ValueError:
            messagebox.showerror("Error", "Please enter valid numbers for Intensity, PPM, and Max Candidates.")
            return
            
        matrix = matrix_var.get()
        sort_by = sort_by_var.get()
        selected_contam = [opt for opt, var in contam_vars.items() if var.get()]
        
        run_btn.config(text="Processing...", state=tk.DISABLED)
        root.update()
        
        for spec in selected_specs:
            print(f"\n========================================")
            print(f"PROCESSING: {os.path.basename(spec)}")
            print(f"========================================")
            # Execute processing logic isolated from GUI
            process_dart_ms(
                spec, db, matrix, min_int, ppm, out_dir, 
                max_candidates=max_cand, 
                contamination_targets=selected_contam, 
                sort_by=sort_by,
                generate_png=export_png_var.get(),
                generate_mol=export_mol_var.get()
            )

        messagebox.showinfo("Success", f"Batch processing complete for {len(selected_specs)} files!\nResults are saved in:\n{out_dir}")
        root.destroy()

    ttk.Label(root, text="DART-MS Auto-Processor", font=("Helvetica", 16, "bold")).grid(row=0, column=0, columnspan=3, pady=(0, 20))
    
    ttk.Label(root, text="Spectrum File(s):").grid(row=1, column=0, sticky="w", pady=5)
    ttk.Entry(root, textvariable=spec_file_var, width=40).grid(row=1, column=1, padx=5)
    ttk.Button(root, text="Browse...", command=browse_spec).grid(row=1, column=2)

    ttk.Label(root, text="Database File:").grid(row=2, column=0, sticky="w", pady=5)
    ttk.Entry(root, textvariable=db_file_var, width=40).grid(row=2, column=1, padx=5)
    ttk.Button(root, text="Browse...", command=browse_db).grid(row=2, column=2)

    ttk.Label(root, text="Output Folder:").grid(row=3, column=0, sticky="w", pady=5)
    ttk.Entry(root, textvariable=out_dir_var, width=40).grid(row=3, column=1, padx=5)
    ttk.Button(root, text="Browse...", command=browse_out_dir).grid(row=3, column=2)

    ttk.Separator(root, orient='horizontal').grid(row=4, column=0, columnspan=3, sticky="ew", pady=15)

    ttk.Label(root, text="Sample Matrix:").grid(row=5, column=0, sticky="w", pady=5)
    matrix_cb = ttk.Combobox(root, textvariable=matrix_var, values=all_matrices, state="readonly", width=18)
    matrix_cb.grid(row=5, column=1, sticky="w", padx=5)

    ttk.Label(root, text="Contamination from:").grid(row=6, column=0, sticky="nw", pady=5)
    contam_frame = ttk.Frame(root)
    contam_frame.grid(row=6, column=1, columnspan=2, sticky="w")
    
    row_idx, col_idx = 0, 0
    for opt in contam_options:
        ttk.Checkbutton(contam_frame, text=opt, variable=contam_vars[opt]).grid(row=row_idx, column=col_idx, sticky="w", padx=5, pady=2)
        col_idx += 1
        if col_idx > 2: 
            col_idx = 0
            row_idx += 1

    ttk.Label(root, text="Min Intensity:").grid(row=7, column=0, sticky="w", pady=5)
    ttk.Entry(root, textvariable=min_int_var, width=21).grid(row=7, column=1, sticky="w", padx=5)

    ttk.Label(root, text="Max Error (ppm):").grid(row=8, column=0, sticky="w", pady=5)
    ttk.Entry(root, textvariable=max_ppm_var, width=21).grid(row=8, column=1, sticky="w", padx=5)

    ttk.Label(root, text="Max Candidates:").grid(row=9, column=0, sticky="w", pady=5)
    ttk.Entry(root, textvariable=max_cand_var, width=21).grid(row=9, column=1, sticky="w", padx=5)

    ttk.Label(root, text="Sort Results By:").grid(row=10, column=0, sticky="w", pady=5)
    sort_cb = ttk.Combobox(root, textvariable=sort_by_var, values=["m/z", "Intensity"], state="readonly", width=18)
    sort_cb.grid(row=10, column=1, sticky="w", padx=5)

    ttk.Label(root, text="Export Files:").grid(row=11, column=0, sticky="w", pady=5)
    export_frame = ttk.Frame(root)
    export_frame.grid(row=11, column=1, columnspan=2, sticky="w")
    
    ttk.Checkbutton(export_frame, text="Images (.png)", variable=export_png_var).pack(side=tk.LEFT, padx=(0, 10))
    ttk.Checkbutton(export_frame, text="Mol Files (.mol)", variable=export_mol_var).pack(side=tk.LEFT)

    run_btn = ttk.Button(root, text="Run Analysis", command=run_analysis)
    run_btn.grid(row=12, column=0, columnspan=3, pady=30, ipadx=20, ipady=10)

    root.mainloop()

if __name__ == '__main__':
    launch_gui()