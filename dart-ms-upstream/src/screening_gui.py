"""
Graphical User Interface for the Cross-Screening Engine.
Provides input validation and user feedback independently from
the pure mathematical operations performed by the engine.
"""

import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os

from screening_engine import run_cross_screening

def launch_screening_gui():
    """
    Constructs the Tkinter interface for the Cross-Screening dashboard.
    Manages user file selections and gracefully handles engine exceptions.
    """
    root = tk.Tk()
    root.title("Cross-Screening Tool")
    root.geometry("550x320")
    root.configure(padx=20, pady=20)
    
    style = ttk.Style()
    style.theme_use('clam')
    
    selected_files = []
    files_var = tk.StringVar(value="No files selected.")
    out_dir_var = tk.StringVar(value=os.getcwd())
    mz_tol_var = tk.StringVar(value="0.001")
    multi_poly_var = tk.BooleanVar(value=False)
    
    def browse_files():
        nonlocal selected_files
        filenames = filedialog.askopenfilenames(title="Select Excel Result Files", filetypes=(("Excel Files", "*.xlsx"), ("All Files", "*.*")))
        if filenames:
            selected_files = list(filenames)
            files_var.set(f"{len(selected_files)} files selected.")

    def browse_out_dir():
        directory = filedialog.askdirectory(title="Select Output Folder")
        if directory: out_dir_var.set(directory)

    def run():
        if not selected_files:
            messagebox.showerror("Error", "Please select at least one Excel file to screen.")
            return
        out_dir = out_dir_var.get()
        if not os.path.exists(out_dir):
            messagebox.showerror("Error", "Please select a valid Output Folder.")
            return
            
        try:
            mz_tol = float(mz_tol_var.get())
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid number for m/z tolerance (e.g. 0.002).")
            return
            
        run_btn.config(text="Processing...", state=tk.DISABLED)
        root.update()
        
        # Engine execution wrapped in try-catch to display potential file-lock errors gracefully
        try:
            success, msg = run_cross_screening(selected_files, out_dir, mz_tol, multi_poly_var.get())
            if success:
                messagebox.showinfo("Success", f"Cross-Screening complete!\nReport saved at:\n{msg}")
            else:
                messagebox.showwarning("Warning", msg)
        except Exception as e:
            messagebox.showerror("Export Error", str(e))
        finally:
            root.destroy()

    ttk.Label(root, text="Cross-Screening Dashboard", font=("Helvetica", 16, "bold")).grid(row=0, column=0, columnspan=3, pady=(0, 20))
    
    ttk.Label(root, text="Result Files:").grid(row=1, column=0, sticky="w", pady=5)
    ttk.Entry(root, textvariable=files_var, width=40, state="readonly").grid(row=1, column=1, padx=5)
    ttk.Button(root, text="Browse...", command=browse_files).grid(row=1, column=2)

    ttk.Label(root, text="Output Folder:").grid(row=2, column=0, sticky="w", pady=5)
    ttk.Entry(root, textvariable=out_dir_var, width=40).grid(row=2, column=1, padx=5)
    ttk.Button(root, text="Browse...", command=browse_out_dir).grid(row=2, column=2)

    ttk.Label(root, text="m/z Grouping Tolerance (Da):").grid(row=3, column=0, sticky="w", pady=5)
    ttk.Entry(root, textvariable=mz_tol_var, width=15).grid(row=3, column=1, sticky="w", padx=5)

    ttk.Checkbutton(root, text="Multi-Polymer Mode", variable=multi_poly_var).grid(row=4, column=0, columnspan=2, sticky="w", pady=10)

    run_btn = ttk.Button(root, text="Run Cross-Screening", command=run)
    run_btn.grid(row=5, column=0, columnspan=3, pady=25, ipadx=20, ipady=10)

    root.mainloop()

if __name__ == '__main__':
    launch_screening_gui()