"""
Main Entry Point for the DART-MS Analytics Suite.
This script provides a unified launcher interface to navigate between
the Auto-Processor and the Cross-Screening Dashboard tools.
"""

import sys
import os
import tkinter as tk
from tkinter import ttk

# Ensure Python looks for modules in the current 'src' directory
# This is crucial for relative imports when launching the script from different working directories.
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

from gui import launch_gui
from screening_gui import launch_screening_gui

def main_menu():
    """
    Initializes and displays the main launcher menu.
    Allows the user to select and launch specific sub-tools of the suite.
    """
    root = tk.Tk()
    root.title("DART-MS Analytics Suite")
    root.geometry("450x300")
    root.configure(padx=30, pady=30)
    
    # Apply a modern styling theme available across all OS
    style = ttk.Style()
    style.theme_use('clam')
    
    def open_processor():
        """Destroys the launcher and opens the Data Processor GUI."""
        root.destroy()
        launch_gui()
        
    def open_screening():
        """Destroys the launcher and opens the Cross-Screening GUI."""
        root.destroy()
        launch_screening_gui()

    # Title and subtitle
    ttk.Label(root, text="DART-MS Analytics Suite", font=("Helvetica", 18, "bold")).pack(pady=(0, 5))
    ttk.Label(root, text="Please select a tool to launch:", font=("Helvetica", 11), foreground="gray").pack(pady=(0, 30))
    
    # Tool 1 Button: Processor
    btn1 = ttk.Button(root, text="1. DART-MS Auto-Processor\n(Process raw spectra against database)", command=open_processor)
    btn1.pack(fill='x', ipady=15, pady=5)
    
    # Tool 2 Button: Cross-Screening
    btn2 = ttk.Button(root, text="2. Cross-Screening Dashboard\n(Compare compounds across multiple samples)", command=open_screening)
    btn2.pack(fill='x', ipady=15, pady=5)

    root.mainloop()

if __name__ == '__main__':
    print("Starting DART-MS Analytics Suite...")
    main_menu()