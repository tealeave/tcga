#!/usr/bin/env python3

import pickle
import os
import re

# Diagnostic script to examine patient ID mapping issues
data_dir = "/dfs6/pub/ddlin/projects/tcga/data/multiomics"

def load_rdata(file_path):
    """Load RData files using rpy2 if available, otherwise return None"""
    try:
        import rpy2.robjects as robjects
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr
        
        # Load the RData file
        robjects.r['load'](file_path)
        
        # Get the loaded objects
        objects = list(robjects.r['ls']())
        results = {}
        
        for obj_name in objects:
            obj = robjects.r[obj_name]
            try:
                # Convert to pandas DataFrame if possible
                with pandas2ri.conversion.localconverter(pandas2ri.converter):
                    if hasattr(obj, 'colnames') and hasattr(obj, 'rownames'):
                        # Matrix-like object
                        import pandas as pd
                        df = pd.DataFrame(obj)
                        df.columns = obj.colnames
                        df.index = obj.rownames
                        results[obj_name] = df
                    else:
                        # Try to convert to pandas
                        results[obj_name] = pandas2ri.ri2py(obj)
            except:
                results[obj_name] = str(obj)
        
        return results
    except ImportError:
        print("rpy2 not available, using basic file inspection")
        return None

def analyze_patient_ids():
    print("=== PATIENT ID MAPPING DIAGNOSTICS ===\n")
    
    files_to_check = [
        "mrna_data.RData",
        "mirna_data.RData", 
        "protein_data.RData",
        "subtype_info.RData"
    ]
    
    diagnostics = {}
    
    for filename in files_to_check:
        filepath = os.path.join(data_dir, filename)
        if os.path.exists(filepath):
            print(f"=== {filename.upper()} ===")
            
            try:
                # Try to load with basic inspection
                import subprocess
                result = subprocess.run(['R', '--vanilla', '--quiet', '-e', 
                    f'load("{filepath}"); print(colnames(get(ls()[1]))[1:10]); print(ls())'], 
                    capture_output=True, text=True)
                
                if result.returncode == 0:
                    output = result.stdout
                    print("R output:", output[:500])
                else:
                    print("Error running R:", result.stderr)
                    
            except Exception as e:
                print(f"Error loading {filename}: {e}")
                
            print()
        else:
            print(f"File not found: {filename}\n")
    
    # Alternative: Check if we can read the files directly
    print("=== ALTERNATIVE ANALYSIS ===")
    
    # Try to understand the structure by checking file sizes and basic patterns
    for filename in files_to_check:
        filepath = os.path.join(data_dir, filename)
        if os.path.exists(filepath):
            size = os.path.getsize(filepath)
            print(f"{filename}: {size/1024/1024:.2f} MB")

if __name__ == "__main__":
    analyze_patient_ids()