# TCGA Breast Cancer Subtype Data Fix

## Problem
The current `get_subtype_info()` function fails with "Cannot find patient ID column in clinical data" because it's using BCR XML clinical data format which has complex nested structure and inconsistent column names.

## Solution Overview
Replace the current approach with `GDCquery_clinic()` which provides clean, pre-processed clinical data with standardized column names. Also fix the patient/sample barcode matching issue.

## Tasks

### ✅ 1. Create tasks/todo.md file following CLAUDE.md workflow
- [x] Set up todo file structure

### 🔄 2. Replace get_subtype_info() function to use GDCquery_clinic() instead of BCR XML
- [ ] Replace lines 413-484 in 06_multiomics_download.R
- [ ] Switch from `GDCquery()` + BCR XML to `GDCquery_clinic()`
- [ ] Use standardized column names: `bcr_patient_barcode` and `paper_BRCA_Subtype_PAM50`
- [ ] Remove complex XML parsing logic

### ⏳ 3. Update create_train_test_sets() function to handle patient/sample barcode mapping
- [ ] Extract patient barcodes from sample barcodes using `substr(..., 1, 12)`
- [ ] Match subtypes at patient level, then map back to samples
- [ ] Handle the patient-to-sample mapping correctly

### ⏳ 4. Test the fix by running the multiomics pipeline
- [ ] Run the updated pipeline
- [ ] Verify subtype data downloads successfully
- [ ] Check patient/sample barcode matching works

### ⏳ 5. Add review section to todo.md with summary of changes
- [ ] Document all changes made
- [ ] Include any lessons learned

## Expected Outcome
- Clinical subtype data will download successfully using the indexed approach
- Patient/sample barcode matching will work correctly
- Pipeline will proceed without the current error