"""Lightweight structural validation without loading the full Excel arrays."""
from pathlib import Path
import zipfile, hashlib, re, sys
root=Path(__file__).resolve().parent.parent
expected={'FTIR_simulated_spectra_dataset_1.xlsx':460,'FTIR_simulated_spectra_dataset_2.xlsx':450}
for file,n in expected.items():
    p=root/'data_synthetic'/file
    with zipfile.ZipFile(p) as z:
        assert z.testzip() is None
        data=z.read('xl/worksheets/sheet1.xml').decode('utf8')
    assert data.count('<row r=')==n+1,(file,'row count')
    assert all(f'VAR_{i}' in data[:150000] for i in range(240,1300)),file
    assert 'SYN_' in data
    print('OK',file,n,'spectra',p.stat().st_size,'bytes')
for name in ('RUN_ALTITUDE.R','RUN_BREED.R','RUN_OUTLIER.R'):
    s=(root/'scripts'/name).read_text()
    assert 'readRDS(' not in s and 'LOFO' not in s.upper(),name
    assert 'data_synthetic' in s and 'R", "FTIR_FUNCTIONS.R' in s,name
print('PASS: structural checks (NOT full R execution)')
