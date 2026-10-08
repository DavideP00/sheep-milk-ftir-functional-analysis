"""Generate independent synthetic spectra in Excel Open XML format. Seed: 20261008.
Does not read original data. Uses only Python stdlib + numpy.
"""
from pathlib import Path
import zipfile,html
import numpy as np

OUT=Path(__file__).resolve().parent.parent/'data_synthetic';OUT.mkdir(exist_ok=True)
rng=np.random.default_rng(20261008)
nu=np.linspace(925.92,5011.54,1060)
peaks=np.array([1050,1160,1545,1650,1745,2850,2920,3300]);width=np.array([70,65,90,100,55,80,80,170]);basis=np.exp(-.5*((nu[:,None]-peaks[None,:])/width[None,:])**2)
mean=.015+.025*(nu-nu.min())/(nu.max()-nu.min())+basis@np.array([.09,.12,.14,.22,.26,.10,.12,.06])
headers=['Azienda','Zona','Matricola']+[f'VAR_{j}' for j in range(240,1300)]

def col(n):
    out=''
    while n:
        n,rem=divmod(n-1,26);out=chr(65+rem)+out
    return out

def textcell(pos,value):
    return f'<c r="{pos}" t="inlineStr"><is><t>{html.escape(str(value))}</t></is></c>'

def save(filename,counts,groups,tag,shift):
    path=OUT/filename
    with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=5,allowZip64=True) as z:
        z.writestr('[Content_Types].xml','''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/></Types>''')
        z.writestr('_rels/.rels','''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>''')
        z.writestr('xl/workbook.xml','''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Synthetic FTIR" sheetId="1" r:id="rId1"/></sheets></workbook>''')
        z.writestr('xl/_rels/workbook.xml.rels','''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/></Relationships>''')
        # stream xml directly inside compressed archive to avoid storing a million cells in memory
        with z.open('xl/worksheets/sheet1.xml','w',force_zip64=True) as f:
            def write(s):f.write(s.encode('utf-8'))
            write('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>')
            write('<row r="1">'+''.join(textcell(f"{col(j)}1",v) for j,v in enumerate(headers,1))+'</row>')
            row=2
            farm_offsets=rng.normal(0,.015,size=(len(counts),len(peaks)))
            for farm,(num,group) in enumerate(zip(counts,groups),1):
                effect={'Pianura':-.010,'Collina':0.,'Montagna':.010}[group]
                weights=rng.normal(0,.009,size=(num,len(peaks)))+farm_offsets[farm-1]+effect*np.array([.4,.2,.9,1,1,.4,.8,.1])
                curves=mean[None,:]+shift+weights@basis.T
                err=rng.normal(0,.003,size=(num,1060))
                curves+=(err+np.roll(err,1,axis=1)+np.roll(err,-1,axis=1))/3+rng.normal(0,.004,size=(num,1))
                for i in range(num):
                    metas=[f'SYN_{tag}_{farm:02d}',group,f'SYN_{tag}_{row-1:04d}']
                    meta=''.join(textcell(f'{col(j)}{row}',v) for j,v in enumerate(metas,1))
                    numcells=''.join(f'<c r="{col(j+4)}{row}"><v>{v:.7f}</v></c>' for j,v in enumerate(curves[i]))
                    write(f'<row r="{row}">{meta}{numcells}</row>');row+=1
            write('</sheetData></worksheet>')
    print(filename,'samples',row-2,'bytes',path.stat().st_size,flush=True)

save('FTIR_simulated_spectra_dataset_1.xlsx',[45,44,44,44,47,47,46,48,48,47],['Pianura']*4+['Collina']*3+['Montagna']*3,'SW',0.)
save('FTIR_simulated_spectra_dataset_2.xlsx',[50]*9,['Pianura']*3+['Collina']*3+['Montagna']*3,'VDB',.012)
