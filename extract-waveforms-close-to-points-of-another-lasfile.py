import sys
import numpy as np
import laspy as lp
import WaveForm as wf
from scipy.spatial import cKDTree as kdtree

basepf = sys.argv[1]
fn_las = basepf + '.las'
fn_wdp = basepf + '.wdp'
fn_las_ext = sys.argv[2]
print('extracting waveforms from %s / %s that are close to points from %s ..' % (fn_las, fn_wdp, fn_las_ext))

print('reading points from %s ..' % fn_las)
fp = lp.file.File(fn_las)
pts = np.c_[fp.x, fp.y, fp.z]
fp.close()

pmin = pts.min(axis = 0)
pts -= pmin

print('reading points from %s ..' % fn_las_ext)
fp = lp.file.File(fn_las_ext)
ext = np.c_[fp.x, fp.y, fp.z] - pmin
fp.close()

print('building kdtree ..')
tr = kdtree(pts)

print('query ..')
d, k = tr.query(ext, n_jobs = -1)
print(d.max(), d.mean(), d.min())

print('create mask ..')
mask = np.ones(pts.shape[0], dtype = 'bool')
mask[k] = False

print('reading waveforms ..')
idx, pts = wf.Read(fn_las, fn_wdp, mask = mask)

print('export to LAS file ..')
wf.ExportLAS('fwf-' + fn_las_ext, pts[:,3], pts)
