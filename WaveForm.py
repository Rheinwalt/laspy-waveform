import sys
import numpy as np
import laspy as lp

def ExportLAS(fn, var, pts, cmap = None):
    '''
    Export point cloud to RGB colored LAS file.

    '''
    import laspy as lp
    if cmap is None:
        from matplotlib.cm import magma_r
        cmap = magma_r

    v = var.astype('float') - var.min()
    v /= v.max()
    rgb = cmap(v)
    rgb = rgb[:, :3]
    rgb *= 65535
    rgb = rgb.astype('uint')
    header = lp.header.Header()
    header.data_format_id = 2
    f = lp.file.File(fn, mode = 'w', header = header)
    f.header.scale = [0.01, 0.01, 0.01]
    f.header.offset = [pts[:,0].min(), pts[:,1].min(), pts[:,2].min()]
    f.x = pts[:, 0]
    f.y = pts[:, 1]
    f.z = pts[:, 2]
    if pts.shape[1] == 4:
        f.intensity = pts[:, 3].astype('int')
    f.set_red(rgb[:, 0])
    f.set_green(rgb[:, 1])
    f.set_blue(rgb[:, 2])
    f.close()

def Read(fn_las, fn_wdp, mask = None):
    '''
    idx, pts = Read(fn_las, fn_wdp, mask = None)

    Reads a LAS file together with its corresponding WDP file.

    Parameters
    ----------
    fn_las : str
        LAS file name.
    fn_wdp : str
        Corresponding WDP file name.
    mask : boolean ndarray, optional
        A boolean array masking those points for which no
        waveform is extracted.

    Returns
    -------
    idx : int ndarry
        Start indices for the individual waveforms.
    pts : float ndarry
        Coordinates and amplitudes for all waveforms.

    Notes
    -----
    We are assuming non-refracted, standard straight-line pulses.
    '''
    fp = lp.file.File(fn_las)

    # read waveform packet descriptors (WPD) from
    # variable length records (VLR)
    wpds = []
    vlrs = fp.header.vlrs
    for vlr in vlrs:
        if 99 < vlr.record_id and vlr.record_id < 355:
            wpds.append(vlr.parsed_body)
    wpds = np.array(wpds)

    # get valid WPDs for points in the LAS file
    # (some points might not have a waveform)
    wf_wpd = fp.wave_packet_desc_index
    valid = (0 < wf_wpd) * (wf_wpd < 256)
    if mask is not None:
        valid *= ~mask
    wf_wpd = wf_wpd[valid] - 1

    # read anchor points, delta vectors, byte offsets for
    # the WDP file and corresponding return point locations
    wf_xyz = np.c_[fp.x, fp.y, fp.z][valid]
    wf_vec = np.c_[fp.x_t, fp.y_t, fp.z_t][valid]
    wf_off = fp.byte_offset_to_waveform_data[valid]
    wf_loc = fp.return_point_waveform_loc[valid]

    # get only unique offsets
    wf_off, i = np.unique(wf_off, return_index = True)
    wf_xyz = wf_xyz[i]
    wf_vec = wf_vec[i]
    wf_loc = wf_loc[i]
    wf_wpd = wf_wpd[i]
    wf_num = len(wf_wpd)

    # pre-compute pulse lengths
    wf_len = np.zeros(wf_num, dtype = 'int')
    for i in range(wf_num):
        wf_len[i] = int(wpds[wf_wpd[i], 2])

    # read amplitudes from waveform data packet (WDP) file
    # and compute sample coordinates
    pts = np.zeros((wf_len.sum(), 4))
    with open(fn_wdp, 'rb') as bf:
        i = 0
        for j, off in enumerate(wf_off):
            bf.seek(off)
            l = wf_len[j]
            b = bf.read(l+l)
            a = np.frombuffer(b, dtype = np.uint8, count = l+l)
            sampling, gain, offset = wpds[wf_wpd[j], 3], wpds[wf_wpd[j], 4], wpds[wf_wpd[j], 5]
            pts[i:i+l, 3] = (a.reshape(l, 2) * np.array([1, 255])).sum(axis = 1) * gain + offset
            pts[i:i+l,:3] = wf_xyz[j,:] + np.outer(wf_loc[j] - np.arange(l) * sampling, wf_vec[j,:])
            i += l

    idx = np.append([0,], np.cumsum(wf_len))
    return (idx, pts)

if __name__ == '__main__':
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as pl

    # example files
    pfix = '100429_152240_2535pt_UTM'
    fn_las = pfix + '.las'
    fn_wdp = pfix + '.wdp'

    # read entire files
    idx, pts = Read(fn_las, fn_wdp)
    ExportLAS('fwf-' + fn_las, pts[:,3], pts)

    # make a 3d plot with waveforms 5, 6 and 7
    pts = pts[idx[5]:idx[8], :]
    x, y, z, a = pts[:,0], pts[:,1], pts[:,2], pts[:,3]
    fg = pl.figure(1, (8, 6))
    ax = fg.add_subplot(111, projection = '3d')
    im = ax.scatter(x, y, z, c = a, cmap = pl.cm.magma_r)
    cb = fg.colorbar(im, ax = ax)
    cb.set_label('Waveform amplitude')
    ax.set_xlabel('UTM X [m]')
    ax.set_ylabel('UTM Y [m]')
    ax.set_zlabel('Elevation [m]')

    pl.tight_layout()
    pl.show()
