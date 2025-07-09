#!/usr/bin/env python3
import argparse
import time
import numpy as np
import h5py
import multiprocessing as mp
from scipy.ndimage import map_coordinates
import os

# Global holder for parameters
global GLOBAL_PARAMS
GLOBAL_PARAMS = {}

# -----------------------------------------------------------------------------
# compute_spectrum.py
# Usage: python compute_spectrum.py --h5 PATH --N 1024 --LN 200 --procs 32
# -----------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute 3D spectra from a series stored in HDF5"
    )
    parser.add_argument('--h5',     type=str,
                        default='/projects/DEIKE/cmartinb/eta/rui_eta_series_kpHs0p3_dt0p2_kp0p157_g9p8.h5',
                        help='Path to input HDF5 file containing dataset "eta"')
    parser.add_argument('--memmap', type=str, default='/tmp/eta_series.dat',
                        help='Path for intermediate memory-mapped file')
    parser.add_argument('--N',      type=int, default=1024,
                        help='Spatial grid size (NxN)')
    parser.add_argument('--LN',     type=int, default=200,
                        help='Number of time samples per subdivision')
    parser.add_argument('--pad',    type=int, default=200,
                        help='Zero-padding in time dimension')
    parser.add_argument('--dt',     type=float, default=0.2,
                        help='Time step between samples (seconds)')
    parser.add_argument('--kp',     type=float, default=0.157,
                        help='Characteristic wavenumber for Hs estimation')
    parser.add_argument('--ak',     type=float, default=0.157,
                        help='Normalization amplitude factor')
    parser.add_argument('--procs',  type=int, default=32,
                        help='Number of parallel worker processes')
    return parser.parse_args()


def setup_memmap(h5_path, memmap_path, dtype=np.float32, chunk=100):
    """Load the HDF5 dataset into a shared memory-mapped file."""
    print("Setting up memmap...")
    with h5py.File(h5_path, 'r') as hf:
        dset = hf['eta']
        shape = dset.shape  # (T, N, N)
        mmap = np.memmap(memmap_path, shape=shape, dtype=dtype, mode='w+')
        for start in range(0, shape[0], chunk):
            end = min(start + chunk, shape[0])
            mmap[start:end] = dset[start:end].astype(dtype)
        mmap.flush()
    print(f"Memmap ready: {memmap_path} with shape {shape}")
    return shape


def init_worker(memmap_path, shape):
    """Initializer for each worker: load memmap and precompute grids."""
    params = GLOBAL_PARAMS
    DataEta = np.memmap(memmap_path, shape=shape, dtype=np.float32, mode='r')
    KP, AK = params['kp'], params['ak']
    LN, N   = params['LN'], params['N']
    PAD_LEN, DT = params['pad'], params['dt']
    hann = np.hanning(LN).astype(np.float32)[:, None, None]
    L0 = 2 * np.pi
    wn = 2 * np.pi * np.fft.fftfreq(N, L0 / N)
    kx = np.fft.fftshift(wn)
    ky = kx.copy()
    theta = np.linspace(0, 2 * np.pi, N, dtype=np.float32)
    k = wn[:N//2]
    k_tile, theta_tile = np.meshgrid(k, theta)
    kxp_tile = k_tile * np.cos(theta_tile)
    kyp_tile = k_tile * np.sin(theta_tile)
    dx = kx[1] - kx[0]
    dy = ky[1] - ky[0]
    col = ((kxp_tile - kx[0]) / dx).ravel()
    row = ((kyp_tile - ky[0]) / dy).ravel()
    interp_coords = np.vstack([row, col])
    dtheta = float(theta[1] - theta[0])
    globals().update({
        'DataEta': DataEta,
        'KP': KP, 'AK': AK,
        'LN': LN, 'N': N,
        'PAD_LEN': PAD_LEN, 'DT': DT,
        'hann': hann,
        'kx': kx, 'ky': ky, 'theta': theta, 'k': k,
        'interp_coords': interp_coords,
        'K_tile': k_tile, 'dtheta': dtheta
    })


def process_subdivision(args_tuple):
    """Compute spectrum for one subdivision and print its time interval."""
    idx, norm_flag = args_tuple
    params = GLOBAL_PARAMS
    start_sec = idx * params['dt']
    end_sec   = (idx + params['LN']) * params['dt']
    print(f"Subdivision {idx}:{idx+params['LN']} -> time [{start_sec:.1f}s - {end_sec:.1f}s]")
    sub_t0 = time.time()
    data = DataEta[idx:idx+params['LN']].astype(np.float32)
    Hs = 4 * KP * np.sqrt(data.var())
    data *= hann
    data = np.pad(data, ((0, PAD_LEN),(0,0),(0,0)), mode='constant')
    spec = np.fft.fftn(data, axes=(0,1,2))
    spec = np.fft.fftshift(spec, axes=(0,1,2))
    P3 = (spec.real**2 + spec.imag**2) / (params['N'] * params['N'])
    amp = np.abs(spec) * np.sqrt(params['LN'] * params['N']**2)
    if norm_flag is None:
        kidx = np.argmin(np.abs(k - 4))
        norm_flag = AK / (amp[0,0,kidx] or 1.0)
    amp *= norm_flag
    omega = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(P3.shape[0], params['dt']))
    Fxy, Fko, Fkt, Ak_arr = [], [], [], []
    for i in range(P3.shape[0]):
        sl = P3[i]
        arr = map_coordinates(sl, interp_coords, order=0).reshape(K_tile.shape)
        Fxy.append(sl)
        Fkt.append(arr)
        Fko.append((arr * K_tile).sum(axis=0) * dtheta)
        Ak_arr.append(np.sqrt((arr**2).sum(axis=0) * dtheta))
    Om, Kg = np.meshgrid(omega, k)
    sub_elapsed = time.time() - sub_t0
    print(f"Finished subdivision {idx}:{idx+params['LN']} in {sub_elapsed:.1f}s")
    return {
        'start_idx': idx,
        'Hs': Hs,
        'norm': norm_flag,
        'omega': omega,
        'Omega': Om,
        'K_grid': Kg,
        'F_xy': Fxy,
        'F_komega': Fko,
        'F_ktheta': Fkt,
        'amplitude': amp,
        'amplitude_k': Ak_arr
    }


def run_spectrum(shape):
    """Run subdivisions in parallel with progress tracking."""
    params = GLOBAL_PARAMS
    init_worker(params['memmap'], shape)
    T = shape[0]
    step = params['LN'] // 2
    idxs = list(range(0, T - params['LN'], step))
    first = process_subdivision((idxs[0], None))
    norm_flag = first['norm']
    tasks = [(i, norm_flag) for i in idxs[1:]]
    results = [first]
    total = len(tasks)
    print(f"Dispatching {total} subdivisions to pool with {params['procs']} workers...")
    with mp.Pool(params['procs'], initializer=init_worker,
                 initargs=(params['memmap'], shape), maxtasksperchild=1) as pool:
        out = pool.map(process_subdivision, tasks, chunksize=1)
    for i, res in enumerate(out, 1):
        print(f"[{i}/{total}] Collected subdivision {res['start_idx']}")
        results.append(res)
    return results


def main():
    args = parse_args()
    GLOBAL_PARAMS.update(vars(args))
    shape = setup_memmap(args.h5, args.memmap)
    t0 = time.time()
    results = run_spectrum(shape)
    print(f"Processing finished in {time.time()-t0:.2f}s across {len(results)} subdivisions.")
    # Memory footprint
    print("Computing memory footprint of F_komega...")
    total_bytes = sum(arr.nbytes for r in results for arr in r['F_komega'])
    print(f"Total F_komega data size: {total_bytes/1e6:.1f} MB")
        # Save only F_komega, omega, and k to HDF5
    h5_out = 'spectrum_results.h5'
    print(f"Saving F_komega, omega, and k to {h5_out}...")
    # Compute k vector once (same across subdivisions)
    N = GLOBAL_PARAMS['N']
    L0 = 2 * np.pi
    wn = 2 * np.pi * np.fft.fftfreq(N, L0 / N)
    k = wn[:N//2]
    # Use omega from first subdivision (same size for all)
    omega = results[0]['omega']
    with h5py.File(h5_out, 'w') as hf:
        hf.create_dataset('k', data=k)
        hf.create_dataset('omega', data=omega)
        grp = hf.create_group('subdivisions')
        for r in results:
            sub = grp.create_group(str(r['start_idx']))
            # Stack F_komega arrays into 2D (omega x k)
            fk2d = np.stack(r['F_komega'])  # shape (n_omega, n_k)
            sub.create_dataset('F_komega', data=fk2d, compression='gzip')
    size_mb = os.path.getsize(h5_out) / 1e6
    print(f"HDF5 file size: {size_mb:.1f} MB")
    print("Done.")

if __name__ == '__main__':
    main()
