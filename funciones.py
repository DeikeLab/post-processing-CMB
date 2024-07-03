"""Compute the Reynolds stress"""

from scipy import *
from prepare import load_object, save_object
from defs import Case, Interface2D
from phase import extract_phase
import numpy as np
import scipy
from prepare import field
from defs import Case, Interface2D
import gc
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

import pandas as pd
import numpy as np
import math as m 
from matplotlib import pyplot as plt
from scipy.signal import butter,filtfilt , hilbert

import sys, os
sys.path.append('/projects/DEIKE/jiarongw/jiarongw-postprocessing/jupyter_notebook/functions/')
sys.path.append('/projects/DEIKE/jiarongw/jiarongw-postprocessing/jupyter_notebook/project_specific/windwave/')
# sys.path.append('/home/jiarong/research/postprocessing/jupyter_notebook/functions/')
from tqdm import tqdm
"""Use customized plotting theme"""
import matplotlib as mpl
plt.style.use('/projects/DEIKE/jiarongw/jiarongw-postprocessing/media/matplotlib/stylelib/pof.mplstyle')



def field (case, time, PRE=False, PLOT=False):
    """Put the field slices together and compute some turbulent statistics.
    Pickle the fileds if it's the first time read-in.
    
    Args:
    case: the case instance providing metadata
    time: the times where we want to compute the field statistics
    
    Returns:
    Write the following fields as case attributes.
    
    TODO: add coordinate transformation.
    
    """
    NGRID = 512; NSLICE = 256; L0 = 2*np.pi
    working_dir = case.path
    case.tstat = time # time that the statistics are computed
    case.eta =[]; case.eta_shift = [] # 1D eta (shift or not)
    case.ux = []; case.uy = []; case.f = [] # 2D ux uy # 2D ux uy
    case.re_stress = []
    """Instead of recording shifted fields, record shifting index
    Then do field_shift = np.roll(field, -idx) in x direction when needed."""
    case.shift_index = []
    case.ux_shift = []; case.uy_shift = []; case.f_shift = [] # 2D ux uy shifted, leave blank
    case.re_stress_shift = [] # 2D re shifted
    case.uxmean = [] # 1D ux mean
    case.ux_center = []; case.uy_center = []; case.f_center = [] # The center non-averaged slice to show turbulence
    case.ux_yzcrest = []; case.uy_yzcrest = []; case.f_yzcrest = [] # The y-z plane slice at crest x
    case.ux_yztrough = []; case.uy_yztrough = []; case.f_yztrough = [] # The y-z plane slice at trough x
    
    for t in tqdm(time):
        # Read in eta (utilize the Interface class)
#         filename = working_dir + 'eta/eta_t%g' % 59
#         snapshot = pd.read_table(filename, delimiter = ',')
#         eta_raw = {'x':np.array(snapshot.x), 'z':np.array(snapshot.z), 'eta':np.array(snapshot.pos)}
#         xarray = np.linspace(-L0/2.,L0/2.,NGRID,endpoint=False)+L0/2**NGRID/2 
#         zarray = np.linspace(-L0/2.,L0/2.,NGRID,endpoint=False)+L0/2**NGRID/2
#         x_tile, z_tile = np.meshgrid(xarray,zarray)
#         eta_tile = griddata((eta_raw['x'].ravel(), eta_raw['z'].ravel()), eta_raw['eta'].ravel(), 
#                                  (x_tile, z_tile), method='nearest')
#         eta_1D = np.average(eta_tile, axis=0)
        """For precursor cases just read the fixed shape"""
        """For moving wave cases just read the fixed shape"""
        if (PRE == True): 
            interface = Interface2D(L0 = 2*np.pi, N = 512, path = case.path, 
                                    pre='eta/eta_pre', t = None, PRUNING=True, filename=case.path+'eta/eta_pre')
            eta_1D = np.average(interface.eta, axis=0)        
        else: 
            interface = Interface2D(L0 = 2*np.pi, N = 512, path = case.path, 
                                    pre='eta/eta_loc_t', t = t, PRUNING=True)
            eta_1D = np.average(interface.eta, axis=0)
        # Filter the data (subtract the mean)
        eta_1D_filtered = butter_lowpass_filter(eta_1D-np.average(eta_1D))
        analytic_signal = hilbert(eta_1D_filtered)
        phase = np.angle(analytic_signal)
        # Shift the velocity field along x axis so that phase starts at 0
        idx = (np.abs(phase - 0)).argmin()
#         eta_1D_shift = np.roll(eta_1D, -idx)
        ux_3D = {'name':'ux', 'value':[]} # axis0 in z, axis1 in x, axis2 in y  (in the code)
        uy_3D = {'name':'uy', 'value':[]}
        f_3D = {'name':'f', 'value':[]}

        # Read in the fields either from pickle or from slice data
        for field in (ux_3D,uy_3D,f_3D):         
            """NOTICE: to accomodate different pickle versions"""
            picklename = working_dir + 'field/' + 'pickle_tiger/' + field['name']+'_t%g' % t +'.pkl'
#             picklename = working_dir + 'field/' + 'pickle_desktop/' + field['name']+'_t%g' % t +'.pkl'
            exists = os.path.exists(picklename)
            # If the pickle is there read in the pickles
            if exists:
                field['value'] = load_object(picklename)
                print('pickle restored!')
            # If no pickle read in from the slice files and pickle dump
            if not exists:
                for sn in range (0, NSLICE-1):
                    filename = working_dir + 'field/'+field['name']+'_t%g_slice%g' % (t,sn)
                    print(filename)
                    snapshot = np.loadtxt(filename, dtype = np.str, delimiter='\t')
                    snapshot.reshape([NGRID,NGRID+1])
                    field['value'].append(snapshot[:,0:NGRID].astype(np.float))
                field['value'] = np.array(field['value'])
                save_object(field['value'], picklename)
                
        # Compute Reynolds stress (not shifted)
        ux_mean = np.tile(np.average(ux_3D['value'], axis=(0,1)), (ux_3D['value'].shape[1], 1)) 
        uy_mean = np.tile(np.average(uy_3D['value'], axis=(0,1)), (uy_3D['value'].shape[1], 1))
        re_stress_3D = (ux_3D['value']-ux_mean)*(uy_3D['value']-uy_mean)*(1-f_3D['value'])
        # Compute wave coherent stress (TODO!)
        
        # Append z direction averaged 2D profile
        case.eta.append(eta_1D)
        case.ux.append(np.average(ux_3D['value'], axis=0))
        case.uy.append(np.average(uy_3D['value'], axis=0))
        case.f.append(np.average(f_3D['value'], axis=0))
        # A few fields for visualization
        # First is the center field
        SAMPLE = int(ux_3D['value'].shape[0]/2)
        print(SAMPLE) # center in z
        case.ux_center.append(ux_3D['value'][SAMPLE].copy())
        case.uy_center.append(uy_3D['value'][SAMPLE].copy())
        case.f_center.append(f_3D['value'][SAMPLE].copy())        
        # Second is at crest and trough (require shifted field)
        # These are only for plotting
        if PLOT == True:
            ux_3D_shift = {'name':'ux_shift', 'value':[]} # axis0 in z, axis1 in x, axis2 in y  (in the code)
            uy_3D_shift = {'name':'uy_shift', 'value':[]}
            f_3D_shift = {'name':'f_shift', 'value':[]}
            for (field,field_shift) in zip((ux_3D,uy_3D,f_3D),(ux_3D_shift,uy_3D_shift,f_3D_shift)):
                field_shift['value'] = np.roll(field['value'], -idx, axis=1)
            case.ux_yzcrest.append(ux_3D_shift['value'][:,256,:].copy()) # z, x, y
            case.uy_yzcrest.append(uy_3D_shift['value'][:,256,:].copy())
            case.f_yzcrest.append(f_3D_shift['value'][:,256,:].copy())
            case.ux_yztrough.append(ux_3D_shift['value'][:,192,:].copy()) # z, x, y
            case.uy_yztrough.append(uy_3D_shift['value'][:,192,:].copy())
            case.f_yztrough.append(f_3D_shift['value'][:,192,:].copy())
        case.shift_index.append(idx)
        case.re_stress.append(np.average(re_stress_3D, axis=0))        
        # Additional 1D profile
        case.uxmean.append(np.average(ux_3D['value'], axis=(0,1)))
        case.yarray = np.linspace(0,case.L0,case.N,endpoint=False)+case.L0/2**case.N/2 # Centered grid for interpolation
        del(ux_3D, uy_3D, f_3D, ux_mean, uy_mean, re_stress_3D)
        gc.collect()
        
def butter_lowpass_filter(data, CUT=3, N=512):
    """A helper function that performs lowpass filtering."""
    T = 1           # Sample Period
    fs = N        # Sample rate, Hz (should be the xarray size)
    cutoff = CUT    # desired cutoff frequency of the filter, Hz
    nyq = 0.5 * fs  # Nyquist Frequency
    order = 4       # sin wave can be approximately represented as quadratic
    n = int(T * fs) # total number of samples
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y
        
    
class Interface():
    '''
    Class for every interface related output. Unstructured grid.   
    '''  
    def __init__(self, L0, N, path, t, PRUNING=True):
        '''
        Input:
            L0, N: The desired output grid number
            working_dir: The case's directory
            t: Time of this eta file.
            PRUNING: If eta is output by multiple processes and have multiple headers
                    (might become obsolete).    
        Output:
            self.eta_tile: 2D interpolated eta
            self.p: 2D interpolated surface pressure
            self.tau: 2D interpolated surface viscous stress
        '''
        self.L0 = L0; self.N = N; self.t = t
        xarray = np.linspace(-self.L0/2.,self.L0/2.,self.N,endpoint=False)+self.L0/2**self.N/2 #Centered grid
        zarray = np.linspace(-self.L0/2.,self.L0/2.,self.N,endpoint=False)+self.L0/2**self.N/2 #size of self.N*self.N
        x_tile, z_tile = np.meshgrid(xarray,zarray)
        filename = path + 'eta/eta_loc_t%g' %self.t
        snapshot = pd.read_table(filename, delimiter = ',')
        if PRUNING:
            snapshot = snapshot[snapshot.x != 'x']
            snapshot = snapshot.astype('float')
            snapshot = snapshot[snapshot.pos < 1 + 0.4/4] # Exclude data over slope 0.4
        # Interpolate over x-z plane
        xdata = np.array(snapshot.x, dtype=float); zdata = np.array(snapshot.z, dtype=float); 
        etadata = np.array(snapshot.pos, dtype=float); 
        del (snapshot); gc.collect()        
        self.eta_tile = griddata((xdata.ravel(), zdata.ravel()), etadata.ravel(), (x_tile, z_tile), method='nearest')

        
from defs import *
def eta_series(case, nframe, tstart, dt=1, PRUNING=True):
        '''
        This function reads in a eta time series and create a Eta object for each time.
        Input:
            nframe: Number of total frames.
            tstart: The starting time.
            dt: Time interval between each read-in.
            PRUNING: If eta is output by multiple processes and have multiple headers
            (might become obsolete).
        Output:
            self.energy_t: energy (std(eta)**2) time series (scalar)
            self.interface_t: time series of Interface object       
        '''       
        case.t = np.zeros(nframe)
        case.energy_t = []
        case.interface_t = []
        for i in tqdm (range (0,nframe)):
            case.t[i] = tstart+i*dt
            interface = Interface(case.L0, self.N, self.path, self.t[i], PRUNING=PRUNING)
            case.interface_t.append(interface)
            case.energy_t.append(np.std(interface.eta_tile)**2)      
        
        case.energy_t = np.array(self.energy_t)
        
        
"""Read 2D surface elevation data"""

def read_eta (case,time):
    # Read in data and some cleaning 
    #path = '/projects/DEIKE/jiarongw/turbulence/curved_fixREtau_boundary_REtau720_BO200_g1_ak0.2_MU16_LEVEL10_emax0.3/'
    
    time_str = "%.2f" % time
    
    # Verify if the number end by '.00'
    if time_str.endswith('.00'):
        time_str = str(int(float(time_str)))  # convert to an integer and then to a string
    
    # Resto del código...
    filename = case.path + 'eta/eta_loc_t%s' % time 
    #filename = case.path + 'eta/eta_loc_t%.2f' % time # String formatting to include two decimals
    filename = case.path + 'eta/eta_loc_t%g' %time
    snapshot = pd.read_table(filename, delimiter = ',')
    snapshot = snapshot[snapshot.x != 'x']
    snapshot =+ snapshot.astype('float')
    snapshot = snapshot[snapshot.pos < 1+ 0.4/4]
    snapshot = snapshot.sort_values(by = ['x']) 

    # Interpolate unstructured grid data onto uniform grid, TF work on uniform grid
    from scipy.interpolate import griddata
    L0 = np.pi*2.; N = 512
    xarray = np.linspace(-L0/2.,L0/2.,N,endpoint=False)+L0/2**N/2 #Centered grid
    zarray = np.linspace(-L0/2.,L0/2.,N,endpoint=False)+L0/2**N/2 #size of self.N*self.N
    x_tile, z_tile = np.meshgrid(xarray,zarray)
    xdata = np.array(snapshot.x)
    zdata = np.array(snapshot.z)
    etadata = np.array(snapshot.pos)
    eta_tile = griddata((xdata.ravel(), zdata.ravel()), etadata.ravel(), (x_tile, z_tile), method='nearest')
    return eta_tile



def pol2cart(rho, phi):
    '''Esta funcion pasa de polares a cartesianas'''
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)



variance =[]
integral = []
polar_integral = []


def spectrum_integration(eta, N,L,CHECK=False):
    """ This function performs azimuthal integration of the 2D spectrum (Notice it's 2D instead of 3D with the frequency as well).
        When in doubt, enable CHECK so that the integration is printed out at each step to make sure that 
        units are consistent and we always recover the variance of the data. """  
    #print('N=%g' %N + ', L=%g' %L)
    varr =  np.var(eta)
    
    if CHECK: print('var', np.var(eta))
    
    variance.append(varr)
    if CHECK:print('mean', np.mean(eta))
    
    T0 = 2*np.pi # L0/k = 2pi/4pi
    deltaf = 1/ 2*np.pi
    wavenumber = 2*m.pi*np.fft.fftfreq(N,L/N)
    kx = np.fft.fftshift(wavenumber); ky = kx
    kx_tile, ky_tile = np.meshgrid(kx,ky)
    theta = 2*np.pi*np.arange(-N/2,N/2)/(N)
    
    k = wavenumber[0:int(N/2)]
    dkx = kx[1] - kx[0]; dky = ky[1] - ky[0]
    dk = k[1]-k[0]; dtheta = theta[1]-theta[0]
    dkx = kx[1]-kx[0]; dky = ky[1]-ky[0]
    
    spectrum = np.fft.fft2(eta)/(N*N)**0.5  # FFT normalization 
    F = np.absolute(spectrum)**2/N**2/(dkx*dky) # Per area normalization
    
    if CHECK: print ('sum F', np.sum(F))
    F_center = np.fft.fftshift(F,axes=(0,1)) # Further normalization by independent variables
    
    k_tile, theta_tile = np.meshgrid(k,theta)
    kxp_tile, kyp_tile = pol2cart(k_tile, theta_tile)
    
    integ = np.sum(F_center)*dkx*dky
    if CHECK:print('integral',np.sum(F_center)*dkx*dky)
    
    integral.append(integ)
    
    F_center_polar = scipy.interpolate.griddata((kx_tile.ravel(),ky_tile.ravel()), F_center.ravel(), (kxp_tile, kyp_tile), method='nearest', fill_value=0)
    F_center_polar_integrated = np.sum(F_center_polar*k_tile, axis=0)*dtheta # Azimuthal integration
    
    int_pol = np.sum(F_center_polar_integrated)*dk
    if CHECK: print ('sum polar integrated', np.sum(F_center_polar_integrated)*dk)
    
    polar_integral.append(int_pol)
    
    return k, F_center, F_center_polar_integrated , F_center_polar , k_tile, kxp_tile, kyp_tile , theta_tile , theta , variance, integral, polar_integral , kx, ky

def calculate_parameters(kpHs, uoc, u, kp):
    # Convert kpHs and uoc to their numeric values
    ak = float(kpHs.replace('p', '.'))
    uoc_val = float(uoc.replace('p', '.'))

    # Determine the value of g based on uoc
    if uoc_val == 0.50:
        g = 1
    elif uoc_val == 0.25:
        g = 4
    elif uoc_val == 0.75:
        g = 0.5
    else:
        raise ValueError("uoc value not recognized. Please use 0p25, 0p50, or 0p75.")
    
    c = u / uoc_val
    omegap = np.sqrt(g * kp)
    Re_water = 2.5 * 10**4
    nu_water = c * 2 * m.pi / (4 * Re_water)

    return ak, c, omegap, nu_water, g, uoc_val