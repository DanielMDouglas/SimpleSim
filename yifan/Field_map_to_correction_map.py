from array import array
import numpy as np
from ROOT import TFile, TH2F, TH3F, TCanvas, TLegend, gStyle, TGaxis
import math

# Binning
xmin = -0.32
xmax = 0.32
ymin = -0.6
ymax = 0.6
zmin = -0.302768
zmax = 0.302768

v_nominal = 1.5136504594138773 #mm/us
E_nominal = 50 #kV/m
temp = 87.5 # Module0 run

n_xbin = 65
n_ybin = 121
n_zbin = 61

xgap = (xmax - xmin) / (n_xbin-1)
ygap = (ymax - ymin) / (n_ybin-1)
zgap = (zmax - zmin) / (n_zbin-1)
tgap = zgap * 1000 / v_nominal # us

x_grid = np.linspace(xmin, xmax, n_xbin)
y_grid = np.linspace(ymin, ymax, n_ybin)
z_grid = np.linspace(zmin, zmax, n_zbin)
if len(z_grid)%2 == 0:
    z_half_len = int(len(z_grid)/2)
else:
    z_half_len = int((len(z_grid)+1)/2)

z_grid_half = z_grid[:z_half_len]

# Define output TH3F's (spatial correction)
th3_Dx_corr = TH3F("th3_Dx_corr", "th3_Dx_corr", int(n_xbin), xmin-xgap*0.5, xmax+xgap*0.5, int(n_ybin), ymin-ygap*0.5, ymax+ygap*0.5, int(n_zbin), zmin-zgap*0.5, zmax+zgap*0.5)
th3_Dy_corr = TH3F("th3_Dy_corr", "th3_Dy_corr", int(n_xbin), xmin-xgap*0.5, xmax+xgap*0.5, int(n_ybin), ymin-ygap*0.5, ymax+ygap*0.5, int(n_zbin), zmin-zgap*0.5, zmax+zgap*0.5)
th3_Dz_corr = TH3F("th3_Dz_corr", "th3_Dz_corr", int(n_xbin), xmin-xgap*0.5, xmax+xgap*0.5, int(n_ybin), ymin-ygap*0.5, ymax+ygap*0.5, int(n_zbin), zmin-zgap*0.5, zmax+zgap*0.5)

# Define input TH3F's (efield)
f_root_read = TFile(f"Efield_TH3_raw.root","READ")
th3_Ex_raw = f_root_read.Get("th3_Ex_raw")
th3_Ey_raw = f_root_read.Get("th3_Ey_raw")
th3_Ez_raw = f_root_read.Get("th3_Ez_raw")

# Get vdrift according to the efield
def drift_speed_helper(params, e, t=87.17):
        
    p1, p2, p3, p4, p5, p6, t0 = params

    vdrift = (1 + p1 * (t - t0) ) * (p3 * e * math.log(1 + abs(p4) / e) + p5 * np.power(e,p6)) + p2 * (t-t0)

    return vdrift


# vdrift mm/us
# Efield kV/cm
def vdrift(e_field):
    
    # for low eField use mobility, but the parametrization is different than the BNL one
    tdiff = temp - 87.302
    eFit = 0.0938163 - 0.0052563 * tdiff - 0.000146981 * np.power(tdiff,2)
    muFit = 5.183987 + 0.01447761 * tdiff - 0.0034972 * np.power(tdiff,2) - 0.0005162374 * np.power(tdiff,3)

    # parameters for drift speed fit
    # p1, p2, p3, p4, p5, p6, t0
    ICARUS_params = np.array([-0.04640231, 0.0171171, 1.881246, 0.9940772, 0.0117183, 4.202141, 105.7491])
    Walkowiak_params = np.array([-0.01481, -0.0075, 0.141, 12.4, 1.627, 0.317, 90.371])

    # for low eField, vdrift model uses mobility * eField
    if e_field < eFit:
        v_drift = muFit * e_field

    # for intermediate eField, vdrift model uses ICARUS parametrization
    elif e_field < 0.619:
        v_drift = drift_speed_helper(ICARUS_params, e_field, temp)

    # for eField between two model ranges
    elif e_field < 0.699:
        v_drift = (0.699 - e_field) / 0.08 * drift_speed_helper(ICARUS_params, e_field, temp)                              + (e_field - 0.619) / 0.08 * drift_speed_helper(Walkowiak_params, e_field, temp)

    # for high eField, vdrift model uses Walkowiak parametrization
    else:
        v_drift = drift_speed_helper(Walkowiak_params, e_field, temp)
                
    return v_drift

# Translation
for ix, x in enumerate(x_grid):
    for iy, y in enumerate(y_grid):        
        for iz, z in enumerate(z_grid_half):
            if iz == 0:
                this_x = x
                this_y = y
                this_z = z
                this_pos = np.array([this_x, this_y, this_z])
                
            corr_vec = this_pos - np.array([x,y,z])
            corr_vec = corr_vec * 100 # m->cm  
               
            th3_Dx_corr.Fill(x,y,z,corr_vec[0])
            th3_Dy_corr.Fill(x,y,z,corr_vec[1])
            th3_Dz_corr.Fill(x,y,z,corr_vec[2])

            if z != 0:
                th3_Dx_corr.Fill(x,y,-z,corr_vec[0])
                th3_Dy_corr.Fill(x,y,-z,corr_vec[1])
                th3_Dz_corr.Fill(x,y,-z,-corr_vec[2])

            # Calculate the next "this_pos"
            Ex = th3_Ex_raw.Interpolate(this_x,this_y,this_z)
            Ey = th3_Ey_raw.Interpolate(this_x,this_y,this_z)
            Ez = th3_Ez_raw.Interpolate(this_x,this_y,this_z)
            E_v = np.array([Ex, Ey, Ez]) #kV/m
            
            if np.linalg.norm(E_v) == 0:
                E_v_unit = np.array([0, 0, 1]) # negative half of the TPC
                this_v = vdrift(E_nominal * 0.01) * E_v_unit
            else:             
                E_v_unit = E_v / np.linalg.norm(E_v)
                this_v = vdrift(np.linalg.norm(E_v) * 0.01) * E_v_unit
                
            
            step = this_v * tgap # mm/us * us -> mm
            step = step * 0.001 # mm -> m
            this_pos = np.array([this_x, this_y, this_z]) + step
            
            # update this_pos to this_x, this_y and this_z
            this_x = this_pos[0]
            this_y = this_pos[1]
            this_z = this_pos[2]
            
# Write in the outputs
f_Dx_root = TFile(f"Spatial_Corr_TH3.root","RECREATE")
th3_Dx_corr.Write("th3_Dx_corr")
th3_Dy_corr.Write("th3_Dy_corr")
th3_Dz_corr.Write("th3_Dz_corr")



