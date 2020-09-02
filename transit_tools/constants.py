#jupiter radius in m, au, R_sun
#solar radius in m, au, R_jup
#year in min, sec

#nomenclature for variables:
#   XYYY_ZZZ
#   X = physical parameter measured
#   YYY = object name
#   ZZZ = Unit that the quantity is in
#   Ex: "X of YYY in ZZZ", "Radius of Jupiter in meters"

# Various radii in meters
Rsolar_m = 696342000      #1 solar radius
Rjup_m = 71492000         #1 Jupiter radius
Rearth_m = 6371008        #1 Earth radius

# Various radii in solar radii
Rjup_sol = 0.102763       #1 Jupiter radius

# Various radii in Jupiter radii
Rearth_jup = 0.0892147    #1 Earth radius

# Various radii in Earth radii
Rsolar_ear = 109.076      #1 solar radius

# Various radii in Astronomical Units
Rsolar_au = 0.00465047    #1 solar radius

# Default stellar params for TLS fitting
default_star_params = {'rstar' : 1.0, 'rlow' : 0.67, 'rhigh' : 2.5,
                       'mstar' : 1.0, 'mlow' : 0.9, 'mhigh' : 1.0,
                       'ab' : [0.4804, 0.1867]}
