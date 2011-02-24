---------------------------------------------
TEST 1 Steady state
----------------------------------------------
1-0-1234 ra = 0   N=16   46x90L26   30 days
1-0-1234 ra = 0   N=32   91x180L26  30 days
1-0-1234 ra = 0   N=64  181x360L26  30 days
1-0-1234 ra = 0   N=128 361x720L26  30 days
1-3-1234 ra = 45  N=64  181x360L26  30 days
1-6-1234 ra = 90  N=64  181x360L26  30 days

NetCDF output: PS,U,V,W,P,R,T,PHIS,Z3,q1,q2,q3,q4 daily

Plot PS pressure field at days 5, 10, 20, 30
Plot T850
Plot total mass, total energy, kinetic energy
Print q4 max & min
Plot q1-q4 integrated masses.  
Plot zonal wind symmetry (2 norms)

---------------------------------------------
TEST 2 Baroclinic wave
----------------------------------------------
2-0-1234 ra = 0   N=16   46x90L26   30 days
2-0-1234 ra = 0   N=32   91x180L26  30 days
2-0-1234 ra = 0   N=64  181x360L26  30 days
2-0-1234 ra = 0   N=128 361x720L26  30 days
2-0-1234 ra = 0   N=128 721x1440L26 15 days
2-3-1234 ra = 45  N=64  181x360L26  30 days
2-6-1234 ra = 90  N=64  181x360L26  30 days

NetCDF output: PS,U,V,W,P,R,T,PHIS,Z3,q1,q2,q3,q4 daily

Plot PS, T850, VORT850 at days 7, 9
Compute: kinetic energy spectra at 700 mb at days 5,20,25,30 (from NetCDF output)
Plot total mass, total energy, kinetic energy time sequence
plot: q1 tracer at 700 mb, 600 mb, 500 mb at days 9 and 15 (equidistant cylindrical)
plot: q2 tracer field near surface at days 9 and 15 (equidistant cylindrical)
plot: q3 tracer field at 850 mb, 500 mb, 300 mb at days 9 and 15 (equidistant cylindrical)
Print q4 max & min
Plot, q1-q4 integrated masses.  

----------------------------------------------
TEST 3 Advection test solid body rotation
----------------------------------------------
3-0-56 ra = 0   N=64   181x360L60   12 days
3-3-56 ra = 45  N=64   181x360L60   12 days
3-6-56 ra = 90  N=64   181x360L60   12 days
3-0-56 ra = 0   N=128  351x720L60   12 days (optional)

Deltaz = 200 m uniform, 60 levels, zbot = 0 m, ztop = 12 km

Plot longitude-height cross section of q5, q6 contours at equator with ra = 0
at 12 days

Plot latitude-longitude cross section of q5 and q6 at z = 4.5 km at 12 days

Output q5, q6 fields in 3D at days 0.75, 1.5, 2.25, 3, 6, 9, 12 days

Compute normalized l1, l2, linf error norms at 12 days

----------------------------------------------
TEST 4 3D Rossby-Haurwitz wave
----------------------------------------------
4-0-0  ra = 0   N=64   181x360L26   30 days

26 levels, zbot = 0 m, ztop must be below 44.307 km

Plot PS,T850,T300,U850,U200, V850, V200, W850,W500,Z500 at days 5, 10, 15, 30

CI = 4 m/s for U,V and .1K for T

CI = 10 mb for PS

Plot time sequence of domain-integrated total energy diff

Try additional simulation without or with reduced diffusion

Plot U,V,T at p = 7.3 mb to see effect of explicit diffusion near model top 
(if applicable)

Check symmetry of wave in N and S hemispheres

Whether and when wave breaks down over 30-day period

----------------------------------------------
TEST 5 Mountain Rossby
----------------------------------------------

5-0-0  N=64   181x360L26   30 days

U700, V700, T700, W700, Z700, T300 

CI 5 m/s for U and V, 3 K for T700, 100 m for Z, 1 K for T300

W (longitude-height at 45 N)

Domain integrated TE

----------------------------------------------
TEST 6 Pure and Inertio-gravity waves
----------------------------------------------
6-0-0  N=64   181x360L20   4 days
6-0-0  N=32    91x180L20   4 days
6-1-0  N=64   181x360L20   4 days
6-1-0  N=32    91x180L20   4 days
6-2-0  N=64   181x360L20   4 days
6-2-0  N=32    91x180L20   4 days
6-3-0  N=64   181x360L20  15 days
6-3-0  N=32    91x180L20  15 days


Notes:

1.  Try with iterative eta routine because upper levels appear to not be in 
geostrophic balance

2.  Try with topo to match eta = 1 surface

3.  
