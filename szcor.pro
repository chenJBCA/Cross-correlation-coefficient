;====================================================================================
;This programme is to calculate the cross-correlation coefficient of two maps 
;using i)polspice ii)ianafast

;cross-correlation coefficient = cross-spectrum/sqrt(auto-spectrum1 X auto-spectrum2)

;T.CHEN 20/02/2018
;====================================================================================
;-----------FUNCTION DEFINATIONS--------------

Function bin, inidata, step    

orinl = n_elements(inidata) ; number of total multipole
binnl = floor(float(orinl)/float(step)) ; number of bins
bindata = fltarr(orinl) 

;average within each bin
for i = 0, binnl-1 do begin
left = i*step
right = (i+1)*step
if right ge orinl then right = orinl-1
bindata[left:right] = mean(inidata[left:right])
endfor

return, bindata
end


;----------------------------------------MAIN PROGRAMME----------------------------------
pro szcor

;-----------------INPUT PARAMETERS/FILES------------------------

;input maps
SZmap = '/scratch/nas_core/tchen/CMB_crosscor/simulation/sz/mlsz/szmap_2048_5arcmin_mK_smo_nocut_mlymlts_pixsize.fits'
galmap = '/scratch/nas_core/tchen/CMB_crosscor/BOSS/density_contrast/rband/BOSSDR13_density_constrast_rband_mag17-22_5arcmin_2048.fits'
oweight = 'overweight.fits'

;healpix option
nside = 2048L ;map nside
lmax = 2*nside ;maximum multipole for power spectrum calculation
beam = 5. ;beam resolution of the two maps in arcmin

;polspice option
path = '/local4/tchen/lib/idlhealpix/lib_idl/PolSpice/src/spice' ;polspice location
apo_angle =20. ;set to half thetamax 
thetamax = 40.; maximum angle for correlation function  integrals to compute power spectrum
anafast = 1  ;Control if want to compare with anafast spectrum; 1-yes, 0-no.

bindata = 1 ;if want to bin the power spectrum; 1-yes, 0-no.
step = 100. ;multipole bin width 

;output filename
autoszcl = 'SZ_spice.fits' ;auto spectrum name for sz map
autogalcl =  'BOSS_spice.fits' ;auto spectrum name for galaxy map
outxcl = 'cross_spice.fits' ;cross spectrum name for szXgalaxy map
outcoeff = 'coeff_spice.fits' ;output cross-correlation coefficient fits file
outps = 'coeff_plot.ps' ;output plot 

anaszcl = 'SZ_anafast.fits'
anagalcl = 'BOSS_anafast.fits'
anaxcl = 'cross_anafast.fits'
anacoeff = 'coeff_anafast.fits'

;--------------------------X-CORRLEATION--------------------------------
ispice, SZmap, outxcl, mapfile2 = galmap,  weightfile1 = oweight, weightfile2 = oweight, nlmax = lmax-1, binpath = path, fwhm1 = beam, fwhm2 = beam, pixelfile = YES, apodizesigma = apo_angle, thetamax = thetamax ; cross-spectrum
ispice, SZmap, autoszcl, weightfile1 = oweight, nlmax = lmax-1, binpath = path, fwhm1 = beam, pixelfile = YES, apodizesigma = apo_angle, thetamax = thetamax ;SZ auto-spectrum
ispice, galmap, autogalcl, weightfile1 = oweight,nlmax = lmax-1, binpath = path, fwhm1 = beam, pixelfile = YES, apodizesigma = apo_angle,  thetamax = thetamax ;galaxy auto-spectrum

read_spice, outxcl, xcl          ;read in cross-spectrum
read_spice, autoszcl, autosz ;read in sz auto-spectrum
read_spice, autogalcl, autogal ;read in galaxy auto-spectrum
xcoeff = xcl/sqrt(autosz*autogal) ;calculate cross-correlation coefficient 
if bindata eq 1 then xcoeff = bin(xcoeff, step) ;bin data


;-----------------------------COMPARE WITH IANAFAST----------------------------
if anafast eq 1 then begin
   read_fits_map, oweight, overweight
   fsky = n_elements(overweight[where(overweight ge 0.5)])/float(n_elements(overweight)) ;fsky
   bl = gaussbeam(beam, lmax-1) ;beam function
   pl = healpixwindow(nside) ;pixel function

   ianafast, SZmap, xclana,  map2_in = galmap,/ring,  maskfile = overweight, nlmax = lmax-1, simul_type = 1 ;cross-spectrum
   xclana = xclana/(fsky*bl^2*pl^2) ;correct fsky, beam, pixel
   ianafast, SZmap, szclana, maskfile = overweight, /ring, nlmax =lmax-1, simul_type = 1 ;SZ auto-spectrum
   szclana = szclana/(fsky*bl^2*pl^2) ;correct fsky, beam, pixel
   ianafast, galmap, galclana, maskfile = overweight, /ring, nlmax =lmax-1, simul_type = 1;galaxy auto-spectrum
   galclana = galclana/(fsky*bl^2*pl^2);correct fsky, beam, pixel
   xcoeffana = xclana/sqrt(szclana*galclana) ;calculate cross-correlation coefficient
   if bindata eq 1 then xcoeffana = bin(xcoeffana, step) ;bin data
endif


;--------------------PLOT AND WRITE OUT------------------------

;write out cl
mwrfits, xcoeff, outcoeff, /create
if anafast eq 1 then begin
mwrfits, szclana, anaszcl, /create
mwrfits, galclana, anagalcl, /create
mwrfits, xclana, anaxcl, /create
mwrfits, xcoeffana, anacoeff, /create
endif

;plot in terminal 

ll = indgen(n_elements(xcoeff))
zero = fltarr(n_elements(ll))

plot, ll, xcoeff, xtitle = 'Multipole l',ytitle  = 'Cross-correlation coefficient', /nodata
oplot, ll, xcoeff, color = 5
oplot, ll, zero
if anafast eq 1 then begin
   oplot, ll, xcoeffana, color = 4
   items = ['ispice', 'ianafast']
   al_legend, items, psym = [0, 0], colors = [5, 4], /bottom, /right
endif

;save plot into ps file
device, decomposed = 0
entry_device = !d.name
set_plot, 'PS'
device, filename = outps, /color


plot, ll, xcoeff,  xtitle = 'Multipole l',ytitle  = 'Cross-correlation coefficient', /nodata
oplot, ll, xcoeff, color = 5
oplot, ll, zero
if anafast eq 1 then begin
   oplot, ll, xcoeffana, color = 4
   items = ['ispice', 'ianafast']
   al_legend, items, psym = [0, 0], colors = [5, 4], /bottom, /right
endif

device, /close_file
set_plot, entry_device




stop
end
