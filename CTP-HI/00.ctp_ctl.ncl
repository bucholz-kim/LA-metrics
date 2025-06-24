begin
;;----------
; Convective Triggering Potential (CTP)
;   CTP     =  integral of curve between the moist adiabat and environmental lapse
;              rate from 100 mb above the ground to 300mb above the ground
;
;  References:  ** Original Method **
;               Findell, K.L., Eltahir, E.A.B. 2003: Atmospheric Controls on 
;               Soil Moisture-Boundary Layer Interactions. Part I: Framework
;               Development. Journal of Hydrometeorology
;
;               Findell, K.L., Eltahir, E.A.B. 2003: Atmospheric Controls on
;               Soil Moisture-Boundary Layer Interactions. Part II: Feedbacks
;               within the Continental United States. Journal of Hydrometeorology
;
;               ** Application and Evaluation **
;               Craig R. Ferguson and Eric F. Wood, 2011: Observed Land–Atmosphere
;               Coupling from Satellite Remote Sensing and Reanalysis. J. Hydrometeor,
;
; Author and Revision History: 
; Original C code       -- Craig Ferguson 
; Converted to F90 code -- Joshua Roundy on Apr 2015
; Modified by           -- A.B. Tawfik on June 2015
; Converted to NCL code -- Sungyoon Kim on Mar 2023
;;----------

;; Read Files
 name = "ctl"
 odir1   ="/glade/scratch/sungyoon/WRF/run/Irrigation/wrfout/ctl_jja/"
 files1  =systemfunc("ls "+ odir1 +"wrfout_d02_*_12*")
 ofile1  = addfiles(files1,"r")

;; Set up dimensions
 xdim = 330 -1
 ydim = 195 -1 
 tdim = 92  -1

;; Read variables
 tk      = wrf_user_getvar(ofile1,"tk",-1)
 plev    = wrf_user_getvar(ofile1,"pressure",-1)
 psfc    = wrf_user_getvar(ofile1,"PSFC",-1)
 psfc    = psfc/100.				; hPa
;printVarSummary(tk)

 geopt    = wrf_user_getvar(ofile1,"geopt",-1)
 gh   = geopt
 gh   = 0.
  
 rh   = plev
 rh   = 100.
 qs   = tk
 qs   = 0.

 qs   = mixhum_ptrh (plev, tk, rh, 1)	; saturation-specific humididy
					; qs kg kg-1

 g0   = 9.807   			; gravity, m s-2
 gh   = geopt/g0			; geopotential height, m


;; Initializaion
 dum  = psfc
 dum  = 0.
 tk100 = dum
 qs100 = dum
 gh100 = dum

 tk150 = dum
 qs150 = dum

 tk300 = dum
 gh300 = dum

 MALR  = dum
 tkmid = dum
 qsmid = dum

;; Above Ground Level (AGL)
 p100 = psfc - 100.
 p150 = psfc - 150.
 p300 = psfc - 300.


do nt = 0 , tdim 

do ny = 0, ydim
do nx = 0, xdim

;; Checking Missing data
 if (.not.all(ismissing(plev(nt,:,ny,nx))))
 if (.not.all(ismissing(p100(nt,ny,nx))))
 if (.not.all(ismissing(p150(nt,ny,nx))))
 if (.not.all(ismissing(p300(nt,ny,nx))))

;; Finding AGL
 plev1d   = ndtooned(plev(nt,:,ny,nx))
 p1001d   = ndtooned(p100(nt,ny,nx))
 p1501d   = ndtooned(p150(nt,ny,nx))
 p3001d   = ndtooned(p300(nt,ny,nx))

 lo100    = max(ind(plev1d - p1001d .ge. 0.)) 
 up100    = min(ind(plev1d - p1001d .lt. 0.)) 

 lo150    = max(ind(plev1d - p1501d .ge. 0.)) 
 up150    = min(ind(plev1d - p1501d .lt. 0.)) 

 lo300    = max(ind(plev1d - p3001d .ge. 0.)) 
 up300    = min(ind(plev1d - p3001d .lt. 0.)) 


;;interpolation
;; AGL 100 hPa
 x_up     = plev(nt,up100,ny,nx) - p100(nt,ny,nx)
 x_lo     = plev(nt,lo100,ny,nx) - p100(nt,ny,nx)

 y_up     = tk(nt,up100,ny,nx)
 y_lo     = tk(nt,lo100,ny,nx)

 tk100(nt,ny,nx) = y_up - ((y_up - y_lo)/(x_up-x_lo))*x_up

 y_up     = qs(nt,up100,ny,nx)
 y_lo     = qs(nt,lo100,ny,nx)

 qs100(nt,ny,nx) = y_up - ((y_up - y_lo)/(x_up-x_lo))*x_up

 y_up     = gh(nt,up100,ny,nx)
 y_lo     = gh(nt,lo100,ny,nx)

 gh100(nt,ny,nx) = y_up - ((y_up - y_lo)/(x_up-x_lo))*x_up

;; AGL 150 hPa
 x_up     = plev(nt,up150,ny,nx) - p150(nt,ny,nx)
 x_lo     = plev(nt,lo150,ny,nx) - p150(nt,ny,nx)

 y_up     = tk(nt,up150,ny,nx)
 y_lo     = tk(nt,lo150,ny,nx)

 tk150(nt,ny,nx) = y_up - ((y_up - y_lo)/(x_up-x_lo))*x_up

 y_up     = qs(nt,up150,ny,nx)
 y_lo     = qs(nt,lo150,ny,nx)

 qs150(nt,ny,nx) = y_up - ((y_up - y_lo)/(x_up-x_lo))*x_up

;; AGL 300 hPa
 x_up     = plev(nt,up300,ny,nx) - p300(nt,ny,nx)
 x_lo     = plev(nt,lo300,ny,nx) - p300(nt,ny,nx)

 y_up     = tk(nt,up300,ny,nx)
 y_lo     = tk(nt,lo300,ny,nx)

 tk300(nt,ny,nx) = y_up - ((y_up - y_lo)/(x_up-x_lo))*x_up

 y_up     = gh(nt,up300,ny,nx)
 y_lo     = gh(nt,lo300,ny,nx)

 gh300(nt,ny,nx) = y_up - ((y_up - y_lo)/(x_up-x_lo))*x_up

 end if
 end if
 end if
 end if

end do 		; nx end
end do		; ny end

end do		; nt end

;; Moisture adiatatic lapse rate (MALR)

 gamma_d  = 9.8 				; dry adiabatic lapse rate, K km-1
 Lv       = 2500000.				; latent heat of vaporization, J kg-1
 Cp       = 1005.7				; specific heat of dry air, J K-1 kg-1
 Rd       = 287. 
 Rv       = 461.5

 tkmid    = 0.5 * (tk100 + tk150)
 qsmid    = 0.5 * (qs100 + qs150)

 MALR     = -gamma_d * ( ( 1. + (Lv   * qsmid)/(   Rd * tkmid))/   \
                         ( 1. + (Lv^2 * qsmid)/(Cp*Rv * tkmid^2)))

 dz       = gh300 - gh100			; thickness of layer, m

;; Making segments for interpolation from 100 to 300 hPa AGL
 nseg = 20

 dz_seg = dz/nseg
 tk_seg = (tk300-tk100)/nseg

 tpar_old = tk100
 tenv_old = tk100

 tpar     = dum
 tenv     = dum
 CTP      = dum

do i = 0,nseg-1

; Getting parcel temperature ,tpar, and env. temperature, tenv
 tpar     = tpar_old + MALR*dz_seg/1000.
 tenv     = tenv_old + tk_seg

 tpar_mid = 0.5 * (tpar + tpar_old)
 tenv_mid = 0.5 * (tenv + tenv_old)

 CTP      = CTP + g0*((tpar_mid - tenv_mid)/tenv_mid)*dz_seg

 tpar_old = tpar
 tenv_old = tenv

end do

printVarSummary(CTP)

;; Maknig binary
 rm = systemfunc("rm -rf ../"+name+".bin")
 fbindirwrite("../"+name+".bin",CTP)

end
