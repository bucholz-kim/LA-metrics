begin
;;----------
; Low-level Humidity Index (HIlow)
;               _       _                _       _
;              |         |              |         |
;   Hi_low  =  | T - T_d |          -   | T - T_d |
;              |_       _|50mb abg      |_       _|150mb abg
;
;              where T_d = dew p;int temperature [K] ;  T = air temperature [K]
;                    abg = above ground;
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
 name    ="ctl"
 odir1   ="/glade/scratch/sungyoon/WRF/run/Irrigation/wrfout/ctl_jja/"
 files1  =systemfunc("ls "+ odir1 + "wrfout_d02_*_12*")
 ofile1  = addfiles(files1,"r")

;; Set up dimensions
 xdim = 330 -1
 ydim = 195 -1 
 tdim = 92  -1

;; Read variables
 tc   = wrf_user_getvar(ofile1,"tc",-1)
 td   = wrf_user_getvar(ofile1,"td",-1)
 plev = wrf_user_getvar(ofile1,"pressure",-1)
 psfc = wrf_user_getvar(ofile1,"PSFC",-1)
 psfc = psfc/100.				; hPa

;; Initializaion
 dum  = psfc
 dum  = 0.
 tc50 = dum
 td50 = dum
 tc150 = dum
 td150 = dum

;; Above Ground Level (AGL)
 p50  = psfc - 50.
 p150 = psfc - 150.

do nt = 0 , tdim 

do ny = 0, ydim
do nx = 0, xdim

;; Checking Missing data
 if (.not.all(ismissing(plev(nt,:,ny,nx))))
 if (.not.all(ismissing(p50(nt,ny,nx))))
 if (.not.all(ismissing(p150(nt,ny,nx))))

;; Finding AGL
 plev1d   = ndtooned(plev(nt,:,ny,nx))
 p501d    = ndtooned(p50(nt,ny,nx))
 p1501d   = ndtooned(p150(nt,ny,nx))

 lo50     = max(ind(plev1d - p501d .ge. 0.)) 
 up50     = min(ind(plev1d - p501d .lt. 0.)) 

 lo150    = max(ind(plev1d - p1501d .ge. 0.)) 
 up150    = min(ind(plev1d - p1501d .lt. 0.)) 
;print(plev1d)
;print(p501d)
;print(p1501d)

;;interpolation
;; AGL 50 hPa
 x_up     = plev(nt,up50,ny,nx) - p50(nt,ny,nx)
 x_lo     = plev(nt,lo50,ny,nx) - p50(nt,ny,nx)

 y_up     = tc(nt,up50,ny,nx)
 y_lo     = tc(nt,lo50,ny,nx)

 tc50(nt,ny,nx) = y_up - ((y_up - y_lo)/(x_up-x_lo))*x_up


 y_up     = td(nt,up50,ny,nx)
 y_lo     = td(nt,lo50,ny,nx)

 td50(nt,ny,nx) = y_up - ((y_up - y_lo)/(x_up-x_lo))*x_up

;; AGL 150 hPa
 x_up     = plev(nt,up150,ny,nx) - p150(nt,ny,nx)
 x_lo     = plev(nt,lo150,ny,nx) - p150(nt,ny,nx)

 y_up     = tc(nt,up150,ny,nx)
 y_lo     = tc(nt,lo150,ny,nx)

 tc150(nt,ny,nx) = y_up - ((y_up - y_lo)/(x_up-x_lo))*x_up

 y_up     = td(nt,up150,ny,nx)
 y_lo     = td(nt,lo150,ny,nx)

 td150(nt,ny,nx) = y_up - ((y_up - y_lo)/(x_up-x_lo))*x_up

 end if
 end if
 end if

end do		; nx end
end do		; ny end
end do		; nt end


; Low-level Humidity Index 
 hi = (tc50 - td50) + (tc150 - td150)

printVarSummary(hi)

;; Maknig binary
 rm = systemfunc("rm -rf ../"+name+".bin")
 fbindirwrite("./"+name+".bin",hi)

end
