;*************************************************
; Buxmann_2015_Fig5_alpha.ncl
;************************************************
;
; Plot 
;
;  in a similar way as Buxmann et al 2015 (Fig 5: alpha case)
;
; Author: Josué Bock
;************************************************
begin
;************************************************
; read in netCDF file
;************************************************

; Buys et al file  
  if (True) then
    outd = "/media/josue/BACKUP2/mafalda-backup/data/Joelle/"
    outdir = outd+"alpha/" ; NOK
    outdir = outd+"alphwlHOBr/" ; panel 1 roughly ok
    ;;outdir = outd+"alphwl1HCHO20"
    ;outdir = outd+"alphwlBr/"
    ;outdir = outd+"alphwl1/"
    ;outdir = outd+"alphaloss/"
    ;outdir = outd+"alphwl1HCHO20_nA/"
    ;outdir = outd+"alphaloss2/"
    ;outdir = outd+"alphaloss1/"
    ;outdir = outd+"alphaloss1HCHO/"
    ;outdir = outd+"alphwl10HCHO20/"
    ;outdir = outd+"alphaloss10/"
  end if

  ;GridFile = addfile(outdir+"grid.nc","r")
  MeteoFile = addfile(outdir+"meteo.nc","r")
  GasFile = addfile(outdir+"gas.nc","r")

; define bounds
; -------------
; plot up to ~600m: end at grid index 62 (remember index 1 is skipped (zero-height layer), index 2 middle of layer = 5m, index 62 = 605m)
  gidx = 1
; start at day=0
  tidx = 0
  eidx = 180

; create a time array, that will be used to redefine time dimension of all arrays
  nday = GasFile->lday(tidx:eidx,0,0,0) ; file is (rec=, noz=1, y=1, x=1) ; units day
  nhour = GasFile->lst(tidx:eidx,0,0,0) ; file is (rec=, noz=1, y=1, x=1) ; units hour
  nmin = GasFile->lmin(tidx:eidx,0,0,0) ; file is (rec=, noz=1, y=1, x=1) ; units min
  ;time = (nday-1)*24. + nhour + nmin/60.
  time = nday*24. + nhour + nmin/60.
  ;time = (nday-2)*24. + nhour + nmin/60.
  time = -20. + (nhour-12)*60. + nmin


; import data
; -----------
  ; rho, to convert units
  rho = MeteoFile->rho(tidx:eidx,gidx,0,0) ; file is (rec=, n=2, y=1, x=1) ; units kg/m3
; convert to m3/mol
  M_air=28.96546e-3 ; kg/mol
  xair = M_air / rho

  ; HOBr
  HOBr = GasFile->HOBr(tidx:eidx,gidx,0,0) ; file is (rec=, n=2, y=1, x=1) ; units mol/m3
; convert to pmol/mol = pptv
  HOBr = HOBr * xair * 1e12

  ; HCl
  HCl = GasFile->HCl(tidx:eidx,gidx,0,0) ; file is (rec=, n=2, y=1, x=1) ; units mol/m3
; convert to pmol/mol = pptv
  HCl = HCl * xair * 1e12

  data1 = new((/3,181/),float)
  data1(0,:) = HOBr
  data1(1,:) = HCl

  ; BrO
  BrO = GasFile->BrO(tidx:eidx,gidx,0,0) ; file is (rec=, n=2, y=1, x=1) ; units mol/m3
; convert to pmol/mol = pptv
  BrO = BrO * xair * 1e12

  ; HBr
  HBr = GasFile->HBr(tidx:eidx,gidx,0,0) ; file is (rec=, n=2, y=1, x=1) ; units mol/m3
; convert to pmol/mol = pptv
  HBr = HBr * xair * 1e12

  data2 = new((/2,181/),float)
  data2(0,:) = BrO
  data2(1,:) = HBr

; PPB for these species
  ; O3
  O3 = GasFile->O3(tidx:eidx,gidx,0,0) ; file is (rec=, n=2, y=1, x=1) ; units mol/m3
; convert to pmol/mol
  O3 = O3 * xair * 1e9

  ; HCHO
  HCHO = GasFile->HCHO(tidx:eidx,gidx,0,0) ; file is (rec=, n=2, y=1, x=1) ; units mol/m3
; convert to pmol/mol
  HCHO = HCHO * xair * 1e9

  data3 = new((/2,181/),float)
  data3(0,:) = O3
  data3(1,:) = HCHO

; Br total: import other Br species
  ; Br2
  Br2 = GasFile->Br2(tidx:eidx,gidx,0,0) ; file is (rec=, n=2, y=1, x=1) ; units mol/m3
; convert to pmol/mol = pptv
  Br2 = Br2 * xair * 1e12

  ; BrCl
  BrCl = GasFile->BrCl(tidx:eidx,gidx,0,0) ; file is (rec=, n=2, y=1, x=1) ; units mol/m3
; convert to pmol/mol = pptv
  BrCl = BrCl * xair * 1e12

  ; BrNO2
  BrNO2 = GasFile->BrNO2(tidx:eidx,gidx,0,0) ; file is (rec=, n=2, y=1, x=1) ; units mol/m3
; convert to pmol/mol = pptv
  BrNO2 = BrNO2 * xair * 1e12

  ; BrNO3
  BrNO3 = GasFile->BrNO3(tidx:eidx,gidx,0,0) ; file is (rec=, n=2, y=1, x=1) ; units mol/m3
; convert to pmol/mol = pptv
  BrNO3 = BrNO3 * xair * 1e12

  ; Br total
  data1(2,:) = HOBr+BrO+HBr+Br2+BrCl+BrNO2+BrNO3



;************************************************
; create plots
;************************************************
  wks = gsn_open_wks("pdf","Buxmann_2015_Fig5_alpha")          ; send graphics to PDF file
  plot = new(3,graphic)

  res                              = True
  res@gsnDraw                      = False
  res@gsnFrame                     = False
  res@gsnMaximize                  = True            ; use full page
  res@gsnPaperOrientation          = "portrait"      ; prevent automatic rotation
  ; adjust plot to get same height and width
  res@vpHeightF           = 0.2
  res@vpWidthF            = 0.75

  ; Axis settings
  ;res@trYMinF              = 0.
  ;res@trYMaxF              = 600.
  res@tiXAxisString     = "Time (min)"
  res@tiYAxisString     = "Mixing ratio (ppt)"
  ; Informations displayed on top
  res@gsnLeftString = ""
  res@gsnRightString = ""
  res@tmYLMinorOn = False
  res@tmYMajorGrid = True
  res@tmYMajorGridLineColor = "Grey70"
  res@tmYMajorGridLineDashPattern = 2
  res@trXMinF     = -20.
  res@trXMaxF     = 150.
  res@tmXBMode = "Manual"
  res@tmXBTickStartF = -20.
  res@tmXBTickEndF = 140.
  res@tmXBTickSpacingF = 20.
  res@xyLineThicknessF = 3

  res@trYMaxF     = 2200.
  res@tmYLMode = "Manual"
  res@tmYLTickStartF = 0.
  res@tmYLTickEndF = 2000.
  res@tmYLTickSpacingF = 500.
  res@xyLineColors = (/"royalblue4","limegreen","royalblue4"/)
  res@xyDashPatterns = (/1,8,0/)
  plot(0) = gsn_csm_xy(wks,time,data1,res)

  res@trYMaxF     = 900.
  res@tmYLTickStartF = 0.
  res@tmYLTickEndF = 800.
  res@tmYLTickSpacingF = 200.
  res@xyLineColors := (/"royalblue4","royalblue4"/)
  res@xyDashPatterns := (/1,2/)
  plot(1) = gsn_csm_xy(wks,time,data2,res)

  res@tiYAxisString     = "Mixing ratio (ppb)"
  res@trYMaxF     = 300.
  res@tmYLTickStartF = 0.
  res@tmYLTickEndF = 300.
  res@tmYLTickSpacingF = 50.
  res@xyLineColors = (/"red2","black"/)
  res@xyDashPatterns = (/2,0/)
  plot(2) = gsn_csm_xy(wks,time,data3,res)


;************************************************
; create panel
;************************************************
  gsn_panel(wks,plot,(/3,1/),False)             ; now draw as one plot

end