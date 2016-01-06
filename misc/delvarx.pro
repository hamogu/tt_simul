;+
; NAME: 
;	DELVARX
; PURPOSE: 
; 	Delete variables for memory management (can call from routines) 
; EXPLANATION:
;	Like intrinsic DELVAR function, but can be used from any calling level
;
; CALLING SEQUENCE:
; 	DELVARX,  a [,b,c,d,e,f,g,h,i,j, /FREE_MEM]
;
; INPUTS: 
;	p0, p1...p9 - variables to delete
;
; OPTIONAL KEYWORD:
;       /FREE_MEM - If set, then free memory associated with pointers 
;                   and objects.   Requires V5.3 or later
; RESTRICTIONS: 
;	Can't use recursively due to EXECUTE function
;
; METHOD: 
;	Uses EXECUTE and TEMPORARY function   
;
; REVISION HISTORY:
;	Copied from the Solar library, written by slf, 25-Feb-1993
;	Added to Astronomy Library,  September 1995
;	Converted to IDL V5.0   W. Landsman   September 1997
;       Modified, 26-Mar-2003, Zarro (EER/GSFC) 26-Mar-2003
;       - added FREE_MEM to free pointer/objects
;	Removed Execute to allow use with IDL virtual machine - replace by brute force - H. M. Guenther 20-Oct-2006
;-

PRO delvarx, p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,free_mem = free_mem
   
;  only delete if defined on inpu (avoids error message)
      defined=n_elements(p0) 
      IF defined GT 0 THEN BEGIN
         if keyword_set(free_mem) then begin
            if !VERSION.RELEASE GE '5.3' then $
                  heap_free,p0
         ENDIF
         p0=0
         dvar=temporary(p0)
      ENDIF
      defined=n_elements(p1) 
      IF defined GT 0 THEN BEGIN
         if keyword_set(free_mem) then begin
            if !VERSION.RELEASE GE '5.3' then $
                  heap_free,p1
         ENDIF
         p1=0
         dvar=temporary(p1)
      ENDIF
      defined=n_elements(p2) 
      IF defined GT 0 THEN BEGIN
         if keyword_set(free_mem) then begin
            if !VERSION.RELEASE GE '5.3' then $
                  heap_free,p2
         ENDIF
         p2=0
         dvar=temporary(p2)
      ENDIF
      defined=n_elements(p3) 
      IF defined GT 0 THEN BEGIN
         if keyword_set(free_mem) then begin
            if !VERSION.RELEASE GE '5.3' then $
                  heap_free,p3
         ENDIF
         p3=0
         dvar=temporary(p3)
      ENDIF
      defined=n_elements(p4) 
      IF defined GT 0 THEN BEGIN
         if keyword_set(free_mem) then begin
            if !VERSION.RELEASE GE '5.3' then $
                  heap_free,p4
         ENDIF
         p4=0
         dvar=temporary(p4)
      ENDIF
      defined=n_elements(p5) 
      IF defined GT 0 THEN BEGIN
         if keyword_set(free_mem) then begin
            if !VERSION.RELEASE GE '5.3' then $
                  heap_free,p5
         ENDIF
         p5=0
         dvar=temporary(p5)
      ENDIF
      defined=n_elements(p6) 
      IF defined GT 0 THEN BEGIN
         if keyword_set(free_mem) then begin
            if !VERSION.RELEASE GE '5.3' then $
                  heap_free,p6
         ENDIF
         p6=0
         dvar=temporary(p6)
      ENDIF
      defined=n_elements(p7) 
      IF defined GT 0 THEN BEGIN
         if keyword_set(free_mem) then begin
            if !VERSION.RELEASE GE '5.3' then $
                  heap_free,p7
         ENDIF
         p70
         dvar=temporary(p7)
      ENDIF
      defined=n_elements(p8) 
      IF defined GT 0 THEN BEGIN
         if keyword_set(free_mem) then begin
            if !VERSION.RELEASE GE '5.3' then $
                  heap_free,p8
         ENDIF
         p8=0
         dvar=temporary(p8)
      ENDIF
      defined=n_elements(p9) 
      IF defined GT 0 THEN BEGIN
         if keyword_set(free_mem) then begin
            if !VERSION.RELEASE GE '5.3' then $
                  heap_free,p9
         ENDIF
         p9=0
         dvar=temporary(p9)
      ENDIF   




   RETURN
END
